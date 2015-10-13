from collections import defaultdict
import datetime
import json
import os
import dateutil.parser
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.matproj.snl import StructureNL
from pymongo import MongoClient
import yaml
from rubicon.utils.snl.egsnl import get_meta_from_structure

DATETIME_HANDLER = lambda obj: obj.isoformat() \
    if isinstance(obj, datetime.datetime) else None
YAML_STYLE = False  # False = YAML is formatted as blocks


class SubmissionMongoAdapterEG(object):
    # This is the user interface to submissions

    def __init__(self, host='localhost', port=27017, db='snl', username=None,
                 password=None):
        self.host = host
        self.port = port
        self.db = db
        self.username = username
        self.password = password

        self.connection = MongoClient(host, port, j=False)
        self.database = self.connection[db]
        if self.username:
            self.database.authenticate(username, password)

        self.jobs = self.database.jobs
        self.reactions = self.database.reactions
        self.id_assigner = self.database.id_assigner
        self.reaction_id_assigner = self.database.reaction_id_assigner

        self._update_indices()

    def _reset(self):
        if "prod" in self.database.name:
            print "PROD database is not supposed to reset, please changed the code to reset"
            exit()
        self._restart_id_assigner_at(1)
        self.jobs.remove()

    def _reset_reactions_collection(self):
        if "prod" in self.database.name:
            print "PROD database is not supposed to reset, please changed the code to reset"
            exit()
        self._restart_reaction_id_assigner_at(1)
        self.reactions.remove()

    def _update_indices(self):
        self.jobs.ensure_index('submission_id', unique=True)
        self.jobs.ensure_index('state')
        self.jobs.ensure_index('submitter_email')

        self.reactions.ensure_index('reaction_id', unique=True)
        self.reactions.ensure_index('all_inchis')
        self.reactions.ensure_index('reactant_inchis')
        self.reactions.ensure_index('product_inchis')
        self.reactions.ensure_index('reactant_submission_ids')
        self.reactions.ensure_index('product_submission_ids')
        self.reactions.ensure_index('submitter_email')

    def _get_next_submission_id(self):
        return self.id_assigner.find_and_modify(
            query={}, update={'$inc': {'next_submission_id': 1}})[
                'next_submission_id']

    def _get_next_reaction_id(self):
        return self.reaction_id_assigner.find_and_modify(
            query={}, update={'$inc': {'next_reaction_id': 1}})[
            'next_reaction_id']

    def _restart_id_assigner_at(self, next_submission_id):
        self.id_assigner.remove()
        self.id_assigner.insert({"next_submission_id": next_submission_id})

    def _restart_reaction_id_assigner_at(self, next_reaction_id):
        self.reaction_id_assigner.remove()
        self.reaction_id_assigner.insert({"next_reaction_id": next_reaction_id})

    def submit_snl(self, snl, submitter_email, parameters=None):
        parameters = parameters if parameters else {}

        d = snl.as_dict()
        d['submitter_email'] = submitter_email
        d['parameters'] = parameters
        d['state'] = 'SUBMITTED'
        d['state_details'] = {}
        d['task_dict'] = {}
        d['submission_id'] = self._get_next_submission_id()
        d['submitted_at'] = datetime.datetime.utcnow().isoformat()
        if 'is_valid' not in d:
            d.update(get_meta_from_structure(snl.structure))

        d.update(snl.structure.as_dict())

        self.jobs.insert(d)
        return d['submission_id']

    def submit_reaction(self, reactant_snls, product_snls, reactant_fragments, product_fragments, submitter_email, parameters=None):
        """
            Submit a reaction. This task will be separated to several single point energy calculations, and submitted
            as individual molecule.

            Args:
                reactant_snls: List of tuple(snl, count, nickname).
                product_snls: List of tuple(snl, count, nickname).
                reactant_fragments: BSSE fragments definition. (BSSEFragment)
                product_fragments: BSSE fragments definition. (BSSEFragment)
                submitter_email: Email.
                parameters: dict of parameter. Expected parameters are 1) method: QChem theoretival method. e.g.
                    B3LYP-XDM/6-31+G*; 2) solvent: implicit solvent in energy calcuation. e.g. THF; ...
        """
        reaction_element_count = defaultdict(lambda: 0)
        for snl, n, nick_name in reactant_snls:
            mol = snl.structure
            for site in mol.sites:
                element = site.specie.symbol
                reaction_element_count[element] += n
        product_element_count = defaultdict(lambda: 0)
        for snl, n, nick_name in product_snls:
            mol = snl.structure
            for site in mol.sites:
                element = site.specie.symbol
                product_element_count[element] += n
        if reaction_element_count != product_element_count:
            raise Exception("Number of atoms is inconsistant in reactant and product")
        reactant_inchis = []
        product_inchis = []
        num_reactants = []
        num_products = []
        reactant_nicknames = []
        product_nicknames = []
        reactant_charges = []
        product_charges = []
        reactant_spin_multiplicities = []
        product_spin_multiplicities = []
        for snl, n, nick_name in reactant_snls:
            mol = snl.structure
            bb = BabelMolAdaptor(mol)
            pbmol = bb.pybel_mol
            inchi = pbmol.write("inchi").strip()
            reactant_inchis.append(inchi)
            reactant_nicknames.append(nick_name)
            num_reactants.append(n)
            reactant_charges.append(mol.charge)
            reactant_spin_multiplicities.append(mol.spin_multiplicity)
        for snl, n, nick_name in product_snls:
            mol = snl.structure
            bb = BabelMolAdaptor(mol)
            pbmol = bb.pybel_mol
            inchi = pbmol.write("inchi").strip()
            product_inchis.append(inchi)
            product_nicknames.append(nick_name)
            num_products.append(n)
            product_charges.append(mol.charge)
            product_spin_multiplicities.append(mol.spin_multiplicity)
        all_inchis = reactant_inchis + product_inchis
        d = dict()
        d['submitter_email'] = submitter_email
        d['parameters'] = parameters
        d['state'] = 'SUBMITTED'
        d['reaction_id'] = self._get_next_reaction_id()
        d['submitted_at'] = datetime.datetime.utcnow().isoformat()
        d["reactant_snls"] = [s[0].as_dict() for s in reactant_snls]
        d["product_snls"] = [s[0].as_dict() for s in product_snls]
        d['all_inchis'] = all_inchis
        d['reactant_inchis'] = reactant_inchis
        d['product_inchis'] = product_inchis
        d['num_reactants'] = num_reactants
        d['num_products'] = num_products
        d['reactant_nicknames'] = reactant_nicknames
        d['product_nicknames'] = product_nicknames
        d['reactant_charges'] = reactant_charges
        d['product_charges'] = product_charges
        d['reactant_spin_multiplicities'] = reactant_spin_multiplicities
        d['product_spin_multiplicities'] = product_spin_multiplicities
        d['reactant_fragments'] = [[frag.to_dict() for frag in specie] for specie in reactant_fragments]
        d['product_fragments'] = [[frag.to_dict() for frag in specie] for specie in product_fragments]
        self.reactions.insert(d)
        dummy_snl = StructureNL.from_dict(d["reactant_snls"][0])
        parameters['reaction_id'] = d['reaction_id']
        self.submit_snl(dummy_snl, submitter_email, parameters)
        return d['reaction_id']

    def resubmit(self, submission_id):
        self.jobs.update(
            {'submission_id': submission_id},
            {'$set': {'state': 'SUBMITTED',
                      'state_details': {},
                      'task_dict': {}}})

    def cancel_submission(self, submission_id):
        # TODO: implement me
        # set state to 'cancelled'
        # in the SubmissionProcessor, detect this state and defuse the FW
        raise NotImplementedError()

    def set_job_state_to_cancel(self, submission_id):
        self.jobs.update({'submission_id': submission_id},
                         {'$set': {'state': 'CANCELLED'}})

    def get_submission_ids_after(self, dt):
        '''
        Return the submission id after a time point
        :param dt: a datetime object after which point you want to get
        :return: a list of IDs
        '''
        job_docs = self.jobs.find({}, fields=["submission_id",
                                              "submitted_at",
                                              "parameters.nick_name"])
        ids = []
        for j in job_docs:
            submission_dt = dateutil.parser.parse(j['submitted_at'])
            if submission_dt > dt:
                ids.append(j['submission_id'])
        return ids

    def get_states(self, crit):
        props = ['state', 'state_details', 'task_dict', 'submission_id',
                 'formula']
        infos = []
        for j in self.jobs.find(crit, dict([(p, 1) for p in props])):
            infos.append(dict([(p, j[p]) for p in props]))
        return infos

    def as_dict(self):
        """
        Note: usernames/passwords are exported as unencrypted Strings!
        """
        d = {'host': self.host, 'port': self.port, 'db': self.db,
             'username': self.username,
             'password': self.password}
        return d

    def update_state(self, submission_id, state, state_details, task_dict):
        self.jobs.find_and_modify({'submission_id': submission_id},
                                  {'$set': {'state': state,
                                            'state_details': state_details,
                                            'task_dict': task_dict}})

    @classmethod
    def from_dict(cls, d):
        return cls(d['host'], d['port'], d['db'], d['username'], d['password'])

    @classmethod
    def auto_load(cls):
        s_dir = os.environ['DB_LOC']
        s_file = os.path.join(s_dir, 'submission_db.yaml')
        return SubmissionMongoAdapterEG.from_file(s_file)

    def to_format(self, f_format='json', **kwargs):
        """
        returns a String representation in the given format
        :param f_format: the format to output to (default json)
        """
        if f_format == 'json':
            return json.dumps(self.as_dict(),
                              default=DATETIME_HANDLER,
                              **kwargs)
        elif f_format == 'yaml':
            # start with the JSON format, and convert to YAML
            return yaml.dump(self.as_dict(), default_flow_style=YAML_STYLE,
                             allow_unicode=True)
        else:
            raise ValueError('Unsupported format {}'.format(f_format))

    @classmethod
    def from_format(cls, f_str, f_format='json'):
        """
        convert from a String representation to its Object
        :param f_str: the String representation
        :param f_format: serialization format of the String (default json)
        """
        if f_format == 'json':
            return cls.from_dict(_reconstitute_dates(json.loads(f_str)))
        elif f_format == 'yaml':
            return cls.from_dict(_reconstitute_dates(yaml.load(f_str)))
        else:
            raise ValueError('Unsupported format {}'.format(f_format))

    def to_file(self, filename, f_format=None, **kwargs):
        """
        Write a serialization of this object to a file
        :param filename: filename to write to
        :param f_format: serialization format, default checks the filename
                         extension
        """
        if f_format is None:
            f_format = filename.split('.')[-1]
        with open(filename, 'w') as f:
            f.write(self.to_format(f_format=f_format, **kwargs))

    @classmethod
    def from_file(cls, filename, f_format=None):
        """
        Load a serialization of this object from a file
        :param filename: filename to read
        :param f_format: serialization format, default (None) checks the
                         filename extension
        """
        if f_format is None:
            f_format = filename.split('.')[-1]
        with open(filename, 'r') as f:
            return cls.from_format(f.read(), f_format=f_format)


def _reconstitute_dates(obj_dict):
    if obj_dict is None:
        return None

    if isinstance(obj_dict, dict):
        return {k: _reconstitute_dates(v) for k, v in obj_dict.items()}

    if isinstance(obj_dict, list):
        return [_reconstitute_dates(v) for v in obj_dict]

    if isinstance(obj_dict, basestring):
        try:
            return datetime.datetime.strptime(obj_dict, "%Y-%m-%dT%H:%M:%S.%f")
        except ValueError:
            pass

    return obj_dict

