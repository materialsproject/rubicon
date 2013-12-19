import datetime
import json
import os
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
        self.id_assigner = self.database.id_assigner

        self._update_indices()

    def _reset(self):
        self._restart_id_assigner_at(1)
        self.jobs.remove()

    def _update_indices(self):
        self.jobs.ensure_index('submission_id', unique=True)
        self.jobs.ensure_index('state')
        self.jobs.ensure_index('submitter_email')

    def _get_next_submission_id(self):
        return self.id_assigner.find_and_modify(
            query={}, update={'$inc': {'next_submission_id': 1}})[
                'next_submission_id']

    def _restart_id_assigner_at(self, next_submission_id):
        self.id_assigner.remove()
        self.id_assigner.insert({"next_submission_id": next_submission_id})

    def submit_snl(self, snl, submitter_email, parameters=None):
        parameters = parameters if parameters else {}

        d = snl.to_dict
        d['submitter_email'] = submitter_email
        d['parameters'] = parameters
        d['state'] = 'SUBMITTED'
        d['state_details'] = {}
        d['task_dict'] = {}
        d['submission_id'] = self._get_next_submission_id()
        d['submitted_at'] = datetime.datetime.utcnow().isoformat()
        if 'is_valid' not in d:
            d.update(get_meta_from_structure(snl.structure))

        d.update(snl.structure.to_dict)

        self.jobs.insert(d)
        return d['submission_id']

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

    def get_states(self, crit):
        props = ['state', 'state_details', 'task_dict', 'submission_id',
                 'formula']
        infos = []
        for j in self.jobs.find(crit, dict([(p, 1) for p in props])):
            infos.append(dict([(p, j[p]) for p in props]))
        return infos

    def to_dict(self):
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
            return json.dumps(self.to_dict(),
                              default=DATETIME_HANDLER,
                              **kwargs)
        elif f_format == 'yaml':
            # start with the JSON format, and convert to YAML
            return yaml.dump(self.to_dict(), default_flow_style=YAML_STYLE,
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

