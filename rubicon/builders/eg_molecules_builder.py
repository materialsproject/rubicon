"""
Build molecules collection
Adapted from Dan Gunter and Wei Chen's vasp materials builder
"""
import copy
import logging
import datetime
from pymongo import ASCENDING
from rubicon.builders import eg_shared

__author__ = "Xiaohui Qu"
__copyright__ = "Copyright 2012-2013, The Electrolyte Genome Project"
__version__ = "1.0"
__maintainer__ = "Xiaohui Qu"
__email__ = "xqu@lbl.gov"
__status__ = "Development"
__date__ = "1/1/14"



_log = logging.getLogger('eg.' + __name__)



class TaskKeys:
    """Keys we need to project from task collection to do
       the work of building the materials collection.
    """
    fields = (
        'task_id', 'snlgroup_id_final', 'inchi_final', 'task_type', 'elements',
        'can', 'smiles', 'charge', 'spin_multiplicity', 'implicit_solvent',
        'user_tags', 'run_tags', 'snl_final', 'task_id', "molecule_final",
        'nelements', 'reduced_cell_formula_abc', 'pretty_formula',
        'pointgroup', 'inchi_root',
        'calculations.scf.energies', 'calculations.scf_pcm.energies',
        'calculations.scf_sm12mk.energies',
        'formula', 'task_id_deprecated', 'svg', 'xyz')


class MoleculesBuilder(eg_shared.ParallelBuilder):
    """Build derived 'molecules' collection.
    """

    # absolute electrode potentials of some metals/molecules
    ref_potentials = {
        'hydrogen': 4.44,
        'magnesium': 2.07,
        'lithium': 1.40}

    def __init__(self, collections, **kwargs):
        """Create new molecules builder.

        Args:
            collections: Set of connected DB collections
                Type: eg_shared.Collections
        """
        eg_shared.ParallelBuilder.__init__(self, **kwargs)
        self._c = collections
        self._c.molecules.remove()
        self.ref_charge = 0
        self.ref_charge_range = (-1, 0, 1)


    def run(self):
        """Run the builder.
        """
        sss = []
        _log.info("Getting distinct root INCHIs")
        inchi_root = self._c.tasks.distinct('inchi_root')
        for ch in self.ref_charge_range:
            self.ref_charge = ch
            map(self.add_item, inchi_root)
            _log.info("Beginning analysis")
            states = self.run_parallel()
            sss.extend(states)
            self._build_indexes()
        return self.combine_status(sss)


    def process_item(self, inchi_root):
        """Create and add material for a given grouping identifer.
        """
        query = {'state': 'successful', 'inchi_root': inchi_root,
                 'task_type': "single point energy"}
        solvents = self._c.tasks.find(query, fields=TaskKeys.fields).distinct(
            "implicit_solvent.solvent_name"
        )
        molecule = dict()
        molecule['charge'] = self.ref_charge
        docs_available = False
        molecule['solvated_properties'] = dict()
        for solvent in solvents:
            query['implicit_solvent.solvent_name'] = solvent
            docs = list(self._c.tasks.find(query, fields=TaskKeys.fields))
            if docs:
                docs_available = True
            d = self.build_molecule_solvated_properties(docs)
            if d and len(d) > 0:
                molecule['solvated_properties'][solvent] = d
        if not docs_available:
            return 1
        if len(molecule['solvated_properties']) == 0:
            return 2
        del query['implicit_solvent.solvent_name']
        d = self.build_molecule_vacuum_properties(copy.deepcopy(query))
        if d and len(d) > 0:
            molecule['vacuum_properties'] = d
        else:
            return 2
        query['charge'] = self.ref_charge
        docs = self._c.tasks.find_one(query, fields=TaskKeys.fields)
        if not docs:
            return 1
        d = self.build_molecule_common_properties(docs)
        if d and len(d) > 0:
            molecule.update(d)
        else:
            return 2
        molecule['created_at'] = datetime.datetime.now()
        molecule['updated_at'] = datetime.datetime.now()
        self._insert_molecule(molecule)
        return 0

    def build_molecule_ipea(self, docs, molecule, solution_phase=True):
        molecule["task_id"] = dict()
        molecule["task_id_deprecated"] = dict()
        molecule["snlgroup_id_final"] = dict()
        molecule["charge"] = dict()
        molecule["spin_multiplicity"] = dict()
        molecule["snl_final"] = dict()
        molecule["molecule"] = dict()
        molecule["xyz"] = dict()
        molecule["inchi"] = dict()
        molecule["can"] = dict()
        molecule["smiles"] = dict()
        for k in docs.keys():
            molecule["task_id"][k] = docs[k]["task_id"]
            molecule["task_id_deprecated"][k] = docs[k]["task_id_deprecated"]
            molecule["snlgroup_id_final"][k] = docs[k]["snlgroup_id_final"]
            molecule["charge"][k] = docs[k]["charge"]
            molecule["spin_multiplicity"][k] = docs[k]["spin_multiplicity"]
            molecule["snl_final"][k] = docs[k]["snl_final"]
            molecule["molecule"][k] = docs[k]["molecule_final"]
            molecule["xyz"][k] = docs[k]["xyz"]
            molecule["inchi"][k] = docs[k]["inchi_final"]
            molecule["can"][k] = docs[k]["can"]
            molecule["smiles"][k] = docs[k]["smiles"]
        if solution_phase:
            # get the solution phase scf key name, scf_pcm, scf_sm12mk, etc.
            scf_all = set(docs["neutral"]["calculations"].keys())
            scf_all.remove('scf')
            scf_name = scf_all.pop()
        else:
            scf_name = 'scf'
        if "cation" in docs:
            molecule["IP"] = \
                docs["cation"]["calculations"][scf_name]["energies"][-1][-1] \
                - \
                docs["neutral"]["calculations"][scf_name]["energies"][-1][-1]
        if "anion" in docs:
            molecule["EA"] = \
                docs["neutral"]["calculations"][scf_name]["energies"][-1][-1] \
                - \
                docs["anion"]["calculations"][scf_name]["energies"][-1][-1]
        molecule['electrode_potentials'] = dict()
        if solution_phase:
            if 'IP' in molecule:
                molecule['electrode_potentials']['oxidation'] = dict()
                for electrode in self.ref_potentials.keys():
                    molecule['electrode_potentials']['oxidation'][electrode] \
                        = -(molecule['IP'] - self.ref_potentials[electrode])
            if 'EA' in molecule:
                molecule['electrode_potentials']['reduction'] = dict()
                for electrode in self.ref_potentials.keys():
                    molecule['electrode_potentials']['reduction'][electrode] \
                        = molecule['EA'] - self.ref_potentials[electrode]
        return scf_name

    def build_molecule_solvated_properties(self, taskdocs):
        docs = dict()
        for td in taskdocs:
            if td["charge"] == self.ref_charge:
                docs["neutral"] = td
            elif td["charge"] == self.ref_charge + 1:
                docs["cation"] = td
            elif td["charge"] == self.ref_charge - 1:
                docs["anion"] = td
        if len(docs) < 2 or ("neutral" not in docs):
            return None
        molecule = dict()
        scf = self.build_molecule_ipea(docs, molecule, solution_phase=True)
        molecule['solvation_energy'] = docs["neutral"]["calculations"]["scf"][
            "energies"][-1][-1] - \
            docs["neutral"]["calculations"][scf]["energies"][-1][-1]
        molecule["implicit_solvent"] = copy.deepcopy(docs['neutral'][
            "implicit_solvent"])
        return molecule

    def build_molecule_vacuum_properties(self, query):
        docs = dict()
        for c, i in zip(["anion", "neutral", "cation"], [-1, 0, 1]):
            query['charge'] = self.ref_charge + i
            taskdocs = self._c.tasks.find_one(query, fields=TaskKeys.fields)
            if not taskdocs:
                return None
            docs[c] = taskdocs
        if len(docs) < 2 or ("neutral" not in docs):
            return None
        molecule = dict()
        self.build_molecule_ipea(docs, molecule, solution_phase=False)
        return molecule

    def build_molecule_common_properties(self, docs):
        """Transforms task document to molecules document.
        """
        molecule = dict()
        molecule["inchi_root"] = docs["inchi_root"]
        molecule["elements"] = copy.deepcopy(docs["elements"])
        molecule["nelements"] = docs["nelements"]
        molecule["user_tags"] = copy.deepcopy(docs["user_tags"])
        molecule["run_tags"] = copy.deepcopy(docs["run_tags"])
        molecule["reduced_cell_formula_abc"] = docs["reduced_cell_formula_abc"]
        molecule["pretty_formula"] = docs["pretty_formula"]
        molecule["formula"] = docs["formula"]
        molecule["pointgroup"] = docs["pointgroup"]
        molecule["svg"] = docs["svg"]
        return molecule

    def _build_indexes(self):
        self._c.molecules.ensure_index(
            [('inchi_root', ASCENDING), ('charge', ASCENDING)], unique=True)
        for key in ['inchi_root', 'charge', 'nelements', 'elements',
                    'reduced_cell_formula', 'pretty_formula']:
            _log.info("Building {} index".format(key))
            self._c.molecules.ensure_index(key)
        _log.info("Building nelements and elements compound index")
        compound_index = [('nelements', ASCENDING), ('elements', ASCENDING)]
        self._c.molecules.ensure_index(compound_index)

    def _insert_molecule(self, doc):
        """All database insertion should be done from this method
        """
        _log.info("Inserting Material with InChI i, ".
            format(i=str(doc['inchi_root'])))
        self._c.molecules.insert(doc)





