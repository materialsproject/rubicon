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
        'pointgroup',
        'calculations.scf.energies', 'calculations.scf_pcm.energies')


class MoleculesBuilder(eg_shared.ParallelBuilder):
    """Build derived 'molecules' collection.
    """

    def __init__(self, collections, **kwargs):
        """Create new molecules builder.

        Args:
            collections: Set of connected DB collections
                Type: eg_shared.Collections
        """
        eg_shared.ParallelBuilder.__init__(self, **kwargs)
        self._c = collections


    def run(self):
        """Run the builder.
        """
        _log.info("Getting distinct INCHI")
        inchi_finals = self._c.tasks.distinct('inchi_final')
        map(self.add_item, inchi_finals)
        _log.info("Beginning analysis")
        states = self.run_parallel()
        self._build_indexes()
        return self.combine_status(states)


    def process_item(self, inchi_finals):
        """Create and add material for a given grouping identifer.
        """
        query = {'state': 'successful', 'inchi_final': inchi_finals,
                 'task_type': "single point energy"}
        docs = list(self._c.tasks.find(query, fields=TaskKeys.fields))
        if not docs:
            return 1
        doc = self.build_molecule(docs)
        if not doc:
            return 2
        doc['created_at'] = datetime.datetime.now()
        doc['updated_at'] = datetime.datetime.now()
        self._insert_molecule(doc)
        return 0

    def build_molecule(self, taskdocs):
        """Transforms task document to molecules document.
        """
        docs = dict()
        for td in taskdocs:
            if td["charge"] == 0:
                docs["neutral"] = td
            elif td["charge"] == 1:
                docs["cation"] = td
            elif td["charge"] == -1:
                docs["anion"] = td
        if len(docs) < 2 or ("neutral" not in docs):
            return None
        molecule = dict()
        molecule["inchi"] = docs["neutral"]["inchi_final"]
        molecule["can"] = docs["neutral"]["can"]
        molecule["smiles"] = docs["neutral"]["smiles"]
        molecule["elements"] = copy.deepcopy(docs["neutral"]["elements"])
        molecule["nelements"] = docs["neutral"]["nelements"]
        molecule["user_tags"] = copy.deepcopy(docs["neutral"]["user_tags"])
        molecule["run_tags"] = copy.deepcopy(docs["neutral"]["run_tags"])
        molecule["reduced_cell_formula_abc"] = docs[
            "neutral"]["reduced_cell_formula_abc"]
        molecule["implicit_solvent"] = copy.deepcopy(docs["neutral"][
            "implicit_solvent"])
        molecule["pretty_formula"] = docs["neutral"]["pretty_formula"]
        molecule["pointgroup"] = docs["neutral"]["pointgroup"]

        molecule["task_id"] = dict()
        molecule["snlgroup_id_final"] = dict()
        molecule["charge"] = dict()
        molecule["spin_multiplicity"] = dict()
        molecule["snl_final"] = dict()
        molecule["molecule"] = dict()
        for k in docs.keys():
            molecule["task_id"][k] = docs[k]["task_id"]
            molecule["snlgroup_id_final"][k] = docs[k]["snlgroup_id_final"]
            molecule["charge"][k] = docs[k]["charge"]
            molecule["spin_multiplicity"][k] = docs[k]["spin_multiplicity"]
            molecule["snl_final"][k] = docs[k]["snl_final"]
            molecule["molecule"][k] = docs[k]["molecule_final"]
        if "cation" in docs:
            molecule["IP"] = dict()
            molecule["IP"]["vacuum"] = \
                docs["cation"]["calculations"]["scf"]["energies"][-1][-1] \
                - \
                docs["neutral"]["calculations"]["scf"]["energies"][-1][-1]
            molecule["IP"]["sol"] = \
                docs["cation"]["calculations"]["scf_pcm"]["energies"][-1][-1] \
                - \
                docs["neutral"]["calculations"]["scf_pcm"]["energies"][-1][-1]
        if "anion" in docs:
            molecule["EA"] = dict()
            molecule["EA"]["vacuum"] = \
                docs["neutral"]["calculations"]["scf"]["energies"][-1][-1] \
                - \
                docs["anion"]["calculations"]["scf"]["energies"][-1][-1]
            molecule["EA"]["sol"] = \
                docs["neutral"]["calculations"]["scf_pcm"]["energies"][-1][-1] \
                - \
                docs["anion"]["calculations"]["scf_pcm"]["energies"][-1][-1]
        return molecule

    def _build_indexes(self):
        self._c.molecules.ensure_index('task_id', unique=True)
        for key in ['snlgroup_id_final', 'inchi_final', 'task_type', 'can',
                    'smiles', 'charge', 'spin_multiplicity', 'nelements',
                    'reduced_cell_formula', 'pretty_formula']:
            _log.info("Building {} index".format(key))
            self._c.molecules.ensure_index(key)
        _log.info("Building nelements and elements compound index")
        compound_index = [('nelements', ASCENDING), ('elements', ASCENDING)]
        self._c.molecules.ensure_index(compound_index)

    def _insert_molecule(self, doc):
        """All database insertion should be done from this method
        """
        _log.info("Inserting Material from task_id i, ".
            format(i=str(doc['task_id'])))
        self._c.molecules.insert(doc)





