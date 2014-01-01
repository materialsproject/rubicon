"""
Build molecules collection
Adapted from Dan Gunter and Wei Chen's vasp materials builder
"""
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
        'user_tags', 'run_tags', 'snl_final', 'task_id',
        'nelements', 'reduced_cell_formula', 'pretty_formula', 'pointgroup',
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
        docs = self._c.tasks.find(query, fields=TaskKeys.fields)
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
        molecule = None
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
        _log.info("Inserting Material from task_id {i}".format(i=doc['task_id']))
        self._c.materials.insert(doc)





