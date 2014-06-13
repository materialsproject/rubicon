"""
Build molecules collection
Adapted from Dan Gunter and Wei Chen's vasp materials builder
"""
import copy
import json
import logging
import datetime
import math
from pymongo import ASCENDING
import sys
import re
from rubicon.builders import eg_shared
from rubicon.submission.submission_mongo_eg import SubmissionMongoAdapterEG
from rubicon.workflows.multistep_ipea_wf import QChemFireWorkCreator

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

    def __init__(self):
        pass

    tasks_fields = (
        'task_id', 'snlgroup_id_final', 'inchi_final', 'task_type', 'elements',
        'can', 'smiles', 'charge', 'spin_multiplicity', 'implicit_solvent',
        'user_tags', 'run_tags', 'snl_final', 'task_id', "molecule_final",
        'nelements', 'reduced_cell_formula_abc', 'pretty_formula',
        'pointgroup', 'inchi_root',
        'calculations',
        'formula', 'task_id_deprecated', 'svg', 'xyz')
    reactions_fields = ('reaction_id', 'num_reactants', 'num_products',  'reactant_nicknames',
                        'product_nicknames', 'reactant_inchis', 'product_inchis',
                        'reactant_submission_ids', 'product_submission_ids', 'all_inchis',
                        'reactant_spin_multiplicities', 'product_spin_multiplicities',
                        'reactant_charges', 'product_charges', "submitter_email")




class ReactionsBuilder(eg_shared.ParallelBuilder):
    """Build derived 'reactions' collection.
    """
    GAS_CONSTANT = 8.3144621 * (0.01036410 * 10E-3)  # eV K^-1 mol^-1
    TEMPERATURE = 298.15 # K

    def __init__(self, collections, **kwargs):
        """Create new molecules builder.

        Args:
            collections: Set of connected DB collections
                Type: eg_shared.Collections
        """
        eg_shared.ParallelBuilder.__init__(self, **kwargs)
        self._c = collections
        self._c.reactions.remove()
        sma = SubmissionMongoAdapterEG.auto_load()
        self.source_reactions = sma.reactions
        logging.basicConfig(level=logging.INFO)
        _log.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        _log.addHandler(sh)

    def run(self):
        """Run the builder.
        """
        _log.info("Getting Reaction Indices")
        reactions = list(self.source_reactions.find({}, fields=TaskKeys.reactions_fields))
        map(self.add_item, reactions)
        _log.info("Beginning analysis")
        states = self.run_parallel()
        return self.combine_status(states)

    def find_reaction_tasks_docs(self, solvent, reaction):
        reactant_freq_docs = []
        reactant_sp_docs = []
        product_freq_docs = []
        product_sp_docs = []
        freq_query_template = {"state": "successful",
                               "task_type": "vibrational frequency",
                               "stationary_type": "minimum"}
        sp_query_template = {"implicit_solvent.solvent_name": solvent,
                             "state": "successful",
                             "task_type": "single point energy"}
        for inchi, charge, spin in zip(reaction["reactant_inchis"],
                                       reaction["reactant_charges"],
                                       reaction["reactant_spin_multiplicities"]):
            freq_query = copy.deepcopy(freq_query_template)
            freq_query["inchi_root"] = inchi
            freq_query["charge"] = charge
            freq_query["spin_multiplicity"] = spin
            freq_doc = self._c.tasks.find_one(freq_query, fields=TaskKeys.tasks_fields)
            if not freq_doc:
                return None
            reactant_freq_docs.append(freq_doc)
            sp_query = copy.deepcopy(sp_query_template)
            sp_query["inchi_root"] = inchi
            sp_query["charge"] = charge
            sp_query["spin_multiplicity"] = spin
            sp_doc = self._c.tasks.find_one(sp_query, fields=TaskKeys.tasks_fields)
            if not sp_doc:
                return None
            reactant_sp_docs.append(sp_doc)

        for inchi, charge, spin in zip(reaction["product_inchis"],
                                       reaction["product_charges"],
                                       reaction["product_spin_multiplicities"]):
            freq_query = copy.deepcopy(freq_query_template)
            freq_query["inchi_root"] = inchi
            freq_query["charge"] = charge
            freq_query["spin_multiplicity"] = spin
            freq_doc = self._c.tasks.find_one(freq_query, fields=TaskKeys.tasks_fields)
            if not freq_doc:
                return None
            product_freq_docs.append(freq_doc)
            sp_query = copy.deepcopy(sp_query_template)
            sp_query["inchi_root"] = inchi
            sp_query["charge"] = charge
            sp_query["spin_multiplicity"] = spin
            sp_doc = self._c.tasks.find_one(sp_query, fields=TaskKeys.tasks_fields)
            if not sp_doc:
                return None
            product_sp_docs.append(sp_doc)

        return [zip(reactant_freq_docs, reactant_sp_docs),
                zip(product_freq_docs, product_sp_docs)]

    def build_reaction_data(self, docs, reaction, solution_phase=True):
        data = dict()
        for side, freq_and_sps, counts in zip(["reactant", "product"],
                                              docs,
                                              [reaction["num_reactants"], reaction["num_products"]]):
            data[side] = []
            for n, freq_and_sp in zip(counts, freq_and_sps):
                specie = dict()
                specie["number"] = n
                specie["task_id"] = dict()
                specie["task_id_deprecated"] = dict()
                specie["snlgroup_id_final"] = freq_and_sp[0]["snlgroup_id_final"]
                specie["charge"] = freq_and_sp[0]["charge"]
                specie["spin_multiplicity"] = freq_and_sp[0]["spin_multiplicity"]
                specie["snl_final"] = freq_and_sp[0]["snl_final"]
                specie["molecule"] = freq_and_sp[0]["molecule_final"]
                specie["xyz"] = freq_and_sp[0]["xyz"]
                specie["inchi"] = freq_and_sp[0]["inchi_final"]
                specie["can"] = freq_and_sp[0]["can"]
                specie["smiles"] = freq_and_sp[0]["smiles"]
                specie["inchi_root"] = freq_and_sp[0]["inchi_root"]
                specie["elements"] = freq_and_sp[0]["elements"]
                specie["nelements"] = freq_and_sp[0]["nelements"]
                specie["user_tags"] = freq_and_sp[0]["user_tags"]
                specie["run_tags"] = freq_and_sp[0]["run_tags"]
                specie["reduced_cell_formula_abc"] = freq_and_sp[0]["reduced_cell_formula_abc"]
                specie["pretty_formula"] = freq_and_sp[0]["pretty_formula"]
                specie["formula"] = freq_and_sp[0]["formula"]
                specie["pointgroup"] = freq_and_sp[0]["pointgroup"]
                specie["svg"] = freq_and_sp[0]["svg"]
                freq_cal_doc = freq_and_sp[0]["calculations"]
                sp_cal_doc = freq_and_sp[1]["calculations"]
                specie["thermo_corrections"] = freq_cal_doc["freq"]["corrections"]
                if solution_phase:
                    # get the solution phase scf key name, scf_pcm, scf_sm12mk, etc.
                    scf_all = set(sp_cal_doc.keys())
                    scf_all.remove('scf')
                    scf_name = scf_all.pop()
                else:
                    scf_name = 'scf'
                specie["scf_energy"] = sp_cal_doc[scf_name]["energies"][-1][-1]
                for task_type, d in zip(["freq", "sp"], freq_and_sp):
                    specie["task_id"][task_type] = d["task_id"]
                    specie["task_id_deprecated"][task_type] = d["task_id_deprecated"]
                data[side].append(specie)
        gibbs_energy = {"reactant": [], "product": []}
        for side in ["reactant", "product"]:
            for specie in data[side]:
                elec_energy = specie["scf_energy"]
                h = specie["thermo_corrections"]["Total Enthalpy"]
                s = specie["thermo_corrections"]["Total Entropy"]
                g = elec_energy + (h - self.TEMPERATURE * s)
                gibbs_energy[side].append(g)
        data["total_gibbs_free_energies"] = gibbs_energy
        reactant_energy = sum([energy*n for energy, n in zip(gibbs_energy["reactant"], reaction["num_reactants"])])
        product_energy = sum([energy*n for energy, n in zip(gibbs_energy["product"], reaction["num_products"])])
        data["delta_g"] = product_energy - reactant_energy
        data["equilibrium_constants"] = math.exp(-data["delta_g"]/(self.GAS_CONSTANT*self.TEMPERATURE))
        return data


    def process_item(self, reaction):
        """Create and add material for a given grouping identifer.
        """

        solvents = self._c.tasks.find({"inchi_root": {"$in": reaction["all_inchis"]},
                                       "state": "successful",
                                       "task_type": {"$in": ["vibrational frequency",
                                                             "single point energy"]}},
                                      fields=TaskKeys.tasks_fields)\
            .distinct("implicit_solvent.solvent_name")
        fe_docs = dict()
        fe_docs["reaction_id"] = reaction["reaction_id"]
        fe_docs["all_inchis"] = reaction["all_inchis"]
        fe_docs["reactant_inchis"] = reaction["reactant_inchis"]
        fe_docs["product_inchis"] = reaction["product_inchis"]
        fe_docs["reactant_nicknames"] = reaction["reactant_nicknames"]
        fe_docs["product_nicknames"] = reaction["product_nicknames"]
        fe_docs["submitter_email"] = reaction["submitter_email"]
        docs_available = False
        fe_docs['solvated_properties'] = dict()
        for solvent in solvents:
            docs = self.find_reaction_tasks_docs(solvent, reaction)
            if docs:
                docs_available = True
            d = self.build_reaction_data(docs, reaction, solution_phase=True)
            if d and len(d) > 0:
                fe_docs['solvated_properties'][solvent] = d
                if "vacuum_properties" not in fe_docs:
                    fe_docs["vacuum_properties"] = self.build_reaction_data(docs, reaction, solution_phase=False)
        if not docs_available or "vacuum_properties" not in fe_docs or fe_docs["vacuum_properties"] is None:
            return 1
        if len(fe_docs['solvated_properties']) == 0:
            return 2

        fe_docs['created_at'] = datetime.datetime.now()
        fe_docs['updated_at'] = datetime.datetime.now()

        self._insert_molecule(fe_docs)
        return 0

    def _build_indexes(self):
        _log.info("Building reaction index")
        self._c.reactions.ensure_index("reaction_id", unique=True)
        self._c.reactions.ensure_index("all_inchis")
        self._c.reactions.ensure_index("reactant_inchis")
        self._c.reactions.ensure_index("product_inchis")
        self._c.reactions.ensure_index("reactant_nicknames")
        self._c.reactions.ensure_index("product_nicknames")
        self._c.reactions.ensure_index("submitter_email")

    def _insert_molecule(self, doc):
        """All database insertion should be done from this method
        """
        _log.info('Inserting Reaction with ID "{i}", '.
                  format(i=str(doc['reaction_id'])))
        self._c.reactions.insert(doc)
