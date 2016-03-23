# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import copy
import datetime
import json
import logging
import math
import os
import sys

from fireworks.core.firework import FireTaskBase, FWAction, Workflow
from fireworks.utilities.fw_serializers import FWSerializable
from pymongo import MongoClient
from six.moves import zip

from pymatgen import Molecule
from pymatgen.analysis.molecule_structure_comparator import \
    MoleculeStructureComparator
from pymatgen.io.qchem import QcInput
from rubicon.borg.hive import DeltaSCFQChemToDbTaskDrone
from rubicon.dupefinders.dupefinder_eg import DupeFinderEG
from rubicon.utils.atomic_charge_mixed_basis_set_generator import \
    AtomicChargeMixedBasisSetGenerator
from rubicon.utils.eg_wf_utils import get_eg_file_loc, \
    get_defuse_causing_qchem_fwid
from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator
from rubicon.utils.snl.egsnl import EGStructureNL
from rubicon.utils.snl.egsnl_mongo import EGSNLMongoAdapter

__author__ = 'xiaohuiqu'


def get_basic_update_specs(fw_spec, d):
    update_specs = {'mol': d["molecule_final"],
                    'egsnl': d["snl_final"],
                    'snlgroup_id': d["snlgroup_id_final"],
                    'inchi_root': fw_spec["inchi_root"]}
    mixed_basis = None
    mixed_aux_basis = None
    if "mixed_basis" in fw_spec:
        mixed_basis = fw_spec["mixed_basis"]
    if "mixed_aux_basis" in fw_spec:
        mixed_aux_basis = fw_spec["mixed_aux_basis"]
    if "_mixed_basis_set_generator" in fw_spec:
        bs_generator_dict = fw_spec["_mixed_basis_set_generator"]
        if isinstance(d["molecule_final"], dict):
            mol = Molecule.from_dict(d["molecule_final"])
        else:
            mol = d["molecule_final"]
        pop_method = None
        if "scf" in d["calculations"]:
            if "nbo" in d["calculations"]["scf"]["charges"]:
                pop_method = "nbo"
            elif "hirshfeld" in d["calculations"]["scf"]["charges"]:
                pop_method = "hirshfeld"
        if pop_method is None:
            raise ValueError(
                "An vacuum single point caculation is require to use mixed basis set generator")
        charges = d["calculations"]["scf"]["charges"][pop_method]
        if isinstance(bs_generator_dict, dict):
            bs_generator = AtomicChargeMixedBasisSetGenerator.from_dict(
                bs_generator_dict)
        else:
            bs_generator = bs_generator_dict
        mixed_basis = bs_generator.get_basis(mol, charges)
    if "_mixed_aux_basis_set_generator" in fw_spec:
        aux_bs_generator_dict = fw_spec["_mixed_aux_basis_set_generator"]
        pop_method = None
        if "scf" in d["calculations"]:
            if "nbo" in d["calculations"]["scf"]["charges"]:
                pop_method = "nbo"
            elif "hirshfeld" in d["calculations"]["scf"]["charges"]:
                pop_method = "hirshfeld"
        if pop_method is None:
            raise ValueError(
                "An vacuum single point caculation is require to use mixed auxiliary basis set generator")
        charges = d["calculations"]["scf"]["charges"][pop_method]
        aux_bs_generator = AtomicChargeMixedBasisSetGenerator(
            aux_bs_generator_dict)
        mixed_aux_basis = aux_bs_generator.get_basis(mol, charges)
    if mixed_basis or mixed_aux_basis:
        update_specs["mixed_basis"] = mixed_basis
    if mixed_aux_basis:
        update_specs["mixed_aux_basis"] = mixed_basis
    return update_specs


def standard_parsing_db_insertion(fw_spec):
    if '_fizzled_parents' in fw_spec and not 'prev_qchem_dir' in fw_spec:
        prev_dir = fw_spec['_fizzled_parents'][0]['launches'][0][
            'launch_dir']
    else:
        prev_dir = fw_spec['prev_qchem_dir']
    db_dir = os.environ['DB_LOC']
    db_path = os.path.join(db_dir, 'tasks_db.json')
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('QChemDrone')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setLevel(getattr(logging, 'INFO'))
    logger.addHandler(sh)
    with open(db_path) as f:
        db_creds = json.load(f)
    drone = DeltaSCFQChemToDbTaskDrone(
        host=db_creds['host'], port=db_creds['port'],
        database=db_creds['database'], user=db_creds['admin_user'],
        password=db_creds['admin_password'],
        collection=db_creds['collection'])
    qcout_path = get_eg_file_loc(os.path.abspath(os.path.join(
        prev_dir, "mol.qcout")))
    t_id, d = drone.assimilate(qcout_path, fw_spec=fw_spec)
    return d, qcout_path, t_id


class QChemGeomOptDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem Geometry Optimization DB Insertion Task"

    def run_task(self, fw_spec):
        d, qcout_path, t_id = standard_parsing_db_insertion(fw_spec)
        update_specs = get_basic_update_specs(fw_spec, d)

        if d["state"] == "successful":
            return FWAction(
                stored_data={'task_id': t_id},
                update_spec=update_specs)
        else:
            if d['state'] == 'rejected' and \
                            d['reject_reason'] == 'structural change':
                defuse_reason = 'structural change'
            else:
                defuse_reason = d.get("errors", "unknown")
            offending_fwid = get_defuse_causing_qchem_fwid(qcout_path)
            return FWAction(
                stored_data={'task_id': t_id},
                update_spec=dict({'defuse_reason': defuse_reason,
                                  'offending_fwid': offending_fwid},
                                 **update_specs),
                defuse_children=True)


class QChemFrequencyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem Frequency DB Insertion Task"
    molecule_perturb_scale = 0.3

    def run_task(self, fw_spec):
        d, qcout_path, t_id = standard_parsing_db_insertion(fw_spec)
        update_specs = get_basic_update_specs(fw_spec, d)

        if d["state"] == "successful":
            if d['stationary_type'] == 'minimum':
                return FWAction(
                    stored_data={'task_id': t_id},
                    update_spec=update_specs)
            else:
                return self.img_freq_action(fw_spec, d, t_id, qcout_path)
        else:
            defuse_reason = d.get("errors", "unknown")
            offending_fwid = get_defuse_causing_qchem_fwid(qcout_path)
            return FWAction(
                stored_data={'task_id': t_id},
                defuse_children=True,
                update_spec=dict({'defuse_reason': defuse_reason,
                                  'offending_fwid': offending_fwid},
                                 **update_specs))

    @staticmethod
    def spawn_opt_freq_wf(mol, molname, mission, additional_user_tags,
                          priority, update_spec, charge,
                          spin_multiplicity, grid, qm_method):
        fw_creator = QChemFireWorkCreator(
            mol=mol, molname=molname, mission=mission,
            additional_user_tags=additional_user_tags, priority=priority,
            update_spec=update_spec)
        geom_fwid_cal, geom_fwid_db = -1, -2
        freq_fwid_cal, freq_fwid_db = -3, -4
        geom_fw_cal, geom_fw_db = fw_creator.geom_fw(
            charge, spin_multiplicity, geom_fwid_cal, geom_fwid_db,
            method=qm_method)
        geom_fw_cal.spec["run_tags"]["task_type_amend"] = "imaginary " \
                                                          "frequency elimination"
        freq_fw_cal, freq_fw_db = fw_creator.freq_fw(
            charge, spin_multiplicity, freq_fwid_cal, freq_fwid_db,
            method=qm_method)
        freq_fw_cal.spec["run_tags"]["task_type_amend"] = "imaginary " \
                                                          "frequency elimination"
        if grid:
            for fw in [geom_fw_cal, geom_fw_db, freq_fw_cal, freq_fw_db]:
                if isinstance(fw.spec["qcinp"], dict):
                    qcinp = QcInput.from_dict(fw.spec["qcinp"])
                else:
                    qcinp = fw.spec["qcinp"]
                for j in qcinp.jobs:
                    j.set_dft_grid(*grid)
                    j.set_integral_threshold(12)
                    if j.params["rem"]["jobtype"] == "opt":
                        j.scale_geom_opt_threshold(0.1, 0.1, 0.1)
                        j.set_geom_max_iterations(100)
                fw.spec["qcinp"] = qcinp.as_dict()
                fw.spec["run_tags"]["grid"] = grid
        wf = Workflow([geom_fw_cal, geom_fw_db, freq_fw_cal, freq_fw_db],
                      links_dict={geom_fwid_db: freq_fwid_cal,
                                  geom_fwid_cal: geom_fwid_db,
                                  freq_fwid_cal: freq_fwid_db})
        return wf

    @staticmethod
    def _check_structure_change(mol1, mol2, fw_spec):
        """
        Check whether structure is changed:

        Return:
            True: structure changed, False: unchanged
        """
        if 'egsnl' not in fw_spec:
            raise ValueError("Can't find initial SNL")
        if isinstance(fw_spec['egsnl'], EGStructureNL):
            egsnl_dict = fw_spec['egsnl'].as_dict()
        else:
            egsnl_dict = fw_spec['egsnl']
        if 'known_bonds' not in egsnl_dict:
            raise ValueError("Can't find known bonds information")
        bonds = egsnl_dict['known_bonds']
        msc = MoleculeStructureComparator(priority_bonds=bonds)
        return not msc.are_equal(mol1, mol2)

    @classmethod
    def perturb_molecule(cls, d, reversed_direction=False):
        old_mol = Molecule.from_dict(d['molecule_final'])
        vib_mode = d['calculations']['freq']['frequencies'][0]["vib_mode"]
        max_dis = max([math.sqrt(sum([x ** 2 for x in mode]))
                       for mode in vib_mode])
        scale = cls.molecule_perturb_scale / max_dis
        normalized_mode = [[x * scale for x in mode]
                           for mode in vib_mode]
        direction = 1.0
        if reversed_direction:
            direction = -1.0
        new_coords = [[c + v * direction for c, v in zip(site.coords, mode)]
                      for site, mode in zip(old_mol.sites, normalized_mode)]
        species = [site.specie.symbol
                   for site in old_mol.sites]
        charge = old_mol.charge
        spin_multiplicity = old_mol.spin_multiplicity
        new_mol = Molecule(species, new_coords, charge=charge,
                           spin_multiplicity=spin_multiplicity)
        return new_mol

    def img_freq_action(self, fw_spec, d, t_id, qcout_path):
        if "img_freq_eli" in d['user_tags']:
            img_freq_eli = copy.deepcopy(d['user_tags']["img_freq_eli"])
            img_freq_eli['current_method_id'] += 1
        else:
            img_freq_eli = {"methods": ["dir_dis_opt", "den_dis_opt",
                                        "alt_den_dis_opt"],
                            "current_method_id": 0}
        update_specs = get_basic_update_specs(fw_spec, d)
        if img_freq_eli['current_method_id'] >= len(img_freq_eli['methods']):
            logging.error("Failed to eliminate imaginary frequency")
            offending_fwid = get_defuse_causing_qchem_fwid(qcout_path)
            return FWAction(
                stored_data={'task_id': t_id},
                defuse_children=True,
                update_spec=dict({'defuse_reason': "imaginary frequency "
                                                   "elimination failed",
                                  'offending_fwid': offending_fwid},
                                 **update_specs))

        new_mol = self.perturb_molecule(d)
        old_mol = Molecule.from_dict(d['molecule_final'])
        structure_changed = self._check_structure_change(
            old_mol, new_mol, fw_spec)
        if structure_changed:
            self.perturb_molecule(d, reversed_direction=True)
            structure_changed = self._check_structure_change(
                old_mol, new_mol, fw_spec)
        if structure_changed:
            offending_fwid = get_defuse_causing_qchem_fwid(qcout_path)
            return FWAction(
                stored_data={'task_id': t_id},
                defuse_children=True,
                update_spec=dict({'perturbed_mol': new_mol.as_dict(),
                                  'defuse_reason': "structural change in imaginary "
                                                   "frequency elimination",
                                  'offending_fwid': offending_fwid},
                                 **update_specs))
        molname = d['user_tags']['molname']
        mission = d['user_tags']['mission']
        additional_user_tags = {"img_freq_eli": img_freq_eli}
        if "initial_charge" in fw_spec["user_tags"]:
            additional_user_tags["initial_charge"] = fw_spec["user_tags"][
                "initial_charge"]
        priority = fw_spec['_priority']
        old_snl = EGStructureNL.from_dict(d['snl_initial'])
        history = old_snl.history
        history.append(
            {'name': 'Electrolyte Genome Project eliminate imaginary '
                     'frequency by perturb molecular geometry',
             'url': 'http://www.materialsproject.org',
             'description': {'task_type': d['task_type'],
                             'task_id': d['task_id'],
                             'max_displacement': self.molecule_perturb_scale},
             'when': datetime.datetime.utcnow()})
        new_snl = EGStructureNL(new_mol, old_snl.authors, old_snl.projects,
                                old_snl.references, old_snl.remarks,
                                old_snl.data, history)

        # enter new SNL into SNL db
        # get the SNL mongo adapter
        sma = EGSNLMongoAdapter.auto_load()
        egsnl, snlgroup_id = sma.add_snl(
            new_snl, snlgroup_guess=d['snlgroup_id_initial'])
        update_specs = {'mol': new_mol.as_dict(),
                        'egsnl': egsnl.as_dict(),
                        'snlgroup_id': fw_spec['snlgroup_id'],
                        'inchi_root': fw_spec['inchi_root']}

        eli_strategy = img_freq_eli["methods"][
            img_freq_eli["current_method_id"]]
        charge = new_mol.charge
        spin_multiplicity = new_mol.spin_multiplicity
        qm_method = fw_spec["qm_method"]
        if eli_strategy == "dir_dis_opt":
            logging.info("Eliminate Imaginary Frequency By Perturbing the "
                         "Structure of Molecule")
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge,
                                        spin_multiplicity,
                                        grid=None, qm_method=qm_method)
        elif eli_strategy == "den_dis_opt":
            logging.info("Eliminate Imaginary Frequency By Perturbing the "
                         "Structure of Molecule, and increase the grid "
                         "density")
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge,
                                        spin_multiplicity,
                                        grid=(128, 302), qm_method=qm_method)
        elif eli_strategy == "alt_den_dis_opt":
            logging.info("Eliminate Imaginary Frequency By Perturbing the "
                         "Structure of Molecule, and increase the grid "
                         "density")
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge,
                                        spin_multiplicity,
                                        grid=(90, 590), qm_method=qm_method)
        else:
            raise Exception("Unknown imaginary frequency fixing method")

        return FWAction(stored_data={'task_id': t_id}, detours=wf,
                        update_spec=update_specs)


def get_bsse_update_specs(fw_spec, d):
    if "super_mol_snlgroup_id" in fw_spec["run_tags"]:
        ghost_atoms = fw_spec["run_tags"].get("ghost_atoms", list())
        from rubicon.workflows.bsse_wf import BSSEFragment
        fragment_type = fw_spec["run_tags"].get("bsse_fragment_type",
                                                BSSEFragment.ISOLATED)
        fragment_key = BasisSetSuperpositionErrorCalculationTask.get_fragment_key(
            ghost_atoms, fragment_type)
        fragment_dict = dict()
        fragment_dict["task_id"] = [d["task_id"]]
        fragment_dict["energy"] = d["calculations"]["scf"]["energies"][-1][-1]
        fragment_dict["ghost_atoms"] = ghost_atoms
        fragment_dict["fragment_type"] = fragment_type
        fragment_dict["super_mol_snlgroup_id"] = fw_spec["run_tags"][
            "super_mol_snlgroup_id"]
        return {fragment_key: fragment_dict}
    else:
        return {}


class QChemSinglePointEnergyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem Single Point Energy DB Insertion Task"

    def run_task(self, fw_spec):
        d, qcout_path, t_id = standard_parsing_db_insertion(fw_spec)
        update_specs = get_basic_update_specs(fw_spec, d)

        if d["state"] == "successful":
            bsse_specs = get_bsse_update_specs(fw_spec, d)
            update_specs.update(bsse_specs)
            return FWAction(
                stored_data={'task_id': t_id},
                update_spec=update_specs)
        else:
            offending_fwid = get_defuse_causing_qchem_fwid(qcout_path)
            return FWAction(
                stored_data={'task_id': t_id},
                defuse_children=True,
                update_spec=dict({'defuse_reason': 'SCF failed',
                                  'offending_fwid': offending_fwid},
                                 **update_specs))


class QChemAIMDDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem Ab Initio Molecule Dynamics DB Insertion Task"

    def run_task(self, fw_spec):
        d, qcout_path, t_id = standard_parsing_db_insertion(fw_spec)
        update_specs = get_basic_update_specs(fw_spec, d)

        if d["state"] == "successful":
            return FWAction(
                stored_data={'task_id': t_id},
                update_spec=update_specs)
        else:
            if d['state'] == 'rejected' and \
                            d['reject_reason'] == 'structural change':
                defuse_reason = 'structural change'
            else:
                defuse_reason = d.get("errors", "unknown")
            offending_fwid = get_defuse_causing_qchem_fwid(qcout_path)
            return FWAction(
                stored_data={'task_id': t_id},
                update_spec=dict({'defuse_reason': defuse_reason,
                                  'offending_fwid': offending_fwid},
                                 **update_specs),
                defuse_children=True)


class BasisSetSuperpositionErrorCalculationTask(FireTaskBase, FWSerializable):
    _fw_name = "BSSE Calculation Task"

    def run_task(self, fw_spec):
        fragment_dicts = fw_spec["fragments"]
        from rubicon.workflows.bsse_wf import BSSEFragment
        if len(fragment_dicts) > 0 and isinstance(fragment_dicts[0], dict):
            fragments = [BSSEFragment.from_dict(d) for d in fragment_dicts]
        else:
            fragments = fragment_dicts
        fragments_def = [d.to_dict() for d in fragments]
        fragments_dict = dict()
        bsse = 0.0
        for frag in fragments:
            if len(frag.ghost_atoms) == 0:
                continue
            fragment_name = self.get_fragment_name(frag.ghost_atoms)
            fragments_dict[fragment_name] = dict()
            ov_fragment_key = self.get_fragment_key(frag.ghost_atoms,
                                                    BSSEFragment.OVERLAPPED)
            fragments_dict[fragment_name][BSSEFragment.OVERLAPPED] = fw_spec[
                ov_fragment_key]
            ov_energy = fragments_dict[fragment_name][BSSEFragment.OVERLAPPED][
                "energy"]
            iso_fragment_key = self.get_fragment_key(frag.ghost_atoms,
                                                     BSSEFragment.ISOLATED)
            fragments_dict[fragment_name][BSSEFragment.ISOLATED] = fw_spec[
                iso_fragment_key]
            fragments_dict[fragment_name]["charge"] = frag.charge
            fragments_dict[fragment_name][
                "spin_multiplicity"] = frag.spin_multiplicity
            iso_energy = fragments_dict[fragment_name][BSSEFragment.ISOLATED][
                "energy"]
            fragments_dict[fragment_name][
                "fragment_bsse"] = ov_energy - iso_energy
            bsse += fragments_dict[fragment_name]["fragment_bsse"]
        result_dict = dict()
        result_dict["fragments_result"] = fragments_dict
        result_dict["fragments_def"] = fragments_def
        result_dict["bsse"] = bsse
        if isinstance(fw_spec["mol"], dict):
            result_dict["mol"] = fw_spec["mol"]
        else:
            result_dict["mol"] = fw_spec["mol"].as_dict()
        if isinstance(fw_spec['egsnl'], EGStructureNL):
            egsnl_dict = fw_spec['egsnl'].as_dict()
        else:
            egsnl_dict = fw_spec['egsnl']
        result_dict["egsnl"] = egsnl_dict
        result_dict["snlgroup_id"] = fw_spec["snlgroup_id"]
        result_dict["user_tags"] = fw_spec["user_tags"]
        result_dict["qm_method"] = fw_spec["qm_method"]
        result_dict["charge"] = fw_spec["charge"]
        result_dict["spin_multiplicity"] = fw_spec["spin_multiplicity"]
        result_dict["inchi_root"] = fw_spec["inchi_root"]
        result_dict['task_type'] = "BSSE Counterpoise Correction"
        result_dict["super_mol_snlgroup_id"] = fw_spec["snlgroup_id"]
        t_id = self._insert_doc(result_dict)
        return FWAction(
            stored_data={'task_id': t_id},
            update_spec={"mol": fw_spec["mol"],
                         "egsnl": fw_spec["egsnl"],
                         "snlgroup_id": fw_spec["snlgroup_id"],
                         "inchi_root": fw_spec["inchi_root"]})

    @classmethod
    def _insert_doc(cls, d, update_duplicates=True):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'tasks_db.json')
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('QChemDrone')
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        logger.addHandler(sh)
        with open(db_path) as f:
            db_creds = json.load(f)
        conn = MongoClient(db_creds['host'], db_creds['port'],
                           connect=False)
        db = conn[db_creds['database']]
        if db_creds['admin_user']:
            db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
        coll = db[db_creds['collection']]

        result = coll.find_one(
            filter={"super_mol_snlgroup_id": d["super_mol_snlgroup_id"],
                    "fragments_def": d["fragments_def"]},
            projection=["task_id", "super_mol_snlgroup_id", "fragments_def"])
        if result is None or update_duplicates:
            d["last_updated"] = datetime.datetime.today()
            if result is None:
                if ("task_id" not in d) or (not d["task_id"]):
                    id_num = db.counter.find_and_modify(
                        query={"_id": "mol_taskid"},
                        update={"$inc": {"c": 1}}
                    )["c"]
                    d["task_id"] = "mol-" + str(id_num)
                    d["task_id_deprecated"] = id_num
                logger.info("Inserting BSSE for snlgroup {} with taskid = {}"
                            .format(d["snlgroup_id"], d["task_id"]))
            elif update_duplicates:
                d["task_id"] = result["task_id"]
                logger.info("Updating BSSE for snlgroup {} with taskid = {}"
                            .format(d["super_mol_snlgroup_id"], d["task_id"]))
            coll.update({"super_mol_snlgroup_id": d["super_mol_snlgroup_id"],
                         "fragments_def": d["fragments_def"]},
                        {"$set": d},
                        upsert=True)
            return d["task_id"]
        else:
            logger.info("Skipping duplicate snlgrup {}".format(
                d["super_mol_snlgroup_id"]))
        conn.close()
        if "task_id" in d:
            return d["task_id"]
        elif "task_id" in result:
            return result["task_id"]
        else:
            return None

    @classmethod
    def get_fragment_name(cls, ghost_atoms):
        return "ga_" + "-".join(
            [str(i) for i in sorted(set(ghost_atoms))]) if len(
            ghost_atoms) > 0 else "ga_none"

    @classmethod
    def get_fragment_key(cls, ghost_atoms, fragment_type):
        return "{}_fragment_{}".format(fragment_type,
                                       cls.get_fragment_name(ghost_atoms))


class CounterpoiseCorrectionGenerationTask(FireTaskBase, FWSerializable):
    _fw_name = "Counterpoise Generation Task"

    def run_task(self, fw_spec):
        molname = fw_spec["user_tags"]["molname"]
        mission = fw_spec["user_tags"]["mission"]
        super_mol_snlgroup_id = fw_spec["snlgroup_id"]
        charge = fw_spec["charge"]
        spin_multiplicity = fw_spec["spin_multiplicity"]
        inchi_root = fw_spec["inchi_root"]
        if isinstance(fw_spec["egsnl"], dict):
            egsnl = fw_spec["egsnl"]
        else:
            egsnl = fw_spec["egsnl"].as_dict()
        qm_method = fw_spec["qm_method"]
        fragment_dicts = fw_spec["fragments"]
        large = fw_spec["large"]
        from rubicon.workflows.bsse_wf import bsse_wf, BSSEFragment
        if len(fragment_dicts) > 0 and isinstance(fragment_dicts[0], dict):
            fragments = [BSSEFragment.from_dict(d) for d in fragment_dicts]
        else:
            fragments = fragment_dicts
        priority = fw_spec.get('_priority', 1)
        dupefinder = fw_spec.get('_dupefinder', DupeFinderEG())
        cc_wf = bsse_wf(super_mol=egsnl, name=molname,
                        super_mol_snlgroup_id=super_mol_snlgroup_id,
                        super_mol_charge=charge,
                        super_mol_spin_multiplicity=spin_multiplicity,
                        super_mol_inchi_root=inchi_root, qm_method=qm_method,
                        fragments=fragments, mission=mission,
                        dupefinder=dupefinder, priority=priority,
                        is_spawnned=True, large=large)

        return FWAction(detours=cc_wf)
