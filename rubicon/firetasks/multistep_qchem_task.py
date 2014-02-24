import copy
import json
import logging
import os
import sys
import math
import datetime

from fireworks.core.firework import FireTaskBase, FWAction, Workflow
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen import Molecule
from pymatgen.analysis.molecule_structure_comparator import MoleculeStructureComparator
from pymatgen.io.qchemio import QcInput

from rubicon.borg.hive import DeltaSCFQChemToDbTaskDrone
from rubicon.utils.snl.egsnl import EGStructureNL
from rubicon.utils.snl.egsnl_mongo import EGSNLMongoAdapter


__author__ = 'xiaohuiqu'


class QChemGeomOptDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem Geometry Optimization DB Insertion Task"

    def run_task(self, fw_spec):
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
        t_id, d = drone.assimilate(os.path.abspath(
            os.path.join(prev_dir, "mol.qcout")), fw_spec=fw_spec)

        if d["state"] == "successful":
            return FWAction(
                stored_data={'task_id': t_id},
                update_spec={
                    'mol': d["molecule_final"],
                    'egsnl': d["snl_final"],
                    'snlgroup_id': d["snlgroup_id_final"],
                    'inchi_root': fw_spec["inchi_root"]})
        else:
            if d['state'] == 'rejected' and \
                    d['reject_reason'] == 'structural change':
                inchi_root = d['snlgroup_id_final']
                defuse_reason = 'structural change'
            else:
                inchi_root = fw_spec['inchi_root']
                defuse_reason = d.get("errors", "unknown")
            return FWAction(
                stored_data={'task_id': t_id},
                update_spec={'mol': d["molecule_final"],
                             'egsnl': d["snl_final"],
                             'snlgroup_id': d["snlgroup_id_final"],
                             'inchi_root': inchi_root,
                             'defuse_reason': defuse_reason,
                             'offending_fw_id': fw_spec['fw_id']},
                defuse_children=True)


class QChemFrequencyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem Frequency DB Insertion Task"
    molecule_perturb_scale = 0.3

    def run_task(self, fw_spec):
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
        t_id, d = drone.assimilate(os.path.abspath(
            os.path.join(prev_dir, "mol.qcout")), fw_spec=fw_spec)

        if d["state"] == "successful":
            if d['stationary_type'] == 'minimum':
                return FWAction(
                    stored_data={'task_id': t_id},
                    update_spec={
                        'mol': d["molecule_final"],
                        'egsnl': d["snl_final"],
                        'snlgroup_id': d["snlgroup_id_final"],
                        'inchi_root': fw_spec["inchi_root"]})
            else:
                return self.img_freq_action(fw_spec, d, t_id)
        else:
            defuse_reason = d.get("errors", "unknown")
            return FWAction(
                stored_data={'task_id': t_id},
                defuse_children=True,
                update_spec={
                    'mol': d["molecule_final"],
                    'egsnl': d["snl_final"],
                    'snlgroup_id': d["snlgroup_id_final"],
                    'inchi_root': fw_spec["inchi_root"],
                    'defuse_reason': defuse_reason,
                    'offending_fw_id': fw_spec['fw_id']})

    @staticmethod
    def spawn_opt_freq_wf(mol, molname, mission, additional_user_tags,
                          priority, update_spec, charge,
                          spin_multiplicity, grid):
        from rubicon.workflows.multistep_ipea_wf \
            import QChemFireWorkCreator
        fw_creator = QChemFireWorkCreator(
            mol=mol, molname=molname, mission=mission,
            additional_user_tags=additional_user_tags, priority=priority,
            update_spec=update_spec)
        geom_fwid_cal, geom_fwid_db = -1, -2
        freq_fwid_cal, freq_fwid_db = -3, -4
        geom_fw_cal, geom_fw_db = fw_creator.geom_fw(charge, spin_multiplicity,
            geom_fwid_cal, geom_fwid_db)
        geom_fw_cal.spec["run_tags"]["task_type_amend"] = "imaginary " \
            "frequency elimination"
        freq_fw_cal, freq_fw_db = fw_creator.freq_fw(charge, spin_multiplicity,
            freq_fwid_cal, freq_fwid_db)
        freq_fw_cal.spec["run_tags"]["task_type_amend"] = "imaginary " \
            "frequency elimination"
        if grid:
            for fw in [geom_fw_cal, geom_fw_db, freq_fw_cal, freq_fw_db]:
                qcinp = QcInput.from_dict(fw.spec["qcinp"])
                for j in qcinp.jobs:
                    j.set_dft_grid(*grid)
                    if j.params["rem"]["jobtype"] == "opt":
                        j.scale_geom_opt_threshold(0.1, 0.1, 0.1)
                        j.set_geom_max_iterations(100)
                fw.spec["qcinp"] = qcinp.to_dict
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
        if 'known_bonds' not in fw_spec['egsnl']:
            raise ValueError("Can't find known bonds information")
        bonds = fw_spec['egsnl']['known_bonds']
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
        new_coords = [[c+v*direction for c, v in zip(site.coords, mode)]
                      for site, mode in zip(old_mol.sites, normalized_mode)]
        species = [site.specie.symbol
                   for site in old_mol.sites]
        charge = old_mol.charge
        spin_multiplicity = old_mol.spin_multiplicity
        new_mol = Molecule(species, new_coords, charge=charge,
                           spin_multiplicity=spin_multiplicity)
        return new_mol

    def img_freq_action(self, fw_spec, d, t_id):
        if "img_freq_eli" in d['user_tags']:
            img_freq_eli = copy.deepcopy(d['user_tags']["img_freq_eli"])
            img_freq_eli['current_method_id'] += 1
        else:
            img_freq_eli = {"methods": ["dir_dis_opt", "den_dis_opt",
                                        "alt_den_dis_opt"],
                            "current_method_id": 0}

        if img_freq_eli['current_method_id'] >= len(img_freq_eli['methods']):
            logging.error("Failed to eliminate imaginary frequency")
            return FWAction(
                stored_data={'task_id': t_id},
                defuse_children=True,
                update_spec={'mol': d["molecule_final"],
                             'egsnl': fw_spec['egsnl'],
                             'snlgroup_id': fw_spec['snlgroup_id'],
                             'inchi_root': fw_spec["inchi_root"],
                             'defuse_reason': "imaginary frequency "
                                              "elimination failed",
                             'offending_fw_id': fw_spec['fw_id']})

        new_mol = self.perturb_molecule(d)
        old_mol = Molecule.from_dict(d['molecule_final'])
        structure_changed = self._check_structure_change(
            old_mol, new_mol, fw_spec)
        if structure_changed:
            self.perturb_molecule(d, reversed_direction=True)
            structure_changed = self._check_structure_change(
                old_mol, new_mol, fw_spec)
        if structure_changed:
            return FWAction(
                stored_data={'task_id': t_id},
                defuse_children=True,
                update_spec={
                    'mol': d["molecule_final"],
                    'perturbed_mol': new_mol.to_dict,
                    'egsnl': fw_spec['egsnl'],
                    'snlgroup_id': fw_spec['snlgroup_id'],
                    'inchi_root': fw_spec["inchi_root"],
                    'defuse_reason': "structural change in imaginary "
                                     "frequency elimination",
                    'offending_fw_id': fw_spec['fw_id']})
        molname = d['user_tags']['molname']
        mission = d['user_tags']['mission']
        additional_user_tags = {"img_freq_eli": img_freq_eli}
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
        update_specs = {'mol': new_mol.to_dict,
                        'egsnl': egsnl.to_dict,
                        'snlgroup_id': fw_spec['snlgroup_id'],
                        'inchi_root': fw_spec['inchi_root']}

        method = img_freq_eli["methods"][img_freq_eli["current_method_id"]]
        charge = new_mol.charge
        spin_multiplicity = new_mol.spin_multiplicity
        if method == "dir_dis_opt":
            logging.info("Eliminate Imaginary Frequency By Perturbing the "
                         "Structure of Molecule")
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge,
                                        spin_multiplicity,
                                        grid=None)
        elif method == "den_dis_opt":
            logging.info("Eliminate Imaginary Frequency By Perturbing the "
                         "Structure of Molecule, and increase the grid "
                         "density")
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge,
                                        spin_multiplicity,
                                        grid=(128, 302))
        elif method == "alt_den_dis_opt":
            logging.info("Eliminate Imaginary Frequency By Perturbing the "
                         "Structure of Molecule, and increase the grid "
                         "density")
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge,
                                        spin_multiplicity,
                                        grid=(90, 590))
        else:
            raise Exception("Unknown imaginary frequency fixing method")

        return FWAction(stored_data={'task_id': t_id}, detours=wf,
                        update_spec=update_specs)


class QChemSinglePointEnergyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem Single Point Energy DB Insertion Task"

    def run_task(self, fw_spec):
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
        t_id, d = drone.assimilate(os.path.abspath(
            os.path.join(prev_dir, "mol.qcout")), fw_spec=fw_spec)

        if d["state"] == "successful":
            return FWAction(
                stored_data={'task_id': t_id},
                update_spec={
                    'mol': d["molecule_final"],
                    'egsnl': d['snl_final'],
                    'snlgroup_id': d['snlgroup_id_final'],
                    'inchi_root': fw_spec["inchi_root"]})
        else:
            return FWAction(
                stored_data={'task_id': t_id},
                defuse_children=True,
                update_spec={
                    'mol': d["molecule_final"],
                    'egsnl': d['snl_final'],
                    'snlgroup_id': d['snlgroup_id_final'],
                    'inchi_root': fw_spec["inchi_root"],
                    'defuse_reason': 'SCF failed',
                    'offending_fw_id': fw_spec['fw_id']})
