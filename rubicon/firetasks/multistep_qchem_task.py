import copy
import json
import logging
import os
import sys
import math

from fireworks.core.firework import FireTaskBase, FWAction, Workflow
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen import Molecule
from pymatgen.io.qchemio import QcBatchInput

from rubicon.borg.hive import DeltaSCFQChemToDbTaskDrone


__author__ = 'xiaohuiqu'


class QChemGeomOptDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem Geometry Optimization DB Insertion Task"

    def run_task(self, fw_spec):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'qchem_calc_db.json')

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
        assi_result = drone.assimilate(os.path.abspath(
            os.path.join(os.getcwd(), "mol.qcout")))

        t_id = None
        d = None
        if assi_result:
            t_id, d = assi_result
        if t_id:
            if d["state"] == "successful":
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': fw_spec['egsnl'],
                                             'snlgroup_id':
                                             fw_spec['snlgroup_id']})
            else:
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': fw_spec['egsnl'],
                                             'snlgroup_id':
                                             fw_spec['snlgroup_id']},
                                defuse_children=True)
        else:
            return FWAction(defuse_children=True,
                            update_spec={'egsnl': fw_spec['egsnl'],
                                         'snlgroup_id': fw_spec['snlgroup_id']})


class QChemFrequencyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem Frequency DB Insertion Task"
    molecule_perturb_scale = 0.3

    def run_task(self, fw_spec):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'qchem_calc_db.json')

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
        assi_result = drone.assimilate(os.path.abspath(
            os.path.join(os.getcwd(), "mol.qcout")))

        t_id = None
        d = None
        if assi_result:
            t_id, d = assi_result
        if t_id:
            if d["state"] == "successful":
                if d['stationary_type'] == 'minimum':
                    return FWAction(stored_data={'task_id': t_id},
                                    update_spec={"mol": d["final_molecule"],
                                                 'egsnl': fw_spec['egsnl'],
                                                 'snlgroup_id':
                                                 fw_spec['snlgroup_id']})
                else:
                    return self.img_freq_action(fw_spec, d, t_id)
            else:
                return FWAction(stored_data={'task_id': t_id},
                                defuse_children=True,
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': fw_spec['egsnl'],
                                             'snlgroup_id':
                                             fw_spec['snlgroup_id']})
        else:
            return FWAction(defuse_children=True,
                            update_spec={'egsnl': fw_spec['egsnl'],
                                         'snlgroup_id': fw_spec['snlgroup_id']})

    @staticmethod
    def spawn_opt_freq_wf(mol, molname, mission, additional_user_tags,
                          priority, update_spec, charge_shift, grid=None):
        from rubicon.workflows.multistep_ipea_wf \
            import QChemFireWorkCreator
        fw_creator = QChemFireWorkCreator(mol, molname, mission,
                                          additional_user_tags, priority,
                                          update_spec)
        geom_fwid, freq_fwid = -1, -2
        geom_fw = fw_creator.geom_fw(charge_shift, geom_fwid)
        freq_fw = fw_creator.freq_fw(charge_shift, freq_fwid)
        if grid:
            for fw in [geom_fw, freq_fw]:
                qcinp = QcBatchInput.from_dict(fw.spec["qcinp"])
                for j in qcinp.jobs:
                    j.set_dft_grids(*grid)
                fw.spec["qcinp"] = qcinp.to_dict
        wf = Workflow([geom_fw, freq_fw],
                      links_dict={geom_fwid: freq_fwid})
        return wf

    @classmethod
    def perturb_molecule(cls, d):
        old_mol = Molecule.from_dict(d['final_molecule'])
        vib_mode = d['calculations']['freq']['frequencies'][0]["vib_mode"]
        max_dis = max([math.sqrt(sum([x ** 2 for x in mode]))
                       for mode in vib_mode])
        scale = cls.molecule_perturb_scale / max_dis
        normalized_mode = [[x * scale for x in mode]
                           for mode in vib_mode]
        new_coords = [[c+v for c, v in zip(site.coords, mode)]
                      for site, mode in zip(old_mol.sites, normalized_mode)]
        species = [site.specie.symbol
                   for site in old_mol.sites]
        charge_shift = d['user_tags']['charge_shift']
        charge = old_mol.charge - charge_shift
        new_mol = Molecule(species, new_coords, charge)
        return new_mol

    def img_freq_action(self, fw_spec, d, t_id):
        if "img_freq_eli" in d['user_tags']:
            img_freq_eli = copy.deepcopy(d['user_tags'])
            img_freq_eli['current_method_id'] += 1
        else:
            img_freq_eli = {"methods": ["dir_dis_opt", "den_dis_opt",
                                        "alt_den_dis_opt"],
                            "current_method_id": 0}

        if img_freq_eli['current_method_id'] >= len(img_freq_eli['methods']):
            return FWAction(stored_data={'task_id': t_id},
                            defuse_children=True,
                            update_spec={"mol": d["final_molecule"],
                                         'egsnl': fw_spec['egsnl'],
                                         'snlgroup_id': fw_spec['snlgroup_id']})

        new_mol = self.perturb_molecule(d)
        molname = d['user_tags']['molname']
        mission = d['user_tags']['mission']
        additional_user_tags = {"img_freq_eli": img_freq_eli}
        priority = fw_spec['_priotity']
        update_specs = {'egsnl': fw_spec['egsnl'],
                        'snlgroup_id': fw_spec['snlgroup_id']}
        charge_shift = d['user_tags']['charge_shift']

        method = img_freq_eli["methods"][img_freq_eli["current_method_id"]]
        if method == "dir_dis_opt":
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge_shift,
                                        grid=None)
        elif method == "den_dis_opt":
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge_shift,
                                        grid=(128, 302))
        elif method == "alt_den_dis_opt":
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge_shift,
                                        grid=(90, 590))
        else:
            raise Exception("Unknown imaginary frequency fixing method")

        return FWAction(stored_data={'task_id': t_id}, detours=wf,
                        update_spec={'egsnl': fw_spec['egsnl'],
                                     'snlgroup_id': fw_spec['snlgroup_id']})


class QChemSinglePointEnergyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem Single Point Energy DB Insertion Task"

    def run_task(self, fw_spec):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'qchem_calc_db.json')

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
        assi_result = drone.assimilate(os.path.abspath(
            os.path.join(os.getcwd(), "mol.qcout")))

        t_id = None
        d = None
        if assi_result:
            t_id, d = assi_result
        if t_id:
            if d["state"] == "successful":
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': fw_spec['egsnl'],
                                             'snlgroup_id':
                                             fw_spec['snlgroup_id']})
            else:
                return FWAction(stored_data={'task_id': t_id},
                                defuse_children=True,
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': fw_spec['egsnl'],
                                             'snlgroup_id':
                                             fw_spec['snlgroup_id']})
        else:
            return FWAction(defuse_children=True,
                            update_spec={'egsnl': fw_spec['egsnl'],
                                         'snlgroup_id':
                                         fw_spec['snlgroup_id']})