import copy
import json
import logging
import os
import sys
import math

from fireworks.core.firework import FireTaskBase, FWAction, Workflow
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen import Molecule
from pymatgen.io.qchemio import QcInput

from rubicon.borg.hive import DeltaSCFQChemToDbTaskDrone
from rubicon.utils.snl.egsnl import EGStructureNL
from rubicon.utils.snl.egsnl_mongo import EGSNLMongoAdapter


__author__ = 'xiaohuiqu'


class QChemGeomOptDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem Geometry Optimization DB Insertion Task"

    def run_task(self, fw_spec):
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
        assi_result = drone.assimilate(os.path.abspath(
            os.path.join(os.getcwd(), "mol.qcout")))

        t_id = None
        d = None
        egsnl = None
        snlgroup_id = None
        if assi_result:
            t_id, d = assi_result
            d['eg_snl'] = fw_spec['egsnl']
            d['snlgroup_id'] = fw_spec['snlgroup_id']
            d['task_type'] = fw_spec['task_type']
            new_s = Molecule.from_dict(d["final_molecule"])
            old_snl = EGStructureNL.from_dict(d['egsnl'])
            history = old_snl.history
            history.append(
                {'name': 'Electrolyte Genome Project structure optimization',
                 'url': 'http://www.materialsproject.org',
                 'description': {'task_type': d['task_type'],
                                 'fw_id': d['fw_id'],
                                 'task_id': d['task_id']}})
            new_snl = EGStructureNL(new_s, old_snl.authors, old_snl.projects,
                                    old_snl.references, old_snl.remarks,
                                    old_snl.data, history)

            # enter new SNL into SNL db
            # get the SNL mongo adapter
            sma = EGSNLMongoAdapter.auto_load()

            # add snl
            egsnl, snlgroup_id = sma.add_snl(new_snl,
                                             snlgroup_guess=d['snlgroup_id'])
            d['snl_final'] = egsnl.to_dict
            d['snlgroup_id_final'] = snlgroup_id
            d['snlgroup_changed'] = (d['snlgroup_id'] !=
                                     d['snlgroup_id_final'])
        if t_id:
            if d["state"] == "successful":
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': egsnl,
                                             'snlgroup_id': snlgroup_id})
            else:
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': egsnl,
                                             'snlgroup_id': snlgroup_id},
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
        fw_creator = QChemFireWorkCreator(
            mol=mol, molname=molname, mission=mission,
            additional_user_tags=additional_user_tags, priority=priority,
            update_spec=update_spec)
        geom_fwid, freq_fwid = -1, -2
        geom_fw = fw_creator.geom_fw(charge_shift, geom_fwid)
        geom_fw.spec["run_tags"]["task_type_amend"] = "imaginary frequency " \
                                                      "elimination"
        freq_fw = fw_creator.freq_fw(charge_shift, freq_fwid)
        freq_fw.spec["run_tags"]["task_type_amend"] = "imaginary frequency " \
                                                      "elimination"
        if grid:
            for fw in [geom_fw, freq_fw]:
                qcinp = QcInput.from_dict(fw.spec["qcinp"])
                for j in qcinp.jobs:
                    j.set_dft_grid(*grid)
                    if j.params["rem"]["jobtype"] == "opt":
                        j.scale_geom_opt_threshold(0.1, 0.1, 0.1)
                        j.set_geom_max_iterations(100)
                fw.spec["qcinp"] = qcinp.to_dict
                fw.spec["run_tags"]["grid"] = grid
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
            logging.error("Failed to eliminate imaginary frequency")
            return FWAction(stored_data={'task_id': t_id},
                            defuse_children=True,
                            update_spec={"mol": d["final_molecule"],
                                         'egsnl': fw_spec['egsnl'],
                                         'snlgroup_id': fw_spec['snlgroup_id']})

        new_mol = self.perturb_molecule(d)
        molname = d['user_tags']['molname']
        mission = d['user_tags']['mission']
        additional_user_tags = {"img_freq_eli": img_freq_eli}
        priority = fw_spec['_priority']
        update_specs = {'egsnl': fw_spec['egsnl'],
                        'snlgroup_id': fw_spec['snlgroup_id']}
        charge_shift = d['user_tags']['charge_shift']

        method = img_freq_eli["methods"][img_freq_eli["current_method_id"]]
        if method == "dir_dis_opt":
            logging.info("Eliminate Imaginary Frequency By Perturbing the "
                         "Structure of Molecule")
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge_shift,
                                        grid=None)
        elif method == "den_dis_opt":
            logging.info("Eliminate Imaginary Frequency By Perturbing the "
                         "Structure of Molecule, and increase the grid "
                         "density")
            wf = self.spawn_opt_freq_wf(new_mol, molname, mission,
                                        additional_user_tags, priority,
                                        update_specs, charge_shift,
                                        grid=(128, 302))
        elif method == "alt_den_dis_opt":
            logging.info("Eliminate Imaginary Frequency By Perturbing the "
                         "Structure of Molecule, and increase the grid "
                         "density")
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