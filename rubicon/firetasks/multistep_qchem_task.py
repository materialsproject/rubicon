import json
import logging
import os
from fireworks.core.firework import FireTaskBase, FWAction, Workflow
from fireworks.utilities.fw_serializers import FWSerializable
import sys
from pymatgen import Molecule
from pymongo import MongoClient
from rubicon.borg.hive import DeltaSCFNwChemToDbTaskDrone, DeltaSCFQChemToDbTaskDrone

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
                                             'snlgroup_id': fw_spec['snlgroup_id']})
            else:
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': fw_spec['egsnl'],
                                             'snlgroup_id': fw_spec['snlgroup_id']},
                                defuse_children=True)
        else:
            return FWAction(defuse_children=True,
                            update_spec={'egsnl': fw_spec['egsnl'],
                                         'snlgroup_id': fw_spec['snlgroup_id']})



class QChemFrequencyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem Frequency DB Insertion Task"

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
                                                 'snlgroup_id': fw_spec['snlgroup_id']})
                else:
                    return self.img_freq_action(fw_spec, d, t_id)
            else:
                return FWAction(stored_data={'task_id': t_id},
                                defuse_children=True,
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': fw_spec['egsnl'],
                                             'snlgroup_id': fw_spec['snlgroup_id']})
        else:
            return FWAction(defuse_children=True,
                            update_spec={'egsnl': fw_spec['egsnl'],
                                         'snlgroup_id': fw_spec['snlgroup_id']})

    def img_freq_action(self, fw_spec, d, t_id):
        if 'img_freq_max_trial' in fw_spec:
            max_fix_time = fw_spec['img_freq_max_trial']
        else:
            max_fix_time = 3
        if 'freq_fix_time' in d['user_tags']:
            fix_time = d['user_tags']['freq_fix_time'] + 1
        else:
            fix_time = 1
            db_dir = os.environ['DB_LOC']
            db_path = os.path.join(db_dir, 'nwchem_calc_db.json')
            with open(db_path) as f:
                db_creds = json.load(f)
            conn = MongoClient(db_creds['host'], db_creds['port'])
            db = conn[db_creds['database']]
            if db_creds['admin_user']:
                db.authenticate(db_creds['admin_user'],
                                db_creds['admin_password'])
            coll = db[db_creds['collection']]
            coll.update({"task_id": t_id},
                        {"$set": {'user_tags.freq_fix_time': 0}}, upsert=True)

        if fix_time > max_fix_time:
            return FWAction(stored_data={'task_id': t_id},
                            defuse_children=True,
                            update_spec={"mol": d["final_molecule"],
                                         'egsnl': fw_spec['egsnl'],
                                         'snlgroup_id': fw_spec['snlgroup_id']})
        else:
            old_mol = Molecule.from_dict(d['final_molecule'])
            vib_mode = d['calculations']['freq']['frequencies'][0][1]
            new_coords = [[c+v for c, v in zip(site.coords, mode)]
                          for site, mode in zip(old_mol.sites, vib_mode)]
            species = [site.specie.symbol
                       for site in old_mol.sites]
            new_mol = Molecule(species, new_coords)
            from rubicon.workflows.multistep_ipea_wf \
                import NWChemFireWorkCreator
            fw_creator = NWChemFireWorkCreator(new_mol,
                             d['user_tags']['molname'],
                             d['user_tags']['mission'],
                             additional_user_tags={'freq_fix_time': fix_time},
                             priority=fw_spec['_priority'],
                             update_spec={'egsnl': fw_spec['egsnl'],
                                          'snlgroup_id': fw_spec['snlgroup_id']},)
            geom_fwid, freq_fwid = -1, -2
            geom_fw = fw_creator.geom_fw(d['user_tags']['charge_shift'], geom_fwid)
            freq_fw = fw_creator.freq_fw(d['user_tags']['charge_shift'], freq_fwid)
            wf = Workflow([geom_fw, freq_fw],
                          links_dict={geom_fwid: freq_fwid})
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
                                             'snlgroup_id': fw_spec['snlgroup_id']})
            else:
                return FWAction(stored_data={'task_id': t_id},
                                defuse_children=True,
                                update_spec={"mol": d["final_molecule"],
                                             'egsnl': fw_spec['egsnl'],
                                             'snlgroup_id': fw_spec['snlgroup_id']})
        else:
            return FWAction(defuse_children=True,
                            update_spec={'egsnl': fw_spec['egsnl'],
                                         'snlgroup_id': fw_spec['snlgroup_id']})

