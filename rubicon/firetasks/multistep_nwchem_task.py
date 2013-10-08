import json
import logging
import os
import datetime
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import FWSerializable
import sys
from pymongo import MongoClient
from rubicon.borg.hive import DeltaSCFNwChemToDbTaskDrone

__author__ = 'xiaohuiqu'


class NWChemGeomOptDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem Geometry Optimization DB Insertion Task"

    def run_task(self, fw_spec):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'nwchem_calc_db.json')

        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('NWChemDrone')
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        logger.addHandler(sh)

        with open(db_path) as f:
            db_creds = json.load(f)
        drone = DeltaSCFNwChemToDbTaskDrone(
            host=db_creds['host'], port=db_creds['port'],
            database=db_creds['database'], user=db_creds['admin_user'],
            password=db_creds['admin_password'],
            collection=db_creds['collection'])
        assi_result = drone.assimilate(os.path.abspath(
                                os.path.join(os.getcwd(), "mol.nwout")))

        t_id = None
        d = None
        if assi_result:
            t_id, d = assi_result
        if t_id:
            if d["state"] == "successful":
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"]})
            else:
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"]},
                                defuse_children=True)
        else:
            return FWAction(defuse_children=True)



class NWChemFrequencyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem Frequency DB Insertion Task"

    def run_task(self, fw_spec):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'nwchem_calc_db.json')

        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('NWChemDrone')
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        logger.addHandler(sh)

        with open(db_path) as f:
            db_creds = json.load(f)
        drone = DeltaSCFNwChemToDbTaskDrone(
            host=db_creds['host'], port=db_creds['port'],
            database=db_creds['database'], user=db_creds['admin_user'],
            password=db_creds['admin_password'],
            collection=db_creds['collection'])
        assi_result = drone.assimilate(os.path.abspath(
                                os.path.join(os.getcwd(), "mol.nwout")))

        t_id = None
        d = None
        if assi_result:
            t_id, d = assi_result
        if t_id:
            if d["state"] == "successful":
                if d['stationary_type'] == 'minimum':
                    return FWAction(stored_data={'task_id': t_id})
                else:
                    return FWAction(stored_data={'task_id': t_id},
                                    defuse_children=True)
            else:
                return FWAction(stored_data={'task_id': t_id},
                                defuse_children=True)
        else:
            return FWAction(defuse_children=True)


class NWChemSinglePointEnergyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem Single Point Energy DB Insertion Task"

    def run_task(self, fw_spec):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'nwchem_calc_db.json')

        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('NWChemDrone')
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        logger.addHandler(sh)

        with open(db_path) as f:
            db_creds = json.load(f)
        drone = DeltaSCFNwChemToDbTaskDrone(
            host=db_creds['host'], port=db_creds['port'],
            database=db_creds['database'], user=db_creds['admin_user'],
            password=db_creds['admin_password'],
            collection=db_creds['collection'])
        assi_result = drone.assimilate(os.path.abspath(
                                os.path.join(os.getcwd(), "mol.nwout")))

        t_id = None
        d = None
        if assi_result:
            t_id, d = assi_result
        if t_id:
            if d["state"] == "successful":
                return FWAction(stored_data={'task_id': t_id})
            else:
                return FWAction(stored_data={'task_id': t_id},
                                defuse_children=True)
        else:
            return FWAction(defuse_children=True)

