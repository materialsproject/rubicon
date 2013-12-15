import json
import logging
import shlex
import os
import socket
import datetime
import sys
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QchemJob
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import FWSerializable

from custodian.custodian import Custodian
from pymatgen.io.qchemio import QcBatchInput
from rubicon.borg.hive import DeltaSCFQChemToDbTaskDrone

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Jun 07, 2013'


DATETIME_HANDLER = lambda obj: obj.isoformat() \
    if isinstance(obj, datetime.datetime) else None


class QChemTask(FireTaskBase, FWSerializable):
    """
    Write QChem task and run QChem
    """

    _fw_name = "QChem Task"

    def run_task(self, fw_spec):
        qcbat = QcBatchInput.from_dict(fw_spec["qcbat"])
        qcbat.write_file("mol.qcinp")
        if 'nid' in socket.gethostname():  # hopper compute nodes
            qc_exe = shlex.split("qchem -np 24")
        elif 'c' in socket.gethostname():  # mendel compute nodes
            qc_exe = shlex.split("qchem -np 8")
        else:
            qc_exe = ["qchem"]

        job = QchemJob(qc_exe)
        handler = QChemErrorHandler()
        c = Custodian(handlers=[handler], jobs=[job])
        c.run()


class QCDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "QChem DB Insertion Task"

    def run_task(self, fw_spec):
        # get the directory containing the db file
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'molecules_db.json')

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
        if assi_result:
            t_id, d = assi_result
        if t_id:
            return FWAction(stored_data={'task_id': t_id})
