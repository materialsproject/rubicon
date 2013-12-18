import json
import logging
import re
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
from pymatgen.core.structure import Molecule
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
        qcbat = QcBatchInput.from_dict(fw_spec["qcinp"])
        if 'mol' in fw_spec:
            mol = Molecule.from_dict(fw_spec["mol"])
            qcbat.jobs[0].mol = mol
        qcbat.write_file("mol.qcinp")
        hopper_name_pattern = re.compile("nid\d+")
        carver_name_pattern = re.compile("c[0-9]{4}-ib")
        if hopper_name_pattern.match(socket.gethostname()):
        # hopper compute nodes
            qc_exe = shlex.split("qchem -np 24")
        elif carver_name_pattern.match(socket.gethostname()):
        # mendel compute nodes
            qc_exe = shlex.split("qchem -np 8")
        else:
            qc_exe = ["qchem"]

        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('QChemDrone')
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        logger.addHandler(sh)

        job = QchemJob(qc_exe, input_file="mol.qcinp", output_file="mol.qcout",
                       qclog_file="mol.qclog")
        handler = QChemErrorHandler()
        c = Custodian(handlers=[handler], jobs=[job], max_errors=50)
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
