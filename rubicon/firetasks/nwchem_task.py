import json
import logging
import shlex
import os
import socket
import datetime
import sys
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.core.fw_config import FWData
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen.io.nwchemio import NwInput

from custodian.custodian import Custodian
from custodian.nwchem.handlers import NwchemErrorHandler
from custodian.nwchem.jobs import NwchemJob
from rubicon.borg.hive import DeltaSCFNwChemToDbTaskDrone

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Jun 07, 2013'


DATETIME_HANDLER = lambda obj: obj.isoformat() if isinstance(obj, datetime.datetime) else None

class NWChemTask(FireTaskBase, FWSerializable):
    """
    Write NWChem task and run NWChem
    """

    _fw_name = "NWChem Task"

    def run_task(self, fw_spec):
        nwi = NwInput.from_dict(fw_spec)
        nwi.write_file('mol.nw')

        fw_data = FWData()

        if 'macqu.dhcp.lbl.gov' == socket.gethostname() \
            or 'MacQu.local' == socket.gethostname() \
            or 'macqu.dhcp.lbnl.us' == socket.gethostname(): # Xiaohui's Laptop
            nwc_exe = ['nwchem']
        elif 'nid' in socket.gethostname():  # hopper compute nodes
            # TODO: can base ncores on FW_submit.script
            if (not fw_data.MULTIPROCESSING) or (fw_data.NODE_LIST is None):
                nwc_exe = shlex.split('aprun -n 24 nwchem')
            else:
                list_str = ','.join(fw_data.NODE_LIST)
                num_str = str(fw_data.SUB_NPROCS)
                nwc_exe = shlex.split('aprun -n ' + num_str +
                                      ' -L ' + list_str + ' nwchem')
        elif 'c' in socket.gethostname():  # mendel compute nodes
            # TODO: can base ncores on FW_submit.script
            if (not fw_data.MULTIPROCESSING) or (fw_data.NODE_LIST is None):
                nwc_exe = shlex.split('mpirun -n 16 nwchem')
            else:
                list_str = ','.join(fw_data.NODE_LIST)
                num_str = str(fw_data.SUB_NPROCS)
                nwc_exe = shlex.split('mpirun -n ' + num_str +
                                      ' --host ' + list_str + ' nwchem')

        job = NwchemJob(nwchem_cmd=nwc_exe)
        handler = NwchemErrorHandler()
        c = Custodian(handlers=[handler], jobs=[job])
        c.run()


class NWDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem DB Insertion Task"

    def run_task(self, fw_spec):


        # get the directory containing the db file
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'molecules_db.json')

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
            t_id = drone.assimilate(os.path.abspath(os.path.join(os.getcwd(), "mol.nwout")))

        if t_id:
            return FWAction(stored_data={'task_id': t_id})

