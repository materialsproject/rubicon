import shlex
import socket
from fireworks.core.firework import FireTaskBase
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen.io.nwchemio import NwInput

from custodian.custodian import Custodian
from custodian.nwchem.handlers import NwchemErrorHandler
from custodian.nwchem.jobs import NwchemJob
import os
from rubicon.borg.hive import DeltaSCFNwChemToDbTaskDrone

from fireworks.core.firework import FWAction
from fireworks.core.fw_config import FWConfig

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Jun 07, 2013'



class NWChemTask(FireTaskBase, FWSerializable):
    """
    Write NWChem task and run NWChem
    """

    _fw_name = "NWChem Task"

    def run_task(self, fw_spec):
        nwi = NwInput.from_dict(fw_spec)
        nwi.write_file('mol.nw')

        fw_conf = FWConfig()

        if 'macqu.dhcp.lbl.gov' == socket.gethostname() \
            or 'MacQu.local' == socket.gethostname(): # Xiaohui's Laptop
            nwc_exe = ['nwchem']
        elif 'nid' in socket.gethostname():  # hopper compute nodes
            # TODO: can base ncores on FW_submit.script
            if (not fw_conf.MULTIPROCESSING) or (fw_conf.NODE_LIST is None):
                nwc_exe = shlex.split('aprun -n 24 nwchem')
            else:
                list_str = ','.join(fw_conf.NODE_LIST)
                num_str = str(24*len(fw_conf.NODE_LIST))
                nwc_exe = shlex.split('aprun -n ' + num_str +
                                      ' -L ' + list_str + ' nwchem')
            print 'running on HOPPER'
        elif 'c' in socket.gethostname():  # mendel compute nodes
            # TODO: can base ncores on FW_submit.script
            if (not fw_conf.MULTIPROCESSING) or (fw_conf.NODE_LIST is None):
                nwc_exe = shlex.split('mpirun -n 16 nwchem')
            else:
                list_str = ','.join(fw_conf.NODE_LIST)
                num_str = str(len(fw_conf.NODE_LIST))
                nwc_exe = shlex.split('mpirun -n ' + num_str +
                                      ' --host ' + list_str + ' nwchem')

        job = NwchemJob(nwchem_cmd=nwc_exe)
        handler = NwchemErrorHandler()
        c = Custodian(handlers=[handler], jobs=[job])
        c.run()
        curdir = os.getcwd()
        drone = DeltaSCFNwChemToDbTaskDrone()
        fwa = FWAction()
        d = drone.assimilate(curdir + "/mol.nwout")
        fwa.stored_data.update({"Result": d})
        return fwa
