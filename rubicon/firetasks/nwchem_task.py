import shlex
import subprocess
from fireworks.core.firework import FireTaskBase
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen.io.nwchemio import NwInput

from custodian.custodian import Custodian
from custodian.nwchem.handlers import NwchemErrorHandler
from custodian.nwchem.jobs import NwchemJob

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

        if 'nid' in socket.gethostname():  # hopper compute nodes
            # TODO: can base ncores on FW_submit.script
            nwc_exe = shlex.split('aprun -n 24 nwchem')
            print 'running on HOPPER'
        elif 'c' in socket.gethostname():  # mendel compute nodes
            # TODO: can base ncores on FW_submit.script
            nwc_exe = shlex.split('mpirun -n 16 nwchem')

        # nwc_exe = shlex.split('aprun -n 24 nwchem nwchem.nw')
        # subprocess.call(nwc_exe)

        job = NwchemJob(nwchem_cmd=nwc_exe)
        handler = NwchemErrorHandler()
        c = Custodian(handlers=[handler], jobs=[job])
        c.run()
