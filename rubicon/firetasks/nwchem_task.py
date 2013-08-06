import shlex
import socket
import subprocess
from fireworks.core.firework import FireTaskBase
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen.io.nwchemio import NwInput

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
        nwi.write_file('nwchem.nw')

        # TODO: replace with a custodian
        if 'nid' in socket.gethostname():  # hopper compute nodes
            # TODO: can base ncores on FW_submit.script
            nwc_exe = shlex.split('aprun -n 24 nwchem nwchem.nw')
            print 'running on HOPPER'
        elif 'c' in socket.gethostname():  # mendel compute nodes
            # TODO: can base ncores on FW_submit.script
            nwc_exe = shlex.split('mpirun -n 16 nwchem nwchem.nw')

        subprocess.call(nwc_exe)