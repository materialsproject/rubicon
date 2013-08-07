import shlex
import subprocess
from fireworks.core.firework import FireTaskBase
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen.io.nwchemio import NwInput
from custodian.custodian import Custodian
from custodian.nwchem import handlers
from custodian.nwchem import jobs


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
        nwc_exe = shlex.split('aprun -n 24 nwchem nwchem.nw')
        #subprocess.call(nwc_exe)
        c= Custodian(NwchemJob("nwchem"),NwchemErrorHandler())
        c.run()