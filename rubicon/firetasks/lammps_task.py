import shlex
import shutil
import subprocess
from fireworks.core.firework import Workflow
from pymatgen import Molecule
from rubicon.gff.lammpsin import DictLammpsInputSet
from rubicon.gff.lamppsio import LmpInput
from pymatgen.packmol.packmol import PackmolRunner
from rubicon.gff.antechamberio import AntechamberRunner

__author__ = 'navnidhirajput'


from fireworks import FireTaskBase, FWAction, explicit_serialize, FireWork, LaunchPad
from custodian import Custodian

@explicit_serialize
class WritelammpsInputTask(FireTaskBase):
    """
    Writes LAMMPS Input files.

    Required params:

        lammps_input_set (str): A string name for the VASP input set. E.g.,
            "MPVaspInputSet" or "MITVaspInputSet".

    Optional params:
        input_set_params (dict): If the input set requires some additional
            parameters, specify them using input_set_params. E.g.,
            {"user_incar_settings": ...}.
    """

    ##required_params = ["molecules"]
    #optional_params = []
    _fw_name = "Lammps Input Writer"


    def run_task(self, fw_spec):
        mols_dict = fw_spec["molecules"]
        mols = [Molecule.from_dict(m) for m in mols_dict]
        ffmol = AntechamberRunner(mols)
        gff_list, top_list = ffmol.get_ff_top_mol(mols,'mol.pdb')
        mols_in_box = PackmolRunner(mols, [{"number":1,"inside box":[0.,0.,0.,40.,40.,40.]},{"number":1},{"number":1}])
        data_lammps=LmpInput(ffmol,mols_in_box)
        data_lammps.write_lammps_data('mol_data.lammps')
        control_lammps = DictLammpsInputSet()
        control_lammps.get_lammps_control('Lammps.json',ensemble='nvt',temp=350)
        control_lammps.write_lampps_control('mol_control.lammps')
        subprocess.check_call(shlex.split("lmp_mac <  mol_control.lammps"))

