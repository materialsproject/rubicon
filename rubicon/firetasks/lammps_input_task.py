import shlex
import subprocess
from monty import logging
from pymatgen import Molecule
try:
    # just a walkaround before the packmol is merged to master branch
    # after packmol is merged to master branch, the try...catch block
    # should be removed
    from pymatgen.packmol.packmol import PackmolRunner
except:
    pass
from rubicon.gff.boxmol import BoxMol
from rubicon.gff.lammpsin import DictLammpsInputSet
from rubicon.gff.lamppsio import LmpInput
from rubicon.gff.antechamberio import AntechamberRunner


__author__ = 'navnidhirajput'


from fireworks import FireTaskBase, explicit_serialize


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

    _fw_name = "Lammps Input Writer"


    def run_task(self, fw_spec):
        mols_dict = fw_spec["molecules"]
        mols = [Molecule.from_dict(m) for m in mols_dict]
        ffmol_list = []
        for mol in mols:
            acr = AntechamberRunner(mol)
            ffmol_list.append(acr.get_ff_top_mol(mol,'mol.pdb'))

        pmr = PackmolRunner(mols, [{"number":1,"inside box":[0.,0.,0.,40.,40.,40.]},{"number":1},{"number":1}])
        mols_coord = pmr.run()
        boxmol= BoxMol.from_packmol(pmr, mols_coord)
        data_lammps=LmpInput(ffmol_list, boxmol)
        data_lammps.write_lammps_data('mol_data.lammps')
        control_lammps = DictLammpsInputSet()
        control_lammps.get_lammps_control('Lammps.json',ensemble='nvt',temp=350)
        control_lammps.write_lampps_control('mol_control.lammps')
        #subprocess.check_call(shlex.split("lmp_mac <  mol_control.lammps"))
        #subprocess.check_call(shlex.split("lmp_hopper <  mol_control.lammps"))
        with open("mol_control.lammps") as f:
            subprocess.check_call(shlex.split("aprun -n 48 lmp_hopper"), stdin=f)


