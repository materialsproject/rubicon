import os
import shlex
import subprocess
from pymatgen import Molecule
try:
    # just a walkaround before the packmol is merged to master branch
    # after packmol is merged to master branch, the try...catch block
    # should be removed
    from pymatgen.packmol.packmol import PackmolRunner
except:
    pass
from rubicon.gff.boxmol import BoxMol
from rubicon.gff.lammps_control import DictLammpsInputSet
from rubicon.gff.lammps_data import LmpInput
from rubicon.gff.antechamberio import AntechamberRunner


__author__ = 'navnidhirajput'


from fireworks import FireTaskBase, explicit_serialize, FWAction


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
        filename = fw_spec['prev_gaussian_freq']
        mols_dict = fw_spec["molecule"]
        mol = mols_dict
        ffmol_list = []
        #for mol in mols:
        acr = AntechamberRunner(mol)
        ffmol_list.append(acr.get_ff_top_mol(mol,filename))
        #ffmol_list.append(acr.get_ff_top_mol(mol,'mol.pdb'))
        pmr = PackmolRunner([mol], [{"number":100,"inside box":[0.,0.,0.,50.,50.,50.]}])
        #pmr = PackmolRunner(mols, [{"number":15,"inside box":[0.,0.,0.,100.,100.,100.]},{"number":30,"inside box":[0.,0.,0.,100.,100.,100.]},{"number":232,"inside box":[0.,0.,0.,100.,100.,100.]}])

        mols_coord = pmr.run()
        boxmol= BoxMol.from_packmol(pmr, mols_coord)

        data_lammps=LmpInput(ffmol_list, boxmol)
        data_lammps.write_lammps_data('mol_data.lammps')
        control_lammps = DictLammpsInputSet()
        control_lammps.get_lammps_control('Lammps.json',ensemble='npt',temp=300)
        control_lammps.write_lampps_control('mol_control.lammps')

        with open("mol_control.lammps") as f:
            subprocess.call(shlex.split("srun -n 48 lmp_edison"), stdin=f)

        prev_lammps_log = os.path.join(os.getcwd(), 'mol.log')
        prev_lammps_trj = os.path.join(os.getcwd(), "mol.lammpstrj")
        prev_lammps_data = os.path.join(os.getcwd(), "mol_data.lammps")

        update_spec = {'prev_lammps_trj': prev_lammps_trj,
                       'prev_lammps_data': prev_lammps_data,
                       'prev_lammps_log': prev_lammps_log}

        return FWAction(update_spec=update_spec)


if __name__ == '__main__':
     test = WritelammpsInputTask()
     test.run_task()




