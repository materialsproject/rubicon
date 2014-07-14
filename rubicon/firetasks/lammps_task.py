import subprocess

from pymatgen import Molecule
try:
    # work around before packmole is merged into master
    from pymatgen.packmol.packmol import PackmolRunner
except:
    pass

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

    ##required_params = ["molecules"]
    #optional_params = []
    _fw_name = "Lammps Input Writer"


    def run_task(self, fw_spec):
        m = Molecule.from_dict(fw_spec["molecules"][0])
        print m
        ffmol = AntechamberRunner(m)
        mols_in_box = PackmolRunner(m, [{"number":1,"inside box":[0.,0.,0.,40.,40.,40.]},{"number":1},{"number":1}])
        data_lammps=LmpInput(ffmol,mols_in_box)
        data_lammps.write_lammps_data('mol_data.lammps')
        control_lammps = DictLammpsInputSet()
        control_lammps.get_lammps_control('Lammps.json',ensemble='nvt',temp=350)
        control_lammps.write_lampps_control('mol_control.lammps')
        subprocess.check_call("lmp_hopper <  mol_control.lammps")
#
#
# @explicit_serialize
# class LamppsCustodianTask(FireTaskBase):
#     """
#     Runs LAMMPS using Custodian.
#
#     Required Params:
#
#     """
#     required_params = []
#
#     optional_params = []
#
#     def run_task(self, fw_spec):
#         subprocess.check_call("lammps < lmps.inp")
#
# @explicit_serialize
# class LammpsAnalyzeTask(FireTaskBase):
#     """
#
#     """
#
#     optional_params = ["vasprun_fname"]
#
#     def run_task(self, fw_spec):
#         pass
