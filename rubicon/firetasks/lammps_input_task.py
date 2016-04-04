# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import shlex
import subprocess

from rubicon.io.lammps.inputs import LammpsAmberData
from rubicon.io.lammps.sets import DictLammpsInputSet_2
from rubicon.io.packmol.packmol import PackmolRunner

from fireworks import FireTaskBase, explicit_serialize, FWAction


__author__ = 'Navnidhi Rajput, Kiran Mathew'



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
        gaussian_file = fw_spec['prev_gaussian_freq']
        mols_dict = fw_spec["molecule"]
        mol = mols_dict
        molecules = [mol]
        param_list = [{"number": 100,
                       "inside box": [-14.82, -14.82, -14.82,
                                      14.82, 14.82, 14.82]}]

        pmr = PackmolRunner(molecules, param_list)
        packed_molecule = pmr.run()
        data_lammps = LammpsAmberData(molecules, param_list["number"],
                                      param_list["inside box"],
                                      packed_molecule, gaussian_file)
        data_lammps.to('mol_data.lammps')
        control_lammps = DictLammpsInputSet_2()
        # control_lammps.get_lammps_control('Lammps.json',ensemble='npt',temp=300)
        control_lammps.get_lammps_control('Lammps.json', ensemble1='npt',
                                          ensemble2='nvt', temp=298)
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
