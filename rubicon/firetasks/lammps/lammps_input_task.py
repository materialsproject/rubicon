# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import shlex
import subprocess

from rubicon.io.lammps.data import LammpsForceFieldData
from rubicon.io.lammps.input import NPTNVTLammpsInput
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
        mol = fw_spec["molecule"]
        molecules = [mol]
        # list of gaussian files, one for each molecule type in molecules
        gaussian_files = [gaussian_file]
        param_list = [{"number": 100,
                       "inside box": [-14.82, -14.82, -14.82,
                                      14.82, 14.82, 14.82]}]

        pmr = PackmolRunner(molecules, param_list)
        packed_molecule = pmr.run()

        # lammps input
        control_filename = 'mol_control.lammps'
        data_filename = 'mol_data.lammps'
        # generate amber force filed data
        lammps_data = LammpsForceFieldData.from_amber(molecules, param_list["number"],
                                                      param_list["inside box"],
                                                      packed_molecule,
                                                      gaussian_files)
        user_lammps_settings = {
                                 "temp": 298,
                                 "fix1": {
                                     "style": "npt",
                                     "ID": "npt",
                                     "Tstart": 298,
                                     "Tstop": 298},
                                 "fix2": {
                                     "style": "nvt",
                                      "ID": "nvt"}
                                }
        lis = NPTNVTLammpsInput(lammps_data=lammps_data,
                                user_lammps_settings=user_lammps_settings)
        lis.write_input(control_filename, data_filename=data_filename)

        with open(control_filename) as f:
            subprocess.call(shlex.split("srun -n 48 lmp_edison"), stdin=f)
        # lammps output
        prev_lammps_log = os.path.join(os.getcwd(), 'mol.log')
        prev_lammps_trj = os.path.join(os.getcwd(), "mol.lammpstrj")
        prev_lammps_data = os.path.join(os.getcwd(), "mol_data.lammps")

        update_spec = {'prev_lammps_trj': prev_lammps_trj,
                       'prev_lammps_data': prev_lammps_data,
                       'prev_lammps_log': prev_lammps_log}

        return FWAction(update_spec=update_spec)
