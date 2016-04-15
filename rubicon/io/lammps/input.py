# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module implements classes for generating Lammps input files.
The input files consist of the main input file with the control parameters and
the data file.
"""

import json
import os
from functools import partial
from collections import OrderedDict

from monty.json import MSONable

from rubicon.io.lammps.data import LammpsData, LammpsForceFieldData


__author__ = 'Kiran Mathew, Navnidhi Rajput'
__email__ = "kmathew@lbl.gov"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class DictLammpsInput(MSONable):
    """
    Implementation of LammpsInputSet that is initialized from a dict
    settings. It is typically used by other LammpsInputSets for
    initialization from json or yaml source files.

    Args:
        name (str): A name for the input set.
        config_dict (dict): The config dictionary to use.
        lammps_data (LammpsData): LammpsData object
        data_filename (str): name of the the lammps data file
        user_lammps_settings (dict): User lammps settings. This allows a user
            to override lammps settings, e.g., setting a different force field
            or bond type.
    """

    def __init__(self, name, config_dict, lammps_data=None, data_filename="in.data",
                 user_lammps_settings={}):
        self.name = name
        self.lines = []
        self.config_dict = config_dict
        self.lammps_data = lammps_data
        self.data_filename = data_filename
        self.config_dict["read_data"] = data_filename
        self.user_lammps_settings = user_lammps_settings
        if self.user_lammps_settings:
            self.config_dict.update(self.user_lammps_settings)

    @staticmethod
    def get_lines_from_dict(in_dict):
        """
        recursively convert the nested dict to an appropriate string
        representation to be used as lammps input file with the control paramters.
        """
        lines = " "
        for k1, v1 in in_dict.items():
            lines = lines + "{} ".format(k1)
            if isinstance(v1,dict):
                lines = lines + DictLammpsInput.get_lines_from_dict(v1)
            elif isinstance(v1, list):
                lines = lines + " ".join([str(x) for x in v1]) + os.linesep
            else:
                lines = lines + " {}{}".format(str(v1), os.linesep)
        return lines

    def __str__(self):
        """
        string representation of the lammps input file with the control parameters
        """
        return self.get_lines_from_dict(self.config_dict)

    def write_input(self, filename):
        """
        get the string representation of the main input file and write it.
        Also writes the data file if the lammps_data attribute is set.

        Args:
            filename (string): name of the input file
        """
        # write the main input file
        with open(filename, 'w') as f:
            f.write(self.__str__())
        # write the data file if present
        if self.lammps_data:
            print("Writing the data to {}".format(self.data_filename))
            self.lammps_data.write_data_file(filename=self.data_filename)

    @staticmethod
    def from_file(name, filename, data_filename=None, is_forcefield=False, **kwargs):
        """
        Read in the input settings from json file as ordereddict. Also reads in the
        datafile if provided.
        Note: with monty.serialization.loadfn the order of paramters in the
        json file is not preserved

        Args:
            filename (string): name of the file with the lamps control paramters
            data_filename (string): name of the data file name
            is_forcefield (bool): whether the data file has forcefield and topology info
                in it

        Returns:
            DictLammpsInput
        """
        with open(filename) as f:
            config_dict = json.load(f, object_pairs_hook=OrderedDict)
        lammps_data = None
        if data_filename:
            if is_forcefield:
                lammps_data = LammpsForceFieldData.from_file(data_filename)
            else:
                lammps_data = LammpsData.from_file(data_filename)
        return DictLammpsInput(name, config_dict, lammps_data, data_filename, **kwargs)

    @property
    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "name": self.name,
                "config_dict": self.config_dict,
                "lammps_data": self.lammps_data,
                "data_filename": self.data_filename,
                "user_lammps_settings": self.user_lammps_settings}

    @classmethod
    def from_dict(cls, d):
        return DictLammpsInput(d["name"], d["config_dict"],
                               lammps_data=d.get("lammps_data"),
                               data_filename=d.get("data_filename"),
                               user_lammps_settings=d.get("user_lammps_settings"))


# NPT
NPTLammpsInput = partial(DictLammpsInput.from_file, "NPT",
                         os.path.join(MODULE_DIR, "data_files", "Lammps_npt.json"))


# NPT followed by NVT
NPTNVTLammpsInput = partial(DictLammpsInput.from_file, "NPT_NVT",
                            os.path.join(MODULE_DIR, "data_files","Lammps_npt_nvt.json"))
