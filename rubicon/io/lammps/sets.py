# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module implements classes for generating Lammps input sets.
The input set consists of the main input file with the control parameters and
the data file.
"""

import json
import os
from functools import partial
from collections import OrderedDict


__author__ = 'Kiran Mathew, Navnidhi Rajput'

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class DictLammpsInputSet(object):
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

    def __init__(self, name, config_dict, lammps_data=None, data_filename=None,
                 user_lammps_settings=None):
        self.name = name
        self.lines = []
        self.config_dict = config_dict
        self.lammps_data = lammps_data
        self.data_filename = data_filename
        self.config_dict["LAMMPSINNPT"]["read_data"] = data_filename
        self.user_lammps_settings = user_lammps_settings
        if user_lammps_settings:
            self.config_dict["LAMMPSINNPT"].update(user_lammps_settings)

    @staticmethod
    def get_lines_from_dict(in_dict):
        """
        recursively convert the nested dict to an appropriate string
        representation to be used as lammps input file.
        """
        lines = " "
        for k1, v1 in in_dict.items():
            lines = lines + "{} ".format(k1)
            if isinstance(v1,dict):
                lines = lines + DictLammpsInputSet.get_lines_from_dict(v1)
            elif isinstance(v1, list):
                lines = lines + " ".join([str(x) for x in v1]) + os.linesep
            else:
                lines = lines + " {}{}".format(str(v1), os.linesep)
        return lines

    def __str__(self):
        return self.get_lines_from_dict(self.config_dict['LAMMPSINNPT'])

    def write_input(self, filename=None, data_filename=None):
        """
        Write the main input file and the data file
        """
        # write the main input file
        with open(filename, 'w') as f:
            f.write(self.__str__())
        # write the data file if present
        if self.lammps_data:
            if data_filename:
                self.lammps_data.write_data_file(filename=data_filename)
            elif self.data_filename:
                self.lammps_data.write_data_file(filename=self.data_filename)
            else:
                dfilename = "lammps_data"
                print("No data filename provided. Writing the data to {"
                      "}".format(dfilename))
                self.lammps_data.write_data_file(filename=dfilename)

    @staticmethod
    def from_file(name, filename, **kwargs):
        """
        Read in the input settings from json file as ordereddict.
        Reads only the main input file from the json file, skips the data file.
        Note: with monty.serialization.loadfn the order of paramters in the
        json file is not preserved
        """
        with open(filename) as f:
            return DictLammpsInputSet(name,
                                      json.load(f,
                                                object_pairs_hook=OrderedDict),
                                      **kwargs)

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "name": self.name,
                "config_dict": self.config_dict,
                "lammps_data": self.lammps_data,
                "data_filename": self.data_filename,
                "user_lammps_settings": self.user_lammps_settings}

    @classmethod
    def from_dict(cls, d):
        return DictLammpsInputSet(d["name"], d["config_dict"],
                                  lammps_data=d.get("lammps_data"),
                                  data_filename=d.get("data_filename"),
                                  user_lammps_settings=d.get("user_lammps_settings"))


# NPT
NPTLammpsInputSet = partial(DictLammpsInputSet.from_file, "NPT",
                            os.path.join(MODULE_DIR, "data_files", "Lammps.json"))


# NPT followed by NVT
NPTNVTLammpsInputSet = partial(DictLammpsInputSet.from_file, "NPT_NVT",
                            os.path.join(MODULE_DIR, "data_files","Lammps_min_npt_nvt.json"))
