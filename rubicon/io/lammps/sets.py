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


__author__ = 'Navnidhi Rajput, Kiran Mathew'

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
            self.config_dict.update(user_lammps_settings)

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

    def write_input(self, filename=None):
        """
        Write the main input file and the data file
        """
        # write the main input file
        with open(filename, 'w') as f:
            f.write(self.__str__())
        # write the data file if present
        if self.lammps_data and self.data_filename:
            self.lammps_data.write_input(filename=self.data_filename)

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


NPTLammpsInputSet = partial(DictLammpsInputSet.from_file, "NPT",
                            os.path.join(MODULE_DIR, "data_files", "Lammps.json"))


NPTNVTLammpsInputSet = partial(DictLammpsInputSet.from_file, "NPT_NVT",
                            os.path.join(MODULE_DIR, "data_files","Lammps_min_npt_nvt.json"))


class DictLammpsInputSet_to_be_replaced(object):
    """
    Concrete implementation of LammpsInputSet that is initialized from a dict
    settings. This allows arbitrary settings to be input. In general,
    this is rarely used directly unless there is a source of settings in JSON
    format (e.g., from a REST interface). It is typically used by other
    JCESRLammpsInputSets for initialization.

    Args:
        name (str): A name fo the input set.
        config_dict (dict): The config dictionary to use.
        user_lammps_settings (dict): User LAMMPS settings. This allows a user
            to override LAMMPS settings, e.g., setting a different force field
             or bond type.
    """

    def __init__(self, name=None, parajson=None, config_dict=None,
                 user_lammps_settings=None):
        self.name = name
        self.lines = []
        self.parajson = parajson
        if user_lammps_settings:
            self.lammps_settings.update(user_lammps_settings)

    def get_lammps_control(self, filename=None, ensemble1="npt",
                           ensemble2="nvt", temp=None):
        jsonfile = open(
            os.path.join(os.path.dirname(__file__), "Lammps_min_npt_nvt.json"))
        self.parajson = json.load(jsonfile, encoding="utf-8")
        self.parajson['LAMMPSINNPT']['temp'] = temp
        self.parajson['LAMMPSINNPT']['fix1']['style'] = ensemble1
        self.parajson['LAMMPSINNPT']['fix1']['ID'] = ensemble1
        self.parajson['LAMMPSINNPT']['fix2']['style'] = ensemble2
        self.parajson['LAMMPSINNPT']['fix2']['ID'] = ensemble2
        self.parajson['LAMMPSINNPT']['fix1']['Tstart'] = temp
        self.parajson['LAMMPSINNPT']['fix1']['Tstop'] = temp
        return self.parajson

    def __str__(self):
        lines = []
        lines.append('log ' + self.parajson['LAMMPSINNPT']['log'])
        lines.append(
            ' # read_restart ' + self.parajson['LAMMPSINNPT']['read_restart'])

        lines.append('units ' + self.parajson['LAMMPSINNPT']['units'])
        lines.append(
            'atom_style ' + self.parajson['LAMMPSINNPT']['atom_style'])
        lines.append('{} {}'.format('boundary ', " ".join(
            (self.parajson['LAMMPSINNPT']['boundary'][:]))))
        lines.append('{} {} {}'.format('pair_style ', "/".join([str(x) for x in
                                                                self.parajson[
                                                                    'LAMMPSINNPT'][
                                                                    'pair_style'][
                                                                    "style"]]),
                                       (self.parajson['LAMMPSINNPT'][
                                            'pair_style']["args"])))
        lines.append('{} {} {}'.format('kspace_style ', 'pppm', (
            self.parajson['LAMMPSINNPT']['kspace_style']["pppm"])))
        lines.append('{} {} {} {} {}'.format('pair_modify ', "tail", (
            self.parajson['LAMMPSINNPT']['pair_modify']["tail"]),
                                             "mix", (
                                                 self.parajson['LAMMPSINNPT'][
                                                     'pair_modify']["mix"])))
        lines.append('{} {}'.format('special_bonds ', (
            self.parajson['LAMMPSINNPT']['special_bonds'])))
        lines.append('{} {}'.format('bond_style ', (
            self.parajson['LAMMPSINNPT']['bond_style'])))
        lines.append('{} {}'.format('angle_style ', (
            self.parajson['LAMMPSINNPT']['angle_style'])))
        lines.append('{} {}'.format('dihedral_style ', (
            self.parajson['LAMMPSINNPT']['dihedral_style'])))
        lines.append('{} {}'.format('improper_style ', (
            self.parajson['LAMMPSINNPT']['improper_style'])))
        lines.append('{} {}'.format('read_data ', (
            self.parajson['LAMMPSINNPT']['read_data'])))
        lines.append('{} {} {}'.format('neighbor ', (
            self.parajson['LAMMPSINNPT']['neighbor']["skin"]),
                                       (
                                           self.parajson['LAMMPSINNPT'][
                                               'neighbor'][
                                               "style"])))
        lines.append('{} {} {} {} {} {} {} {} {} {} {}'.format('neigh_modify ',
                                                               " delay", (
                                                                   self.parajson[
                                                                       'LAMMPSINNPT'][
                                                                       'neigh_modify'][
                                                                       "delay"]),
                                                               "every", (
                                                                   self.parajson[
                                                                       'LAMMPSINNPT'][
                                                                       'neigh_modify'][
                                                                       "every"]),
                                                               "check", (
                                                                   self.parajson[
                                                                       'LAMMPSINNPT'][
                                                                       'neigh_modify'][
                                                                       "check"]),
                                                               "page", (
                                                                   self.parajson[
                                                                       'LAMMPSINNPT'][
                                                                       'neigh_modify'][
                                                                       "page"]),
                                                               "one", (
                                                                   self.parajson[
                                                                       'LAMMPSINNPT'][
                                                                       'neigh_modify'][
                                                                       "one"])))
        lines.append('{} {}'.format('timestep ', (
            self.parajson['LAMMPSINNPT']['timestep']["dt"])))
        lines.append('{} {} {} {} {}'.format('minimize ', (
            self.parajson['LAMMPSINNPT']['minimize']["etol"]),
                                             (self.parajson['LAMMPSINNPT'][
                                                  'minimize']["ftol"]),
                                             (self.parajson['LAMMPSINNPT'][
                                                  'minimize']["maxiter"]),
                                             (self.parajson['LAMMPSINNPT'][
                                                  'minimize']["maxeval"])))
        lines.append('{} {} {} {} {} {} {}'.format('velocity ',
                                                   (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'velocity1'][
                                                           "group-id"]),
                                                   (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'velocity1'][
                                                           "style"]),
                                                   (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'temp']),
                                                   (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'velocity1'][
                                                           "seed"]),
                                                   "units", (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'velocity1'][
                                                           "units value"])))
        lines.append('{} {} {} {} {} {}'.format('velocity ',
                                                (self.parajson['LAMMPSINNPT'][
                                                     'velocity2']["group-id"]),
                                                (self.parajson['LAMMPSINNPT'][
                                                     'velocity2']["style"]),
                                                (self.parajson['LAMMPSINNPT'][
                                                     'velocity2']["args"]),
                                                "units", (
                                                    self.parajson[
                                                        'LAMMPSINNPT'][
                                                        'velocity2'][
                                                        "units value"])))
        lines.append('{} {} {} {} {} {} {}'.format('dump ',
                                                   (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'dump']["ID"]),
                                                   (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'dump'][
                                                           "group-id"]),
                                                   (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'dump']["style"]),
                                                   (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'dump']["N"]),
                                                   (
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'dump']["file"]),
                                                   " ".join((
                                                       self.parajson[
                                                           'LAMMPSINNPT'][
                                                           'dump'][
                                                           "possible attributes"][
                                                       :]))))
        lines.append('{} {} {}'.format('thermo_style ', (
            self.parajson['LAMMPSINNPT']['thermo_style']["style"]),
                                       " ".join((self.parajson['LAMMPSINNPT'][
                                                     'thermo_style'][
                                                     "possible attributes"][
                                                 :]))))
        lines.append('{} {}'.format('thermo ', (
            self.parajson['LAMMPSINNPT']['thermo']["N"])))

        lines.append('{} {} {} {} {} {} {} {} {} {} {} {}'.format('fix ',
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "ID"]),
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "group-id"]),
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "style"]),
                                                                  "temp",
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "Tstart"]),
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "Tstop"]),
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "Tdamp"]),
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "iso"]),
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "Pstart"]),
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "Pstop"]),
                                                                  (
                                                                      self.parajson[
                                                                          'LAMMPSINNPT'][
                                                                          'fix1'][
                                                                          "Pdamp"])))
        lines.append(
            '{} {}'.format('run ',
                           (self.parajson['LAMMPSINNPT']['run1']["N"])))

        lines.append(
            '{} {}'.format('unfix ',
                           (self.parajson['LAMMPSINNPT']['unfix1']["N"])))

        lines.append('{} {} {} {} {} {} {} {}'.format('fix ',
                                                      (
                                                          self.parajson[
                                                              'LAMMPSINNPT'][
                                                              'fix2'][
                                                              "ID"]),
                                                      (
                                                          self.parajson[
                                                              'LAMMPSINNPT'][
                                                              'fix2'][
                                                              "group-id"]),
                                                      (
                                                          self.parajson[
                                                              'LAMMPSINNPT'][
                                                              'fix2'][
                                                              "style"]),
                                                      "temp",
                                                      (
                                                          self.parajson[
                                                              'LAMMPSINNPT'][
                                                              'fix2'][
                                                              "Tstart"]),
                                                      (
                                                          self.parajson[
                                                              'LAMMPSINNPT'][
                                                              'fix2'][
                                                              "Tstop"]),
                                                      (
                                                          self.parajson[
                                                              'LAMMPSINNPT'][
                                                              'fix2'][
                                                              "Tdamp"])))
        lines.append(
            '{} {}'.format('run ',
                           (self.parajson['LAMMPSINNPT']['run2']["N"])))

        lines.append(
            '{} {}'.format('unfix ',
                           (self.parajson['LAMMPSINNPT']['unfix2']["N"])))

        # elif self.parajson['LAMMPSINNPT']['fix1']["style"] == "nvt":
        #     lines.append('{} {} {} {} {} {} {} {}'.format('fix1 ',
        #                                                   (self.parajson[
        #                                                        'LAMMPSINNPT'][
        #                                                        'fix1']["ID"]),
        #                                                   (self.parajson[
        #                                                        'LAMMPSINNPT'][
        #                                                        'fix1'][
        #                                                        "group-id"]),
        #                                                   (self.parajson[
        #                                                        'LAMMPSINNPT'][
        #                                                        'fix1']["style"]),
        #                                                   "temp",
        #                                                   (self.parajson[
        #                                                        'LAMMPSINNPT'][
        #                                                        'fix1'][
        #                                                        "Tstart"]),
        #                                                   (self.parajson[
        #                                                        'LAMMPSINNPT'][
        #                                                        'fix1']["Tstop"]),
        #                                                   (self.parajson[
        #                                                        'LAMMPSINNPT'][
        #                                                        'fix1'][
        #                                                        "Tdamp"])))
        # lines.append('{} {} {} {}'.format('restart ',
        #                                   self.parajson['LAMMPSINNPT'][
        #                                       "restart"]['N'],
        #                                   self.parajson['LAMMPSINNPT'][
        #                                       "restart"]['file1'],
        #                                   self.parajson['LAMMPSINNPT'][
        #                                       "restart"]['file2']))
        # lines.append(
        #     '{} {}'.format('run ', (self.parajson['LAMMPSINNPT']['run']["N"])))
        lines.append('{} {}'.format('write_restart ',
                                    self.parajson['LAMMPSINNPT'][
                                        'write_restart']["file"]))

        return '\n'.join(lines)

    def write_lampps_control(self, filename=None):
        with open(filename, 'w') as f:
            f.write(self.__str__())

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "parajson": self.parajson}

    @classmethod
    def from_dict(cls, d):
        return DictLammpsInputSet_to_be_replaced(parajson=d["parajson"])
