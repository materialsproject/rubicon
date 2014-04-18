import abc
from functools import partial
import json
import os
from pymatgen.serializers.json_coders import MSONable

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

__author__ = 'navnidhirajput'


class DictLammpsInputSet():

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
    def __init__(self,name=None, parajson=None,config_dict=None,user_lammps_settings=None):

        self.name=name
        self.lines =[]
        self.parajson=parajson
        if user_lammps_settings:
            self.lammps_settings.update(user_lammps_settings)

    def get_lammpsin(self,filename=None,ensemble=None,temp=None):


        jsonfile=open(filename)
        self.parajson=json.load(jsonfile, encoding="utf-8")
        self.parajson['LAMMPSINNPT']['temp'] = temp
        self.parajson['LAMMPSINNPT']['fix']['style'] = ensemble
        self.parajson['LAMMPSINNPT']['fix']['ID'] = ensemble
        self.parajson['LAMMPSINNPT']['fix']['Tstart'] = temp
        self.parajson['LAMMPSINNPT']['fix']['Tstop'] = temp
        return self.parajson

    def __str__(self):
        lines=[]
        lines.append('log '+ self.parajson['LAMMPSINNPT']['log'])
        lines.append('read_start '+ self.parajson['LAMMPSINNPT']['read_start'])
        lines.append('units '+ self.parajson['LAMMPSINNPT']['units'])
        lines.append('atom_style '+ self.parajson['LAMMPSINNPT']['atom_style'])
        lines.append('{} {}'.format('boundary ', " ".join((self.parajson['LAMMPSINNPT']['boundary'][:]))))
        lines.append('{} {} {}'.format('pair_style ', "/".join([str(x) for x in self.parajson['LAMMPSINNPT']['pair_style']["style"]]),
                                       (self.parajson['LAMMPSINNPT']['pair_style']["args"])))
        lines.append('{} {}'.format('kspace_style ', (self.parajson['LAMMPSINNPT']['kspace_style']["pppm"])))
        lines.append('{} {} {} {} {}'.format('pair_modify ', "tail",(self.parajson['LAMMPSINNPT']['pair_modify']["tail"]),
                                    "mix",(self.parajson['LAMMPSINNPT']['pair_modify']["mix"])))
        lines.append('{} {}'.format('special_bonds ', (self.parajson['LAMMPSINNPT']['special_bonds'])))
        lines.append('{} {}'.format('bond_style ', (self.parajson['LAMMPSINNPT']['bond_style'])))
        lines.append('{} {}'.format('dihedral_style ', (self.parajson['LAMMPSINNPT']['dihedral_style'])))
        lines.append('{} {}'.format('read_data ', (self.parajson['LAMMPSINNPT']['read_data'])))
        lines.append('{} {} {}'.format('neighbor ', (self.parajson['LAMMPSINNPT']['neighbor']["skin"]),
                                    (self.parajson['LAMMPSINNPT']['neighbor']["style"])))
        lines.append('{} {} {} {} {} {} {} {} {} {} {}'.format('neigh_modify ',
                   " delay", (self.parajson['LAMMPSINNPT']['neigh_modify']["delay"]),
                   "every",(self.parajson['LAMMPSINNPT']['neigh_modify']["every"]),
                   "check",(self.parajson['LAMMPSINNPT']['neigh_modify']["check"]),
                   "page",(self.parajson['LAMMPSINNPT']['neigh_modify']["page"]),
                   "one",(self.parajson['LAMMPSINNPT']['neigh_modify']["one"])))
        lines.append('{} {}'.format('timestep ', (self.parajson['LAMMPSINNPT']['timestep']["dt"])))
        lines.append('{} {} {} {} {}'.format('minimize ',(self.parajson['LAMMPSINNPT']['minimize']["etol"]),
                                             (self.parajson['LAMMPSINNPT']['minimize']["ftol"]),
                                             (self.parajson['LAMMPSINNPT']['minimize']["maxiter"]),
                                             (self.parajson['LAMMPSINNPT']['minimize']["maxeval"])))
        lines.append('{} {} {} {} {} {} {}'.format('velocity ',
                   (self.parajson['LAMMPSINNPT']['velocity1']["group-id"]),
                   (self.parajson['LAMMPSINNPT']['velocity1']["style"]),
                   (self.parajson['LAMMPSINNPT']['velocity1']["temp"]),
                   (self.parajson['LAMMPSINNPT']['velocity1']["seed"]),
                   "units",(self.parajson['LAMMPSINNPT']['velocity1']["units value"])))
        lines.append('{} {} {} {} {} {}'.format('velocity ',
                   (self.parajson['LAMMPSINNPT']['velocity2']["group-id"]),
                   (self.parajson['LAMMPSINNPT']['velocity2']["style"]),
                   (self.parajson['LAMMPSINNPT']['velocity2']["args"]),
                   "units",(self.parajson['LAMMPSINNPT']['velocity2']["units value"])))
        lines.append('{} {} {} {} {} {} {}'.format('dump ',
                   (self.parajson['LAMMPSINNPT']['dump']["ID"]),
                   (self.parajson['LAMMPSINNPT']['dump']["group-id"]),
                   (self.parajson['LAMMPSINNPT']['dump']["style"]),
                   (self.parajson['LAMMPSINNPT']['dump']["N"]),
                   (self.parajson['LAMMPSINNPT']['dump']["file"]),
                   " ".join((self.parajson['LAMMPSINNPT']['dump']["possible attributes"][:]))))
        lines.append('{} {} {}'.format('thermo_style ',(self.parajson['LAMMPSINNPT']['thermo_style']["style"]),
                                    " ".join((self.parajson['LAMMPSINNPT']['thermo_style']["possible attributes"][:]))))
        lines.append('{} {}'.format('thermo ', (self.parajson['LAMMPSINNPT']['thermo']["N"])))

        if self.parajson['LAMMPSINNPT']['fix']["style"] == "npt":

            lines.append('{} {} {} {} {} {} {} {} {} {} {} {}'.format('fix ',
                       (self.parajson['LAMMPSINNPT']['fix']["ID"]),
                       (self.parajson['LAMMPSINNPT']['fix']["group-id"]),
                       (self.parajson['LAMMPSINNPT']['fix']["style"]),"temp",
                       (self.parajson['LAMMPSINNPT']['fix']["Tstart"]),
                       (self.parajson['LAMMPSINNPT']['fix']["Tstop"]),
                       (self.parajson['LAMMPSINNPT']['fix']["Tdamp"]),
                       (self.parajson['LAMMPSINNPT']['fix']["iso"]),
                       (self.parajson['LAMMPSINNPT']['fix']["Pstart"]),
                       (self.parajson['LAMMPSINNPT']['fix']["Pstop"]),
                       (self.parajson['LAMMPSINNPT']['fix']["Pdamp"])))
        elif self.parajson['LAMMPSINNPT']['fix']["style"] == "nvt":
                        lines.append('{} {} {} {} {} {} {} {}'.format('fix ',
                       (self.parajson['LAMMPSINNPT']['fix']["ID"]),
                       (self.parajson['LAMMPSINNPT']['fix']["group-id"]),
                       (self.parajson['LAMMPSINNPT']['fix']["style"]),"temp",
                       (self.parajson['LAMMPSINNPT']['fix']["Tstart"]),
                       (self.parajson['LAMMPSINNPT']['fix']["Tstop"]),
                       (self.parajson['LAMMPSINNPT']['fix']["Tdamp"])))
        lines.append('{} {} {} {}'.format('restart ',self.parajson['LAMMPSINNPT']["restart"]['N'],
                                          self.parajson['LAMMPSINNPT']["restart"]['file1'],
                                          self.parajson['LAMMPSINNPT']["restart"]['file2']))
        lines.append('{} {}'.format('run ', (self.parajson['LAMMPSINNPT']['run']["N"])))
        lines.append('{} {}'.format('write_restart ', self.parajson['LAMMPSINNPT']['write_restart']["file"]))


        return '\n'.join(lines)

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "parajson": self.parajson}


    @classmethod
    def from_dict(cls, d):
        return DictLammpsInputSet(parajson=d["parajson"])



