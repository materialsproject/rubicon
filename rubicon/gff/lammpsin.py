__author__ = 'navnidhirajput'

class DictLammpsInputSet():

    def __init__(self,name,config_dict,user_lammps_settings=None):
        self.name = name
        self.lammps_settings = config_dict

        if user_lammps_settings:
            self.lammps_settings.update(user_lammps_settings)


