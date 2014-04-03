__author__ = 'navnidhirajput'

class DictLammpsInputSet():

    def __init__(self,user_lammps_settings=None):
        self.lines =[]

        if user_lammps_settings:
            self.lammps_settings.update(user_lammps_settings)

    def get_lammps_config(self):
        lines=[]
        lines.append("#read lammps data file, write a restart file")

        return '\n'.join(lines)


class JCESRLammpsInputSet(DictLammpsInputSet):

    def __init__(self):
        pass

