

__author__ = 'navnidhirajput'


""" A wrapper for Antechamber which generates force field files
    for a specified molecule using pdb file as input
"""

import os
import subprocess
import copy
from pymatgen.core.structure import Molecule
from pymatgen import write_mol
import shlex
from GFF import GFF
import glob



class Antechamber():
    """
    A wrapper for Antechamber software

    """



    def __init__(self,molecule=None):
        """
        Args:
        molecule: The input molecule. Default is None
        """
        self.molecule=molecule


    def convert_to_pdb(self,molecule,filename):

        """
        generate pdb file for a given molecule
        """

        write_mol(molecule,filename)
        self.filename=filename





    def run_antechamber(self,filename=None):
        """
        generate and run antechamber command for specified pdb file

        Args:
            filename = pdb file of the molecule
        """

        command=('antechamber -i ' +filename +' -fi pdb -o ' +filename[:-4]+
             " -fo charmm")
        return_cmd=subprocess.call(shlex.split(command))
        self.molname = filename.split('.')[0]

        return return_cmd


    def parse_output(self):

        """
        create an object of forcefield_type_molecule
        and call the function
        """

        ffm=GFF()
        ffm.read_forcefield_para(self.molname+".prm")
        ffm.read_mass(self.molname+".rtf")
        return ffm


    def clean_files(self):
        for filename in glob.glob('*.AC'):
            os.remove(filename)
        for filename in glob.glob('*.AC0'):
            os.remove(filename)
        for filename in glob.glob('*.INF'):
            os.remove(filename)
        for filename in glob.glob('*.prm'):
            os.remove(filename)
        for filename in glob.glob('*.rtf'):
            os.remove(filename)
        for filename in glob.glob('*.inp'):
            os.remove(filename)
        for filename in glob.glob('*.pdb'):
            os.remove(filename)






