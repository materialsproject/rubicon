

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
from gff import GFF
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
        self.topbondFF=dict()
        self.topangleFF=dict()
        self.topdihedralFF=dict()
        self.topimdihedralFF=dict()


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

    def run_parmchk(self,filename=None):
        """
        run parmchk using ANTECHAMBER_AC.AC file

        Args:
            filename = pdb file of the molecule
        """

        command_parmchk=('parmchk -i '+filename +' -f ac -o mol.frcmod')

        return_cmd=subprocess.call(shlex.split(command_parmchk))
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
        for filename in glob.glob('*.frcmod'):
            os.remove(filename)


    def get_FF_bonds(self,bonds,top_bond,atom_gaff):

        self.gaff_info=[]
        for keys, values in bonds.iteritems():
                self.gaff_info=[keys,values]
        for item in top_bond:
            d1=item[0]+'-'+ item[1]
            a1,a2 = atom_gaff[item[0]],atom_gaff[item[1]]
            if a1+'-'+a2 in bonds:
                self.topbondFF[d1]=(a1+'-'+a2,bonds[a1+'-'+a2])
            else:
                self.topbondFF[d1]=(a2+'-'+a1,bonds[a2+'-'+a1])
                self.num_bond_types = len(set(self.topbondFF.keys()))


    def get_FF_angles(self,angles,top_angle,atom_gaff):

        self.gaff_info=[]
        for keys, values in angles.iteritems():
                self.gaff_info=[keys,values]
        for item in top_angle:
            d1=item[0]+'-'+ item[1]+'-'+item[2]
            a1,a2,a3 = atom_gaff[item[0]],atom_gaff[item[1]],atom_gaff[item[2]]
            if a1+'-'+a2+'-'+a3 in angles:
                self.topangleFF[d1]=(a1+'-'+a2+'-'+a3,angles[a1+'-'+a2+'-'+a3])

            elif reversed(a1 +'-'+a2+'-'+a3):
                self.topangleFF[d1]=(a1+'-'+a2+'-'+a3,angles[a1+'-'+a2+'-'+a3])
            self.num_ang_types = len(set(self.topangleFF.keys()))


    def get_FF_dihedrals(self,dihedrals,top_dihedral,atom_gaff):

        self.gaff_info=[]
        for keys, values in dihedrals.iteritems():
                self.gaff_info=[keys,values]
        for item in top_dihedral:
            d1=item[0]+'-'+ item[1]+'-'+item[2]+ '-'+item[3]
            a1,a2,a3,a4 = atom_gaff[item[0]],atom_gaff[item[1]],atom_gaff[item[2]],atom_gaff[item[3]]
            if 'X-'+a2+'-'+a3+'-'+'X' in dihedrals:
                self.topdihedralFF[d1]=('X-'+a2+'-'+a3+'-'+'X',dihedrals['X-'+a2+'-'+a3+'-'+'X'])

            elif reversed(a1+'-'+a2+'-'+a3+'-'+ a4):
                self.topdihedralFF[d1]=('X-'+a2+'-'+a3+'-'+'X',dihedrals['X-'+a2+'-'+a3+'-'+'X'])
            self.num_dih_types = len(set(self.topdihedralFF.keys()))


    def get_FF_imdihedrals(self,imdihedrals,top_imdihedral,atom_gaff):

        self.gaff_info=[]
        for keys, values in imdihedrals.iteritems():
                self.gaff_info=[keys,values]
        for item in top_imdihedral:
            d1=item[0]+'-'+ item[1]+'-'+item[2]+ '-'+item[3]
            a1,a2,a3,a4 = atom_gaff[item[0]],atom_gaff[item[1]],atom_gaff[item[2]],atom_gaff[item[3]]
            if a1+'-'+a2+'-'+a3+'-'+ a4 in imdihedrals:
                self.topimdihedralFF[d1]=(a1+'-'+a2+'-'+a3+'-'+ a4,imdihedrals[a1+'-'+a2+'-'+a3+'-'+ a4])

            elif reversed(a1+'-'+a2+'-'+a3+'-'+ a4):
                self.topimdihedralFF[d1]=(a1+'-'+a2+'-'+a3+'-'+ a4,imdihedrals[a1+'-'+a2+'-'+a3+'-'+ a4])
            self.num_imdih_types = len(set(self.topimdihedralFF.keys()))







