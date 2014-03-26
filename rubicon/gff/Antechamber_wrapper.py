

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
from monty.io import ScratchDir
import tempfile



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
        #scratch = tempfile.gettempdir()
        #with ScratchDir(scratch,copy_from_current_on_enter=True) as d:

            #os.chdir(d)
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

        command_parmchk=('parmchk -i '+ filename +'-f ac -o mol.frcmod -a Y')

        return_cmd=subprocess.call(shlex.split(command_parmchk))
        return return_cmd


    def parse_output(self):

        """
        create an object of forcefield_type_molecule
        and call the function
        """

        ffm=GFF()
        ffm.read_forcefield_para(self.molname+".frcmod")
        return ffm


    def clean_files(self):

        """
        clean Antechamber files
        """
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


    def get_FF_bonds(self,gff_bonds,top_bond,atom_gaff):

        self.gaff_info=[]
        for keys, values in gff_bonds.iteritems():
                self.gaff_info=[keys,values]
        for item in top_bond:
            d1=item[0]+' '+ item[1]
            a1,a2 = atom_gaff[item[0]],atom_gaff[item[1]]
            if (str(a1),str(a2)) in gff_bonds:
                self.topbondFF[d1]=((str(a1),str(a2)),gff_bonds[(str(a1),str(a2))])
            elif (str(a2),str(a1)) in gff_bonds:
                bond_type=tuple(sorted((str(a1),str(a2))))
                self.topbondFF[d1]=((str(a1),str(a2)),gff_bonds[bond_type])
                self.num_bond_types = len(set(self.topbondFF.keys()))


    def get_FF_angles(self,gff_angles,top_angle,atom_gaff):

        self.gaff_info=[]
        for keys, values in gff_angles.iteritems():
                self.gaff_info=[keys,values]
        for item in top_angle:
            d1=item[0]+' '+ item[1]+' '+item[2]
            a1,a2,a3 = atom_gaff[item[0]],atom_gaff[item[1]],atom_gaff[item[2]]
            if (str(a1),str(a2),str(a3)) in gff_angles:
                self.topangleFF[d1]=(((str(a1),str(a2),str(a3))),gff_angles[(str(a1),str(a2),str(a3))])
            elif (str(a3),str(a2),str(a1)) in gff_angles:
                angle_type = tuple(sorted((str(a1),str(a2),str(a3))))
                self.topangleFF[d1]=(((str(a1),str(a2),str(a3))),gff_angles[angle_type])
            self.num_ang_types = len(set(self.topangleFF.keys()))



    def get_FF_dihedrals(self,gff_dihedrals,top_dihedral,atom_gaff):

        self.gaff_info=[]
        for keys, values in gff_dihedrals.iteritems():
                self.gaff_info=[keys,values]

        for item in top_dihedral:
            d1=item[0]+' '+ item[1]+' '+item[2]+ ' '+item[3]
            a1,a2,a3,a4 = atom_gaff[item[0]],atom_gaff[item[1]],atom_gaff[item[2]],atom_gaff[item[3]]
            dihedral_label=(a1,a2,a3,a4)
            if dihedral_label[0]>dihedral_label[3]:
                        dihedral_label=tuple(reversed(list(dihedral_label)))
            if dihedral_label in gff_dihedrals:
                self.topdihedralFF[d1]=(dihedral_label,gff_dihedrals[dihedral_label])
            self.num_dih_types = len(set(self.topdihedralFF.keys()))
        print self.topdihedralFF


    def get_FF_imdihedrals(self,gff_imdihedrals,top_imdihedral,atom_gaff):

        self.gaff_info=[]
        for keys, values in gff_imdihedrals.iteritems():
                self.gaff_info=[keys,values]
        for item in top_imdihedral:
            d1=item[0]+' '+ item[1]+' '+item[2]+ ' '+item[3]
            a1,a2,a3,a4 = atom_gaff[item[0]],atom_gaff[item[1]],atom_gaff[item[2]],atom_gaff[item[3]]

            if (str(a1),str(a2),str(a3),str(a4)) in gff_imdihedrals:
                self.topimdihedralFF[d1]=((str(a1),str(a2),str(a3),str(a4)),gff_imdihedrals[(str(a1),str(a2),str(a3),str(a4))])

            elif (str(a4),str(a3),str(a2),str(a1)) in gff_imdihedrals:
                self.topimdihedralFF[d1]=((str(a1),str(a2),str(a3),str(a4)),gff_imdihedrals[(str(a1),str(a2),str(a3),str(a4))])
            self.num_imdih_types = len(set(self.topimdihedralFF.keys()))









