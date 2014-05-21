import copy
from rubicon.gff.atomtyper import AtomTyper

__author__ = 'navnidhirajput'

""" A wrapper for AntechamberRunner which generates force field files
    for a specified molecule using pdb file as input
"""

import subprocess
from pymatgen import write_mol
import shlex
from gff import Gff
from monty.io import ScratchDir
import tempfile


from rubicon.gff.topology import TopMol

ANTECHAMBER_DEBUG = False


class AntechamberRunner(AtomTyper):
    """
    A wrapper for AntechamberRunner software

    """

    def __init__(self, filename=None):
        self.filename = filename

    def _convert_to_pdb(self, molecule, filename=None):
        """
        generate pdb file for a given molecule
        """

        write_mol(molecule, filename)


    def _run_parmchk(self, filename=None):
        """
        run parmchk using ANTECHAMBER_AC.AC file

        Args:
            filename = pdb file of the molecule
        """
        command_parmchk = (
        'parmchk -i ' + filename + ' -f ac -o mol.frcmod -a Y')
        return_cmd = subprocess.call(shlex.split(command_parmchk))
        return return_cmd


    def get_ffmol(self,mols,filename=None):
        """
        generate and run antechamber command for specified pdb file

        Args:
            filename = pdb file of the molecule
        """


        scratch = tempfile.gettempdir()
        return_cmd = None

        with ScratchDir(scratch,
                        copy_to_current_on_exit=ANTECHAMBER_DEBUG) as d:
            gff_list = []
            top_list=[]

            for mol in mols:
                self._convert_to_pdb(mol, 'mol.pdb')
                command = (
                'antechamber -i ' + filename + ' -fi pdb -o ' + filename[:-4] +
                " -fo charmm")

                return_cmd = subprocess.call(shlex.split(command))
                self.molname = filename.split('.')[0]
                self._run_parmchk('ANTECHAMBER_AC.AC')
                top = TopMol.from_file('mol.rtf')
                my_gff = Gff.from_forcefield_para('mol.frcmod')
                my_gff.read_atom_index('ANTECHAMBER_AC.AC')
                gff_list.append(my_gff)
                top_list.append(top)
            self.gff_list=gff_list
            self.top_list=top_list


            return  gff_list, top_list

                #gff = self._parse_output()
                #self.read_atom_index('ANTECHAMBER_AC.AC')
#                self.read_atomType('ANTECHAMBER_AC.AC')

                #my_gff.read_forcefield_para('mol.frcmod')




                #print "+++++++",my_gff.atom_index_gaff

#                self.read_atomType('ANTECHAMBER_AC.AC')
                #my_gff = self._parse_output()
                #self._get_ff_bonds(my_gff.bonds, top.bonds, self.atom_gaff)
                #self._get_ff_angles(my_gff.angles, top.angles,
                #                    self.atom_gaff)
               # top._get_ff_dihedrals(my_gff.dihedrals, top.dihedrals,
                #                       gff.atom_gaff)
                #self._get_ff_imdihedrals(my_gff.imdihedrals, top.imdihedrals,
                #                         self.atom_gaff)
               # print "MYGFFBONDS",my_gff.bonds


                #gff_list.append(copy.deepcopy(my_gff))
               # print '*', type(my_gff)
               # print '**', gff_list

               # print '***', [x.bonds for x in gff_list]
                #print "GFFLIST0",gff_list[0].bonds,gff_list


            #print "GFFLIST0",gff_list[0].bonds
            #print "TOPLIST0",top_list[2].bonds






    def _parse_output(self):

        """
        create an object of forcefield_type_molecule
        and call the function
        """
        ffm = Gff()
        ffm.read_forcefield_para(self.molname + ".frcmod")
        return ffm











