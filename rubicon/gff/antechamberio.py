

__author__ = 'navnidhirajput'

""" A wrapper for AntechamberRunner which generates force field files
    for a specified molecule using pdb file as input
"""

import subprocess
from pymatgen import write_mol, Molecule
import shlex
from gff import Gff
from monty.io import ScratchDir
import tempfile
from pymatgen.packmol.packmol import PackmolRunner

from rubicon.gff.topology import TopMol

ANTECHAMBER_DEBUG = False


class AntechamberRunner():
    """
    A wrapper for AntechamberRunner software

    """

    def __init__(self, molecule=None):
        """
        Args:
        molecule: The input molecule. Default is None
        """
        self.molecule = molecule
        self.topbondff = dict()
        self.topangleff = dict()
        self.topdihedralff = dict()
        self.topimdihedralff = dict()
        self.atom_index=dict()
        self.atom_index_gaff=dict()
        self.atom_gaff=dict()


    def read_atom_index(self,filename=None):

        with open(filename) as f:

            for line in f.readlines():
                token = line.split()
                if token[0]=='ATOM':
                    index=int(token[1])
                    atom_name=token[2]
                    atom_gaff=token[9]
                    self.atom_index[index]=atom_name

                    self.atom_index_gaff[index]=atom_gaff


    def read_atomType(self,filename=None):

        with open(filename) as f:

            for line in f.readlines():
                token = line.split()
                if token[0]=='ATOM':
                    atom_name=token[2]
                    gaff_name=token[-1]
                    self.atom_gaff[atom_name]=gaff_name
            self.atom_gaff.update(self.atom_gaff)
        self.num_types = len(set(self.atom_gaff.values()))


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


    def _run_antechamber(self, filename=None, mols=[]):
        """
        generate and run antechamber command for specified pdb file

        Args:
            filename = pdb file of the molecule
        """
        scratch = tempfile.gettempdir()
        return_cmd = None
        my_gff = None


        with ScratchDir(scratch,
                        copy_to_current_on_exit=ANTECHAMBER_DEBUG) as d:
            gff_list = []

            for mol in mols:

                self._convert_to_pdb(mol, 'mol.pdb')
                command = (
                'antechamber -i ' + filename + ' -fi pdb -o ' + filename[:-4] +
                " -fo charmm")

                return_cmd = subprocess.call(shlex.split(command))
                self.molname = filename.split('.')[0]
                self._run_parmchk('ANTECHAMBER_AC.AC')
                gff = self._parse_output()
                self.read_atom_index('ANTECHAMBER_AC.AC')
                self.read_atomType('ANTECHAMBER_AC.AC')
                top = TopMol.from_file('mol.rtf')
                my_gff = Gff()
                my_gff.read_forcefield_para('mol.frcmod')
                self.read_atomType('ANTECHAMBER_AC.AC')
                self._get_ff_bonds(my_gff.bonds, top.bonds, self.atom_gaff)
                self._get_ff_angles(my_gff.angles, top.angles,
                                    self.atom_gaff)
                self._get_ff_dihedrals(my_gff.dihedrals, top.dihedrals,
                                       self.atom_gaff)
                self._get_ff_imdihedrals(my_gff.imdihedrals, top.imdihedrals,
                                         self.atom_gaff)
                gff_list.append(gff)


            return my_gff, gff_list,top




    def _parse_output(self):

        """
        create an object of forcefield_type_molecule
        and call the function
        """
        ffm = Gff()
        ffm.read_forcefield_para(self.molname + ".frcmod")
        return ffm

    def _get_ff_bonds(self, gff_bonds, top_bond, atom_gaff):

        self.gaff_info = []
        for keys, values in gff_bonds.iteritems():
            self.gaff_info = [keys, values]
        for item in top_bond:
            d1 = item[0] + ' ' + item[1]
            a1, a2 = atom_gaff[item[0]], atom_gaff[item[1]]
            if (str(a1), str(a2)) in gff_bonds:
                self.topbondff[d1] = (
                (str(a1), str(a2)), gff_bonds[(str(a1), str(a2))])
            elif (str(a2), str(a1)) in gff_bonds:
                bond_type = tuple(sorted((str(a1), str(a2))))
                self.topbondff[d1] = ((str(a1), str(a2)), gff_bonds[bond_type])
                self.num_bond_types = len(set(self.topbondff.keys()))


    def _get_ff_angles(self, gff_angles, top_angle, atom_gaff):

        self.gaff_info = []
        for keys, values in gff_angles.iteritems():
            self.gaff_info = [keys, values]
        for item in top_angle:
            d1 = item[0] + ' ' + item[1] + ' ' + item[2]
            a1, a2, a3 = atom_gaff[item[0]], atom_gaff[item[1]], atom_gaff[
                item[2]]
            if (str(a1), str(a2), str(a3)) in gff_angles:
                self.topangleff[d1] = (((str(a1), str(a2), str(a3))),
                                       gff_angles[(str(a1), str(a2), str(a3))])
            elif (str(a3), str(a2), str(a1)) in gff_angles:
                angle_type = tuple(sorted((str(a1), str(a2), str(a3))))
                self.topangleff[d1] = (
                ((str(a1), str(a2), str(a3))), gff_angles[angle_type])
            self.num_ang_types = len(set(self.topangleff.keys()))


    def _get_ff_dihedrals(self, gff_dihedrals, top_dihedral, atom_gaff):

        self.gaff_info = []
        for keys, values in gff_dihedrals.iteritems():
            self.gaff_info = [keys, values]

        for item in top_dihedral:
            d1 = item[0] + ' ' + item[1] + ' ' + item[2] + ' ' + item[3]
            a1, a2, a3, a4 = atom_gaff[item[0]], atom_gaff[item[1]], atom_gaff[
                item[2]], atom_gaff[item[3]]
            dihedral_label = (a1, a2, a3, a4)
            if dihedral_label[0] > dihedral_label[3]:
                dihedral_label = tuple(reversed(list(dihedral_label)))
            if dihedral_label in gff_dihedrals:
                self.topdihedralff[d1] = (
                dihedral_label, gff_dihedrals[dihedral_label])
            self.num_dih_types = len(set(self.topdihedralff.keys()))


    def _get_ff_imdihedrals(self, gff_imdihedrals, top_imdihedral, atom_gaff):

        self.gaff_info = []
        for keys, values in gff_imdihedrals.iteritems():
            self.gaff_info = [keys, values]

        for item in top_imdihedral:
            d1 = item[0] + ' ' + item[1] + ' ' + item[2] + ' ' + item[3]
            a1, a2, a3, a4 = atom_gaff[item[0]], atom_gaff[item[1]], atom_gaff[
                item[2]], atom_gaff[item[3]]
            imdihedral_label = (a1, a2, a3, a4)
            if imdihedral_label[0] > imdihedral_label[3]:
                imdihedral_label = tuple(reversed(list(imdihedral_label)))
            if imdihedral_label in gff_imdihedrals:
                self.topimdihedralff[d1] = (
                imdihedral_label, gff_imdihedrals[imdihedral_label])
            self.num_imdih_types = len(set(self.topimdihedralff.keys()))

    @staticmethod
    def run(mols):

        gff_list = []
        for mol in mols:
            my_ant = AntechamberRunner(mol)
            my_ant._convert_to_pdb(mol, 'mol.pdb')
            my_ant._run_antechamber('mol.pdb')
            my_ant._run_parmchk('ANTECHAMBER_AC.AC')

            gff = my_ant._parse_output()
            my_ant.read_atom_index('ANTECHAMBER_AC.AC')
            my_ant.read_atomType('ANTECHAMBER_AC.AC')

            pmr = PackmolRunner([mol, mol], [
                {"number": 1, "inside box": [0., 0., 0., 40., 40., 40.]},
                {"number": 2}])
            top = TopMol.from_file('mol.rtf')

            my_gff = Gff()
            my_gff.read_forcefield_para('mol.frcmod')
            my_ant.read_atomType('ANTECHAMBER_AC.AC')
            my_ant._get_ff_bonds(my_gff.bonds, top.bonds, my_ant.atom_gaff)
            my_ant._get_ff_angles(my_gff.angles, top.angles,
                                  my_ant.atom_gaff)
            my_ant._get_ff_dihedrals(my_gff.dihedrals, top.dihedrals,
                                     my_ant.atom_gaff)
            my_ant._get_ff_imdihedrals(my_gff.imdihedrals, top.imdihedrals,
                                       my_ant.atom_gaff)
            gff_list.append(gff)
            #my_lampps = LmpInput()
            #my_lampps.set_coeff(my_gff, top, pmr, my_ant)
            #my_lampps.set_atom(pmr, my_gff,my_ant)
            #my_lampps.set_bonds(pmr, my_gff, top,my_ant)
            #my_lampps.set_angles(pmr, my_gff, top,my_ant)
            #my_lampps.set_dihedrals(pmr, my_gff, top, my_ant)
            #my_lampps.set_imdihedrals(pmr, my_gff, top,my_ant)
            return my_gff, my_ant








