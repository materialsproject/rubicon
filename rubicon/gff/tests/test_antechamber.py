from unittest import TestCase
from pymatgen import Molecule
from rubicon.gff.antechamber_wrapper import AntechamberRunner
from rubicon.gff.gff import GFF_library, Gff
from rubicon.gff.topology import TopMol, AC

__author__ = 'navnidhirajput'

coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]

mol = Molecule(["C", "H", "H", "H", "H"], coords)


class TestAntechamber(TestCase):



    def test_get_FF_bonds(self):
        my_ant = AntechamberRunner(mol)
        my_gff = Gff()
        top = TopMol.from_file('mol.rtf')
        atom_gaff = AC()
        my_ant._get_ff_bonds(my_gff.bonds, top.bonds, atom_gaff.atom_gaff)
        ans_bond={'C-H1': ('c3-hc', (337.3, 1.092)), 'C-H2': ('c3-hc', (337.3, 1.092)), 'C-H3': ('c3-hc', (337.3, 1.092)), 'C-H': ('c3-hc', (337.3, 1.092))}
        self.assertEquals(ans_bond,my_ant.topbondff)

    def test_get_FF_angles(self):
        my_ant = AntechamberRunner(mol)
        my_gff = Gff()
        top = TopMol.from_file('mol.rtf')
        atom_gaff = AC()
        my_ant._get_ff_angles(my_gff.angles, top.angles, atom_gaff.atom_gaff)
        ans_angle={'H1-C-H3': ('hc-c3-hc', (39.43, 108.35)), 'H1-C-H2': ('hc-c3-hc', (39.43, 108.35)), 'H2-C-H3': ('hc-c3-hc', (39.43, 108.35)), 'H-C-H3': ('hc-c3-hc', (39.43, 108.35)), 'H-C-H2': ('hc-c3-hc', (39.43, 108.35)), 'H-C-H1': ('hc-c3-hc', (39.43, 108.35))}
        self.assertEquals(ans_angle,my_ant.topangleff)

    def test_get_FF_dihedrals(self):
        my_ant = AntechamberRunner(mol)
        my_gff = Gff()
        top = TopMol.from_file('mol.rtf')
        atom_gaff = AC()
        my_ant._get_ff_dihedrals(my_gff.dihedrals, top.dihedrals, atom_gaff.atom_gaff)
        ans_dihedral={}
        self.assertEquals(ans_dihedral,my_ant.topdihedralff)

    def test_get_FF_imdihedrals(self):
        my_ant = AntechamberRunner(mol)
        my_gff = Gff()
        top = TopMol.from_file('mol.rtf')
        atom_gaff = AC()
        my_ant._get_ff_imdihedrals(my_gff.imdihedrals, top.imdihedrals, atom_gaff.atom_gaff)
        ans_imdihedral={}
        self.assertEquals(ans_imdihedral,my_ant.topimdihedralff)