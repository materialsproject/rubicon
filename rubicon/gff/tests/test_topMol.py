from unittest import TestCase
import unittest
from rubicon.gff.topology import TopMol

__author__ = 'navnidhirajput'


class TestTopMol(TestCase):
    def test_from_file(self):
        top = TopMol.from_file('mol.rtf')
        ans_bonds = [['C', 'H'], ['C', 'H1'], ['C', 'O'], ['C', 'C2'], ['O', 'C1'], ['C1', 'O1'], ['C1', 'O2'], ['O2', 'C2'], ['C2', 'H2'], ['C2', 'H3']]
        ans_angles=[['C', 'O', 'C1'], ['C', 'C2', 'O2'], ['C', 'C2', 'H2'], ['C', 'C2', 'H3'], ['H', 'C', 'H1'], ['H', 'C', 'O'], ['H', 'C', 'C2'], ['H1', 'C', 'O'], ['H1', 'C', 'C2'], ['O', 'C', 'C2'], ['O', 'C1', 'O1'], ['O', 'C1', 'O2'], ['C1', 'O2', 'C2'], ['O1', 'C1', 'O2'], ['O2', 'C2', 'H2'], ['O2', 'C2', 'H3'], ['H2', 'C2', 'H3']]
        ans_dihedrals=[['H', 'C', 'O', 'C1'], ['H1', 'C', 'O', 'C1'], ['C2', 'C', 'O', 'C1'], ['H', 'C', 'C2', 'O2'], ['H1', 'C', 'C2', 'O2'], ['O', 'C', 'C2', 'O2'], ['H', 'C', 'C2', 'H2'], ['H1', 'C', 'C2', 'H2'], ['O', 'C', 'C2', 'H2'], ['H', 'C', 'C2', 'H3'], ['H1', 'C', 'C2', 'H3'], ['O', 'C', 'C2', 'H3'], ['C', 'O', 'C1', 'O1'], ['C', 'O', 'C1', 'O2'], ['O', 'C1', 'O2', 'C2'], ['O1', 'C1', 'O2', 'C2'], ['C1', 'O2', 'C2', 'C'], ['C1', 'O2', 'C2', 'H2'], ['C1', 'O2', 'C2', 'H3']]
        ans_imdihedrals=[['O1', 'O', 'C1', 'O2']]
        self.assertEquals(ans_bonds,top.bonds)
        self.assertEquals(ans_angles,top.angles)
        self.assertEquals(ans_dihedrals,top.dihedrals)
        self.assertEquals(ans_imdihedrals,top.imdihedrals)


if __name__ == '__main__':
    unittest.main()