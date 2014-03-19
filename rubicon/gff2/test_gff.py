from unittest import TestCase
import unittest
from gff import GFF

__author__ = 'navnidhirajput'


class TestGFF(TestCase):


     def test_from_dict(self):

        my_gff = GFF()
        my_gff.read_forcefield_para('mol.prm')
        my_gff.read_forcefield_para('mol.frcmod')
        my_gff.read_mass('mol.prm')

        d1 = my_gff.to_dict
        qc2 = my_gff.from_dict(d1)
        d2 = qc2.to_dict
        self.assertEqual(d1, d2)



if __name__ == '__main__':
    unittest.main()