from unittest import TestCase
import unittest
from GFF import GFF

__author__ = 'navnidhirajput'


class TestGFF(TestCase):
    #def test_from_antechamber_read_file(self):
    #    self.fail()

    #def test_to_dict(self):
    #    self.fail()

    def test_from_dict(self):
        my_ff=GFF()
        my_ff.read_forcefield_para('mol1.prm')
        my_ff.read_mass('mol1.rtf')

        d1 = my_ff.to_dict
        qc2 = my_ff.from_dict(d1)
        d2 = qc2.to_dict
        self.assertEqual(d1, d2)
        print my_ff.to_dict


if __name__ == '__main__':
    unittest.main()