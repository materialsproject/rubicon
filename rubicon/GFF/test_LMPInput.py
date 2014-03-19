from unittest import TestCase
from rubicon.GFF.gff import GFF
from rubicon.GFF.lamppsio import LMPInput
from rubicon.GFF.topology import TopMol

__author__ = 'navnidhirajput'


class TestLMPInput(TestCase):
    def test_write_data_file(self):
        self.fail()

    def test_get_data_file(self):
        self.fail()

    def test_set_atom(self):
        self.fail()

    def test_set_coeff(self):
        my_lampps=LMPInput()
        my_gff = GFF()
        top = TopMol.from_file('mol.rtf')

        lampps_data=my_lampps.set_coeff(my_gff,top)
        ans='''
        LAMMPS Data File


1 atom type
1 bond type
1 angle type
0 dihedral type
0 improper dihedral type


4 atoms
6 bonds
0 angles
0 dihedrals


Masses


1 12.01
2 1.008


Pair Coeffs


1 1.908 0.1094
2 1.487 0.0157


Bond Coeffs


1 337.3 1.092


Angle Coeffs


1 39.43 108.35


Dihedral Coeffs




Imp Dihedral Coeffs
'''

        self.assertEquals(ans,lampps_data)

    def test_set_atom(self):
        self.fail()

    def test_set_dimension(self):
        self.fail()