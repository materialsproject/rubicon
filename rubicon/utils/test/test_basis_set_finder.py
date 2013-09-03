from unittest import TestCase
import unittest
from rubicon.utils.basis_set_finder import NWChemBasisSetFinder

__author__ = 'xiaohuiqu'


class TestNWChemBasisSetFinder(TestCase):
    def setUp(self):
        self.finder = NWChemBasisSetFinder()

    def test_basis_set_exist(self):
        self.assertEqual(self.finder.basis_set_exist('Fe', '6-311++G**'), False)
        self.assertEqual(self.finder.basis_set_exist('F', '6-311++G**'), True)
        self.assertEqual(self.finder.basis_set_exist('Fe', '6-31G*'), True)
        self.assertEqual(self.finder.basis_set_exist('Uuu', 'CRENBL ECP'), True)
        self.assertEqual(self.finder.basis_set_exist('Uuu', '6-31G*'), False)

    def test_get_first_available_basis_set(self):
        self.assertEqual(self.finder.get_first_available_basis_set('Fe', ['6-311++G**', '6-31G*']), '6-31G*')
        self.assertEqual(self.finder.get_first_available_basis_set('F', ['6-311++G**', '6-31G*']), '6-311++G**')
        self.assertEqual(self.finder.get_first_available_basis_set('Rn', ['6-311++G**', '6-31G*']), None)

if __name__ == '__main__':
    unittest.main()