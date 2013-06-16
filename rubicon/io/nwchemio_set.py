#!/usr/bin/env python

"""
TODO: Modify module doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "6/15/13"


from pymatgen.io.nwchemio import NwTask, NwInput


class JCESRDeltaSCFInputSet(object):
    """
    Current works only for DFT based functionals.
    """

    def __init__(self):
        self.functional = "b3lyp"
        self.geom_opt_bset = "aug-cc-pVDZ"
        self.scf_bset = "aug-cc-pVTZ"

    def get_nwinput(self, mol):
        tasks = [
            NwTask.dft_task(mol, operation="optimize", xc=self.functional,
                            basis_set=self.geom_opt_bset),
            NwTask.dft_task(mol, operation="freq", xc=self.functional,
                            basis_set=self.geom_opt_bset),
            NwTask.dft_task(mol, operation="energy", xc=self.functional,
                            basis_set=self.scf_bset),
            NwTask.dft_task(mol, charge=mol.charge + 1, operation="energy",
                            xc=self.functional, basis_set=self.scf_bset),
            NwTask.dft_task(mol, charge=mol.charge - 1, operation="energy",
                            xc=self.functional, basis_set=self.scf_bset)
        ]
        return NwInput(mol, tasks)

    def write_input(self, mol, filename):
        nwi = self.get_nwinput(mol)
        nwi.write_file(filename)


import unittest

from pymatgen import Molecule


class JCESRDeltaSCFInputSetTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)

    def test_get_nwinput(self):
        jis = JCESRDeltaSCFInputSet()
        nwi = jis.get_nwinput(self.mol)
        self.assertEqual(nwi.tasks[0].theory, "dft")
        self.assertEqual(nwi.tasks[0].theory_directives["xc"], "b3lyp")
        self.assertEqual(nwi.tasks[0].basis_set["C"], "aug-cc-pVDZ")
        self.assertEqual(nwi.tasks[-1].basis_set["C"], "aug-cc-pVTZ")

if __name__ == "__main__":
    unittest.main()