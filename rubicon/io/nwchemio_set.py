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
    Currently works only for DFT based functionals.
    """

    def __init__(self):
        pass

    def get_nwinput(self, mol):
        functional = "b3lyp"
        geom_opt_bset = {}
        scf_bset = {}
        for el in mol.composition.elements:
            geom_opt_bset[el.symbol] = "6-31+g*" if el.Z <= 18 else "6-311G*"
            scf_bset[el.symbol] = "6-31++g**" if el.Z <= 18 else "6-311G*"

        if len(mol) > 1:
            #Insert opt and freq job for molecules with more than one atom.
            tasks = [
                NwTask.dft_task(mol, operation="optimize", xc=functional,
                                basis_set=geom_opt_bset),
                NwTask.dft_task(mol, operation="freq", xc=functional,
                                basis_set=geom_opt_bset)
            ]
        else:
            tasks = []

        tasks.extend([
            NwTask.dft_task(mol, operation="energy", xc=functional,
                            basis_set=scf_bset),
            NwTask.dft_task(mol, charge=mol.charge + 1, operation="energy",
                            xc=functional, basis_set=scf_bset),
            NwTask.dft_task(mol, charge=mol.charge - 1, operation="energy",
                            xc=functional, basis_set=scf_bset)
        ])

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
        self.assertEqual(nwi.tasks[0].basis_set["C"], "6-31+g*")
        self.assertEqual(nwi.tasks[-1].basis_set["C"], "6-311g*")

        atom = Molecule(["C"], [[0, 0, 0]])
        nwi = jis.get_nwinput(atom)
        self.assertEqual(len(nwi.tasks), 4)

if __name__ == "__main__":
    unittest.main()