# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import json
import os
import unittest
from unittest import TestCase

from pymatgen.io.smart import read_mol
from rubicon.utils.snl.egsnl import get_meta_from_structure

__author__ = 'xiaohuiqu'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files")


class TestEGStructureNL(TestCase):
    def test_get_meta_from_structure(self):
        mol = read_mol(os.path.join(test_dir, "thiophene1.xyz"))
        meta = get_meta_from_structure(mol)
        json.dumps(meta, indent=4, sort_keys=True)
        meta_ref_text = '''
{
    "anonymized_formula": "AB3C3D5",
    "chemsystem": "C-Cl-H-S",
    "elements": ["C", "Cl", "H", "S"],
    "formula": "H3 C5 S1 Cl3",
    "inchi": "InChI=1S/C5H3Cl3S/c6-5(7,8)4-1-2-9-3-4/h1-3H",
    "is_valid": true,
    "known_bonds": [[0, 1], [1, 2], [1, 10], [1, 11], [2, 3], [2, 5], [3, 4],
                    [3, 9], [5, 6], [5, 7], [7, 8], [7, 9]],
    "nelements": 4,
    "nsites": 12,
    "reduced_cell_formula": "H3C5SCl3",
    "reduced_cell_formula_abc": "C5 Cl3 H3 S1"
}'''
        meta_ref = json.loads(meta_ref_text)
        self.assertEqual(meta, meta_ref)


if __name__ == "__main__":
    unittest.main()
