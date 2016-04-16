# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import unittest
import os

from rubicon.io.lammps.input import DictLammpsInput


__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestLammpsData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.lammps_input = DictLammpsInput.from_file("NVT",
                                                     os.path.join(module_dir,"test_files", "NVT.json"),
                                                     data_file=os.path.join(module_dir,"test_files", "nvt.data"),
                                                     is_forcefield=True)

    def test_string_rep(self):
        self.lammps_input.config_dict["read_data"] = "nvt.data"
        with open(os.path.join(module_dir, "test_files", "nvt.inp")) as f:
            self.assertEqual(str(self.lammps_input), "".join(f.readlines()))

if __name__ == "__main__":
    unittest.main()
