from unittest import TestCase
from pymatgen import Molecule
from rubicon.io.mopacio.mopacio import MopTask

__author__ = 'xiaohuiqu'


coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]
mol = Molecule(["C", "H", "H", "H", "Cl"], coords)


class TestMopTask(TestCase):

    def test_str(self):
        mop = MopTask(mol, charge=0, jobtype="opt", title="first test methane",
                      optional_params={"PRECISE": None, "CYCLES": 500})
        ans = """PM7 EF CHARGE=0 PRECISE XYZ CYCLES=500
first test methane

 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000

"""
        self.assertEqual(ans, str(mop))
        mop = MopTask(mol, charge=-1, jobtype="opt", sqm_method="pm6-dh2",
                      title="first test methane. Use a long title to exceed the limit"
                            "of one line purposly. Doest it work with three line long "
                            "text? I don't know, just try make up more text here. Oh,"
                            "I even can't figure out how write two lines title")
        ans = """PM6-DH2 EF CHARGE=-1 XYZ
first test methane. Use a long title to exceed the limitof one line purposly.
Doest it work with three line long text? I don't know, just try make up more
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000

"""
        self.assertEqual(ans, str(mop))