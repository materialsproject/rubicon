import os
from unittest import TestCase
from pymatgen import Molecule
from rubicon.io.mopacio.mopacio import MopTask, MopOutput

__author__ = 'xiaohuiqu'


test_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                 "test_files", "mopac"))

coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]
mol = Molecule(["C", "H", "H", "H", "Cl"], coords)


class TestMopTask(TestCase):

    def elementary_io_verify(self, moptask):
        self.to_and_from_dict_verify(moptask)
        self.from_string_verify(str(moptask), moptask.as_dict())

    def to_and_from_dict_verify(self, moptask):
        """
        Helper function. This function should be called in each specific test.
        """
        d1 = moptask.as_dict()
        mop2 = MopTask.from_dict(d1)
        d2 = mop2.as_dict()
        self.assertEqual(d1, d2)

    def from_string_verify(self, contents, ref_dict):
        moptask = MopTask.from_string(contents)
        d2 = moptask.as_dict()
        self.assertEqual(ref_dict, d2)

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
        self.elementary_io_verify(mop)
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
        self.elementary_io_verify(mop)

    def test_use_precise(self):
        mop = MopTask(mol, charge=0, jobtype="opt", title="first test methane")
        self.assertFalse("PRECISE" in mop.keywords)
        mop.use_precise()
        self.assertTrue("PRECISE" in mop.keywords)
        mop.use_precise(False)
        self.assertFalse("PRECISE" in mop.keywords)
        self.elementary_io_verify(mop)

    def test_use_bfgs(self):
        mop = MopTask(mol, charge=0, jobtype="opt", title="first test methane")
        self.assertTrue("EF" in mop.keywords)
        mop.use_bfgs()
        self.assertTrue("BFGS" in mop.keywords)
        self.assertFalse("EF" in mop.keywords)


class TestMopOutput(TestCase):

    def test_successful_message(self):
        moo = MopOutput(os.path.join(test_dir, "ch3cl_ef.out"))
        self.assertFalse(moo.data["has_error"])
        self.assertEqual(moo.data["jobtype"], "OPT")
        self.assertTrue(moo.data["gracefully_terminated"])

    def test_energies(self):
        moo = MopOutput(os.path.join(test_dir, "ch3cl_ef.out"))
        self.assertAlmostEqual(moo.data["energies"][0][1], -6.34798632251, 3)
        self.assertAlmostEqual(moo.data["energies"][1][1], -1685.76391, 3)

    def test_parse_structures(self):
        moo = MopOutput(os.path.join(test_dir, "ch3cl_ef.out"))
        ans_mol1 = """Molecule Summary (H9 C4 N1 O3)
Reduced Formula: H9C4NO3
Charge = 0, Spin Mult = 1
Sites (17)
1 C    -1.500500     1.510500    -0.753900
2 C    -2.618000     0.880700     0.073000
3 C    -3.572900     1.988900     0.611400
4 H    -4.132700     2.465100    -0.227600
5 O    -2.857300     3.085300     1.127500
6 H    -2.369800     2.782000     1.884000
7 C    -4.534800     1.408600     1.645800
8 H    -5.165800     0.626600     1.204800
9 H    -5.200800     2.185700     2.043700
10 H    -4.003400     0.963000     2.498200
11 N    -3.390800    -0.115500    -0.704200
12 O    -1.498900     1.782200    -1.940900
13 O    -0.358100     1.767200    -0.079200
14 H    -2.169800     0.329600     0.939000
15 H    -2.775300    -0.817100    -1.055900
16 H    -3.863300     0.327900    -1.464000
17 H     0.271500     2.195900    -0.650900"""
        self.assertEqual(ans_mol1, str(moo.data["molecules"][0]))
        ans_mol2 = """Molecule Summary (H9 C4 N1 O3)
Reduced Formula: H9C4NO3
Charge = 0, Spin Mult = 1
Sites (17)
1 C    -1.517682     1.490044    -0.707280
2 C    -2.693150     0.817317    -0.013378
3 C    -3.626038     1.922600     0.582224
4 H    -4.351657     2.287764    -0.183361
5 O    -2.880302     3.079022     0.858364
6 H    -2.121236     2.897494     1.451907
7 C    -4.328070     1.419830     1.836240
8 H    -4.915870     0.516629     1.620984
9 H    -5.016872     2.181694     2.224643
10 H    -3.625366     1.180514     2.640740
11 N    -3.473778    -0.006085    -0.940194
12 O    -1.354202     1.649287    -1.885729
13 O    -0.581079     1.866540     0.201170
14 H    -2.306531     0.161269     0.812776
15 H    -2.960140    -0.811422    -1.264447
16 H    -3.793316     0.519019    -1.743723
17 H     0.183750     2.356539    -0.192487"""
        self.assertEqual(ans_mol2, str(moo.data["molecules"][1]))

    def test_geom_opt_failed(self):
        moo = MopOutput(os.path.join(test_dir, "quino_salt_geom_failed.out"))
        self.assertTrue(moo.data["has_error"])
        self.assertTrue(moo.data["errors"], ['Geometry optimization failed'])

    def test_scf_failed(self):
        moo = MopOutput(os.path.join(test_dir, "quino_salt_scf_failed.out"))
        self.assertTrue(moo.data["has_error"])
        self.assertTrue(moo.data["errors"], ['Bad SCF convergence'])
