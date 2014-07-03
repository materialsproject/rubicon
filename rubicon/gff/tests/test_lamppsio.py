import os
from unittest import TestCase
import unittest
from pymatgen import Molecule
from pymatgen.packmol.packmol import PackmolRunner
from rubicon.gff.antechamberio import AntechamberRunner
from rubicon.gff.lamppsio import LmpInput

__author__ = 'navnidhirajput'


coords_n1c=[[4.522,   8.999,   5.512],
           [6.666,   9.383,   5.971],
           [7.522,   9.623,   6.578],
           [6.567,   9.114,   4.645],
           [7.322,   9.072,   3.879],
           [5.388,   9.312,   6.496],
           [5.231,   8.881,   4.373],
           [5.041,   9.544,   7.906],
           [5.944,   9.842,   8.432],
           [4.305,  10.344,   7.985],
           [4.653,   8.630,   8.355],
           [4.699,   8.519,   3.03],
           [3.676,   8.890,   2.982],
           [5.285,   9.084,   2.312],
           [4.774,   7.019,   2.765],
           [4.386,   6.818,   1.765],
           [5.802,   6.657,   2.807],
           [4.176,   6.448,   3.480],
           [3.050,   8.855,   5.663],
           [2.657,   8.107,   4.974],
           [2.796,   8.543,   6.676],
           [2.542,   9.803,   5.463]]

coords_pc=[[2.516,   1.762,   7.803],
           [3.382,   2.859,   7.178],
           [4.673,   2.240,   7.066],
           [3.476,   3.751,   7.797],
           [3.037,   3.136,   6.177],
           [4.526,   0.887,   7.121],
           [3.246,   0.565,   7.439],
           [5.398,   0.103,   6.918],
           [1.090,   1.684,   7.302],
           [2.539,   1.829,  8.896],
           [1.069,   1.564,   6.216],
           [0.552,   2.599,   7.564],
           [0.568,   0.839,   7.755]]



coords_tfn=[[3.734,   8.436,  10.848],
           [5.713,   9.869,   7.928],
           [4.712,   9.816,   8.816],
           [4.981,   8.437,  10.088],
           [6.222,   8.807,  10.764],
           [4.651,  11.005,  9.445],
           [5.009,   7.118,   9.161],
           [5.914,   5.626,   7.385],
           [3.564,   9.647,   8.144],
           [6.318,   6.404,   8.552],
           [7.545,   7.196,   8.484],
           [7.774,   4.365,   9.475],
           [5.668,   4.219,   9.986],
           [6.692,   5.074,   9.850],
           [6.947,   5.614,  11.049]]


tfn = Molecule(["O", "F", "C", "S", "O", "F", "N", "O", "F", "S", "O", "F", "F", "C", "F"], coords_tfn)
n1c = Molecule(["C", "C", "H", "C", "H", "N", "N", "C", "H", "H", "H", "C", "H", "H", "C", "H", "H", "H", "C", "H", "H", "H"], coords_n1c)
pc = Molecule(["C", "C", "O", "H", "H", "C", "O", "O", "C", "H", "H", "H", "H"], coords_pc)


class TestLmpInput(TestCase):

    def tearDown(self):

        super(TestLmpInput, self).tearDown()

    @classmethod
    def setUpClass(cls):
        mols=[tfn,n1c,pc]
        cls.ffmol = AntechamberRunner(mols)
        cls.gff_list, top_list = cls.ffmol.get_ff_top_mol(mols,'mol.pdb')
        cls.mols_in_box = PackmolRunner(mols, [{"number":1,"inside box":[0.,0.,0.,40.,40.,40.]},{"number":1},{"number":1}])
        super(TestLmpInput, cls).setUpClass()

    def test__set_gff_types(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        gff_types = data_lammps._set_gff_types(self.ffmol, self.mols_in_box)

        ans="""LAMMPS Data File

# 1 mol 1 molecule
# 1 mol 2 molecule
# 1 mol 3 molecule


14 atom type
23 bond type
41 angle type
56 dihedral type
6 improper type
"""
        self.assertEqual(ans,gff_types)

    def test__set_top_types(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        top_types = data_lammps._set_top_types(self.ffmol, self.mols_in_box)


        ans="""50 atoms
49 bonds
87 angles
114 dihedrals
6 impropers
"""
        self.assertEqual(ans,top_types)

    def test__set_box_dimensions(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        box_dimension = data_lammps._set_box_dimensions(self.mols_in_box)


        ans="""0.0 40.0 xlo  xhi
0.0 40.0 ylo  yhi
0.0 40.0 zlo  zhi
"""
        self.assertEqual(ans,box_dimension)

    def test__set_masses(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        masses = data_lammps._set_masses(self.ffmol, self.mols_in_box)


        ans="""Masses

1 32.06 # sy 1
2 19.0 # f 1
3 14.01 # ne 1
4 16.0 # o 1
5 32.06 # s6 1
6 12.01 # c3 1
7 12.01 # cc 2
8 14.01 # na 2
9 1.008 # h1 2
10 1.008 # h4 2
11 12.01 # cd 2
12 1.008 # hc 2
13 12.01 # c 3
14 16.0 # os 3

"""
        self.assertEqual(ans,masses)



    def test__set_pair_coeffs(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        pair_coeff = data_lammps._set_pair_coeffs(self.ffmol, self.mols_in_box)

        ans="""Pair Coeffs

1 0.25 3.5636 # sy 1
2 0.061 3.11815 # f 1
3 0.17 3.2500032 # ne 1
4 0.21 2.95992616 # o 1
5 0.25 3.5636 # s6 1
6 0.1094 3.3996744 # c3 1
7 0.086 3.3996744 # cc 2
8 0.17 3.2500032 # na 2
9 0.0157 2.4713566 # h1 2
10 0.015 2.5105562 # h4 2
11 0.086 3.3996744 # cd 2
12 0.0157 2.6495366 # hc 2
13 0.086 3.3996744 # c 3
14 0.17 3.00001666 # os 3

"""
        self.assertEqual(ans,pair_coeff)


    def test__set_bond_coeffs(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        bond_coeff = data_lammps._set_bond_coeffs(self.ffmol, self.mols_in_box)


        ans="""Bond Coeffs

1 257.1  1.752 # ne sy 1
2 444.6  1.551 # ne s6 1
3 363.8  1.344 # c3 f 1
4 541.1  1.436 # o s6 1
5 254.0  1.774 # c3 s6 1
6 248.9  1.782 # c3 sy 1
7 493.0  1.466 # o sy 1
8 334.7  1.456 # c3 na 2
9 350.1  1.083 # cc h4 2
10 337.3  1.499 # c3 cc 2
11 438.8  1.371 # cc na 2
12 303.1  1.535 # c3 c3 2
13 337.3  1.092 # c3 hc 2
14 438.8  1.371 # cd na 2
15 350.1  1.083 # cd h4 2
16 335.9  1.093 # c3 h1 2
17 504.0  1.371 # cc cd 2
18 411.3  1.343 # c os 3
19 301.5  1.439 # c3 os 3
20 337.3  1.092 # c3 hc 3
21 303.1  1.535 # c3 c3 3
22 335.9  1.093 # c3 h1 3
23 648.0  1.214 # c o 3

"""
        self.assertEqual(ans,bond_coeff)

    def test__set_angle_coeffs(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        angle_coeff = data_lammps._set_angle_coeffs(self.ffmol, self.mols_in_box)

        ans="""Angle Coeffs

1 43.1 107.06 # ne sy o 1
2 81.22 109.67 # f c3 s6 1
3 46.66 119.73 # o s6 o 1
4 9.72 119.18 # s6 ne sy 1
5 39.68 112.95 # c3 s6 ne 1
6 41.28 108.48 # c3 sy o 1
7 81.22 109.67 # f c3 sy 1
8 41.66 108.32 # c3 s6 o 1
9 39.42 103.12 # c3 sy ne 1
10 71.26 107.16 # f c3 f 1
11 45.3 121.88 # o sy o 1
12 45.1 116.41 # ne s6 o 1
13 49.9 109.45 # h1 c3 na 2
14 62.56 125.09 # c3 na cc 2
15 72.91 109.42 # cc cd na 2
16 65.44 122.99 # c3 cc na 2
17 72.91 109.42 # cd cc na 2
18 46.37 110.05 # c3 c3 hc 2
19 73.65 109.33 # na cc na 2
20 50.22 119.66 # h4 cc na 2
21 47.19 129.11 # cc cd h4 2
22 39.18 109.55 # h1 c3 h1 2
23 47.19 129.11 # cd cc h4 2
24 46.36 110.07 # c3 c3 h1 2
25 68.94 109.9 # cc na cc 2
26 63.88 128.01 # cc na cd 2
27 47.2 110.86 # cc c3 hc 2
28 39.43 108.35 # hc c3 hc 2
29 65.73 112.81 # c3 c3 na 2
30 62.56 125.09 # c3 na cd 2
31 50.22 119.66 # h4 cd na 2
32 63.63 115.14 # c os c3 3
33 46.37 110.05 # c3 c3 hc 3
34 46.36 110.07 # c3 c3 h1 3
35 50.84 108.82 # h1 c3 os 3
36 63.21 110.63 # c3 c3 c3 3
37 67.78 108.42 # c3 c3 os 3
38 39.43 108.35 # hc c3 hc 3
39 75.93 123.33 # o c os 3
40 76.45 111.38 # os c os 3
41 39.18 109.55 # h1 c3 h1 3

"""

        self.assertEqual(ans,angle_coeff)

    def test__set_dihedral_coeffs(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        dihedral_coeff = data_lammps._set_dihedral_coeffs(self.ffmol, self.mols_in_box)


        ans="""Dihedral Coeffs

1  6.667  2  180.0 # o s6 ne sy 1
2  0.144  3  0.0 # f c3 s6 o 1
3  6.8  1  180.0 # o sy ne s6 1
4  0.5  3  180.0 # o sy ne s6 1
5  6.667  2  180.0 # c3 s6 ne sy 1
6  6.8  1  180.0 # c3 sy ne s6 1
7  0.5  3  180.0 # c3 sy ne s6 1
8  0.144  3  0.0 # f c3 s6 ne 1
9  0.144  3  0.0 # f c3 sy o 1
10  0.144  3  0.0 # f c3 sy ne 1
11  4.0  2  180.0 # h4 cc cd na 2
12  4.0  2  180.0 # h4 cc cd h4 2
13  4.0  2  180.0 # na cc cd na 2
14  0.156  3  0.0 # h1 c3 c3 hc 2
15  0.0  2  0.0 # cd na c3 h1 2
16  1.7  2  180.0 # cd na cc na 2
17  0.0  3  0.0 # hc c3 cc na 2
18  4.0  2  180.0 # h4 cd cc h4 2
19  4.0  2  180.0 # na cd cc na 2
20  1.7  2  180.0 # c3 na cc c3 2
21  1.7  2  180.0 # c3 na cc h4 2
22  1.7  2  180.0 # c3 na cd h4 2
23  1.7  2  180.0 # cc na cc h4 2
24  1.7  2  180.0 # c3 na cd cc 2
25  1.7  2  180.0 # cc na cd h4 2
26  1.7  2  180.0 # c3 cc na c3 2
27  1.7  2  180.0 # cc cd na cc 2
28  1.7  2  180.0 # c3 na cc cd 2
29  1.7  2  180.0 # cc na cd cc 2
30  4.0  2  180.0 # h4 cd cc na 2
31  1.7  2  180.0 # c3 na cc na 2
32  1.7  2  180.0 # cc na cc cd 2
33  0.156  3  0.0 # hc c3 c3 na 2
34  0.0  2  0.0 # c3 c3 na cc 2
35  0.0  2  0.0 # c3 c3 na cd 2
36  1.7  2  180.0 # c3 cc na cd 2
37  1.7  2  180.0 # cc na cc na 2
38  1.7  2  180.0 # c3 cc na cc 2
39  0.0  2  0.0 # cc na c3 h1 2
40  0.25  1  0.0 # hc c3 c3 os 3
41  0.0  3  0.0 # hc c3 c3 os 3
42  2.7  2  180.0 # c3 os c os 3
43  0.25  1  0.0 # h1 c3 c3 os 3
44  0.0  3  0.0 # h1 c3 c3 os 3
45  0.156  3  0.0 # h1 c3 c3 h1 3
46  0.156  3  0.0 # h1 c3 c3 hc 3
47  0.8  1  180.0 # c os c3 c3 3
48  0.383  3  0.0 # c os c3 c3 3
49  1.4  1  180.0 # c3 os c o 3
50  2.7  2  180.0 # c3 os c o 3
51  0.156  3  0.0 # c3 c3 c3 h1 3
52  0.16  3  0.0 # c3 c3 c3 hc 3
53  0.144  3  0.0 # os c3 c3 os 3
54  1.175  2  0.0 # os c3 c3 os 3
55  0.383  3  0.0 # c os c3 h1 3
56  0.156  3  0.0 # c3 c3 c3 os 3

"""


        self.assertEqual(ans,dihedral_coeff)


    def test__set_improper_coeffs(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        imdihedrals_coeff = data_lammps._set_improper_coeffs(self.ffmol, self.mols_in_box)

        ans="""Imp Dihedral Coeffs

1 1.1  -1.0  2.0 # c3 cc na cd 2
2 1.1  -1.0  2.0 # cc h4 cd na 2
3 1.1  -1.0  2.0 # cd h4 cc na 2
4 1.1  -1.0  2.0 # c3 cc na cc 2
5 1.1  -1.0  2.0 # c3 na cc na 2
6 10.5  -1.0  2.0 # o os c os 3

"""

        self.assertEqual(ans,imdihedrals_coeff)

    def test__set_atom(self):
        data_lammps=LmpInput(self.mols_in_box,self.ffmol)
        atoms = data_lammps._set_atom(self.ffmol, self.mols_in_box)
        ans="""Atoms

1  1  4  3.841  38.685  36.79 # 1  o O
2  1  2  1.427  38.189  39.692 # 1  f F
3  1  6  2.269  38.874  38.907 # 1  c3 C
4  1  5  2.826  37.843  37.418 # 1  s6 S
5  1  4  1.59  37.562  36.691 # 1  o O1
6  1  2  1.631  39.979  38.478 # 1  f F1
7  1  3  3.568  36.597  38.123 # 1  ne N
8  1  4  3.699  34.562  39.552 # 1  o O2
9  1  2  3.313  39.264  39.652 # 1  f F2
10  1  1  2.911  35.172  38.485 # 1  sy S1
11  1  4  1.452  35.104  38.556 # 1  o O3
12  1  2  2.877  32.876  37.123 # 1  f F3
13  1  2  4.682  34.042  36.8 # 1  f F4
14  1  6  3.354  34.126  36.968 # 1  c3 C1
15  1  2  2.824  34.635  35.847 # 1  f F5
16  2  7  2.604  5.594  7.056 # 2  cc C
17  2  7  3.435  7.464  6.181 # 2  cc C1
18  2  10  3.46  8.506  5.909 # 2  h4 H
19  2  11  4.333  6.471  5.966 # 2  cd C2
20  2  10  5.291  6.484  5.476 # 2  h4 H1
21  2  8  2.364  6.905  6.855 # 2  na N
22  2  8  3.804  5.315  6.511 # 2  na N1
23  2  6  1.16  7.634  7.279 # 2  c3 C3
24  2  9  1.236  8.655  6.912 # 2  h1 H2
25  2  9  0.27  7.17  6.853 # 2  h1 H3
26  2  9  1.092  7.651  8.367 # 2  h1 H4
27  2  6  4.495  3.996  6.514 # 2  c3 C4
28  2  9  3.72  3.23  6.511 # 2  h1 H5
29  2  9  5.022  3.93  5.568 # 2  h1 H6
30  2  6  5.441  3.843  7.701 # 2  c3 C5
31  2  12  5.927  2.868  7.641 # 2  hc H7
32  2  12  6.218  4.608  7.692 # 2  hc H8
33  2  12  4.911  3.898  8.655 # 2  hc H9
34  2  6  1.684  4.635  7.722 # 2  c3 C6
35  2  12  2.241  3.862  8.253 # 2  hc H10
36  2  12  1.051  5.148  8.445 # 2  hc H11
37  2  12  1.035  4.144  6.991 # 2  hc H12
38  3  6  1.008  1.371  3.592 # 3  c3 C
39  3  6  1.397  2.451  4.605 # 3  c3 C1
40  3  14  0.77  2.009  5.819 # 3  os O
41  3  9  1.017  3.442  4.359 # 3  h1 H
42  3  9  2.479  2.496  4.765 # 3  h1 H1
43  3  13  0.49  0.679  5.718 # 3  c C2
44  3  14  0.726  0.243  4.455 # 3  os O1
45  3  4  0.095  0.001  6.613 # 3  o O2
46  3  6  2.071  1.008  2.578 # 3  c3 C3
47  3  9  0.07  1.633  3.091 # 3  h1 H2
48  3  12  2.992  0.697  3.078 # 3  hc H3
49  3  12  2.292  1.87  1.943 # 3  hc H4
50  3  12  1.731  0.191  1.938 # 3  hc H5

"""

        self.assertEqual(ans,atoms)

    def test__set_bonds(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        data_lammps._set_atom(self.ffmol, self.mols_in_box)
        bonds = data_lammps._set_bonds(self.ffmol,self.mols_in_box)

        ans="""Bonds

1 4 1 4 # 1 O S
2 3 2 3 # 1 F C
3 5 3 4 # 1 C S
4 3 3 6 # 1 C F1
5 3 3 9 # 1 C F2
6 4 4 5 # 1 S O1
7 2 4 7 # 1 S N
8 1 7 10 # 1 N S1
9 7 8 10 # 1 O2 S1
10 7 10 11 # 1 S1 O3
11 6 10 14 # 1 S1 C1
12 3 12 14 # 1 F3 C1
13 3 13 14 # 1 F4 C1
14 3 14 15 # 1 C1 F5
15 11 16 21 # 2 C N
16 11 16 22 # 2 C N1
17 10 16 34 # 2 C C6
18 9 17 18 # 2 C1 H
19 17 17 19 # 2 C1 C2
20 11 17 21 # 2 C1 N
21 15 19 20 # 2 C2 H1
22 14 19 22 # 2 C2 N1
23 8 21 23 # 2 N C3
24 8 22 27 # 2 N1 C4
25 16 23 24 # 2 C3 H2
26 16 23 25 # 2 C3 H3
27 16 23 26 # 2 C3 H4
28 16 27 28 # 2 C4 H5
29 16 27 29 # 2 C4 H6
30 12 27 30 # 2 C4 C5
31 13 30 31 # 2 C5 H7
32 13 30 32 # 2 C5 H8
33 13 30 33 # 2 C5 H9
34 13 34 35 # 2 C6 H10
35 13 34 36 # 2 C6 H11
36 13 34 37 # 2 C6 H12
37 21 38 39 # 3 C C1
38 19 38 44 # 3 C O1
39 21 38 46 # 3 C C3
40 22 38 47 # 3 C H2
41 19 39 40 # 3 C1 O
42 22 39 41 # 3 C1 H
43 22 39 42 # 3 C1 H1
44 18 40 43 # 3 O C2
45 18 43 44 # 3 C2 O1
46 23 43 45 # 3 C2 O2
47 20 46 48 # 3 C3 H3
48 20 46 49 # 3 C3 H4
49 20 46 50 # 3 C3 H5

"""

        self.assertEqual(ans,bonds)


    def test__set_angles(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        data_lammps._set_atom(self.ffmol, self.mols_in_box)
        angles = data_lammps._set_angles(self.ffmol, self.mols_in_box)

        ans="""Angles

1  8  1  4  3  #  1  O  S  C
2  3  1  4  5  #  1  O  S  O1
3  12  1  4  7  #  1  O  S  N
4  2  2  3  4  #  1  F  C  S
5  10  2  3  6  #  1  F  C  F1
6  10  2  3  9  #  1  F  C  F2
7  8  3  4  5  #  1  C  S  O1
8  5  3  4  7  #  1  C  S  N
9  2  4  3  6  #  1  S  C  F1
10  2  4  3  9  #  1  S  C  F2
11  4  4  7  10  #  1  S  N  S1
12  12  5  4  7  #  1  O1  S  N
13  10  6  3  9  #  1  F1  C  F2
14  1  7  10  8  #  1  N  S1  O2
15  1  7  10  11  #  1  N  S1  O3
16  9  7  10  14  #  1  N  S1  C1
17  11  8  10  11  #  1  O2  S1  O3
18  6  8  10  14  #  1  O2  S1  C1
19  7  10  14  12  #  1  S1  C1  F3
20  7  10  14  13  #  1  S1  C1  F4
21  7  10  14  15  #  1  S1  C1  F5
22  6  11  10  14  #  1  O3  S1  C1
23  10  12  14  13  #  1  F3  C1  F4
24  10  12  14  15  #  1  F3  C1  F5
25  10  13  14  15  #  1  F4  C1  F5
26  25  16  21  17  #  2  C  N  C1
27  14  16  21  23  #  2  C  N  C3
28  26  16  22  19  #  2  C  N1  C2
29  14  16  22  27  #  2  C  N1  C4
30  27  16  34  35  #  2  C  C6  H10
31  27  16  34  36  #  2  C  C6  H11
32  27  16  34  37  #  2  C  C6  H12
33  21  17  19  20  #  2  C1  C2  H1
34  15  17  19  22  #  2  C1  C2  N1
35  14  17  21  23  #  2  C1  N  C3
36  23  18  17  19  #  2  H  C1  C2
37  20  18  17  21  #  2  H  C1  N
38  17  19  17  21  #  2  C2  C1  N
39  30  19  22  27  #  2  C2  N1  C4
40  31  20  19  22  #  2  H1  C2  N1
41  19  21  16  22  #  2  N  C  N1
42  16  21  16  34  #  2  N  C  C6
43  13  21  23  24  #  2  N  C3  H2
44  13  21  23  25  #  2  N  C3  H3
45  13  21  23  26  #  2  N  C3  H4
46  16  22  16  34  #  2  N1  C  C6
47  13  22  27  28  #  2  N1  C4  H5
48  13  22  27  29  #  2  N1  C4  H6
49  29  22  27  30  #  2  N1  C4  C5
50  22  24  23  25  #  2  H2  C3  H3
51  22  24  23  26  #  2  H2  C3  H4
52  22  25  23  26  #  2  H3  C3  H4
53  18  27  30  31  #  2  C4  C5  H7
54  18  27  30  32  #  2  C4  C5  H8
55  18  27  30  33  #  2  C4  C5  H9
56  22  28  27  29  #  2  H5  C4  H6
57  24  28  27  30  #  2  H5  C4  C5
58  24  29  27  30  #  2  H6  C4  C5
59  28  31  30  32  #  2  H7  C5  H8
60  28  31  30  33  #  2  H7  C5  H9
61  28  32  30  33  #  2  H8  C5  H9
62  28  35  34  36  #  2  H10  C6  H11
63  28  35  34  37  #  2  H10  C6  H12
64  28  36  34  37  #  2  H11  C6  H12
65  37  38  39  40  #  3  C  C1  O
66  34  38  39  41  #  3  C  C1  H
67  34  38  39  42  #  3  C  C1  H1
68  32  38  44  43  #  3  C  O1  C2
69  33  38  46  48  #  3  C  C3  H3
70  33  38  46  49  #  3  C  C3  H4
71  33  38  46  50  #  3  C  C3  H5
72  37  39  38  44  #  3  C1  C  O1
73  36  39  38  46  #  3  C1  C  C3
74  34  39  38  47  #  3  C1  C  H2
75  32  39  40  43  #  3  C1  O  C2
76  35  40  39  41  #  3  O  C1  H
77  35  40  39  42  #  3  O  C1  H1
78  40  40  43  44  #  3  O  C2  O1
79  39  40  43  45  #  3  O  C2  O2
80  41  41  39  42  #  3  H  C1  H1
81  37  44  38  46  #  3  O1  C  C3
82  35  44  38  47  #  3  O1  C  H2
83  39  44  43  45  #  3  O1  C2  O2
84  34  46  38  47  #  3  C3  C  H2
85  38  48  46  49  #  3  H3  C3  H4
86  38  48  46  50  #  3  H3  C3  H5
87  38  49  46  50  #  3  H4  C3  H5

"""
        self.assertEqual(ans,angles)


    def test__set_dihedrals(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        data_lammps._set_atom(self.ffmol, self.mols_in_box)
        dihedrals = data_lammps._set_dihedrals(self.ffmol, self.mols_in_box)

        ans="""Dihedrals

1  2  2  3  4  5  #  1  F1  C  S  O
2  6  9  3  4  7  #  1  F1  C  S  O
3  2  9  3  4  1  #  1  F1  C  S  O
4  5  4  7  10  14  #  1  F1  C  S  O
5  5  4  7  10  14  #  1  F1  C  S  O
6  2  6  3  4  1  #  1  F1  C  S  O
7  6  6  3  4  7  #  1  F1  C  S  O
8  4  3  4  7  10  #  1  F1  C  S  O
9  2  6  3  4  5  #  1  F1  C  S  O
10  1  5  4  7  10  #  1  F1  C  S  O
11  8  7  10  14  12  #  1  F1  C  S  O
12  8  7  10  14  13  #  1  F1  C  S  O
13  8  7  10  14  15  #  1  F1  C  S  O
14  7  8  10  14  13  #  1  F1  C  S  O
15  7  8  10  14  15  #  1  F1  C  S  O
16  7  8  10  14  12  #  1  F1  C  S  O
17  1  1  4  7  10  #  1  F1  C  S  O
18  2  2  3  4  1  #  1  F1  C  S  O
19  6  2  3  4  7  #  1  F1  C  S  O
20  2  9  3  4  5  #  1  F1  C  S  O
21  7  11  10  14  12  #  1  F1  C  S  O
22  7  11  10  14  15  #  1  F1  C  S  O
23  7  11  10  14  13  #  1  F1  C  S  O
24  3  4  7  10  11  #  1  F1  C  S  O
25  3  4  7  10  11  #  1  F1  C  S  O
26  3  4  7  10  8  #  1  F1  C  S  O
27  3  4  7  10  8  #  1  F1  C  S  O
28  37  17  21  23  26  #  2  N1  C  N  C3
29  37  17  21  23  24  #  2  N1  C  N  C3
30  37  17  21  23  25  #  2  N1  C  N  C3
31  11  21  17  19  22  #  2  N1  C  N  C3
32  32  16  22  27  30  #  2  N1  C  N  C3
33  28  21  17  19  20  #  2  N1  C  N  C3
34  24  34  16  22  27  #  2  N1  C  N  C3
35  23  20  19  22  16  #  2  N1  C  N  C3
36  37  16  21  23  24  #  2  N1  C  N  C3
37  37  16  21  23  25  #  2  N1  C  N  C3
38  30  19  17  21  16  #  2  N1  C  N  C3
39  37  16  21  23  26  #  2  N1  C  N  C3
40  29  21  16  22  27  #  2  N1  C  N  C3
41  13  19  22  27  29  #  2  N1  C  N  C3
42  13  19  22  27  28  #  2  N1  C  N  C3
43  21  18  17  21  16  #  2  N1  C  N  C3
44  37  16  22  27  28  #  2  N1  C  N  C3
45  15  22  16  34  37  #  2  N1  C  N  C3
46  15  22  16  34  36  #  2  N1  C  N  C3
47  15  22  16  34  35  #  2  N1  C  N  C3
48  20  20  19  22  27  #  2  N1  C  N  C3
49  12  28  27  30  33  #  2  N1  C  N  C3
50  12  28  27  30  32  #  2  N1  C  N  C3
51  12  28  27  30  31  #  2  N1  C  N  C3
52  15  21  16  34  36  #  2  N1  C  N  C3
53  14  21  16  22  19  #  2  N1  C  N  C3
54  34  34  16  22  19  #  2  N1  C  N  C3
55  15  21  16  34  35  #  2  N1  C  N  C3
56  15  21  16  34  37  #  2  N1  C  N  C3
57  37  16  22  27  29  #  2  N1  C  N  C3
58  25  17  19  22  16  #  2  N1  C  N  C3
59  10  18  17  19  20  #  2  N1  C  N  C3
60  9  18  17  19  22  #  2  N1  C  N  C3
61  29  22  16  21  23  #  2  N1  C  N  C3
62  35  22  16  21  17  #  2  N1  C  N  C3
63  24  34  16  21  23  #  2  N1  C  N  C3
64  36  34  16  21  17  #  2  N1  C  N  C3
65  31  22  27  30  31  #  2  N1  C  N  C3
66  22  17  19  22  27  #  2  N1  C  N  C3
67  31  22  27  30  33  #  2  N1  C  N  C3
68  31  22  27  30  32  #  2  N1  C  N  C3
69  26  19  17  21  23  #  2  N1  C  N  C3
70  12  29  27  30  32  #  2  N1  C  N  C3
71  12  29  27  30  33  #  2  N1  C  N  C3
72  12  29  27  30  31  #  2  N1  C  N  C3
73  19  18  17  21  23  #  2  N1  C  N  C3
74  33  19  22  27  30  #  2  N1  C  N  C3
75  48  47  38  44  43  #  3  O1  C  C1  H
76  38  44  38  46  50  #  3  O1  C  C1  H
77  38  44  38  46  50  #  3  O1  C  C1  H
78  38  44  38  46  49  #  3  O1  C  C1  H
79  38  44  38  46  49  #  3  O1  C  C1  H
80  38  44  38  46  48  #  3  O1  C  C1  H
81  38  44  38  46  48  #  3  O1  C  C1  H
82  43  46  38  44  43  #  3  O1  C  C1  H
83  43  46  38  44  43  #  3  O1  C  C1  H
84  39  39  40  43  44  #  3  O1  C  C1  H
85  46  39  38  46  48  #  3  O1  C  C1  H
86  42  47  38  46  49  #  3  O1  C  C1  H
87  44  39  40  43  45  #  3  O1  C  C1  H
88  44  39  40  43  45  #  3  O1  C  C1  H
89  46  39  38  46  50  #  3  O1  C  C1  H
90  46  39  38  46  49  #  3  O1  C  C1  H
91  43  38  39  40  43  #  3  O1  C  C1  H
92  43  38  39  40  43  #  3  O1  C  C1  H
93  41  47  38  39  42  #  3  O1  C  C1  H
94  43  39  38  44  43  #  3  O1  C  C1  H
95  43  39  38  44  43  #  3  O1  C  C1  H
96  47  44  38  39  40  #  3  O1  C  C1  H
97  47  44  38  39  40  #  3  O1  C  C1  H
98  40  44  38  39  41  #  3  O1  C  C1  H
99  40  44  38  39  41  #  3  O1  C  C1  H
100  42  47  38  46  50  #  3  O1  C  C1  H
101  39  40  43  44  38  #  3  O1  C  C1  H
102  40  47  38  39  40  #  3  O1  C  C1  H
103  40  47  38  39  40  #  3  O1  C  C1  H
104  40  44  38  39  42  #  3  O1  C  C1  H
105  40  44  38  39  42  #  3  O1  C  C1  H
106  41  47  38  39  41  #  3  O1  C  C1  H
107  49  46  38  39  40  #  3  O1  C  C1  H
108  44  45  43  44  38  #  3  O1  C  C1  H
109  44  45  43  44  38  #  3  O1  C  C1  H
110  45  46  38  39  41  #  3  O1  C  C1  H
111  48  42  39  40  43  #  3  O1  C  C1  H
112  42  47  38  46  48  #  3  O1  C  C1  H
113  48  41  39  40  43  #  3  O1  C  C1  H
114  45  46  38  39  42  #  3  O1  C  C1  H

"""
        self.assertEqual(ans,dihedrals)

    def test__set_imdihedrals(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        data_lammps._set_atom(self.ffmol, self.mols_in_box)
        imdihedrals = data_lammps._set_imdihedrals(self.ffmol, self.mols_in_box)

        ans="""Impropers

1  5  34  22  16  21  #  2  C6  N1  C  N
2  3  19  18  17  21  #  2  C2  H  C1  N
3  2  17  20  19  22  #  2  C1  H1  C2  N1
4  4  23  17  21  16  #  2  C3  C1  N  C
5  1  27  16  22  19  #  2  C4  C  N1  C2
6  6  45  44  43  40  #  3  O2  O1  C2  O

"""
        self.assertEqual(ans,imdihedrals)

    def test_get_lammps_string(self):
        data_lammps=LmpInput(self.ffmol,self.mols_in_box)
        data_lammps._set_atom(self.ffmol, self.mols_in_box)
        lammps_data = data_lammps.get_lammps_string()


        ans="""LAMMPS Data File

# 1 mol 1 molecule
# 1 mol 2 molecule
# 1 mol 3 molecule


14 atom type
23 bond type
41 angle type
56 dihedral type
6 improper type

50 atoms
49 bonds
87 angles
114 dihedrals
6 impropers

0.0 40.0 xlo  xhi
0.0 40.0 ylo  yhi
0.0 40.0 zlo  zhi

Masses

1 32.06 # sy 1
2 19.0 # f 1
3 14.01 # ne 1
4 16.0 # o 1
5 32.06 # s6 1
6 12.01 # c3 1
7 12.01 # cc 2
8 14.01 # na 2
9 1.008 # h1 2
10 1.008 # h4 2
11 12.01 # cd 2
12 1.008 # hc 2
13 12.01 # c 3
14 16.0 # os 3


Pair Coeffs

1 0.25 3.5636 # sy 1
2 0.061 3.11815 # f 1
3 0.17 3.2500032 # ne 1
4 0.21 2.95992616 # o 1
5 0.25 3.5636 # s6 1
6 0.1094 3.3996744 # c3 1
7 0.086 3.3996744 # cc 2
8 0.17 3.2500032 # na 2
9 0.0157 2.4713566 # h1 2
10 0.015 2.5105562 # h4 2
11 0.086 3.3996744 # cd 2
12 0.0157 2.6495366 # hc 2
13 0.086 3.3996744 # c 3
14 0.17 3.00001666 # os 3


Bond Coeffs

1 257.1  1.752 # ne sy 1
2 444.6  1.551 # ne s6 1
3 363.8  1.344 # c3 f 1
4 541.1  1.436 # o s6 1
5 254.0  1.774 # c3 s6 1
6 248.9  1.782 # c3 sy 1
7 493.0  1.466 # o sy 1
8 334.7  1.456 # c3 na 2
9 350.1  1.083 # cc h4 2
10 337.3  1.499 # c3 cc 2
11 438.8  1.371 # cc na 2
12 303.1  1.535 # c3 c3 2
13 337.3  1.092 # c3 hc 2
14 438.8  1.371 # cd na 2
15 350.1  1.083 # cd h4 2
16 335.9  1.093 # c3 h1 2
17 504.0  1.371 # cc cd 2
18 411.3  1.343 # c os 3
19 301.5  1.439 # c3 os 3
20 337.3  1.092 # c3 hc 3
21 303.1  1.535 # c3 c3 3
22 335.9  1.093 # c3 h1 3
23 648.0  1.214 # c o 3


Angle Coeffs

1 43.1 107.06 # ne sy o 1
2 81.22 109.67 # f c3 s6 1
3 46.66 119.73 # o s6 o 1
4 9.72 119.18 # s6 ne sy 1
5 39.68 112.95 # c3 s6 ne 1
6 41.28 108.48 # c3 sy o 1
7 81.22 109.67 # f c3 sy 1
8 41.66 108.32 # c3 s6 o 1
9 39.42 103.12 # c3 sy ne 1
10 71.26 107.16 # f c3 f 1
11 45.3 121.88 # o sy o 1
12 45.1 116.41 # ne s6 o 1
13 49.9 109.45 # h1 c3 na 2
14 62.56 125.09 # c3 na cc 2
15 72.91 109.42 # cc cd na 2
16 65.44 122.99 # c3 cc na 2
17 72.91 109.42 # cd cc na 2
18 46.37 110.05 # c3 c3 hc 2
19 73.65 109.33 # na cc na 2
20 50.22 119.66 # h4 cc na 2
21 47.19 129.11 # cc cd h4 2
22 39.18 109.55 # h1 c3 h1 2
23 47.19 129.11 # cd cc h4 2
24 46.36 110.07 # c3 c3 h1 2
25 68.94 109.9 # cc na cc 2
26 63.88 128.01 # cc na cd 2
27 47.2 110.86 # cc c3 hc 2
28 39.43 108.35 # hc c3 hc 2
29 65.73 112.81 # c3 c3 na 2
30 62.56 125.09 # c3 na cd 2
31 50.22 119.66 # h4 cd na 2
32 63.63 115.14 # c os c3 3
33 46.37 110.05 # c3 c3 hc 3
34 46.36 110.07 # c3 c3 h1 3
35 50.84 108.82 # h1 c3 os 3
36 63.21 110.63 # c3 c3 c3 3
37 67.78 108.42 # c3 c3 os 3
38 39.43 108.35 # hc c3 hc 3
39 75.93 123.33 # o c os 3
40 76.45 111.38 # os c os 3
41 39.18 109.55 # h1 c3 h1 3


Dihedral Coeffs

1  6.667  2  180.0 # o s6 ne sy 1
2  0.144  3  0.0 # f c3 s6 o 1
3  6.8  1  180.0 # o sy ne s6 1
4  0.5  3  180.0 # o sy ne s6 1
5  6.667  2  180.0 # c3 s6 ne sy 1
6  6.8  1  180.0 # c3 sy ne s6 1
7  0.5  3  180.0 # c3 sy ne s6 1
8  0.144  3  0.0 # f c3 s6 ne 1
9  0.144  3  0.0 # f c3 sy o 1
10  0.144  3  0.0 # f c3 sy ne 1
11  4.0  2  180.0 # h4 cc cd na 2
12  4.0  2  180.0 # h4 cc cd h4 2
13  4.0  2  180.0 # na cc cd na 2
14  0.156  3  0.0 # h1 c3 c3 hc 2
15  0.0  2  0.0 # cd na c3 h1 2
16  1.7  2  180.0 # cd na cc na 2
17  0.0  3  0.0 # hc c3 cc na 2
18  4.0  2  180.0 # h4 cd cc h4 2
19  4.0  2  180.0 # na cd cc na 2
20  1.7  2  180.0 # c3 na cc c3 2
21  1.7  2  180.0 # c3 na cc h4 2
22  1.7  2  180.0 # c3 na cd h4 2
23  1.7  2  180.0 # cc na cc h4 2
24  1.7  2  180.0 # c3 na cd cc 2
25  1.7  2  180.0 # cc na cd h4 2
26  1.7  2  180.0 # c3 cc na c3 2
27  1.7  2  180.0 # cc cd na cc 2
28  1.7  2  180.0 # c3 na cc cd 2
29  1.7  2  180.0 # cc na cd cc 2
30  4.0  2  180.0 # h4 cd cc na 2
31  1.7  2  180.0 # c3 na cc na 2
32  1.7  2  180.0 # cc na cc cd 2
33  0.156  3  0.0 # hc c3 c3 na 2
34  0.0  2  0.0 # c3 c3 na cc 2
35  0.0  2  0.0 # c3 c3 na cd 2
36  1.7  2  180.0 # c3 cc na cd 2
37  1.7  2  180.0 # cc na cc na 2
38  1.7  2  180.0 # c3 cc na cc 2
39  0.0  2  0.0 # cc na c3 h1 2
40  0.25  1  0.0 # hc c3 c3 os 3
41  0.0  3  0.0 # hc c3 c3 os 3
42  2.7  2  180.0 # c3 os c os 3
43  0.25  1  0.0 # h1 c3 c3 os 3
44  0.0  3  0.0 # h1 c3 c3 os 3
45  0.156  3  0.0 # h1 c3 c3 h1 3
46  0.156  3  0.0 # h1 c3 c3 hc 3
47  0.8  1  180.0 # c os c3 c3 3
48  0.383  3  0.0 # c os c3 c3 3
49  1.4  1  180.0 # c3 os c o 3
50  2.7  2  180.0 # c3 os c o 3
51  0.156  3  0.0 # c3 c3 c3 h1 3
52  0.16  3  0.0 # c3 c3 c3 hc 3
53  0.144  3  0.0 # os c3 c3 os 3
54  1.175  2  0.0 # os c3 c3 os 3
55  0.383  3  0.0 # c os c3 h1 3
56  0.156  3  0.0 # c3 c3 c3 os 3


Imp Dihedral Coeffs

1 1.1  -1.0  2.0 # c3 cc na cd 2
2 1.1  -1.0  2.0 # cc h4 cd na 2
3 1.1  -1.0  2.0 # cd h4 cc na 2
4 1.1  -1.0  2.0 # c3 cc na cc 2
5 1.1  -1.0  2.0 # c3 na cc na 2
6 10.5  -1.0  2.0 # o os c os 3


Atoms

1  1  4  3.841  38.685  36.79 # 1  o O
2  1  2  1.427  38.189  39.692 # 1  f F
3  1  6  2.269  38.874  38.907 # 1  c3 C
4  1  5  2.826  37.843  37.418 # 1  s6 S
5  1  4  1.59  37.562  36.691 # 1  o O1
6  1  2  1.631  39.979  38.478 # 1  f F1
7  1  3  3.568  36.597  38.123 # 1  ne N
8  1  4  3.699  34.562  39.552 # 1  o O2
9  1  2  3.313  39.264  39.652 # 1  f F2
10  1  1  2.911  35.172  38.485 # 1  sy S1
11  1  4  1.452  35.104  38.556 # 1  o O3
12  1  2  2.877  32.876  37.123 # 1  f F3
13  1  2  4.682  34.042  36.8 # 1  f F4
14  1  6  3.354  34.126  36.968 # 1  c3 C1
15  1  2  2.824  34.635  35.847 # 1  f F5
16  2  7  2.604  5.594  7.056 # 2  cc C
17  2  7  3.435  7.464  6.181 # 2  cc C1
18  2  10  3.46  8.506  5.909 # 2  h4 H
19  2  11  4.333  6.471  5.966 # 2  cd C2
20  2  10  5.291  6.484  5.476 # 2  h4 H1
21  2  8  2.364  6.905  6.855 # 2  na N
22  2  8  3.804  5.315  6.511 # 2  na N1
23  2  6  1.16  7.634  7.279 # 2  c3 C3
24  2  9  1.236  8.655  6.912 # 2  h1 H2
25  2  9  0.27  7.17  6.853 # 2  h1 H3
26  2  9  1.092  7.651  8.367 # 2  h1 H4
27  2  6  4.495  3.996  6.514 # 2  c3 C4
28  2  9  3.72  3.23  6.511 # 2  h1 H5
29  2  9  5.022  3.93  5.568 # 2  h1 H6
30  2  6  5.441  3.843  7.701 # 2  c3 C5
31  2  12  5.927  2.868  7.641 # 2  hc H7
32  2  12  6.218  4.608  7.692 # 2  hc H8
33  2  12  4.911  3.898  8.655 # 2  hc H9
34  2  6  1.684  4.635  7.722 # 2  c3 C6
35  2  12  2.241  3.862  8.253 # 2  hc H10
36  2  12  1.051  5.148  8.445 # 2  hc H11
37  2  12  1.035  4.144  6.991 # 2  hc H12
38  3  6  1.008  1.371  3.592 # 3  c3 C
39  3  6  1.397  2.451  4.605 # 3  c3 C1
40  3  14  0.77  2.009  5.819 # 3  os O
41  3  9  1.017  3.442  4.359 # 3  h1 H
42  3  9  2.479  2.496  4.765 # 3  h1 H1
43  3  13  0.49  0.679  5.718 # 3  c C2
44  3  14  0.726  0.243  4.455 # 3  os O1
45  3  4  0.095  0.001  6.613 # 3  o O2
46  3  6  2.071  1.008  2.578 # 3  c3 C3
47  3  9  0.07  1.633  3.091 # 3  h1 H2
48  3  12  2.992  0.697  3.078 # 3  hc H3
49  3  12  2.292  1.87  1.943 # 3  hc H4
50  3  12  1.731  0.191  1.938 # 3  hc H5


Bonds

1 4 1 4 # 1 O S
2 3 2 3 # 1 F C
3 5 3 4 # 1 C S
4 3 3 6 # 1 C F1
5 3 3 9 # 1 C F2
6 4 4 5 # 1 S O1
7 2 4 7 # 1 S N
8 1 7 10 # 1 N S1
9 7 8 10 # 1 O2 S1
10 7 10 11 # 1 S1 O3
11 6 10 14 # 1 S1 C1
12 3 12 14 # 1 F3 C1
13 3 13 14 # 1 F4 C1
14 3 14 15 # 1 C1 F5
15 11 16 21 # 2 C N
16 11 16 22 # 2 C N1
17 10 16 34 # 2 C C6
18 9 17 18 # 2 C1 H
19 17 17 19 # 2 C1 C2
20 11 17 21 # 2 C1 N
21 15 19 20 # 2 C2 H1
22 14 19 22 # 2 C2 N1
23 8 21 23 # 2 N C3
24 8 22 27 # 2 N1 C4
25 16 23 24 # 2 C3 H2
26 16 23 25 # 2 C3 H3
27 16 23 26 # 2 C3 H4
28 16 27 28 # 2 C4 H5
29 16 27 29 # 2 C4 H6
30 12 27 30 # 2 C4 C5
31 13 30 31 # 2 C5 H7
32 13 30 32 # 2 C5 H8
33 13 30 33 # 2 C5 H9
34 13 34 35 # 2 C6 H10
35 13 34 36 # 2 C6 H11
36 13 34 37 # 2 C6 H12
37 21 38 39 # 3 C C1
38 19 38 44 # 3 C O1
39 21 38 46 # 3 C C3
40 22 38 47 # 3 C H2
41 19 39 40 # 3 C1 O
42 22 39 41 # 3 C1 H
43 22 39 42 # 3 C1 H1
44 18 40 43 # 3 O C2
45 18 43 44 # 3 C2 O1
46 23 43 45 # 3 C2 O2
47 20 46 48 # 3 C3 H3
48 20 46 49 # 3 C3 H4
49 20 46 50 # 3 C3 H5


Angles

1  8  1  4  3  #  1  O  S  C
2  3  1  4  5  #  1  O  S  O1
3  12  1  4  7  #  1  O  S  N
4  2  2  3  4  #  1  F  C  S
5  10  2  3  6  #  1  F  C  F1
6  10  2  3  9  #  1  F  C  F2
7  8  3  4  5  #  1  C  S  O1
8  5  3  4  7  #  1  C  S  N
9  2  4  3  6  #  1  S  C  F1
10  2  4  3  9  #  1  S  C  F2
11  4  4  7  10  #  1  S  N  S1
12  12  5  4  7  #  1  O1  S  N
13  10  6  3  9  #  1  F1  C  F2
14  1  7  10  8  #  1  N  S1  O2
15  1  7  10  11  #  1  N  S1  O3
16  9  7  10  14  #  1  N  S1  C1
17  11  8  10  11  #  1  O2  S1  O3
18  6  8  10  14  #  1  O2  S1  C1
19  7  10  14  12  #  1  S1  C1  F3
20  7  10  14  13  #  1  S1  C1  F4
21  7  10  14  15  #  1  S1  C1  F5
22  6  11  10  14  #  1  O3  S1  C1
23  10  12  14  13  #  1  F3  C1  F4
24  10  12  14  15  #  1  F3  C1  F5
25  10  13  14  15  #  1  F4  C1  F5
26  25  16  21  17  #  2  C  N  C1
27  14  16  21  23  #  2  C  N  C3
28  26  16  22  19  #  2  C  N1  C2
29  14  16  22  27  #  2  C  N1  C4
30  27  16  34  35  #  2  C  C6  H10
31  27  16  34  36  #  2  C  C6  H11
32  27  16  34  37  #  2  C  C6  H12
33  21  17  19  20  #  2  C1  C2  H1
34  15  17  19  22  #  2  C1  C2  N1
35  14  17  21  23  #  2  C1  N  C3
36  23  18  17  19  #  2  H  C1  C2
37  20  18  17  21  #  2  H  C1  N
38  17  19  17  21  #  2  C2  C1  N
39  30  19  22  27  #  2  C2  N1  C4
40  31  20  19  22  #  2  H1  C2  N1
41  19  21  16  22  #  2  N  C  N1
42  16  21  16  34  #  2  N  C  C6
43  13  21  23  24  #  2  N  C3  H2
44  13  21  23  25  #  2  N  C3  H3
45  13  21  23  26  #  2  N  C3  H4
46  16  22  16  34  #  2  N1  C  C6
47  13  22  27  28  #  2  N1  C4  H5
48  13  22  27  29  #  2  N1  C4  H6
49  29  22  27  30  #  2  N1  C4  C5
50  22  24  23  25  #  2  H2  C3  H3
51  22  24  23  26  #  2  H2  C3  H4
52  22  25  23  26  #  2  H3  C3  H4
53  18  27  30  31  #  2  C4  C5  H7
54  18  27  30  32  #  2  C4  C5  H8
55  18  27  30  33  #  2  C4  C5  H9
56  22  28  27  29  #  2  H5  C4  H6
57  24  28  27  30  #  2  H5  C4  C5
58  24  29  27  30  #  2  H6  C4  C5
59  28  31  30  32  #  2  H7  C5  H8
60  28  31  30  33  #  2  H7  C5  H9
61  28  32  30  33  #  2  H8  C5  H9
62  28  35  34  36  #  2  H10  C6  H11
63  28  35  34  37  #  2  H10  C6  H12
64  28  36  34  37  #  2  H11  C6  H12
65  37  38  39  40  #  3  C  C1  O
66  34  38  39  41  #  3  C  C1  H
67  34  38  39  42  #  3  C  C1  H1
68  32  38  44  43  #  3  C  O1  C2
69  33  38  46  48  #  3  C  C3  H3
70  33  38  46  49  #  3  C  C3  H4
71  33  38  46  50  #  3  C  C3  H5
72  37  39  38  44  #  3  C1  C  O1
73  36  39  38  46  #  3  C1  C  C3
74  34  39  38  47  #  3  C1  C  H2
75  32  39  40  43  #  3  C1  O  C2
76  35  40  39  41  #  3  O  C1  H
77  35  40  39  42  #  3  O  C1  H1
78  40  40  43  44  #  3  O  C2  O1
79  39  40  43  45  #  3  O  C2  O2
80  41  41  39  42  #  3  H  C1  H1
81  37  44  38  46  #  3  O1  C  C3
82  35  44  38  47  #  3  O1  C  H2
83  39  44  43  45  #  3  O1  C2  O2
84  34  46  38  47  #  3  C3  C  H2
85  38  48  46  49  #  3  H3  C3  H4
86  38  48  46  50  #  3  H3  C3  H5
87  38  49  46  50  #  3  H4  C3  H5


Dihedrals

1  2  2  3  4  5  #  1  F1  C  S  O
2  6  9  3  4  7  #  1  F1  C  S  O
3  2  9  3  4  1  #  1  F1  C  S  O
4  5  4  7  10  14  #  1  F1  C  S  O
5  5  4  7  10  14  #  1  F1  C  S  O
6  2  6  3  4  1  #  1  F1  C  S  O
7  6  6  3  4  7  #  1  F1  C  S  O
8  4  3  4  7  10  #  1  F1  C  S  O
9  2  6  3  4  5  #  1  F1  C  S  O
10  1  5  4  7  10  #  1  F1  C  S  O
11  8  7  10  14  12  #  1  F1  C  S  O
12  8  7  10  14  13  #  1  F1  C  S  O
13  8  7  10  14  15  #  1  F1  C  S  O
14  7  8  10  14  13  #  1  F1  C  S  O
15  7  8  10  14  15  #  1  F1  C  S  O
16  7  8  10  14  12  #  1  F1  C  S  O
17  1  1  4  7  10  #  1  F1  C  S  O
18  2  2  3  4  1  #  1  F1  C  S  O
19  6  2  3  4  7  #  1  F1  C  S  O
20  2  9  3  4  5  #  1  F1  C  S  O
21  7  11  10  14  12  #  1  F1  C  S  O
22  7  11  10  14  15  #  1  F1  C  S  O
23  7  11  10  14  13  #  1  F1  C  S  O
24  3  4  7  10  11  #  1  F1  C  S  O
25  3  4  7  10  11  #  1  F1  C  S  O
26  3  4  7  10  8  #  1  F1  C  S  O
27  3  4  7  10  8  #  1  F1  C  S  O
28  37  17  21  23  26  #  2  N1  C  N  C3
29  37  17  21  23  24  #  2  N1  C  N  C3
30  37  17  21  23  25  #  2  N1  C  N  C3
31  11  21  17  19  22  #  2  N1  C  N  C3
32  32  16  22  27  30  #  2  N1  C  N  C3
33  28  21  17  19  20  #  2  N1  C  N  C3
34  24  34  16  22  27  #  2  N1  C  N  C3
35  23  20  19  22  16  #  2  N1  C  N  C3
36  37  16  21  23  24  #  2  N1  C  N  C3
37  37  16  21  23  25  #  2  N1  C  N  C3
38  30  19  17  21  16  #  2  N1  C  N  C3
39  37  16  21  23  26  #  2  N1  C  N  C3
40  29  21  16  22  27  #  2  N1  C  N  C3
41  13  19  22  27  29  #  2  N1  C  N  C3
42  13  19  22  27  28  #  2  N1  C  N  C3
43  21  18  17  21  16  #  2  N1  C  N  C3
44  37  16  22  27  28  #  2  N1  C  N  C3
45  15  22  16  34  37  #  2  N1  C  N  C3
46  15  22  16  34  36  #  2  N1  C  N  C3
47  15  22  16  34  35  #  2  N1  C  N  C3
48  20  20  19  22  27  #  2  N1  C  N  C3
49  12  28  27  30  33  #  2  N1  C  N  C3
50  12  28  27  30  32  #  2  N1  C  N  C3
51  12  28  27  30  31  #  2  N1  C  N  C3
52  15  21  16  34  36  #  2  N1  C  N  C3
53  14  21  16  22  19  #  2  N1  C  N  C3
54  34  34  16  22  19  #  2  N1  C  N  C3
55  15  21  16  34  35  #  2  N1  C  N  C3
56  15  21  16  34  37  #  2  N1  C  N  C3
57  37  16  22  27  29  #  2  N1  C  N  C3
58  25  17  19  22  16  #  2  N1  C  N  C3
59  10  18  17  19  20  #  2  N1  C  N  C3
60  9  18  17  19  22  #  2  N1  C  N  C3
61  29  22  16  21  23  #  2  N1  C  N  C3
62  35  22  16  21  17  #  2  N1  C  N  C3
63  24  34  16  21  23  #  2  N1  C  N  C3
64  36  34  16  21  17  #  2  N1  C  N  C3
65  31  22  27  30  31  #  2  N1  C  N  C3
66  22  17  19  22  27  #  2  N1  C  N  C3
67  31  22  27  30  33  #  2  N1  C  N  C3
68  31  22  27  30  32  #  2  N1  C  N  C3
69  26  19  17  21  23  #  2  N1  C  N  C3
70  12  29  27  30  32  #  2  N1  C  N  C3
71  12  29  27  30  33  #  2  N1  C  N  C3
72  12  29  27  30  31  #  2  N1  C  N  C3
73  19  18  17  21  23  #  2  N1  C  N  C3
74  33  19  22  27  30  #  2  N1  C  N  C3
75  48  47  38  44  43  #  3  O1  C  C1  H
76  38  44  38  46  50  #  3  O1  C  C1  H
77  38  44  38  46  50  #  3  O1  C  C1  H
78  38  44  38  46  49  #  3  O1  C  C1  H
79  38  44  38  46  49  #  3  O1  C  C1  H
80  38  44  38  46  48  #  3  O1  C  C1  H
81  38  44  38  46  48  #  3  O1  C  C1  H
82  43  46  38  44  43  #  3  O1  C  C1  H
83  43  46  38  44  43  #  3  O1  C  C1  H
84  39  39  40  43  44  #  3  O1  C  C1  H
85  46  39  38  46  48  #  3  O1  C  C1  H
86  42  47  38  46  49  #  3  O1  C  C1  H
87  44  39  40  43  45  #  3  O1  C  C1  H
88  44  39  40  43  45  #  3  O1  C  C1  H
89  46  39  38  46  50  #  3  O1  C  C1  H
90  46  39  38  46  49  #  3  O1  C  C1  H
91  43  38  39  40  43  #  3  O1  C  C1  H
92  43  38  39  40  43  #  3  O1  C  C1  H
93  41  47  38  39  42  #  3  O1  C  C1  H
94  43  39  38  44  43  #  3  O1  C  C1  H
95  43  39  38  44  43  #  3  O1  C  C1  H
96  47  44  38  39  40  #  3  O1  C  C1  H
97  47  44  38  39  40  #  3  O1  C  C1  H
98  40  44  38  39  41  #  3  O1  C  C1  H
99  40  44  38  39  41  #  3  O1  C  C1  H
100  42  47  38  46  50  #  3  O1  C  C1  H
101  39  40  43  44  38  #  3  O1  C  C1  H
102  40  47  38  39  40  #  3  O1  C  C1  H
103  40  47  38  39  40  #  3  O1  C  C1  H
104  40  44  38  39  42  #  3  O1  C  C1  H
105  40  44  38  39  42  #  3  O1  C  C1  H
106  41  47  38  39  41  #  3  O1  C  C1  H
107  49  46  38  39  40  #  3  O1  C  C1  H
108  44  45  43  44  38  #  3  O1  C  C1  H
109  44  45  43  44  38  #  3  O1  C  C1  H
110  45  46  38  39  41  #  3  O1  C  C1  H
111  48  42  39  40  43  #  3  O1  C  C1  H
112  42  47  38  46  48  #  3  O1  C  C1  H
113  48  41  39  40  43  #  3  O1  C  C1  H
114  45  46  38  39  42  #  3  O1  C  C1  H


Impropers

1  5  34  22  16  21  #  2  C6  N1  C  N
2  3  19  18  17  21  #  2  C2  H  C1  N
3  2  17  20  19  22  #  2  C1  H1  C2  N1
4  4  23  17  21  16  #  2  C3  C1  N  C
5  1  27  16  22  19  #  2  C4  C  N1  C2
6  6  45  44  43  40  #  3  O2  O1  C2  O

"""
        self.assertEqual(ans,lammps_data)





if __name__ == '__main__':
    unittest.main()
