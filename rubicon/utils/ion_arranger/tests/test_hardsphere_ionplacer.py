import os
from unittest import TestCase
import math

import openbabel as ob
from pymatgen.core.structure import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.qchemio import QcOutput

from rubicon.utils.ion_arranger.ion_arranger import IonPlacer


__author__ = 'xiaohuiqu'


test_dir = os.path.join(os.path.dirname(__file__),
                        "../../../../test_files")


class TestHardSphereIonPlacer(TestCase):
    def setUp(self):
        obconv = ob.OBConversion()
        obconv.SetInFormat("xyz")
        mol_file = os.path.join(test_dir, "2AcetoxyQ.xyz")
        obmol = ob.OBMol()
        obconv.ReadFile(obmol, mol_file)
        self.obmol_template = obmol
        pymagen_sodium = Molecule(species=['Na'],
                                  coords=[[999.0, 999.0, 999.0]],
                                  charge=1)
        # noinspection PyProtectedMember
        sodium_obmol = BabelMolAdaptor(pymagen_sodium)._obmol
        pymagen_chloride = Molecule(species=['Cl'],
                                    coords=[[-200.0, -200.0, -200.0]],
                                    charge=-1)
        # noinspection PyProtectedMember
        chloride_obmol = BabelMolAdaptor(pymagen_chloride)._obmol
        acetoxyq_anion_qcout = QcOutput(
            os.path.join(test_dir, "acetoxyq_anion.out"))
        pymatgen_acetoxyq = acetoxyq_anion_qcout.data[0]["molecules"][-1]
        # noinspection PyProtectedMember
        acetoxyq_obmol = BabelMolAdaptor(pymatgen_acetoxyq)._obmol
        acetoxyq_charges = acetoxyq_anion_qcout.data[0]["charges"]['chelpg']
        self.acetoxyq_nacl_placer = IonPlacer(
            acetoxyq_obmol, acetoxyq_charges, sodium_obmol, [1.0],
            chloride_obmol, [-1.0])
        tfsi_qcout = QcOutput(os.path.join(test_dir, "tfsi.qcout"))
        pymatgen_tfsi = tfsi_qcout.data[0]["molecules"][-1]
        # noinspection PyProtectedMember
        tfsi_obmol = BabelMolAdaptor(pymatgen_tfsi)._obmol
        tfsi_charges = tfsi_qcout.data[0]["charges"]['chelpg']
        self.acetoxyq_natfsi_placer = IonPlacer(
            acetoxyq_obmol, acetoxyq_charges, sodium_obmol, [1.0],
            tfsi_obmol, tfsi_charges)

    def get_copy_of_mol(self):
        obmol = ob.OBMol(self.obmol_template)
        return obmol

    def test_normalize_molecule(self):
        obmol = self.get_copy_of_mol()
        coords, radius, elements = IonPlacer.normalize_molecule(
            obmol, 2.0, 0.9)
        x, y, z = 0.0, 0.0, 0.0
        for c in coords:
            x += c[0]
            y += c[1]
            z += c[2]
        self.assertAlmostEqual(x, 0.0, 5)
        self.assertAlmostEqual(y, 0.0, 5)
        self.assertAlmostEqual(z, 0.0, 5)
        ref_coords = [[-6.53199, -0.59843, 5e-05], [-4.81611, -2.56038, 1e-05],
                      [-2.1887, -2.03609, 0.0], [-1.34938, 0.53786, 6e-05],
                      [-3.15885, 2.52026, 0.0001], [-5.70081, 1.9539, 9e-05],
                      [-8.54594, -1.0004, 6e-05], [-5.41052, -4.52467, -2e-05],
                      [-2.4763, 4.45601, 0.00013],
                      [-7.08962, 3.46651, 0.00012],
                      [2.75804, -0.81316, -2e-05], [1.91223, -3.37868, -6e-05],
                      [1.15959, 1.10448, 8e-05], [-0.4971, -3.98645, -5e-05],
                      [3.29338, -4.89733, -5e-05], [5.55185, -0.24909, 0.0001],
                      [7.06607, -1.99415, 0.00029],
                      [6.36898, 2.48666, -0.00013],
                      [5.61277, 3.47056, -1.65816],
                      [5.61191, 3.47122, 1.65707],
                      [8.43049, 2.57137, 0.00032]]
        for c1, c2, in zip(coords, ref_coords):
            for x1, x2 in zip(c1, c2):
                self.assertAlmostEqual(x1, x2, 3)
        ref_radius = [2.759, 2.759, 2.759, 2.759, 2.759, 2.759, 1.17163,
                      1.17163, 1.17163, 1.17163, 2.759, 2.759, 2.68341,
                      2.68341, 1.17163, 2.759, 2.49444, 2.759, 1.17163,
                      1.17163, 1.17163]
        for r1, r2 in zip(radius, ref_radius):
            self.assertAlmostEqual(r1, r2, 3)
        ref_elements = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'C',
                        'C', 'N', 'N', 'H', 'C', 'O', 'C', 'H', 'H', 'H']
        self.assertEqual(elements, ref_elements)

    def test_rotate(self):
        obmol1 = self.get_copy_of_mol()
        IonPlacer.rotate(obmol1, math.pi/2, math.pi/4)
        coords1 = IonPlacer.get_mol_coords(obmol1)
        ref_coords1 = [[-4.85126, -2e-05, 4.66393], [-5.02525, 3e-05, 2.06331],
                       [-2.79667, 3e-05, 0.57618], [-0.38313, -3e-05, 1.80275],
                       [-0.26084, -6e-05, 4.484], [-2.45876, -5e-05, 5.88096],
                       [-6.55957, -2e-05, 5.80377], [-6.83453, 6e-05, 1.09466],
                       [1.59057, -9e-05, 5.37015], [-2.37122, -8e-05, 7.93257],
                       [1.56594, 5e-05, -2.05695],
                       [-0.84624, 0.0001, -3.27297],
                       [1.79164, -5e-05, 0.42929], [-2.97964, 9e-05, -1.99908],
                       [-0.94345, 9e-05, -5.32344],
                       [3.94032, -7e-05, -3.63362],
                       [3.7771, -0.00026, -5.93828],
                       [6.45259, 0.00016, -2.27695],
                       [6.61359, 1.65819, -1.0465],
                       [6.61345, -1.65704, -1.04543],
                       [7.9702, -0.00028, -3.67476]]
        for c1, c2, in zip(coords1, ref_coords1):
            for x1, x2 in zip(c1, c2):
                self.assertAlmostEqual(x1, x2, 3)
        obmol2 = self.get_copy_of_mol()
        IonPlacer.rotate(obmol2, math.pi/3, -math.pi * 0.6)
        coords2 = IonPlacer.get_mol_coords(obmol2)
        ref_coords2 = [[2.18824, -0.06625, -6.3635],
                       [3.27396, -1.04718, -4.20654],
                       [2.03023, -0.78504, -1.84804],
                       [-0.34917, 0.50189, -1.73864],
                       [-1.42281, 1.49306, -3.99007],
                       [-0.17082, 1.20988, -6.25605],
                       [3.14166, -0.26724, -8.17131],
                       [5.07553, -2.0293, -4.24618],
                       [-3.2281, 2.4609, -3.85898],
                       [-0.98752, 1.96617, -7.98169],
                       [-0.50564, -0.17355, 2.52932],
                       [1.86882, -1.45628, 2.41148],
                       [-1.59118, 0.78518, 0.49589],
                       [3.11391, -1.76017, 0.28273],
                       [2.69282, -2.2156, 4.13146],
                       [-1.83362, 0.10838, 5.03542],
                       [-0.86434, -0.76431, 6.94251],
                       [-4.33929, 1.47645, 5.08046],
                       [-4.12755, 3.4043, 4.35413],
                       [-5.70431, 0.53356, 3.8409],
                       [-5.04632, 1.51843, 7.01833]]
        for c1, c2, in zip(coords2, ref_coords2):
            for x1, x2 in zip(c1, c2):
                self.assertAlmostEqual(x1, x2, 3)
        obmol3 = self.get_copy_of_mol()
        obmol4 = self.get_copy_of_mol()
        obconv = ob.OBConversion()
        obconv.SetOutFormat("xyz")
        xyz3 = obconv.WriteString(obmol3)
        xyz4 = obconv.WriteString(obmol4)
        self.assertEqual(xyz3, xyz4)
        IonPlacer.rotate(obmol4, math.pi/3, math.pi/4)
        xyz3 = obconv.WriteString(obmol3)
        xyz4 = obconv.WriteString(obmol4)
        self.assertNotEqual(xyz3, xyz4)

    def test_pair_energy(self):
        coords1 = [[1.0, 0.0, 0.0]]
        charges1 = [1.0]
        radius1 = [0.4]
        coords2 = [[0.0, 0.0, 0.0]]
        charges2 = [-1.0]
        radius2 = [0.4]
        energy = IonPlacer.pair_energy(coords1, charges1, radius1,
                                                 coords2, charges2, radius2)
        self.assertEqual(energy, -1.0)
        radius1, radius2 = [0.8], [0.8]
        energy = IonPlacer.pair_energy(coords1, charges1, radius1,
                                                 coords2, charges2, radius2)
        self.assertGreater(energy, 1000.0)

    def test_set_bound(self):
        ref_nacl_upper_bound = [27.62783, 27.62783, 27.62783, 27.62783,
                                27.62783, 27.62783]
        ref_nacl_lower_bound = [-27.62783, -27.62783, -27.62783, -27.62783,
                                -27.62783, -27.62783]
        for x1, x2 in zip(self.acetoxyq_nacl_placer.bounder.upper_bound,
                          ref_nacl_upper_bound):
            self.assertAlmostEqual(x1, x2, 3)
        for x1, x2 in zip(self.acetoxyq_nacl_placer.bounder.lower_bound,
                          ref_nacl_lower_bound):
            self.assertAlmostEqual(x1, x2, 3)
        ref_natfsi_upper_bound = [41.62831, 41.62831, 41.62831, 41.62831,
                                  41.62831, 41.62831, 3.14159, 3.14159]
        ref_natfsi_lower_bound = [-41.62831, -41.62831, -41.62831, -41.62831,
                                  -41.62831, -41.62831, 0.0, -3.14159]
        for x1, x2 in zip(self.acetoxyq_natfsi_placer.bounder.upper_bound,
                          ref_natfsi_upper_bound):
            self.assertAlmostEqual(x1, x2, 3)
        for x1, x2 in zip(self.acetoxyq_natfsi_placer.bounder.lower_bound,
                          ref_natfsi_lower_bound):
            self.assertAlmostEqual(x1, x2, 3)

    def test_calc_energy(self):
        c = [0.0] * 6
        cation_coords, anion_coords = self.acetoxyq_nacl_placer\
            .decode_solution(c)
        energy = self.acetoxyq_nacl_placer.calc_energy(cation_coords,
                                                       anion_coords)
        self.assertGreater(energy, 1000.0)
        c = [20.0, 0.0, 0.0, 23.0, 0.0, 0.0]
        cation_coords, anion_coords = self.acetoxyq_nacl_placer\
            .decode_solution(c)
        energy = self.acetoxyq_nacl_placer.calc_energy(cation_coords,
                                                       anion_coords)
        self.assertAlmostEqual(energy, -0.12330098044441913, 5)
        c = [0.0] * 8
        cation_coords, anion_coords = self.acetoxyq_natfsi_placer.\
            decode_solution(c)
        energy = self.acetoxyq_natfsi_placer.calc_energy(cation_coords,
                                                         anion_coords)
        self.assertGreater(energy, 1000.0)
        c = [25.0, 1.0, 0.0, 20.0, 0.0, 0.0, 1.0, -2.0]
        cation_coords, anion_coords = self.acetoxyq_natfsi_placer.\
            decode_solution(c)
        energy = self.acetoxyq_natfsi_placer.calc_energy(cation_coords,
                                                         anion_coords)
        self.assertAlmostEqual(energy, -0.046013805864322584)

    def test_gravitational_energy(self):
        c = [20.0, 0.0, 0.0, 10.0, 0.0, 0.0]
        cation_coords, anion_coords = self.acetoxyq_nacl_placer\
            .decode_solution(c)
        cation_g_energy = self.acetoxyq_nacl_placer.gravitational_energy(
            cation_coords[0], self.acetoxyq_nacl_placer.cation_radius)
        anion_g_energy = self.acetoxyq_nacl_placer.gravitational_energy(
            anion_coords[0], self.acetoxyq_nacl_placer.anion_radius)
        self.assertAlmostEqual(cation_g_energy, 0.026525587104203413, 3)
        self.assertAlmostEqual(anion_g_energy, 0.005372882799983831, 3)
        c = [-20.0, 0.0, 0.0, 10.0, 0.0, 0.0, 1.0, 2.0]
        cation_coords, anion_coords = self.acetoxyq_natfsi_placer\
            .decode_solution(c)
        cation_g_energy = self.acetoxyq_natfsi_placer.gravitational_energy(
            cation_coords[0], self.acetoxyq_natfsi_placer.cation_radius)
        anion_g_energy = self.acetoxyq_natfsi_placer.gravitational_energy(
            anion_coords[0], self.acetoxyq_natfsi_placer.anion_radius)
        self.assertAlmostEqual(cation_g_energy, 0.026730174838683875, 3)
        self.assertAlmostEqual(anion_g_energy, 0.0031244355633821674, 3)