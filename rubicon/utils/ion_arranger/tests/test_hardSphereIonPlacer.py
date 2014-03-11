import os
from unittest import TestCase
import unittest
import openbabel as ob
from rubicon.utils.ion_arranger.ion_arranger import HardSphereIonPlacer

__author__ = 'xiaohuiqu'


test_dir = os.path.join(os.path.dirname(__file__),
                        "../../../../test_files")


class TestHardSphereIonPlacer(TestCase):
    def test_normalize_molecule(self):
        obconv = ob.OBConversion()
        obconv.SetInFormat("xyz")
        mol_file = os.path.join(test_dir, "2AcetoxyQ.xyz")
        obmol = ob.OBMol()
        obconv.ReadFile(obmol, mol_file)
        coords, radius, elements = HardSphereIonPlacer.normalize_molecule(
            obmol, 2.0)
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
        ref_radius = [1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 0.62, 0.62, 0.62,
                      0.62, 1.46, 1.46, 1.42, 1.42, 0.62, 1.46, 1.32, 1.46,
                      0.62, 0.62, 0.62]
        for r1, r2 in zip(radius, ref_radius):
            self.assertAlmostEqual(r1, r2)
        print elements




if __name__ == '__main__':
    unittest.main()