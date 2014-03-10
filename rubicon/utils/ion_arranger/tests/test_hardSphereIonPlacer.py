import os
from unittest import TestCase
import unittest
import openbabel as ob
from rubicon.utils.ion_arranger.ion_arranger import HardSphereIonPlacer

__author__ = 'xiaohuiqu'


test_dir = os.path.join(os.path.dirname(__file__),
                        "../../..//workflows/test_mols")


class TestHardSphereIonPlacer(TestCase):
    def test_normalize_molecule(self):
        obconv = ob.OBConversion()
        obconv.SetInFormat("xyz")
        mol_file = os.path.join(test_dir, "2AcetoxyQ.xyz")
        obmol = ob.OBMol()
        obconv.ReadFile(obmol, mol_file)
        coords, radius = HardSphereIonPlacer.normalize_molecule(obmol, 2.0)
        x, y, z = 0.0, 0.0, 0.0
        for c in coords:
            x += c[0]
            y += c[1]
            z += c[2]
        self.assertAlmostEqual(x, 0.0, 5)
        self.assertAlmostEqual(y, 0.0, 5)
        self.assertAlmostEqual(z, 0.0, 5)
        ref_coords = [[-3.45658, -0.31668, 3e-05], [-2.54857, -1.35489, 1e-05],
                      [-1.15821, -1.07745, 0.0], [-0.71406, 0.28462, 3e-05],
                      [-1.67159, 1.33366, 5e-05], [-3.01674, 1.03396, 5e-05],
                      [-4.52231, -0.52939, 3e-05], [-2.86312, -2.39435, -1e-05],
                      [-1.3104, 2.35802, 7e-05], [-3.75166, 1.8344, 6e-05],
                      [1.45949, -0.43031, -1e-05], [1.01191, -1.78792, -3e-05],
                      [0.61363, 0.58446, 4e-05], [-0.26305, -2.10954, -3e-05],
                      [1.74278, -2.59155, -3e-05], [2.93791, -0.13181, 5e-05],
                      [3.73921, -1.05526, 0.00015], [3.37032, 1.31588, -7e-05],
                      [2.97015, 1.83654, -0.87746], [2.9697, 1.83689, 0.87689],
                      [4.46122, 1.36071, 0.00017]]
        for c1, c2, in zip(coords, ref_coords):
            for x1, x2 in zip(c1, c2):
                self.assertAlmostEqual(x1, x2, 3)
        ref_radius = [1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 0.62, 0.62, 0.62,
                      0.62, 1.46, 1.46, 1.42, 1.42, 0.62, 1.46, 1.32, 1.46,
                      0.62, 0.62, 0.62]
        for r1, r2 in zip(radius, ref_radius):
            self.assertAlmostEqual(r1, r2)




if __name__ == '__main__':
    unittest.main()