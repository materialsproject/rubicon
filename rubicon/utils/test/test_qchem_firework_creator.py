from unittest import TestCase
from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator

__author__ = 'xiaohuiqu'


class TestQChemFireWorkCreator(TestCase):
    def test_get_dielectric_constant(self):
        solvent = "ethylenecarbonate"
        dielectric, probe_radius, solvent_name = QChemFireWorkCreator.get_dielectric_constant(solvent)
        self.assertEqual(dielectric, 90.0)
        self.assertEqual(probe_radius, 2.38)
        self.assertEqual(solvent_name, "ethylenecarbonate")

        solvent = {"components": {"ethylenecarbonate": 3.0,
                                  "ethylmethylcarbonate": 7.0},
                   "metrics": "volume"}
        dielectric, probe_radius, solvent_name = QChemFireWorkCreator.get_dielectric_constant(solvent)
        self.assertAlmostEqual(dielectric, 28.998908694815448, 3)
        self.assertAlmostEqual(probe_radius, 2.6208686663466936, 3)
        self.assertEqual(solvent_name, "30.00% ethylenecarbonate, 70.00% ethylmethylcarbonate in volume")
