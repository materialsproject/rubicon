from unittest import TestCase
from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator

__author__ = 'xiaohuiqu'



class TestQChemFireWorkCreator(TestCase):
    def test_get_dielectric_constant(self):
        solvent = "ethylenecarbonate"
        dielectric, probe_radius = QChemFireWorkCreator.get_dielectric_constant(solvent)
        self.assertEqual(dielectric, 90.0)
        self.assertEqual(probe_radius, 2.38)