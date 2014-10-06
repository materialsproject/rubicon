import os
from unittest import TestCase
from pymatgen.io.qchemio import QcOutput
from rubicon.utils.atomic_charge_mixed_basis_set_generator import AtomicChargeMixedBasisSetGenerator

__author__ = 'xiaohuiqu'

test_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                           "..", "..", "..", "test_files"))


class TestAtomicChargeMixedBasisSetGenerator(TestCase):

    def setUp(self):
        self.generator = AtomicChargeMixedBasisSetGenerator(charge_threshold=-0.5)

    def test_get_basis(self):
        filename = os.path.join(test_dir, "quinoxaline_anion.qcout")
        qcout = QcOutput(filename)
        charges = qcout.data[0]["charges"]["nbo"]
        mol = qcout.data[0]["molecules"][-1]
        basis = self.generator.get_basis(mol, charges)
        ans = [('C', '6-31G*'), ('C', '6-31G*'), ('C', '6-31G*'), ('C', '6-31G*'),
               ('C', '6-31G*'), ('C', '6-31G*'), ('H', '6-31G*'), ('H', '6-31G*'),
               ('H', '6-31G*'), ('H', '6-31G*'), ('C', '6-31G*'), ('C', '6-31G*'),
               ('H', '6-31G*'), ('N', '6-31+G*'), ('N', '6-31+G*'), ('H', '6-31G*')]
        self.maxDiff = None
        self.assertEqual(basis, [(b[0], b[1].lower()) for b in ans])

    def test_to_and_from_dict(self):
        g1 = AtomicChargeMixedBasisSetGenerator(
            charge_threshold=-0.2, normal_basis_set="cc-PVDZ", diffuse_basis_set="aug-cc-PVDZ")
        g2 = AtomicChargeMixedBasisSetGenerator.from_dict(g1.to_dict())
        self.assertEqual(g2.charge_threshold, -0.2)
        self.assertEqual(g2.normal_basis_set, "cc-PVDZ".lower())
        self.assertEqual(g2.diffuse_basis_set, "aug-cc-PVDZ".lower())
