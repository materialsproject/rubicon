import os
from unittest import TestCase
from pymatgen import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.qchemio import QcOutput
from rubicon.utils.ion_arranger.hard_sphere_energy_evaluators import AtomicRadiusUtils, LargestContactGapEnergyEvaluator
from rubicon.utils.ion_arranger.ion_arranger import IonPlacer

__author__ = 'xiaohuiqu'

test_dir = os.path.join(os.path.dirname(__file__),
                        "../../../../test_files")

class TestLargestContactGapEnergyEvaluator(TestCase):
    def setUp(self):
        pymagen_sodium = Molecule(species=['Na'],
                                  coords=[[999.0, 999.0, 999.0]],
                                  charge=1)
        # noinspection PyProtectedMember
        sodium_obmol = BabelMolAdaptor(pymagen_sodium)._obmol
        acetoxyq_anion_qcout = QcOutput(
            os.path.join(test_dir, "acetoxyq_anion.out"))
        pymatgen_acetoxyq = acetoxyq_anion_qcout.data[0]["molecules"][-1]
        acetoxyq_obmol = BabelMolAdaptor(pymatgen_acetoxyq)._obmol
        tfsi_qcout = QcOutput(os.path.join(test_dir, "tfsi.qcout"))
        pymatgen_tfsi = tfsi_qcout.data[0]["molecules"][-1]
        # noinspection PyProtectedMember
        tfsi_obmol = BabelMolAdaptor(pymatgen_tfsi)._obmol
        self.acetoxyq_natfsi_placer = IonPlacer(
            acetoxyq_obmol, sodium_obmol, tfsi_obmol, 1, 1, None)
        rad_util = AtomicRadiusUtils(covalent_radius_scale=3.0, metal_radius_scale=1.5)
        mol_radius = rad_util.get_radius(acetoxyq_obmol)
        cation_radius = rad_util.get_radius(sodium_obmol)
        anion_radius = rad_util.get_radius(tfsi_obmol)
        mol_coords = IonPlacer.normalize_molecule(acetoxyq_obmol)
        bounder = IonPlacer.get_bounder(mol_coords, sodium_obmol, tfsi_obmol, 1, 1)
        max_cap = max(bounder.upper_bound)
        self.evaluator = LargestContactGapEnergyEvaluator(mol_coords, mol_radius, cation_radius, anion_radius, )

    def test_calc_energy(self):
        self.fail()