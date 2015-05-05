import os
from unittest import TestCase
from pymatgen import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.qchemio import QcOutput
from rubicon.utils.ion_arranger.hard_sphere_energy_evaluators import AtomicRadiusUtils, ContactDetector
from rubicon.utils.ion_arranger.ion_arranger import IonPlacer
import numpy as np

__author__ = 'xiaohuiqu'


test_dir = os.path.join(os.path.dirname(__file__),
                        "../../../../test_files")

class TestContactDetector(TestCase):
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
        fragments = [sodium_obmol, tfsi_obmol]
        nums_fragments = [1, 1]
        self.acetoxyq_natfsi_placer = IonPlacer(
            acetoxyq_obmol, fragments, nums_fragments, None)
        rad_util = AtomicRadiusUtils(covalent_radius_scale=3.0, metal_radius_scale=1.5)
        mol_radius = rad_util.get_radius(acetoxyq_obmol)
        cation_radius = rad_util.get_radius(sodium_obmol)
        anion_radius = rad_util.get_radius(tfsi_obmol)
        mol_coords = IonPlacer.normalize_molecule(acetoxyq_obmol)
        fragments_atom_radius = [cation_radius, anion_radius]
        nums_fragments = [1, 1]
        self.detector = ContactDetector(mol_coords, mol_radius, fragments_atom_radius, nums_fragments)

    def test_get_contact_matrix(self):
        c = [-40.0, 0.0, 0.0, 20.0, 0.0, 0.0, 1.0, 2.0]
        fragments_coords = self.acetoxyq_natfsi_placer.decode_solution(c)
        contact_matrix = self.detector._get_contact_matrix(fragments_coords)
        ans = np.zeros((3, 3), int)
        self.assertEqual(str(contact_matrix), str(ans))
        c = [-20.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 2.0]
        fragments_coords = self.acetoxyq_natfsi_placer.decode_solution(c)
        contact_matrix = self.detector._get_contact_matrix(fragments_coords)
        ans = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], int)
        self.assertEqual(str(contact_matrix), str(ans))
        c = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0]
        fragments_coords = self.acetoxyq_natfsi_placer.decode_solution(c)
        contact_matrix = self.detector._get_contact_matrix(fragments_coords)
        ans = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], int)
        self.assertEqual(str(contact_matrix), str(ans))

    def test_get_distance_matrix(self):
        start = np.zeros((5, 5), int)
        start[0, 1] = 1
        start[1, 2] = 1
        start[1, 4] = 1
        start[3, 4] = 1
        start[1, 0] = 1
        start[2, 1] = 1
        start[4, 1] = 1
        start[4, 3] = 1
        floydAPSP = self.detector._get_distance_matrix(start)
        self.assertEqual((5, 5), floydAPSP.shape)
        self.assertEquals(1, floydAPSP[0, 1])
        self.assertEquals(2, floydAPSP[0, 2])
        self.assertEquals(3, floydAPSP[0, 3])
        self.assertEquals(2, floydAPSP[0, 4])
        self.assertEquals(1, floydAPSP[1, 2])
        self.assertEquals(2, floydAPSP[1, 3])
        self.assertEquals(1, floydAPSP[1, 4])
        self.assertEquals(3, floydAPSP[2, 3])
        self.assertEquals(2, floydAPSP[2, 4])
        self.assertEquals(1, floydAPSP[3, 4])

    def test_is_contact(self):
        c = [-20.0, 0.0, 0.0, 10.0, 0.0, 0.0, 1.0, 2.0]
        fragments_coords = self.acetoxyq_natfsi_placer.decode_solution(c)
        is_contact = self.detector.is_contact(fragments_coords)
        self.assertFalse(is_contact)
        c = [-20.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 2.0]
        fragments_coords = self.acetoxyq_natfsi_placer.decode_solution(c)
        is_contact = self.detector.is_contact(fragments_coords)
        self.assertFalse(is_contact)
        c = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0]
        fragments_coords = self.acetoxyq_natfsi_placer.decode_solution(c)
        is_contact = self.detector.is_contact(fragments_coords)
        self.assertTrue(is_contact)

