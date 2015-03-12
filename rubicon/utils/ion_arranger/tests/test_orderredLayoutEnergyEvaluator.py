from unittest import TestCase
from rubicon.utils.ion_arranger.hard_sphere_energy_evaluators import OrderredLayoutEnergyEvaluator

__author__ = 'xiaohuiqu'


class TestOrderredLayoutEnergyEvaluator(TestCase):
    def test_calc_energy(self):
        pass

    def test_get_frag_ranks(self):
        frag_coords = [[[3.0, 0.0, 0.0]], [[0.0, 0.0, -.5]], [[0.0, -.0, 0.0]]]
        ranks = OrderredLayoutEnergyEvaluator._get_frag_ranks(frag_coords)
        self.assertEqual(ranks, [3, 1, 2])
        frag_coords = [[[3.0, 0.0, 0.0]], [[0.0, 0.0, 15.0]], [[0.0, 2.0, 0.0]], [[0.0, 4.0, 0.0]], [[0.0, 0.0, 3.0]]]
        ranks = OrderredLayoutEnergyEvaluator._get_frag_ranks(frag_coords)
        self.assertEqual(ranks, [5, 2, 3, 4, 1])

    def test_spearsman_rank_coefficient(self):
        rank_y = [3, 1, 2]
        spearsamn = OrderredLayoutEnergyEvaluator._spearsman_rank_coefficient(rank_y)
        self.assertEqual(spearsamn, -0.5)
        rank_y = [2, 3, 6, 5, 1, 9, 4, 7, 10, 8]
        spearsamn = OrderredLayoutEnergyEvaluator._spearsman_rank_coefficient(rank_y)
        self.assertAlmostEqual(spearsamn, 0.68484848, places=5)
        rank_y = [1, 6, 8, 7, 10, 9, 3, 5, 2, 4]
        spearsamn = OrderredLayoutEnergyEvaluator._spearsman_rank_coefficient(rank_y)
        self.assertAlmostEqual(spearsamn, -0.17575757, places=5)