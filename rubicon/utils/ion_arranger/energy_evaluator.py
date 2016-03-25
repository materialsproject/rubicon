# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from abc import ABCMeta, abstractmethod

import six

__author__ = 'xiaohuiqu'


class EnergyEvaluator(six.with_metaclass(ABCMeta, object)):
    """
    Public interface of energy evaluator, all the energy evaluators
    in salt alignment should inherits this interface.
    This interface define the calc_energy() method that will be called
    by the salt alignment module.
    """

    def __init__(self, mol_coords):
        self.mol_coords = mol_coords

    @abstractmethod
    def calc_energy(self, fragments_coords):
        pass

    @abstractmethod
    def taboo_current_position(self):
        pass
