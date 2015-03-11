from abc import ABCMeta, abstractmethod

__author__ = 'xiaohuiqu'


class EnergyEvaluator(object):
    """
    Public interface of energy evaluator, all the energy evaluators
    in salt alignment should inherits this interface.
    This interface define the calc_energy() method that will be called
    by the salt alignment module.
    """
    __metaclass__ = ABCMeta

    def __init__(self, mol_coords):
        self.mol_coords = mol_coords


    @abstractmethod
    def calc_energy(self, fragments_coords):
        pass

    @abstractmethod
    def taboo_current_position(self):
        pass