# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from fireworks.utilities.fw_serializers import FWSerializable
from six.moves import zip

from pymatgen.core.structure import Molecule

__author__ = 'xiaohuiqu'


class AtomicChargeMixedBasisSetGenerator(FWSerializable):
    """
    Apply mixed basis set scheme to molecule to reduce computational cost.
    Only atom with significantly negative charge will get the diffuse function.

    Args:
        charge_threshold (int): the threshold of atomic charge to consider the
            atom as significantly negatively charged.
        normal_basis_set (str): normal basis set without diffuse function.
        diffuse_basis_set (str): diffuse function augmented basis set

    """

    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "charge_threshold": self.charge_threshold,
                "normal_basis_set": self.normal_basis_set,
                "diffuse_basis_set": self.diffuse_basis_set}

    @classmethod
    def from_dict(cls, m_dict):
        return AtomicChargeMixedBasisSetGenerator(
            charge_threshold=m_dict["charge_threshold"],
            normal_basis_set=m_dict["normal_basis_set"],
            diffuse_basis_set=m_dict["diffuse_basis_set"])

    def __init__(self, charge_threshold=-0.5, normal_basis_set="6-31G*",
                 diffuse_basis_set="6-31+G*"):
        self.charge_threshold = charge_threshold
        self.normal_basis_set = normal_basis_set.lower()
        self.diffuse_basis_set = diffuse_basis_set.lower()

    def get_basis(self, mol, charges):
        if not isinstance(mol, Molecule):
            raise ValueError("Only pymatgen Molecule object is accepted")
        if len(mol) != len(charges):
            raise ValueError("Number of atoms and charges must match")
        elements = [site.species_string for site in mol.sites]
        basis = [
            self.normal_basis_set if c >= self.charge_threshold else self.diffuse_basis_set
            for c in charges]
        return list(zip(elements, basis))
