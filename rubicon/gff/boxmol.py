from pymatgen import Molecule

__author__ = 'navnidhirajput'
class BoxMol:
    def __init__(self):
        self.mols = None
        self.no_molecules = None
        self.box_dimension = None
        self.param_list = None
        self.mols_in_box = None

    @classmethod
    def from_packmol(cls, pmr, mols_in_box):
        bm = cls()
        bm.mols = pmr.mols
        bm.param_list = pmr.param_list
        bm.mols_in_box = mols_in_box
        return bm


    @classmethod
    def from_YongRunner(cls, pmr, mols_in_box):
        bm = cls()
        bm.mols = pmr.mols
        bm.param_list = pmr.param_list
        bm.mols_in_box = mols_in_box
        return bm

    @classmethod
    def form_coords(cls, species, coords):
        bm = cls()
        mols = [Molecule(s,c) for s, c in species, coords]
        bm.mols = mols
