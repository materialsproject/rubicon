from pymatgen import Molecule

__author__ = 'navnidhirajput'
class BoxMol:
    def __init__(self, mols, num_mols, box_dims, mol_coords):
        self.mols = mols
        self.box_dims = box_dims
        self.num_mols = num_mols
        self.mols_coords = mol_coords

    @classmethod
    def from_packmol(cls, pmr, mols_coords):
        bm = cls(pmr.mols,
                 [p["number"] for p in pmr.param_list],
                 [p['inside box'] for p in pmr.param_list],
                 mols_coords)
        return bm

    @classmethod
    def from_YongRunner(cls, pmr, mols_in_box):
        bm = cls()
        bm.mols = pmr.mols
        bm.param_list = pmr.param_list
        bm.mols_in_box = mols_in_box
        return bm

    @classmethod
    def from_coords(cls, species, coords):
        bm = cls()
        mols = [Molecule(s,c) for s, c in species, coords]
        bm.mols = mols
