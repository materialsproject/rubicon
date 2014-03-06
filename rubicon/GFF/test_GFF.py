__author__ = 'navnidhirajput'

from GFF import GFF
import os
from pymatgen.core.structure import Molecule


my_ff=GFF()
#os.chdir (change the working directory)
print my_ff.read_mass('mol1.rtf')






coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]
mol=Molecule(["C", "H", "H", "H", "H"], coords)
'''
test_ATA=AtomTypeAssigner(mol)
print test_ATA.get_atom_types()
'''




