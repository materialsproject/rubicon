#!/usr/bin/env python



from __future__ import print_function

import os

from rubicon.workflows.mol_to_wf import mol_to_wf_nwchem

from pymatgen.analysis.molecule_matcher import InchiMolAtomMapper
from pymatgen.core.structure import Molecule, FunctionalGroups
from pymatgen.io.babel import BabelMolAdaptor

"""
Defining a molecule to be substituted.
coords =[[0.000000, 0.000000, 0.000000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
            [0.000000, 0.000000,   1.540000],
            [0.000000,  -1.010000,  1.900000],
            [-0.870000,    0.510000,     1.900000],
            [0.870000,   0.510000,    1.900000]]
mol = Molecule(["C","H", "H", "H","C","H","H","H"], coords)
"""

"""
Read the xyz coordinate of a molecule from a file and convert into molecule
"""
f = open("thiane")
species_list = []
coords_list = []
for line in f:
    line_splitted = line.split()
    species_list.append(line_splitted[0])
    line_splitted.pop(0)
    line_splitted_float = []
    for coord in line_splitted:
        coord = float(coord)
        line_splitted_float.append(coord)
    coords_list.append(line_splitted_float)
mol = Molecule(species_list, coords_list)

"""
Find equivalent non-H sites then store a list of non-H atom without equivalent sites (label_list).
"""
mapper = InchiMolAtomMapper()
labelInfo = mapper._inchi_labels(BabelMolAdaptor(mol)._obmol)
# print(labelInfo,len(labelInfo))
label_list = list(labelInfo[0])
for equivalent_group in labelInfo[1]:
    for atomLabel in range(0, len(equivalent_group)):
        if atomLabel > 0:
            label_list.remove(labelInfo[0][equivalent_group[atomLabel] - 1])
print("Atomic numbers of non-H non-equivalent sites are: ", label_list)

"""
Find a list of sites that can be substituted using functional groups. These sites are hydrogen that
directly bound to either C or N atoms.
"""
substitute_sitelist = []
for i in label_list:
    if mol[i - 1].species_string in ["C", "N"]:
        """
        If site is carbon or nitrogen, find a hydrogen that directly bind to the site (<1.3 in
        distance) and store its index in "substitute_sitelist" for future substitution.
        Note there might be more than one H bound to C or N but only one is substituted.
        """
        for j, site_neighbor in enumerate(mol):
            if mol[i - 1].distance(
                    site_neighbor) < 1.3 and site_neighbor.species_string in [
                "H"]:
                substitute_sitelist += [j]
                break
print("Indices (in molecule) of atoms to be substituted are: ", end='')
print(substitute_sitelist)

"""
Do the actual substitution, first-order force-field structure clean-up
and store derivatives in dictionary "derivatives" (name as key and molecule object as value).
"""

derivatives = {}
for index in substitute_sitelist:
    for func_grp in FunctionalGroups():
        name = str(index) + "_" + func_grp
        mol2 = mol.copy()
        mol2.substitute(index, func_grp)
        a = BabelMolAdaptor(mol2)
        a.localopt()
        derivatives[name] = a.pymatgen_mol

print("Total of ", end='')
print(len(derivatives.keys()), end='')
print(" derivatives are created")

'''
Mol to WF
'''
module_dir = os.path.dirname(os.path.abspath(__file__))

for name in derivatives:
    print(name)
    wf = mol_to_wf_nwchem(derivatives[name], name)
    wf.to_file(os.path.join(module_dir, 'thiane_wfs', name + '.yaml'))
