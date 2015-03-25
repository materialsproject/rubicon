from __future__ import unicode_literals

import argparse
import copy
from pymatgen.io.qchemio import QcInput

__author__ = 'xiaohuiqu'

def assign_basis_set_by_range(qcinp, basis_def_file):
    with open(basis_def_file) as f:
        basis_def_texts = f.readlines()
    basis_def_token_strs = [line.split() for line in basis_def_texts]
    basis_def_token_values = [[int(t[0]), int(t[1]), t[2]]
                              for t in basis_def_token_strs]
    for i, t in enumerate(basis_def_token_values):
        if t[1] < t[0]:
            print "The current line is:"
            print basis_def_texts[i]
            raise ValueError("The atom order specifications must be in ascending order inside a line")
        if i > 0:
            prev_t = basis_def_token_values[i-1]
            if t[0] - prev_t[1] != 1:
                print "The previous line is:"
                print basis_def_texts[i-1]
                print "The current line is:"
                print basis_def_texts[i]
                raise ValueError("The index of first atom must be exactly 1 larger than the previous line")
    mol = qcinp.jobs[0].mol
    if basis_def_token_values[-1][1] != len(mol):
        raise ValueError("The provided numbers atom in the basis and QChem input file must be consistent")
    elements = [site.species_string for site in mol.sites]
    basis = []
    for t in basis_def_token_values:
        natoms = t[1] - t[0] + 1
        basis.extend([t[2]] * natoms)
    elements_and_basis = zip(elements, basis)
    qcinp_with_basis = copy.deepcopy(qcinp)
    for j in qcinp_with_basis.jobs:
        j.set_basis_set(elements_and_basis)
    return qcinp_with_basis



def main():
    parser = argparse.ArgumentParser(
        description="Assign mixed basis set by list of ranges")
    parser.add_argument("-i", "--input", dest="input", type=str,
                        required=True,
                        help="the QChem input filename")
    parser.add_argument("-b", "--basis", dest="basis", type=str,
                        required=True,
                        help="The file contain the list of basis sets by four columns: 1) start atom index (1 based);"
                             "2) end  atom index; 3) basis set; 4) comments")
    parser.add_argument("-o", "--output", dest="output", type=str,
                        required=True,
                        help="the QChem input filename with the mixed basis set")
    options = parser.parse_args()

    qcinp_no_basis = QcInput.from_file(options.input)
    mol_with_basis = assign_basis_set_by_range(qcinp=qcinp_no_basis, basis_def_file=options.basis)
    mol_with_basis.write_file(options.output)


if __name__ == "__main__":
    main()