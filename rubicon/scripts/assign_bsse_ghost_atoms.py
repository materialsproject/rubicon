# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import argparse
import copy
import itertools

from six.moves import range
from six.moves import zip

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem import QcInput
from rubicon.workflows.bsse_wf import BSSEFragment, get_sub_mol

__author__ = 'xiaohuiqu'


def get_bsse_fragment_files(qcinp_file_name, bsse_fragments):
    parent_qcinp = QcInput.from_file(qcinp_file_name)
    parent_mol = parent_qcinp.jobs[0].mol
    overlapped_qcinps = []
    for bsse_frag in bsse_fragments:
        frag_qcinp_dict = parent_qcinp.as_dict()
        for j in frag_qcinp_dict["jobs"]:
            j['ghost_atoms'] = bsse_frag.ghost_atoms
            j['charge'] = bsse_frag.charge
            j['spin_multiplicity'] = bsse_frag.spin_multiplicity
        frag_qcinp = QcInput.from_dict(frag_qcinp_dict)
        overlapped_qcinps.append(frag_qcinp)
    isolated_qcinps = []
    for bsse_frag in bsse_fragments:
        frag_mol = get_sub_mol(parent_mol, bsse_frag)
        frag_qcinp = copy.deepcopy(parent_qcinp)
        for j in frag_qcinp.jobs:
            if isinstance(j.mol, Molecule):
                j.mol = frag_mol
                j.charge = bsse_frag.charge
                j.spin_multiplicity = bsse_frag.spin_multiplicity
            if 'basis' in j.params:
                parent_basis_def = j.params['basis']
                frag_basis_def = [bs for i, bs in enumerate(parent_basis_def)
                                  if i not in bsse_frag.ghost_atoms]
                j.set_basis_set(frag_basis_def)
        isolated_qcinps.append(frag_qcinp)
    return overlapped_qcinps, isolated_qcinps


def parse_fragments_definition(frag_def_file_name, qcinp_file_name):
    with open(frag_def_file_name) as f:
        frag_def_text = f.readlines()
    parent_qcinp = QcInput.from_file(qcinp_file_name)
    parent_mol = parent_qcinp.jobs[0].mol
    frag_def_token_str = [line.split()[:3] for line in frag_def_text]
    fragment_atoms = []
    bsse_fragments = []
    for i, token in enumerate(frag_def_token_str):
        frag_atoms = []
        for t in token[0].split(','):
            if '-' in t:
                if t.count('-') != 1:
                    print("current line is \'{}\'".format(frag_def_text[i]))
                    raise ValueError("can't understand range \'{}\'".format(t))
                start_index, end_index = [int(x) for x in t.split('-')]
                frag_atoms.extend(list(range(start_index, end_index + 1)))
            else:
                frag_atoms.append(int(t))
        frag_atoms = [x - 1 for x in frag_atoms]
        fragment_atoms.append(frag_atoms)
        ghost_atoms = BSSEFragment.get_host_atoms(frag_atoms, parent_mol)
        charge = int(token[1])
        spin = int(token[2])
        frag = BSSEFragment(charge, spin, ghost_atoms)
        bsse_fragments.append(frag)
    all_atom_from_fragments = itertools.chain(*fragment_atoms)
    if set(all_atom_from_fragments) != set(range(len(parent_mol))):
        print("The current fragments is:")
        for i, frag in enumerate(fragment_atoms):
            print("Fragment {}:".format(i), ', '.join([str(x) for x in frag]))
        raise ValueError("all the fragment should form the complete molecule")
    return bsse_fragments


def main():
    parser = argparse.ArgumentParser(
        description="Assign mixed basis set by list of ranges")
    parser.add_argument("-i", "--input", dest="input", type=str,
                        required=True,
                        help="the QChem input filename")
    parser.add_argument("--fragments", dest="fragments", type=str,
                        required=True,
                        help="The file contains the fragments by four columns: "
                             "1) ranges or single atoms separated by comma (index starts from 1, no spaces);"
                             "2) charge; 3) spin multiplicity; 4) comments")
    parser.add_argument("--overlapped", dest="overlapped", type=str, nargs='+',
                        required=True,
                        help="the list of QChem input files to write isolated fragments")
    parser.add_argument("--isolated", dest="isolated", type=str, nargs='+',
                        required=True,
                        help="the list of QChem input files to write isolated fragments")
    options = parser.parse_args()

    bsse_fragments = parse_fragments_definition(options.fragments,
                                                options.input)
    overlapped_qcinps, isolated_qcinps = get_bsse_fragment_files(options.input,
                                                                 bsse_fragments)
    if len(bsse_fragments) != len(options.overlapped):
        raise ValueError(
            "Please specify {} names of QChem input file for overlapped "
            "fragments".format(len(bsse_fragments)))
    if len(bsse_fragments) != len(options.isolated):
        raise ValueError(
            "Please specify {} names of QChem input file for isolated "
            "fragments".format(len(bsse_fragments)))
    for qcinp, filename in list(zip(overlapped_qcinps, options.overlapped)) + \
            list(zip(isolated_qcinps, options.isolated)):
        qcinp.write_file(filename)


if __name__ == "__main__":
    main()
