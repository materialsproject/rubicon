import json

__author__ = 'xiaohuiqu'

import copy
import logging
import re
import shlex
import os
import socket
from custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QchemJob
import math
from pymatgen import Molecule
from pymatgen.io.qchem import QcInput, QcOutput
import sys

__author__ = 'xiaohuiqu'



def perturb_molecule(old_mol, vib_mode, reversed_direction=False, perturb_scale=0.3):
    max_dis = max([math.sqrt(sum([x ** 2 for x in mode]))
                   for mode in vib_mode])
    scale = perturb_scale / max_dis
    normalized_mode = [[x * scale for x in mode]
                       for mode in vib_mode]
    direction = 1.0
    if reversed_direction:
        direction = -1.0
    new_coords = [[c+v*direction for c, v in zip(site.coords, mode)]
                  for site, mode in zip(old_mol.sites, normalized_mode)]
    species = [site.specie.symbol
               for site in old_mol.sites]
    charge = old_mol.charge
    spin_multiplicity = old_mol.spin_multiplicity
    new_mol = Molecule(species, new_coords, charge=charge,
                       spin_multiplicity=spin_multiplicity)
    return new_mol

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Run A QChem Job for a QChem Input File")
    parser.add_argument("-i", "--input", dest="input", type=str,
                        required=True,
                        help="the QChem output file with imaginary frequency")
    parser.add_argument("-o", "--output", dest="output", type=str,
                        required=True,
                        help="the QChem input file with perturbed geometry")
    parser.add_argument("-s", "--scale", dest="scale", type=float,
                        default=0.3,
                        help="the scale factor to perturb molecule")
    parser.add_argument("-r", "--reverse", dest="reverse", type=bool,
                        action="store_true",
                        help="the scale factor to perturb molecule")
    parser.add_argument("-v", "--verbose", dest="verbose", type=bool,
                        action="store_true",
                        help="print parameters")
    options = parser.parse_args()
    qcout = QcOutput(options.input)
    charge = None
    spin_multiplicity = None
    for d in qcout.data:
        j = d["input"]
        if j.charge is not None:
            charge = j.charge
        if j.spin_multiplicity is not None:
            spin_multiplicity = j.spin_multiplicity
    if qcout.data[-1]['frequencies'][0]["frequency"] < -0.00:
        old_mol = qcout.data[-1]["molecules"][-1]
        vib_mode = qcout.data[-1]['frequencies'][0]["vib_mode"]
        if options.verbose:
            if options.reverse:
                direction_text = "User forward direction"
            else:
                direction_text = "User reversed direction"
            print("{} with scale factor {}".format(direction_text, options.scale))
        new_mol = perturb_molecule(old_mol, vib_mode,
                                   reversed_direction=options.reverse,
                                   perturb_scale=options.scale)
        qctask_freq = qcout.data[-1]["input"]
        qctask_freq.mol = "read"
        qctask_freq.charge = charge
        qctask_freq.spin_multiplicity = spin_multiplicity
        qctask_opt = copy.deepcopy(qctask_freq)
        qctask_opt.params["rem"]["jobtype"] = "opt"
        qctask_opt.params["rem"].pop("scf_guess", None)
        qctask_opt.mol = new_mol
        qcinp = QcInput([qctask_opt, qctask_freq])
        qcinp.write_file(options.output)
    else:
        raise ValueError("Must have an imaginary frequency to perturb the molecule")


if __name__ == '__main__':
    main()

