# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import copy

from pymatgen.io.qchem import QcInput, QcOutput
from rubicon.scripts.run_qchem_standalone import perturb_molecule

__author__ = 'xiaohuiqu'


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
    parser.add_argument("-r", "--reverse", dest="reverse",
                        action="store_true",
                        help="use reversed direction to perturb molecule")
    parser.add_argument("-v", "--verbose", dest="verbose",
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
            print(
                "{} with scale factor {}".format(direction_text,
                                                 options.scale))
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
        raise ValueError(
            "Must have an imaginary frequency to perturb the molecule")


if __name__ == '__main__':
    main()
