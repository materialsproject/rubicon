import argparse

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem import QcInput

__author__ = 'xiaohuiqu'


def main():
    parser = argparse.ArgumentParser(
        description="Replace the atom coordinates of the first job in QChem input file with the coordinates from"
                    "an XYZ file")
    parser.add_argument("-i", "--input", dest="input", type=str,
                        required=True,
                        help="the QChem input filename")
    parser.add_argument("-c", "--coords", dest="coords", type=str,
                        required=True,
                        help="The XYZ file contains the new coords")
    parser.add_argument("-o", "--output", dest="output", type=str,
                        required=True,
                        help="the QChem input filename with the coordinates from the XYZ file")
    options = parser.parse_args()
    qcinp = QcInput.from_file(options.input)
    charge, spin = qcinp.jobs[0].charge, qcinp.jobs[0].spin_multiplicity
    new_mol = Molecule.from_file(options.coords)
    if charge is not None:
        new_mol.set_charge_and_spin(charge, spin)
    qcinp.jobs[0].mol = new_mol
    qcinp.write_file(options.output)
    print(
    "created new QChem input file {new_file} using {old_inp} as an template and filled with coordinates " \
    "from {coord_file}".format(old_inp=options.input,
                               coord_file=options.coords,
                               new_file=options.output))


if __name__ == "__main__":
    main()
