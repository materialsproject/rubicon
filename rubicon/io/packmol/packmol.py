# coding: utf-8

from __future__ import division, print_function, absolute_import, \
    unicode_literals

from io import open
import os
import tempfile
from subprocess import Popen, PIPE

import numpy as np
import six

try:
    import pybel as pb
except ImportError:
    pb = None

from monty.os.path import which
from monty.dev import requires
from monty.tempfile import ScratchDir

from pymatgen import Molecule
from pymatgen.io.babel import BabelMolAdaptor

__author__ = 'Michael Humbert, Kiran Mathew'
__email__ = 'mhumbert@nd.edu, kmathew@lbl.gov'


class PackmolRunner(object):
    """
    Create MD simulation box using packmol.
    """

    @requires(which(b'packmol'),
              "PackmolRunner requires the executable 'packmol' to be in "
              "the path. Please download packmol from "
              "https://github.com/leandromartinez98/packmol "
              "and follow the instructions in the README to compile. "
              "Don't forget to add the packmol binary to your path")
    def __init__(self, mols, param_list, input_file="pack.inp",
                 tolerance=2.0, filetype="xyz",
                 control_params={"maxit": 20, "nloop": 600},
                 auto_box=True, output_file="packed.xyz"):
        """
        Create PackMolRunner

        Args:
              mols:
                   list of Molecules to pack
              input_file:
                        name of the packmol input file
              tolerance:
                        min distance between the atoms
              filetype:
                       input/output structure file type
              control_params:
                           packmol control parameters dictionary. Basically
                           all parameters other than structure/atoms
              param_list:
                    list of parameters containing dicts for each molecule
              auto_box:
                    put the molecule assembly in a box
              output_file:
                    output file name. The extension will be adjusted
                    according to the filetype
        """
        self.mols = mols
        self.param_list = param_list
        self.input_file = input_file
        self.boxit = auto_box
        self.control_params = control_params
        if not self.control_params.get("tolerance"):
            self.control_params["tolerance"] = tolerance
        if not self.control_params.get("filetype"):
            self.control_params["filetype"] = filetype
        if not self.control_params.get("output"):
            self.control_params["output"] = \
                "{}.{}".format(output_file.split(".")[0],
                               self.control_params["filetype"])
        if self.boxit:
            self._set_box()

    def _format_param_val(self, param_val):
        """
        Internal method to format values in the packmol parameter dictionaries

        Args:
              param_val:
                   Some object to turn into String

        Returns:
               string representation of the object
        """
        if isinstance(param_val, list):
            return ' '.join(str(x) for x in param_val)
        else:
            return str(param_val)

    def _set_box(self):
        """
        Box the molecule assembly
        """
        volume = 0.0
        for idx, mol in enumerate(self.mols):
            lx, ly, lz = np.max(mol.cart_coords, 0) - np.min(mol.cart_coords,
                                                             0)
            lx += 2.0
            ly += 2.0
            lz += 2.0
            length = max(lx, ly, lz)
            volume += length ** 3.0 * float(self.param_list[idx]['number'])
        length = volume ** (1.0 / 3.0)
        for idx, mol in enumerate(self.mols):
            self.param_list[idx]['inside box'] = '0.0 0.0 0.0 {} {} {}'.format(
                length, length, length)

    def run(self):
        """
        Runs packmol

        Returns:
                Molecule object
        """
        scratch = tempfile.gettempdir()
        with ScratchDir(scratch, copy_to_current_on_exit=True) as d:
            with open(os.path.join(d, self.input_file), 'wt',
                      encoding="utf-8") as inp:
                for k, v in six.iteritems(self.control_params):
                    inp.write('{} {}\n'.format(k, self._format_param_val(v)))
                for idx, mol in enumerate(self.mols):
                    a = BabelMolAdaptor(mol)
                    pm = pb.Molecule(a.openbabel_mol)
                    pm.write(self.control_params["filetype"],
                             filename=os.path.join(d, '{}.{}'.format(idx,
                                                                     self.control_params[
                                                                         "filetype"])).encode(
                                 "ascii"),
                             overwrite=True)
                    inp.write("\n")
                    inp.write(
                        "structure {}.{}\n".format(os.path.join(d, str(idx)),
                                                   self.control_params[
                                                       "filetype"]))
                    for k, v in six.iteritems(self.param_list[idx]):
                        inp.write(
                            '  {} {}\n'.format(k, self._format_param_val(v)))
                    inp.write('end structure\n')
            proc = Popen(['packmol'],
                         stdin=open(os.path.join(d, self.input_file), 'r'),
                         stdout=PIPE)
            (stdout, stderr) = proc.communicate()
            output_file = os.path.join(d, self.control_params["output"])
            if os.path.isfile(output_file):
                packed_mol = BabelMolAdaptor.from_file(output_file)
                print("packed molecule written to {}".format(
                    self.control_params["output"]))
                return packed_mol.pymatgen_mol
            else:
                print("Packmol execution failed")
                print(stdout, stderr)
                return None


if __name__ == '__main__':
    ethanol_coords = [[0.00720, -0.56870, 0.00000],
                      [-1.28540, 0.24990, 0.00000],
                      [1.13040, 0.31470, 0.00000],
                      [0.03920, -1.19720, 0.89000],
                      [0.03920, -1.19720, -0.89000],
                      [-1.31750, 0.87840, 0.89000],
                      [-1.31750, 0.87840, -0.89000],
                      [-2.14220, -0.42390, -0.00000],
                      [1.98570, -0.13650, -0.00000]]
    ethanol = Molecule(["C", "C", "O", "H", "H", "H", "H", "H", "H"],
                       ethanol_coords)
    water_coords = [[9.626, 6.787, 12.673],
                    [9.626, 8.420, 12.673],
                    [10.203, 7.604, 12.673]]
    water = Molecule(["H", "H", "O"], water_coords)
    pmr = PackmolRunner([ethanol, water],
                        [{"number": 1, "fixed": [0, 0, 0, 0, 0, 0],
                          "centerofmass": ""},
                         {"number": 15, "inside sphere": [0, 0, 0, 5]}],
                        input_file="packmol_input.inp", tolerance=2.0,
                        filetype="xyz",
                        control_params={"nloop": 1000},
                        auto_box=False, output_file="cocktail.xyz")
    s = pmr.run()
