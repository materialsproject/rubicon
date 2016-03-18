#!/usr/bin/env python
import os
import tempfile
from subprocess import Popen, PIPE

import numpy as np

try:
    import pybel as pb
except ImportError:
    pb = None

from monty.os.path import which
from monty.dev import requires
from monty.tempfile import ScratchDir

from pymatgen import Molecule
from pymatgen.io.babel import BabelMolAdaptor


class PackmolRunner(object):
    """
    Create MD simulation box using packmol.
    """
    @requires(which('packmol'),
              "PackmolRunner requires the executable 'packmol' to be in "
              "the path. Please download packmol from "
              "https://github.com/leandromartinez98/packmol "
              "and follow the instructions in the README to compile. "
              "Don't forget to add the packmol binary to your path")
    def __init__(self, mols, param_list):
        """
        Create PackMolRunner

        Args:
              mols:
                   list of Molecules to pack
              param_list:
                    list of parameters containing dicts for each molecule
        """
        self.mols = mols
        self.param_list = param_list

    def _format_packmol_str(self, some_obj):
        """
        Internal method to format packmol strings

        Args:
              some_obj:
                   Some object to turn into String

        Returns:
               string representation of the object
        """
        if isinstance(some_obj, list):
            return ' '.join(str(x) for x in some_obj)
        else:
            return some_obj

    def _get_auto_boxsize(self):
        """
        TODO: Put docs here!!
        :param idx:
        :raise NotImplementedError:
        """
        volume = 0.0
        for idx, mol in enumerate(self.mols):
            lx, ly, lz = np.max(mol.cart_coords, 0) - np.min(mol.cart_coords,0)
            lx += 2.0
            ly += 2.0
            lz += 2.0
            length = max(lx, ly, lz)
            volume += length ** (3.0) * float(self.param_list[idx]['number'])
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
            # convert mols to pdb files
            for idx, mol in enumerate(self.mols):
                a = BabelMolAdaptor(mol)
                pm = pb.Molecule(a.openbabel_mol)
                pm.write("pdb", filename=os.path.join(d, '{}.pdb'.format(idx)),
                         overwrite=True)
            # TODO: also check if user specified outside box, etc.
            # Do not use auto mode if user specified any type of box
            if 'inside box' not in self.param_list[idx]:
                self._get_auto_boxsize()
            with open(os.path.join(d, 'pack.inp'), 'w') as inp:
                # create packmol control file
                inp.write('tolerance 2.0\n')
                inp.write('filetype pdb\n')
                inp.write('output {}\n'.format(os.path.join(d, "box.pdb")))
                for idx, mol in enumerate(self.mols):
                    inp.write('\n')
                    inp.write(
                        'structure {}.pdb\n'.format(os.path.join(d, str(idx))))

                    for k, v in self.param_list[idx].iteritems():
                        inp.write(
                            '  {} {}\n'.format(k, self._format_packmol_str(v)))
                    inp.write('end structure\n')
            proc = Popen(['packmol'],
                         stdin=open(os.path.join(d, 'pack.inp'), 'r'),
                         stdout=PIPE)
            (stdout, stderr) = proc.communicate()
            a = BabelMolAdaptor.from_file(os.path.join(d, "box.pdb"), "pdb")
            return a.pymatgen_mol


if __name__ == '__main__':
    coords = [[0.000000, 0.000000, 0.000000],
              [0.000000, 0.000000, 1.089000],
              [1.026719, 0.000000, -0.363000],
              [-0.513360, -0.889165, -0.363000],
              [-0.513360, 0.889165, -0.363000]]
    mol = Molecule(["C", "H", "H", "H", "H"], coords)
    #    pmr = PackmolRunner([mol, mol], [{"number":4,"inside box":[0.,0.,0.,40.,40.,40.]}, {"number":5, "inside box":[0.,0.,0.,40.,40.,40.]}])
    pmr = PackmolRunner([mol, mol], [
        {"number": 4, "inside box": [0., 0., 0., 40., 40., 40.]},
        {"number": 5}])
    s = pmr.run()
    print s
