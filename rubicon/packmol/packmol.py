#!/usr/bin/env python
import tempfile
from monty.io import ScratchDir

from pymatgen import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
import pybel as pb
from subprocess import Popen, PIPE


class PackmolRunner(object):
    """
    Create MD simulation box using packmol.
    """

    def __init__(self, mols, param_list):
        """
        Create PackMolRunner
        :param mols: [Molecule] - list of Molecules to pack
        :param param_list: [{}, {}] - array of parameters containing dicts for each Structure
        """
        self.mols = mols
        self.param_list = param_list

    def _format_packmol_str(self, some_obj):
        """
        Internal method to format packmol strings
        :param some_obj: Some object to turn into String
        :return:
        """
        if isinstance(some_obj,list):
            return ' '.join(str(x) for x in some_obj)
        else:
            return some_obj

    def run(self):
        """
        Runs packmol
        :return: a Molecule object
        """

        scratch = tempfile.gettempdir()
        with ScratchDir(scratch) as d:
            # convert mols to pdb files
            for idx, mol in enumerate(self.mols):
                a = BabelMolAdaptor(mol)
                pm = pb.Molecule(a.openbabel_mol)
                pm.write("pdb", filename='{}.pdb'.format(idx), overwrite=True)

            with open('pack.inp', 'w') as inp:
                # create packmol control file
                inp.write('tolerance 2.0\n')
                inp.write('filetype pdb\n')
                inp.write('output box.pdb\n')
                for idx, mol in enumerate(self.mols):
                    inp.write('\n')
                    inp.write('structure {}.pdb\n'.format(idx))
                    for k, v in self.param_list[idx].iteritems():
                        inp.write('  {} {}\n'.format(k, self._format_packmol_str(v)))
                    inp.write('end structure\n')

            proc = Popen(['./packmol'], stdin=open('pack.inp', 'r'),stdout=PIPE)
            (stdout, stderr) = proc.communicate()

        a = BabelMolAdaptor.from_file("box.pdb", "pdb")
        return a.pymatgen_mol


if __name__ == '__main__':
    coords = [[0.000000, 0.000000, 0.000000],
           [0.000000, 0.000000, 1.089000],
           [1.026719, 0.000000, -0.363000],
           [-0.513360, -0.889165, -0.363000],
           [-0.513360, 0.889165, -0.363000]]
    mol = Molecule(["C", "H", "H", "H", "H"], coords) 
    pmr = PackmolRunner([mol, mol], [{"number":4,"inside box":[0.,0.,0.,40.,40.,40.]}, {"number":5,"inside box":[0.,0.,0.,40.,40.,40.]}])
    s = pmr.run()
    print s
