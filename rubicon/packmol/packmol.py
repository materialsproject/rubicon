#!/usr/bin/env python
import tempfile
from monty.io import ScratchDir

from pymatgen import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
import pybel as pb
from subprocess import Popen, PIPE
import os

class PackmolRunner(object):
    """
    Create MD simulation box using packmol.
    """

    def __init__(self,mols,param_list):
        self.mols = mols
        self.param_list = param_list

    def _format_packmol_str(self,v):
        if isinstance(v,list): 
            return ' '.join(str(x) for x in v)
        else:
            return v

    def run(self):

        scratch = tempfile.gettempdir()
        with ScratchDir(scratch) as d:
            # convert mols to pdb files
            for idx, mol in enumerate(self.mols):
                a = BabelMolAdaptor(mol)
                pm = pb.Molecule(a.openbabel_mol)
                pm.write("pdb", filename='{}.pdb'.format(idx), overwrite=True)

            with open('pack.inp', 'w') as input:
                # create packmol control file
                input.write('tolerance 2.0\n')
                input.write('filetype pdb\n')
                input.write('output box.pdb\n')
                for idx, mol in enumerate(self.mols):
                    input.write('\n')
                    input.write('structure {}.pdb\n'.format(idx))
                    for k, v in self.param_list[idx].iteritems():
                        input.write('  {} {}\n'.format(k, self._format_packmol_str(v)))
                    input.write('end structure\n')

            proc = Popen(['./packmol'], stdin=open('pack.inp', 'r'),stdout=PIPE)
            (stdout, stderr) = proc.communicate()

        print stdout
        print proc.returncode, 'IS THE RETURNCODE'
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
