#!/usr/bin/env python

from pymatgen import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
import pybel as pb
from subprocess import Popen, PIPE

class PackmolRunner(object):
    """
    Create MD simulation box using packmol.
    """

    def __init__(self,mols,param_list):
        self.mols = mols
        self.param_list = param_list

    def mols2pdb(self):
        for idx, mol in enumerate(self.mols):
            a = BabelMolAdaptor(mol)
            pm = pb.Molecule(a.openbabel_mol)
            pm.write("pdb",filename="{}.pdb".format(idx),overwrite=True)
        
    def list2str(self,v):
        if isinstance(v,list): 
            return ' '.join(str(x) for x in v)
        else:
            return v

    def generate_packmol_inp(self):
        output=open('pack.inp','w')
        output.write('tolerance 2.0\n')
        output.write('filetype pdb\n')
        output.write('output box.pdb\n')
        for idx,mol in enumerate(self.mols):
            output.write('\n')
            output.write('structure {}.pdb\n'.format(idx))
            for k, v in self.param_list[idx].iteritems():
                output.write('  {} {}\n'.format(k, self.list2str(v)))
            output.write('end structure\n')

    def run_packmol(self):
        infile=open('pack.inp','r')
#        proc = Popen(['./packmol','< pack.inp'],stdin=infile,stdout=None)
        proc = Popen(['./packmol'],stdin=infile,stdout=PIPE)


    def pdb2mol(self):
        a = BabelMolAdaptor.from_file("box.pdb", "pdb")
        return a.pymatgen_mol

    def cleanfiles(self):
        pass


if __name__ == '__main__':
    coords = [[0.000000, 0.000000, 0.000000],
           [0.000000, 0.000000, 1.089000],
           [1.026719, 0.000000, -0.363000],
           [-0.513360, -0.889165, -0.363000],
           [-0.513360, 0.889165, -0.363000]]
    mol = Molecule(["C", "H", "H", "H", "H"], coords) 
    pmr = PackmolRunner([mol, mol], [{"number":4,"inside box":[0.,0.,0.,40.,40.,40.]}, {"number":5,"inside box":[0.,0.,0.,40.,40.,40.]}])
    pmr.mols2pdb()
    pmr.generate_packmol_inp()
    pmr.run_packmol()
#    pmr.pdb2mol()




