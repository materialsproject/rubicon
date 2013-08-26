__author__ = 'xiaohuiqu'

from fireworks.core.launchpad import LaunchPad
from pymatgen.io.gaussianio import GaussianInput
from rubicon.io.nwchemio_set import JCESRDeltaSCFInputSet
from fireworks.core.firework import FireWork
from rubicon.firetasks.nwchem_task import NWChemTask
from fireworks.core.firework import Workflow
from rubicon.workflows.mol_to_wf import mol_to_wf_nwchem

lp = LaunchPad()

with open("../rubicon/testset/g2-97_cart_neut.txt") as f:
    txt = f.read()

for tok in txt.split("--link1--"):
    gau = GaussianInput.from_string(tok)
    gau_mol = gau.molecule
    if len(gau_mol) < 6:
        wf = mol_to_wf_nwchem(mol=gau_mol, name=gau.title.split()[0])
        lp.add_wf(wf)
