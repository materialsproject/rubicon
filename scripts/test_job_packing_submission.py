__author__ = 'xiaohuiqu'

from fireworks.core.launchpad import LaunchPad
from pymatgen.io.gaussianio import GaussianInput
from rubicon.io.nwchemio_set import JCESRDeltaSCFInputSet
from fireworks.core.firework import FireWork
from rubicon.firetasks.nwchem_task import NWChemTask
from fireworks.core.firework import Workflow

lp = LaunchPad()

with open("../rubicon/testset/g2-97_cart_neut.txt") as f:
    txt = f.read()

for tok in txt.split("--link1--"):
    gau = GaussianInput.from_string(tok)
    gau_mol = gau.molecule
    if len(gau_mol) < 6:
        deltaSCFInputWriter = JCESRDeltaSCFInputSet()
        nwi = deltaSCFInputWriter.get_nwinput(gau_mol)
        fw = FireWork(tasks=[NWChemTask()], spec=nwi.to_dict, name=gau.title.split()[0])
        print fw
        wf = Workflow.from_FireWork(fw)
        print wf
        lp.add_wf(Workflow.from_FireWork(fw))
