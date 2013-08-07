from fireworks.core.launchpad import LaunchPad
from pymatgen import Element, Molecule
from pymatgen import Element, Molecule
from pymatgen.io.gaussianio import GaussianInput
import glob
from rubicon.workflows.mol_to_wf import mol_to_wf, mol_to_wf_general, mol_to_wf_nwchem

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Apr 02, 2013'


def parse_file(filename):
    with open(filename) as f:
        txt = f.read()
    toks = txt.split("--link1--")
    for t in toks:

        gau = GaussianInput.from_string(t.strip())
        yield gau


fireworks_file = 'launchpad_jcesr_PROD.yaml'

lp = LaunchPad.from_file(fireworks_file)
lp.reset('', require_password=False)
all_mol = []
for f in glob.glob("G3/g305*"):
    for gi in parse_file(f):
        clean_sites = []
        for site in gi.molecule:
            if Element.is_valid_symbol(site.specie.symbol):
                clean_sites.append(site)

        mol = Molecule([site.species_and_occu for site in clean_sites],
                   [site.coords for site in clean_sites],
                   charge=gi.charge, spin_multiplicity=gi.spin_multiplicity,
                   validate_proximity=False,
                   site_properties=None)

        wf = mol_to_wf_nwchem(mol, gi.title.split(' ')[0])
        lp.add_wf(wf)