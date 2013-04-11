import os
from fireworks.core.firework import FireWork, Workflow
from pymatgen.io.xyzio import XYZ
from rubicon.firetasks.gaussian_task import GaussianTask
from pymatgen import Molecule

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Apr 10, 2013'


def mol_to_wf(mol, name):
    spec = {'molecule': mol.to_dict,
            'charge': 0,
            'spin_multiplicity': 1,
            'title': 'first test job',
            'functional': 'B3LYP',
            'basis_set': '6-31+G(d)',
            'route_parameters': {'Opt': '', 'SCF': 'Tight'},
            'input_parameters': None,
            'link0_parameters': {'%mem': '100MW', '%chk': 'molecule'},
            '_category': 'Molecules',
            'name': name}

    fw = FireWork([GaussianTask()], spec)

    return Workflow.from_FireWork(fw)


if __name__ == '__main__':
    """
    coords = [[0.000000, 0.000000, 0.000000],
              [0.000000, 0.000000, 1.089000],
              [1.026719, 0.000000, -0.363000],
              [-0.513360, -0.889165, -0.363000],
              [-0.513360, 0.889165, -0.363000]]
    mol = Molecule(["C", "H", "H", "H", "H"], coords)

    wf = mol_to_wf(mol)

    wf.to_file('CH4.yaml')
    """

    module_dir = os.path.dirname(os.path.abspath(__file__))
    for f in os.listdir(os.path.join(module_dir, 'test_mols')):
        if '.xyz' in f:
            mol = XYZ.from_file(os.path.join(module_dir, 'test_mols', f)).molecule
            print mol
            mol_name = f.split('.')[0]
            wf = mol_to_wf(mol, mol_name)
            wf.to_file(os.path.join(module_dir, 'test_wfs', mol_name+'.yaml'))