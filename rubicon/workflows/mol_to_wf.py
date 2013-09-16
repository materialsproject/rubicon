import os
from fireworks.core.firework import FireWork, Workflow
from pymatgen.io.nwchemio import NwTask, NwInput
from pymatgen.io.xyzio import XYZ
from rubicon.firetasks.gaussian_task import GaussianTask
from rubicon.firetasks.nwchem_task import NWChemTask, NWDBInsertionTask


__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Apr 10, 2013'

'''
The workflow decided during the initial meeting between Lei, Shyue, Anubhav, Kristin (~Apr. 10, 2013):

* Use B3LYP functional
* 6-31+G(d) for geometry
* g-311+G(2d,p) for static (or do tests?)
* Vertical IE/EA (delta-SCF)
* Then relax those for doing adiabitic IE/EA
* Vertical and adiabatic (both) for now
* ASMD (PCM) ; test acetone, water, etc
'''

def mol_to_wf(mol, name):
    spec = {'molecule': mol.to_dict,
            'charge': mol.charge,
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


def mol_to_wf_general(gi, name):
    spec = {'molecule': gi.molecule.to_dict,
            'title': gi.title,
            'charge': gi.charge,
            'spin_multiplicity': gi.spin_multiplicity,
            'functional': gi.functional,
            'basis_set': gi.basis_set,
            'route_parameters': gi.route_parameters,
            'input_parameters': gi.input_parameters,
            'link0_parameters': gi.link0_parameters,
            'name': 'gaussian'}

    fw = FireWork([GaussianTask()], spec, name=name)

    return Workflow.from_FireWork(fw, name=name)


def mol_to_wf_nwchem(mol, name):

    tasks = [
        NwTask.dft_task(mol, operation="optimize", xc="b3lyp",
                        basis_set="6-31++G*",theory_directives={"iterations": 300}),
        NwTask.dft_task(mol, operation="freq", xc="b3lyp",
                        basis_set="6-31++G*",theory_directives={"iterations": 300}),
        NwTask.dft_task(mol, operation="energy", xc="b3lyp",
                        basis_set="6-311++G**",theory_directives={"iterations": 300}),
        NwTask.dft_task(mol, charge=mol.charge + 1, operation="energy",
                        xc="b3lyp", basis_set="6-311++G**",theory_directives={"iterations": 300}),
        NwTask.dft_task(mol, charge=mol.charge - 1, operation="energy",
                        xc="b3lyp", basis_set="6-311++G**",theory_directives={"iterations": 300}),
        NwTask.dft_task(mol, charge=mol.charge + 1, operation="energy",
                        xc="b3lyp", basis_set="6-311++G**",alternate_directives={'cosmo':"cosmo"},theory_directives={"iterations": 300}),
        NwTask.dft_task(mol, charge=mol.charge - 1, operation="energy",
                        xc="b3lyp", basis_set="6-311++G**",alternate_directives={'cosmo':"cosmo"},theory_directives={"iterations": 300}),
        NwTask.dft_task(mol, charge=mol.charge, operation="energy",
                        xc="b3lyp", basis_set="6-311++G**",alternate_directives={'cosmo':"cosmo"},theory_directives={"iterations": 300}),
        NwTask.esp_task(mol, charge=mol.charge + 1, operation="",
                         basis_set="6-311++G**"),
        NwTask.esp_task(mol, charge=mol.charge - 1, operation="",
                         basis_set="6-311++G**"),
        NwTask.esp_task(mol, charge=mol.charge, operation="",
                         basis_set="6-311++G**")
    ]

    nwi = NwInput(mol, tasks)
    fw = FireWork([NWChemTask(), NWDBInsertionTask()], spec=nwi.to_dict, name=name)

    return Workflow.from_FireWork(fw, name=name)





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
