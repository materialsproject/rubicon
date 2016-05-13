# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import glob
import os

from rubicon.firetasks.lammps.lammps_properties_task_to_be_replaced import ParselammpsProperties
from rubicon.firetasks.lammps.old_tasks.lammps_output_task_to_be_replaced import WritelammpsOutputTask

import rubicon
from fireworks import Workflow, LaunchPad
from fireworks.core.firework import Firework
from pymatgen import Molecule
from rubicon.firetasks.lammps.old_tasks.input_tasks import WritelammpsInputFromGaussian

__author__ = 'navnidhirajput'

if __name__ == '__main__':
    task1 = WritelammpsInputFromGaussian()
    task2 = WritelammpsOutputTask()
    task3 = ParselammpsProperties()

    coords = []
    sp = []
    solvent_molecules_path = os.path.join(rubicon.__path__[0],
                                          'workflows/test_mols/solvent_molecules')
    moleculelist = glob.glob(solvent_molecules_path + '/*.pdb')
    for filename in moleculelist:
        mol = Molecule.from_file(filename)
        for site in mol:
            coords.append([c for c in site.coords])
            sp.append(site.specie.symbol)
        mol2 = Molecule(sp, coords, site_properties={
            "mol_name": [filename[48:-4]] * len(coords)})
        fw1 = Firework([task1], name='Run Lammps',
                       spec={"molecules": [mol2.as_dict()]}, fw_id=1)
        fw2 = Firework([task2], name='Lammps Log Parsing', fw_id=2)
        fw3 = Firework([task3], name='Lammps Properties Parser', fw_id=3)
        depen = {1: [2, 3]}
        wf = Workflow([fw1, fw2, fw3], name="LAMMPS", links_dict=depen)

        lp = LaunchPad.auto_load()
        lp.add_wf(wf)
