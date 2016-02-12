from rubicon.firetasks.lammps_input_task import WritelammpsInputTask
from rubicon.firetasks.multistep_gauss_task import \
    GaussianGeomOptDBInsertionTask, GaussianFreqESPDBInsertionTask
from rubicon.gaussian.gaussian_input_task import WritegaussianGeoTask, \
    WritegaussianFreqESPTask

__author__ = 'navnidhirajput'

from fireworks import  Workflow, LaunchPad
from fireworks.core.firework import Firework


from pymatgen import write_mol, Molecule
import glob
from pymatgen.io import gaussian

import glob
from pymatgen import Molecule
import os


if __name__ == '__main__':

    task_geo = WritegaussianGeoTask()
    task_geo_dbinsert = GaussianGeomOptDBInsertionTask()
    task_freq_esp = WritegaussianFreqESPTask()
    task_freq_esp_dbinsert = GaussianFreqESPDBInsertionTask()
    task_lammps_inp = WritelammpsInputTask()
    coords = []
    sp = []
    moleculelist = glob.glob("/Users/navnidhirajput/Dropbox/solvent_molecules/*")
    for filename in moleculelist:
        mol = Molecule.from_file(filename)
        file_name = os.path.basename(filename)
        fw1 = Firework([task_geo],name = 'Gaussian geometry optimization', spec= {"molecule":mol, "mol_name": os.path.splitext(file_name)[0], "charge": 0,"spin_multiplicity":1}, fw_id=1)
        fw2 = Firework([task_geo_dbinsert],name='Gaussian Geometry DB insertion', spec= {"molecule":mol, "mol_name": os.path.splitext(file_name)[0], "charge": 0,"spin_multiplicity":1}, fw_id=2)
        fw3 = Firework([task_freq_esp],name='Gaussian Frequency and ESP', spec= {"mol_name": os.path.splitext(file_name)[0], "charge": 0,"spin_multiplicity":1}, fw_id=3)
        #fw4 = Firework([task_freq_esp_dbinsert],name='Gaussian Frequency and ESP DB insertion', spec= {"molecule":mol, "mol_name": os.path.splitext(file_name)[0], "charge": 0,"spin_multiplicity":1}, fw_id=4)
        #fw4 = Firework([task_lammps_inp],name = 'Run Lammps', spec={"molecule":mol, "mol_name": os.path.splitext(file_name)[0]}, fw_id=4)
        fw4 = Firework([task_lammps_inp],name = 'Run Lammps', spec={"molecule":mol, "mol_name": [mol]*len(coords)}, fw_id=4)


        depen = {1:2, 2:3,3:4}
        wf = Workflow([fw1,fw2,fw3, fw4], name="GAUSSIAN", links_dict=depen)

        #depen = {1:2, 2:3}
        #wf = Workflow([fw1,fw2,fw3], name="GAUSSIAN", links_dict=depen)

        lp = LaunchPad.auto_load()
        lp.add_wf(wf)

