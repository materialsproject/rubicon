from rubicon.firetasks.lammps_input_task import WritelammpsInputTask
from rubicon.firetasks.lammps_output_task import WritelammpsOutputTask
from rubicon.firetasks.lammps_properties_task import ParselammpsProperties
from rubicon.firetasks.multistep_gauss_task import \
    GaussianGeomOptDBInsertionTask, GaussianFreqESPDBInsertionTask
from rubicon.gaussian.gaussian_input_task import WritegaussianGeoTask, \
    WritegaussianFreqESPTask

__author__ = 'navnidhirajput'

from fireworks import  Workflow, LaunchPad
from fireworks.core.firework import Firework


from pymatgen import write_mol, Molecule
import glob
from pymatgen.io import gaussianio

import glob
from pymatgen import Molecule
import os


if __name__ == '__main__':

    task_geo = WritegaussianGeoTask()
    task_geo_dbinsert = GaussianGeomOptDBInsertionTask()
    task_freq_esp = WritegaussianFreqESPTask()
    task_freq_esp_dbinsert = GaussianFreqESPDBInsertionTask()
    task_lammps_inp = WritelammpsInputTask()
    task_lammps_log_dbinsert= WritelammpsOutputTask()
    task_lammps_prop_dbinsert = ParselammpsProperties()

    coords = []
    sp = []
    moleculelist = glob.glob("/Users/navnidhirajput/Dropbox/solvent_molecules/*")
    for filename in moleculelist:
        mol = Molecule.from_file(filename)
        file_name = os.path.basename(filename)
        print "filename before", filename
        mol_with_site_prop = Molecule(mol.species, mol.cart_coords,site_properties={"mol_name":[os.path.splitext(file_name)[0]]*len(mol.cart_coords)})
        fw1 = Firework([task_geo],name = 'Gaussian geometry optimization', spec= {"molecule":mol, "mol_name": os.path.splitext(file_name)[0], "charge": 0,"spin_multiplicity":1}, fw_id=1)
        print "filename after", filename
        fw2 = Firework([task_geo_dbinsert],name='Gaussian Geometry DB insertion', spec= {"molecule":mol, "mol_name": os.path.splitext(file_name)[0], "charge": 0,"spin_multiplicity":1}, fw_id=2)
        fw3 = Firework([task_freq_esp],name='Gaussian Frequency and ESP', spec= {"mol_name": os.path.splitext(file_name)[0], "charge": 0,"spin_multiplicity":1}, fw_id=3)
        fw4 = Firework([task_lammps_inp],name = 'Run Lammps', spec={"molecule":mol_with_site_prop}, fw_id=4)
        fw5 = Firework([task_lammps_log_dbinsert],name='Lammps Log Parsing', fw_id=5)
        fw6 = Firework([task_lammps_prop_dbinsert],name='Lammps Properties Parser', fw_id=6)

        depen = {1:2, 2:3,3:4,4:[5,6]}
        wf = Workflow([fw1,fw2,fw3, fw4,fw5], name="LAMMPS", links_dict=depen)

        # depen = {1:2, 2:3,3:4}
        # wf = Workflow([fw1,fw2,fw3, fw4], name="GAUSSIAN", links_dict=depen)


        lp = LaunchPad.auto_load()
        lp.add_wf(wf)

