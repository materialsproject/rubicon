from unittest import TestCase
from rubicon.gff.lammpsin import DictLammpsInputSet

__author__ = 'navnidhirajput'


class TestDictLammpsInputSet(TestCase):
    def test_get_lammpsin(self):
        mylammpsin=DictLammpsInputSet()
        control_file=mylammpsin.get_lammps_control('Lammps.json',ensemble="npt",temp=450)
        ans='''log mol.log
read_start mol.restart
units real
atom_style full
boundary  p p p
pair_style  lj/cut/coul/long 12
kspace_style  0.0001
pair_modify  tail yes mix arithmetic
special_bonds  amber
bond_style  harmonic
dihedral_style  charmm
read_data  data.mol
neighbor  2.0 bin
neigh_modify   delay 0 every 1 check yes page 1000000 one 20000
timestep  1.0
minimize  0.0001 1e-06 10000 10000
velocity  all create 298 314159265 units box
velocity  all zero linear units box
dump  DUMP all custom 2000 mol.lammpstrj id type x y z vx vy vz mol
thermo_style  custom step vol temp press ke pe etotal enthalpy evdwl ecoul epair ebond eangle edihed eimp emol elong etail lx ly lz xy xz yz
thermo  1000
fix  npt all npt temp 450 450 100.0 iso 1.0 1.0 100.0
restart  5000 restart.mol.1 restart.mol.2
run  2000000
write_restart  restart.mol
'''

        self.assertEquals(ans,control_file)

    def test_from_dict(self):
        mylammpsin=DictLammpsInputSet()
        d1=mylammpsin.to_dict
        qc2 = mylammpsin.from_dict(d1)
        d2=qc2.to_dict
        self.assertEqual(d1, d2)