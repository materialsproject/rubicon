__author__ = 'navnidhirajput'

from pymatgen.core.bonds import CovalentBond

from pymatgen.core.sites import Site

from topology import TopBond


site1 = Site("C", [0, 0, 0])
site2 = Site("H", [0, 0.7, 0.6])

my_test=CovalentBond(site1,site2)


my_test1=TopBond(site1,site2)


print my_test1.sites()


