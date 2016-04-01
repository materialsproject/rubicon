# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from pymatgen import Molecule

__author__ = 'navnidhirajput'


class BoxMol:
    def __init__(self, mols, num_mols, box_dims, mol_coords):
        self.mols = mols
        self.box_dims = box_dims
        self.num_mols = num_mols
        self.mols_coords = mol_coords

    @classmethod
    def from_packmol(cls, pmr, mols_coords):
        bm = cls(pmr.mols,
                 [p["number"] for p in pmr.param_list],
                 [p['inside box'] for p in pmr.param_list],
                 mols_coords)
        return bm

    @classmethod
    def from_YongRunner(cls, pmr, mols_in_box):
        bm = cls()
        bm.mols = pmr.mols
        bm.param_list = pmr.param_list
        bm.mols_in_box = mols_in_box
        return bm

    @classmethod
    def from_coords(cls, species, coords):
        bm = cls()
        mols = [Molecule(*sc) for sc in zip(species, coords)]
        bm.mols = mols


class TimeData(object):
    """
    Uses a lammps trajectory file and log file to determine the
    timestep size and the trajectory print frequency

    """

    def dt(self, logfilename):
        dt = None
        logfile = open(logfilename)
        foundtimestep = False
        while foundtimestep == False:
            inline = logfile.readline()
            inline = inline.split()
            if len(inline) > 0:
                if inline[0] == 'timestep':
                    dt = float(inline[1])
                    foundtimestep = True
        logfile.close()
        return dt

    def jump(self, trjfilename):
        trjfile = open(trjfilename)
        trjfile.readline()
        t1 = trjfile.readline()
        t1 = int(t1)
        trjfile.readline()
        n = int(trjfile.readline())
        for i in range(0, n + 6):
            trjfile.readline()
        t2 = int(trjfile.readline())
        tsjump = t2 - t1
        return tsjump