# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module computes various properties that can extracted from lammps data.
TODO: finish implementation
"""

import numpy as np
import scipy.integrate as sp_integrate


__author__ = "Kiran Mathew"


class TransportProperties(object):
    def __init__(self, lammpsrun):
        self.lammpsrun = lammpsrun

    def get_integrated_correlation(self, array):
        auto_corr_full = np.correlate(array, array,mode="full")
        auto_corr = auto_corr_full[auto_corr_full.size / 2:]
        time = self.lammpsrun.traj_timesteps * self.lammpsrun.timestep
        return sp_integrate.simps(auto_corr, time)

    @property
    def electrical_conductivity(self):
        mol_current = self.lammpsrun.current
        kappa = [self.get_integrated_correlation(mol_current[:,dim]) for dim
                 in range(3)]
        return kappa

    @property
    def diffusivity(self):
        pass

    @property
    def nernst_einstein_conductivity(self):
        pass

    @property
    def viscosity(self, skip):
        if not self.lammpsrun.lammpslog.get('pxy'):
            print("no pressure data")
            raise KeyError
        else:
            nu = [self.get_integrated_correlation(np.array(
                self.lammpsrun.lammpslog[comp][skip:]))
                          for comp in ['pxy', 'pxz', 'pyz', 'pxx', 'pyy',
                                       'pzz'] ]
            return nu

class GeometricProperties(object):
    def __init__(self, lammpsrun):
        self.lammpsrun = lammpsrun

    def radial_distribution(self):
        pass

    def coordination_number(self):
        pass