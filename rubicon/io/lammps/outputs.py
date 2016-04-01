# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module implements classes for processing Lammps output files.
"""

import re
import copy
from io import open

from scipy.integrate import cumtrapz

import numpy as np

from monty.json import MSONable


__author__ = 'Navnidhi Rajput, Michael Humbert, Kiran Mathew'


class LammpsLog(MSONable):
    """
    Parser for LAMMPS log file (parse function).
    Saves the output properties (log file) in the form of a dictionary (LOG)
    with the key being the LAMMPS output property (see 'thermo_style custom'
    command in the LAMMPS documentation).
    For example, LOG['temp'] will return the temperature data array in the log file.
    """

    def __init__(self, llog, avgs=None):
        """
        Args:
            llog (dict):
                Dictionary LOG has all the output property data as numpy 1D arrays with the property name as the key
            avgs:
                Dictionary of averages, will be generated automatically if unspecified
        """
        self.llog = llog
        if avgs:
            self.avgs = avgs
        else:
            self.avgs = {}
            # calculate the average
            for key in self.llog.keys():
                self.avgs[str(key)] = np.mean(self.llog[key])

    @classmethod
    def from_file(cls, filename):
        """
        Parses the log file. 
        """
        md = 0  # To avoid reading the minimization data steps
        header = 0
        footer_blank_line = 0
        llog = {}

        with open(filename, 'r') as logfile:
            total_lines = len(logfile.readlines())
            logfile.seek(0)

            for line in logfile:

                # timestep
                time = re.search('timestep\s+([0-9]+)', line)
                if time:
                    timestep = float(time.group(1))
                    llog['timestep'] = timestep

                # total steps of MD
                steps = re.search('run\s+([0-9]+)', line)
                if steps:
                    md_step = float(steps.group(1))
                    md = 1

                # save freq to log
                thermo = re.search('thermo\s+([0-9]+)', line)
                if thermo:
                    log_save_freq = float(thermo.group(1))

                # log format
                format = re.search('thermo_style.+', line)
                if format:
                    data_format = format.group().split()[2:]

                if all(isinstance(x, float) for x in
                       list(_list2float(line.split()))) and md == 1: break

                header += 1

            # note: we are starting from the "break" above
            for line in logfile:
                if line == '\n':
                    footer_blank_line += 1
            print(int(md_step / log_save_freq))

            if total_lines >= header + md_step / log_save_freq:
                rawdata = np.genfromtxt(fname=filename, dtype=float,
                                        skip_header=header, skip_footer=int(
                        total_lines - header - md_step / log_save_freq - 1) - footer_blank_line)

            else:
                rawdata = np.genfromtxt(fname=filename, dtype=float,
                                        skip_header=header, skip_footer=1)

            for column, property in enumerate(data_format):
                llog[property] = rawdata[:, column]

            return LammpsLog(llog)

    @property
    def timestep(self):
        return self.llog['timestep']

    @property
    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "llog": self.llog,
                "avgs": self.avgs}

    @classmethod
    def from_dict(cls, d):
        return LammpsLog(d['llog'], d['avgs'])


class LammpsRun(object):
    def __init__(self, data_file, trajectory_file, log_file=None):
        self.data_file = data_file
        self.trajectory_file = trajectory_file
        if log_file:
            self.log_file = log_file
            self.llog = LammpsLog.from_file(log_file)

    def natoms(self):
        """
        number of atoms
        """
        n = None
        datfile = open(self.data_file)
        foundnumatoms = False
        for j in range(0, 3):
            datfile.readline()
        while foundnumatoms == False:
            line = datfile.readline()
            line = line.split()
            if len(line) >= 2:
                if line[1] == 'atoms':
                    n = int(line[0])
                    foundnumatoms = True
        datfile.close()
        return n

    def get_mol_charges(self, n):
        """
        get molecule charges
        """
        datfile = open(self.data_file)
        for j in range(0, 4):
            datfile.readline()
        atomcharges = [0 for x in range(0, n)]
        mol = [0 for x in range(0, n)]
        foundatoms = False
        readingcharges = True
        while not foundatoms:
            line = datfile.readline()
            line = line.split()
            if len(line) > 0:
                if line[0] == 'Atoms':
                    foundatoms = True
                    datfile.readline()
        while readingcharges:
            line = datfile.readline()
            line = line.split()
            if len(line) >= 6:
                atomcharges[int(line[0]) - 1] = float(line[3])
                mol[int(line[0]) - 1] = int(line[1])
            else:
                readingcharges = False
        nummol = max(mol)
        molcharges = [0 for x in range(0, nummol)]
        for atom in range(0, n):
            molcharges[mol[atom] - 1] += atomcharges[atom]
        datfile.close()
        return molcharges, atomcharges, n

    def get_mol_charge_dict(self, molcharges, moltypel, moltype):
        """
        molecule charges as dict
        """
        molcharge = {}
        for molecules in range(0, len(moltypel)):
            molcharge[moltypel[molecules]] = molcharges[
                moltype.index(molecules)]
        return molcharge

    def get_mols(self):
        # determines molecule types and number of each molecule type
        # also creates a list of molecule type of each molecule
        datfile = open(self.data_file)
        datfile.readline()
        datfile.readline()
        nummoltype = []
        moltypel = []
        moltype = []
        readingmolecules = True
        while readingmolecules == True:
            line = datfile.readline()
            line = line.split()
            if len(line) == 4:
                nummoltype.append(int(line[1]))
                moltypel.append(line[2])
            else:
                readingmolecules = False
        for i in range(0, len(moltypel)):
            for j in range(0, nummoltype[i]):
                moltype.append(int(i))
        datfile.close()
        return nummoltype, moltypel, moltype

    @property
    def timestep(self):
        return self.llog.timestep

    def jump(self, trajectory_file=None):
        """
        trajectory print frequency
        """
        if trajectory_file:
            traj_file = trajectory_file
        else:
            traj_file = self.trajectory_file
        if isinstance(traj_file, list):
            trjfile = open(traj_file[0])
        else:
            trjfile = open(traj_file)
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


class CoordinationNumber(object):
    """
    This program finds the first three local minima
    and finds the coordination number integrating until
    there. Na-H20 represents the coordination number
    for water around sodium.
    """

    def compute(self, output, nummoltype, moltypel, V):
        output['Coordination_Number'] = {}
        output['Coordination_Number'][
            'units'] = 'Minima in angstroms, Coordination numbers in Angstroms'
        output['Coordination_Number'][
            'explanation'] = 'This program finds the first three local minima ' \
                             'and finds the coordination number integrating until ' \
                             'there. Na-H20 represents the coordination number ' \
                             'for water around sodium.'
        pairlist = list(output['RDF'].keys())
        pairlist.remove('units')
        pairlist.remove('distance')
        r = output['RDF']['distance']
        for i in range(0, len(pairlist)):
            g = output['RDF'][pairlist[i]]
            split = pairlist[i].split('-')
            mol1 = split[0]
            mol2 = split[1]
            (minima, index) = self.get_minima(g, r)
            output['Coordination_Number'][
                '{0} around {1}'.format(mol1, mol2)] = {}
            integral = self.integrate(g, r, nummoltype, moltypel, V, mol1)
            output['Coordination_Number']['{0} around {1}'.format(mol1, mol2)][
                'Cumulative_Integral'] = copy.deepcopy(integral)
            output['Coordination_Number']['{0} around {1}'.format(mol1, mol2)][
                'Minima'] = minima
            coord = []
            for j in range(0, len(minima)):
                coord.append(integral[index[j]])
            output['Coordination_Number']['{0} around {1}'.format(mol1, mol2)][
                'Coordination_Numbers'] = coord
            if mol2 != mol1:
                output['Coordination_Number'][
                    '{0} around {1}'.format(mol2, mol1)] = {}
                integral = self.integrate(g, r, nummoltype, moltypel, V, mol2)
                output['Coordination_Number'][
                    '{0} around {1}'.format(mol2, mol1)][
                    'Cumulative_Integral'] = copy.deepcopy(integral)
                output['Coordination_Number'][
                    '{0} around {1}'.format(mol2, mol1)]['Minima'] = minima
                coord = []
                for j in range(0, len(minima)):
                    coord.append(integral[index[j]])
                output['Coordination_Number'][
                    '{0} around {1}'.format(mol2, mol1)][
                    'Coordination_Numbers'] = coord
        return output

    def get_minima(self, g, r):
        """
        return first 3 minima
        """
        foundpositive = False
        minima = []
        index = []
        i = 0
        while not foundpositive:
            if g[i] > 1:
                foundpositive = True
                i += 1
            else:
                i += 1
        while len(minima) < 3 and i < len(g) - 2:
            if g[i - 1] > g[i] and g[i + 1] > g[i]:
                minima.append(r[i])
                index.append(i)
            i += 1
        return minima, index

    def integrate(self, g, r, nummoltype, moltypel, V, mol):
        integrallist = []
        for i in range(0, len(g)):
            integrallist.append(
                g[i] * nummoltype[moltypel.index(mol)] / V * 4 * np.pi * r[
                    i] ** 2)
        integral = cumtrapz(integrallist, x=r)
        integral = integral.tolist()
        return integral


def _list2float(seq):
    for x in seq:
        try:
            yield float(x)
        except ValueError:
            yield x


if __name__ == '__main__':
    filename = 'visc.log'
    log = LammpsLog.from_file(filename)
    # log.viscosity(100001)
