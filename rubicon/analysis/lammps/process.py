# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines classes that extracts information such as the atomic
charges, coordination muber, molecule data and time data from the lammps
output files
"""

import copy

import numpy
from scipy.integrate import cumtrapz

__author__ = "Michael Humbert, Kiran Mathew"


class AtomicCharges(object):
    def findnumatoms(self, datfilename):
        n = None
        datfile = open(datfilename)
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

    def getmolcharges(self, datfilename, n):
        datfile = open(datfilename)
        for j in range(0, 4):
            datfile.readline()
        atomcharges = [0 for x in range(0, n)]
        mol = [0 for x in range(0, n)]
        foundatoms = False
        readingcharges = True

        while foundatoms == False:
            line = datfile.readline()
            line = line.split()

            if len(line) > 0:
                if line[0] == 'Atoms':
                    foundatoms = True
                    datfile.readline()

        while readingcharges == True:
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

    def molchargedict(self, molcharges, moltypel, moltype):
        molcharge = {}
        for molecules in range(0, len(moltypel)):
            molcharge[moltypel[molecules]] = molcharges[
                moltype.index(molecules)]
        return molcharge


class CoordinationNumber(object):
    def calccoordinationnumber(self, output, nummoltype, moltypel, V):
        output['Coordination_Number'] = {}
        output['Coordination_Number'][
            'units'] = 'Minima in angstroms, Coordination numbers in Angstroms'
        output['Coordination_Number'][
            'explanation'] = 'This program finds the first three local minima and finds the coordination number integrating until there. Na-H20 represents the coordination number for water around sodium.'
        pairlist = list(output['RDF'].keys())
        pairlist.remove('units')
        pairlist.remove('distance')
        r = output['RDF']['distance']
        for i in range(0, len(pairlist)):
            g = output['RDF'][pairlist[i]]
            split = pairlist[i].split('-')
            mol1 = split[0]
            mol2 = split[1]
            (minima, index) = self.findfirst3minima(g, r)
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

    def findfirst3minima(self, g, r):
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
                g[i] * nummoltype[moltypel.index(mol)] / V * 4 * numpy.pi * r[
                    i] ** 2)
        integral = cumtrapz(integrallist, x=r)
        integral = integral.tolist()
        return integral


class MoleculeData(object):
    """

             Determines molecule types and number of each molecule type and
             creates a list of molecule type of each molecule

             Requires the following comments in the lammps data file starting
            at the third line

            # "number" "molecule" molecules

            where "number" is the number of that molecule type and
            "molecule" is a name for that molecule

            Do not include blank lines in between the molecule types

    """

    def getmoltype(self, datfilename):
        # determines molecule types and number of each molecule type
        # also creates a list of molecule type of each molecule
        datfile = open(datfilename)
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


class TimeData(object):
    """
            Uses a lammps trajectory file and log file to determine the
            timestep size and the trajectory print frequency

    """

    def getdt(self, logfilename):
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

    def getjump(self, trjfilename):
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