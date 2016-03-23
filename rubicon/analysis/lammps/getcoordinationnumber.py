# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import copy

import numpy
from scipy.integrate import cumtrapz
from six.moves import range

__author__ = "mhumbert"


class getcoordinationnumber(object):
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
