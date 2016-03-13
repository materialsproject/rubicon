# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:07:34 2015

@author: mhumbert
"""

import numpy as np
# import matplotlib.pyplot as plt
# import os
import copy


class COMradialdistribution:
    def runradial(self, datfilename, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2,
                  Lz2, output, nummoltype, moltypel, moltype, firststep=1):
        (maxr, binsize, numbins, count, g) = self.setgparam(Lx2, Ly2, Lz2,
                                                            firststep,
                                                            moltypel)
        while count < len(comx):
            (count) = self.radialdistribution(g, len(comx[1]), moltype, comx,
                                              comy, comz, Lx, Ly, Lz, binsize,
                                              numbins, maxr, count)
        (radiuslist) = self.radialnormalization(numbins, binsize, Lx, Ly, Lz,
                                                nummoltype, count, g,
                                                firststep)
        self.append_dict(radiuslist, g, output, moltypel)
        return output

    def setgparam(self, Lx2, Ly2, Lz2, firststep, moltypel):
        # uses side lengths to set the maximum radius for box and number of bins
        # also sets the first line using data on firststep and number of atoms
        maxr = min(Lx2, Ly2, Lz2)
        binsize = 0.1
        numbins = int(np.ceil(maxr / binsize))
        count = firststep
        g = np.zeros((len(moltypel), len(moltypel), numbins))
        return (maxr, binsize, numbins, count, g)

    def radialdistribution(self, g, nummol, moltype, comx, comy, comz, Lx, Ly,
                           Lz, binsize, numbins, maxr, count):
        # implements a FORTRAN code to calculate the number of molecules within each shell

        for molecule1 in range(0, nummol - 1):
            for molecule2 in range(molecule1 + 1, nummol):
                dx = comx[count][molecule1] - comx[count][molecule2]
                dy = comy[count][molecule1] - comy[count][molecule2]
                dz = comz[count][molecule1] - comz[count][molecule2]

                dx -= Lx * round(dx / Lx)
                dy -= Ly * round(dy / Ly)
                dz -= Lz * round(dz / Lz)

                r2 = dx ** 2 + dy ** 2 + dz ** 2
                r = np.sqrt(r2)

                if r <= maxr:
                    g[moltype[molecule1]][moltype[molecule2]][
                        int(round(r / (binsize)) - 1)] += 1
                    g[moltype[molecule2]][moltype[molecule1]][
                        int(round(r / (binsize)) - 1)] += 1

        count += 1
        return (count)

    def radialnormalization(self, numbins, binsize, Lx, Ly, Lz, nummol, count,
                            g, firststep):
        # normalizes g to box density
        radiuslist = (np.arange(numbins) + 1) * binsize
        for i in range(0, len(g)):
            for j in range(0, len(g)):
                g[i][j] *= Lx * Ly * Lz / nummol[i] / nummol[j] / 4 / np.pi / (
                                                                              radiuslist) ** 2 / binsize / (
                           count - firststep)
        return (radiuslist)

    def append_dict(self, radiuslist, g, output, moltypel):
        for i in range(0, len(moltypel)):
            for j in range(i, len(moltypel)):
                output['RDF']['{0}-{1}'.format(moltypel[i],
                                               moltypel[j])] = copy.deepcopy(
                    g[i][j].tolist())
        if 'distance' not in output['RDF'].keys():
            output['RDF']['distance'] = copy.deepcopy(radiuslist.tolist())
