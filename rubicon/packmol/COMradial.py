# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:07:34 2015

@author: mhumbert
"""

import numpy as np

try:
    import comradial
except:
    import comradial2 as comradial


class COMradialdistribution:
    '''
            Calculates the center of mass radial distribution function for a
            given pair of molecules in a system
            
            Outputs are stored in a dictionary called output to later be stored
           in JSON format
           
   '''

    def runradial(self, datfilename, mol1, mol2, comx, comy, comz, Lx, Ly, Lz,
                  Lx2, Ly2, Lz2, output, nummoltype, moltypel, moltype,
                  firststep=1):
        (maxr, binsize, numbins, count, g) = self.setgparam(Lx2, Ly2, Lz2,
                                                            firststep)
        (nummol1, nummol2, mol1, mol2) = self.getnummol(moltypel, nummoltype,
                                                        mol1, mol2)
        while count < len(comx):
            (count) = self.radialdistribution(g, mol1, mol2, len(comx[1]),
                                              moltype, comx, comy, comz, Lx,
                                              Ly, Lz, binsize, numbins, maxr,
                                              count)
        (radiuslist) = self.radialnormalization(numbins, binsize, Lx, Ly, Lz,
                                                nummol1, nummol2, count, g,
                                                firststep)
        self.append_dict(radiuslist, g, output, mol1, mol2, moltypel)
        return output

    def setgparam(self, Lx2, Ly2, Lz2, firststep):
        # uses side lengths to set the maximum radius for box and number of bins
        # also sets the first line using data on firststep and number of atoms
        maxr = min(Lx2, Ly2, Lz2)
        binsize = 0.1
        numbins = int(np.ceil(maxr / binsize))
        count = firststep
        g = np.zeros(numbins)
        return (maxr, binsize, numbins, count, g)

    def getnummol(self, moltypel, nummoltype, mol1, mol2):
        # returns number of each mpolecule type and converts the molecule type to an integer
        nummol1 = nummoltype[moltypel.index(mol1)]
        nummol2 = nummoltype[moltypel.index(mol2)]
        mol1 = int(moltypel.index(mol1))
        mol2 = int(moltypel.index(mol2))
        return (nummol1, nummol2, mol1, mol2)

    def radialdistribution(self, g, mol1, mol2, nummol, moltype, comx, comy,
                           comz, Lx, Ly, Lz, binsize, numbins, maxr, count):
        # implements a FORTRAN code to calculate the number of molecules within each shell

        g1 = comradial.comradial(mol1, mol2, nummol, moltype, comx[count],
                                 comy[count], comz[count], Lx, Ly, Lz, binsize,
                                 numbins, maxr)
        g += g1

        count += 1
        return (count)

    def radialnormalization(self, numbins, binsize, Lx, Ly, Lz, nummol1,
                            nummol2, count, g, firststep):
        # normalizes g to box density
        radiuslist = (np.arange(numbins) + 1) * binsize
        g *= Lx * Ly * Lz / nummol1 / nummol2 / 4 / np.pi / (
                                                            radiuslist) ** 2 / binsize / (
             count - firststep)
        return (radiuslist)

    def append_dict(self, radiuslist, g, output, mol1, mol2, moltypel):
        output['RDF']['{0}-{1}'.format(moltypel[mol1], moltypel[mol2])] = g
        output['RDF']['distance'] = radiuslist
