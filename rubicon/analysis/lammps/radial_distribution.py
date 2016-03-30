# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import copy

from six.moves import range

import linecache
import os

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
    print("Install matplotlib")

from rubicon.analysis.lammps._md_analyzer import comradial, siteradial

__author__ = "Michael Humbert, Kiran Mathew"


class RadialDistribution(object):
    """
    Calculates the center of mass radial distribution function for a
    given pair of molecules in a system
    Outputs are stored in a dictionary called output to later be stored
    in JSON format
    """

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
        return maxr, binsize, numbins, count, g

    def getnummol(self, moltypel, nummoltype, mol1, mol2):
        # returns number of each mpolecule type and converts the molecule type to an integer
        nummol1 = nummoltype[moltypel.index(mol1)]
        nummol2 = nummoltype[moltypel.index(mol2)]
        mol1 = int(moltypel.index(mol1))
        mol2 = int(moltypel.index(mol2))
        return nummol1, nummol2, mol1, mol2

    def radialdistribution(self, g, mol1, mol2, nummol, moltype, comx, comy,
                           comz, Lx, Ly, Lz, binsize, numbins, maxr, count):
        # uses FORTRAN code to calculate the number of molecules within each
        # shell
        g1 = comradial.comradial(mol1, mol2, nummol, moltype, comx[count],
                                 comy[count], comz[count], Lx, Ly, Lz, binsize,
                                 numbins, maxr)
        g += g1
        count += 1
        return count

    def radialnormalization(self, numbins, binsize, Lx, Ly, Lz, nummol1,
                            nummol2, count, g, firststep):
        # normalizes g to box density
        radiuslist = (np.arange(numbins) + 1) * binsize
        g *= Lx * Ly * Lz / nummol1 / nummol2 / 4 / np.pi / (
                                                                radiuslist) ** 2 / binsize / (
                 count - firststep)
        return radiuslist

    def append_dict(self, radiuslist, g, output, mol1, mol2, moltypel):
        output['RDF']['{0}-{1}'.format(moltypel[mol1], moltypel[mol2])] = g
        output['RDF']['distance'] = radiuslist

class SiteRadialDistribution(object):
    def runatomradial(self, trjfilename, datfilename, moltype1, atype1, anum1,
                      moltype2, atype2, anum2, firststep=1):
        (num_lines, n, num_timesteps, count) = self.getnum(trjfilename)
        (Lx, Lx2, Ly, Ly2, Lz, Lz2) = self.getdimensions(trjfilename)
        (maxr, binsize, numbins, line) = self.setgparam(Lx2, Ly2, Lz2, n,
                                                        firststep)
        (g, x, y, z, mol, atype, aid) = self.createarrays(numbins, n)
        (xcol, ycol, zcol, molcol, typecol, idcol) = self.getcolumns(
            trjfilename)
        (linenum, nummoltype, moltypel, moltype) = self.getmoltype(datfilename)
        (atomcount1, atomid1, atomcount2, atomid2, nummol1, nummol2, aindex1,
         aindex2) = self.getnumatominmol(trjfilename, moltype1, atype1, anum1,
                                         moltype2, atype2, anum2, n, line,
                                         idcol, molcol, typecol, moltype, aid,
                                         mol, atype, nummoltype, moltypel)
        while line < num_lines:
            (x, y, z, mol, line) = self.readdata(trjfilename, n, line, idcol,
                                                 xcol, ycol, zcol, molcol, x,
                                                 y, z, mol)
            (aindex1, aindex2) = self.getindeces(count, nummol1, nummol2,
                                                 atomid1, atomid2, atomcount1,
                                                 atomcount2, aindex1, aindex2)
            (g, count) = self.radialdistribution(x, y, z, aindex1, aindex2,
                                                 mol, Lx, Ly, Lz, binsize,
                                                 numbins, maxr, g, count)
            print('timestep ' + str(count + firststep - 1) + ' of ' + str(
                num_timesteps - 1) + ' finished')
        (radiuslist, g) = self.radialnormalization(numbins, binsize, Lx, Ly,
                                                   Lz, nummol1, nummol2, count,
                                                   g)
        self.plot(radiuslist, g)
        self.writetofile(numbins, radiuslist, g)

    def getnum(self, trjfilename):
        # uses the trjectory file and returns the number of lines and the number of atoms
        num_lines = int(sum(1 for line in open(trjfilename)))
        n = int(linecache.getline(trjfilename, 4))
        num_timesteps = int(num_lines / (n + 9))
        count = 0
        return num_lines, n, num_timesteps, count

    def getdimensions(self, trjfilename):
        # uses trjectory file to get the length of box sides
        xbounds = linecache.getline(trjfilename, 6)
        xbounds = xbounds.split()
        ybounds = linecache.getline(trjfilename, 7)
        ybounds = ybounds.split()
        zbounds = linecache.getline(trjfilename, 8)
        zbounds = zbounds.split()
        Lx = float(xbounds[1]) - float(xbounds[0])
        Lx2 = Lx / 2
        Ly = float(ybounds[1]) - float(ybounds[0])
        Ly2 = Ly / 2
        Lz = float(zbounds[1]) - float(zbounds[0])
        Lz2 = Lz / 2
        return Lx, Lx2, Ly, Ly2, Lz, Lz2

    def setgparam(self, Lx2, Ly2, Lz2, n, firststep):
        # uses side lengths to set the maximum radius for box and number of bins
        # also sets the first line using data on firststep and number of atoms
        maxr = min(Lx2, Ly2, Lz2)
        binsize = 0.1
        numbins = int(np.ceil(maxr / binsize))
        line = 10 + int(firststep) * (n + 9)
        return maxr, binsize, numbins, line

    def createarrays(self, numbins, n):
        # creates numpy arrays for data reading
        g = np.zeros(numbins)
        x = np.zeros(n)
        y = np.zeros(n)
        z = np.zeros(n)
        mol = np.zeros(n)
        atype = np.zeros(n)
        aid = np.zeros(n)
        aid = aid.tolist()
        return g, x, y, z, mol, atype, aid

    def getcolumns(self, trjfilename):
        # defines the columns each data type is in in the trjectory file
        inline = linecache.getline(trjfilename, 9)
        inline = inline.split()
        inline.remove('ITEM:')
        inline.remove('ATOMS')
        xcol = inline.index('x')
        ycol = inline.index('y')
        zcol = inline.index('z')
        molcol = inline.index('mol')
        typecol = inline.index('type')
        idcol = inline.index('id')
        return xcol, ycol, zcol, molcol, typecol, idcol

    def getmoltype(self, datfilename):
        # determines molecule types and number of each molecule type
        # also creates a list of molecule type of each molecule
        linenum = 3
        nummoltype = []
        moltypel = []
        moltype = []
        readingmolecules = True
        while readingmolecules == True:
            line = linecache.getline(datfilename, linenum)
            line = line.split()
            if len(line) == 4:
                nummoltype.append(int(line[1]))
                moltypel.append(line[2])
                linenum += 1

            else:
                readingmolecules = False

        for i in range(0, len(moltypel)):
            for j in range(0, nummoltype[i]):
                moltype.append(int(i))

        return linenum, nummoltype, moltypel, moltype

    def getnumatominmol(self, trjfilename, moltype1, atype1, anum1, moltype2,
                        atype2, anum2, n, line, idcol, molcol, typecol,
                        moltype, aid, mol, atype, nummoltype, moltypel):
        for a in range(0, n):
            inline = linecache.getline(trjfilename, line + a)
            inline = inline.split()
            aid[a] = int(inline[idcol])
            mol[a] = int(inline[molcol])
            atype[a] = int(inline[typecol])

        molnum = moltype.index(moltypel.index(moltype1)) + 1
        molaid = []
        molaidtype = []
        typeid = []
        for atom in range(0, n):
            if mol[atom] == molnum:
                molaid.append(aid[atom])
                molaidtype.append(atype[atom])
        for atom in range(0, len(molaid)):
            if molaidtype[atom] == atype1:
                typeid.append(molaid[atom])
        typeid.sort()
        atomcount1 = len(molaid)
        atomid1 = typeid[anum1 - 1]

        molnum = moltype.index(moltypel.index(moltype2)) + 1
        molaid = []
        molaidtype = []
        typeid = []
        for atom in range(0, n):
            if mol[atom] == molnum:
                molaid.append(aid[atom])
                molaidtype.append(atype[atom])
        for atom in range(0, len(molaid)):
            if molaidtype[atom] == atype2:
                typeid.append(molaid[atom])
        typeid.sort()
        atomcount2 = len(molaid)
        atomid2 = typeid[anum1 - 1]
        nummol1 = nummoltype[moltypel.index(moltype1)]
        nummol2 = nummoltype[moltypel.index(moltype2)]
        aindex1 = []
        aindex2 = []

        return (
            atomcount1, atomid1, atomcount2, atomid2, nummol1, nummol2,
            aindex1,
            aindex2)

    def readdata(self, trjfilename, n, line, idcol, xcol, ycol, zcol, molcol,
                 x, y, z, mol):
        for a in range(0, n):
            inline = linecache.getline(trjfilename, line + a)
            inline = inline.split()
            x[int(inline[idcol]) - 1] = float(inline[xcol])
            y[int(inline[idcol]) - 1] = float(inline[ycol])
            z[int(inline[idcol]) - 1] = float(inline[zcol])
            mol[int(inline[idcol]) - 1] = float(inline[molcol])
        line += 9 + n
        return x, y, z, mol, line

    def getindeces(self, count, nummol1, nummol2, atomid1, atomid2, atomcount1,
                   atomcount2, aindex1, aindex2):
        if count == 0:
            for mol1 in range(0, nummol1):
                aindex1.append(atomid1 + mol1 * atomcount1)
            for mol2 in range(0, nummol2):
                aindex2.append(int(atomid2 + mol2 * atomcount2))

        return aindex1, aindex2

    def radialdistribution(self, x, y, z, aindex1, aindex2, mol, Lx, Ly, Lz,
                           binsize, numbins, maxr, g, count):
        g1 = siteradial.siteradial(x, y, z, aindex1, aindex2, mol, Lx, Ly, Lz,
                                   binsize, numbins, maxr)
        g += g1
        count += 1

        return g, count

    def radialnormalization(self, numbins, binsize, Lx, Ly, Lz, nummol1,
                            nummol2, count, g):
        # normalizes g to box density
        radiuslist = (np.arange(numbins) + 1) * binsize
        g *= Lx * Ly * Lz / nummol1 / nummol2 / 4 / np.pi / (
                                                                radiuslist) ** 2 / binsize / count

        return radiuslist, g

    def plot(self, radiuslist, g):
        # plots radial distribution function
        plt.figure()
        plt.plot(radiuslist, g)
        plt.xlabel('radius')
        plt.ylabel('g(r)')
        plt.title('Pair Distribution Function')
        plt.savefig('Pairdist2.png')
        plt.show()

    def writetofile(self, numbins, radiuslist, g):
        try:
            os.remove("agr2.dat")
        except OSError:
            pass

        gfile = open("agr2.dat", "a")
        for radius in range(0, numbins):
            gfile.write(
                str(radiuslist[radius]) + "     " + str(g[radius]) + "\n")
        gfile.close()


class RadialDistributionPure:
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
        return maxr, binsize, numbins, count, g

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
                        int(round(r / binsize) - 1)] += 1
                    g[moltype[molecule2]][moltype[molecule1]][
                        int(round(r / binsize) - 1)] += 1

        count += 1
        return count

    def radialnormalization(self, numbins, binsize, Lx, Ly, Lz, nummol, count,
                            g, firststep):
        # normalizes g to box density
        radiuslist = (np.arange(numbins) + 1) * binsize
        for i in range(0, len(g)):
            for j in range(0, len(g)):
                g[i][j] *= Lx * Ly * Lz / nummol[i] / nummol[j] / 4 / np.pi / (
                                                                                  radiuslist) ** 2 / binsize / (
                               count - firststep)
        return radiuslist

    def append_dict(self, radiuslist, g, output, moltypel):
        for i in range(0, len(moltypel)):
            for j in range(i, len(moltypel)):
                if not all([v == 0 for v in g[i][j]]):
                    output['RDF']['{0}-{1}'.format(moltypel[i],
                                                   moltypel[
                                                       j])] = copy.deepcopy(
                        g[i][j].tolist())
        if 'distance' not in list(output['RDF'].keys()):
            output['RDF']['distance'] = copy.deepcopy(radiuslist.tolist())