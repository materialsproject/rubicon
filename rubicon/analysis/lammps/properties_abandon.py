# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import linecache
import os

from matplotlib import pyplot as plt
from scipy import stats

"""
This module defines classes that uses the processed lammps data to compute
properties such as viscosity, center of mass, conductivities
"""

import copy
from io import open
from multiprocessing import Pool

import numpy as np

from scipy.integrate import cumtrapz
from scipy.optimize import curve_fit

from six.moves import range
from six.moves import zip

from rubicon.io.lammps.output import LammpsLog, LammpsRun
from rubicon.analysis.lammps._md_analyzer import calccom as calccomf, comradial, \
    siteradial

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
    print("Install matplotlib")


__author__ = "Michael Humbert, A Sharma"


def autocorrelate(a):
    b = np.concatenate((a, np.zeros(len(a))), axis=1)
    c = np.fft.ifft(np.fft.fft(b) * np.conjugate(np.fft.fft(b))).real
    d = c[:len(c) / 2]
    d = d / (np.array(list(range(len(a)))) + 1)[::-1]
    return d


class Viscosity(object):
    """
    All the properties are evaluated based on the *properties input.
    """

    def __init__(self, log_filename):
        self.llog = LammpsLog(log_filename)

    def compute(self, cutoff, ncores=4):

        """
        cutoff: initial lines ignored during the calculation
        output: a file named viscosity_parallel.txt,
        which saves the correlation function and its integration which is teh viscosity in cP
        """
        p = Pool(ncores)
        a1 = self.llog['pxy'][cutoff:]
        a2 = self.llog['pxz'][cutoff:]
        a3 = self.llog['pyz'][cutoff:]
        a4 = self.llog['pxx'][cutoff:] - self.llog['pyy'][cutoff:]
        a5 = self.llog['pyy'][cutoff:] - self.llog['pzz'][cutoff:]
        a6 = self.llog['pxx'][cutoff:] - self.llog['pzz'][cutoff:]
        array_array = [a1, a2, a3, a4, a5, a6]
        pv = p.map(autocorrelate, array_array)
        pcorr = (pv[0] + pv[1] + pv[2]) / 6 + (pv[3] + pv[4] + pv[5]) / 24
        self.pcorr = pcorr
        temp = np.mean(self.llog['temp'][cutoff:])
        self.viscosity = (cumtrapz(pcorr,
                                          self.llog['step'][:len(pcorr)])) * \
                    self.llog['timestep'] * 10 ** -15 * 1000 * 101325. ** 2 * \
                    self.llog['vol'][-1] * 10 ** -30 / (1.38 * 10 ** -23 * temp)
        p.close()

    def to(self, filename='viscosity.txt'):
        """
        write to file
        """
        with open(filename, 'w') as output:
            output.write('#Time (fs), '
                         'Average Pressure Correlation (atm^2), '
                         'Viscosity (cp)\n')
            for line in zip(np.array(self.llog['step'][:len(pcorr) - 1]) *
                                    self.llog['timestep'], self.pcorr, self.viscosity):
                output.write(' '.join(str(x) for x in line) + '\n')


class CenterOfMass(object):
    """
    Calculates the center of mass for all molecules in a system from
    a lammps trajectory file and a lammps data file

    Requires the following comments in the lammps data file starting
    at the third line

    # "number" "molecule" molecules

    where "number" is the number of that molecule type and
    "molecule" is a name for that molecule

    Do not include blank lines in between the molecule types
    """
    def calcCOM(self, trjfilename, datfilename):
        (num_lines, n, num_timesteps, count, line) = self.getnum(
            trjfilename)
        (Lx, Lx2, Ly, Ly2, Lz, Lz2) = self.getdimensions(trjfilename[0])
        (x, y, z, mol, atype) = self.createarrays(n)
        (xcol, ycol, zcol, molcol, typecol) = self.getcolumns(trjfilename[0])
        atommass = self.getmass(datfilename)
        for i in range(0, len(trjfilename)):
            trjfile = open(trjfilename[i])
            while line[i] < num_lines[i]:
                (x, y, z, mol, atype, line) = self.readdata(trjfile, n, line,
                                                            x, y, z, mol,
                                                            atype, xcol, ycol,
                                                            zcol, molcol,
                                                            typecol, i)
                if count == 0:
                    (nummol, comx, comy, comz, molmass) = self.comprep(mol, n,
                                                                       atype,
                                                                       atommass,
                                                                       num_timesteps)
                (comx, comy, comz, count) = self.calccom(comx, comy, comz, x,
                                                         y, z, mol, atype,
                                                         atommass, molmass, Lx,
                                                         Ly, Lz, Lx2, Ly2, Lz2,
                                                         n, count, nummol)
            trjfile.close()
        return comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2

    def getnum(self, trjfiles):
        # uses the trjectory file and returns the number of lines and the number of atoms
        trjfilename = trjfiles[0]
        trjfile = open(trjfilename)
        for i in range(0, 3):
            trjfile.readline()
        n = int(trjfile.readline())
        trjfile.close()
        num_timesteps = 1
        num_lines = []
        for i in range(0, len(trjfiles)):
            num_lines.append(int(sum(1 for line in open(trjfiles[i]))))
            num_timesteps += int(num_lines[i] / (n + 9)) - 1
        line = [10 for x in trjfiles]
        for j in range(1, len(trjfiles)):
            line[j] += n + 9
        count = 0
        return num_lines, n, num_timesteps, count, line

    def getdimensions(self, trjfilename):
        # uses trjectory file to get the length of box sides
        trjfile = open(trjfilename)
        for i in range(0, 5):
            trjfile.readline()
        xbounds = trjfile.readline()
        xbounds = xbounds.split()
        ybounds = trjfile.readline()
        ybounds = ybounds.split()
        zbounds = trjfile.readline()
        zbounds = zbounds.split()
        Lx = float(xbounds[1]) - float(xbounds[0])
        Lx2 = Lx / 2
        Ly = float(ybounds[1]) - float(ybounds[0])
        Ly2 = Ly / 2
        Lz = float(zbounds[1]) - float(zbounds[0])
        Lz2 = Lz / 2
        trjfile.close()
        return Lx, Lx2, Ly, Ly2, Lz, Lz2

    def createarrays(self, n):
        # creates numpy arrays for data reading
        x = np.zeros(n)
        y = np.zeros(n)
        z = np.zeros(n)
        mol = np.zeros(n)
        atype = np.zeros(n)
        return x, y, z, mol, atype

    def getcolumns(self, trjfilename):
        # defines the columns each data type is in in the trjectory file
        trjfile = open(trjfilename)
        for j in range(0, 8):
            trjfile.readline()
        inline = trjfile.readline()
        inline = inline.split()
        inline.remove('ITEM:')
        inline.remove('ATOMS')
        xcol = inline.index('x')
        ycol = inline.index('y')
        zcol = inline.index('z')
        molcol = inline.index('mol')
        typecol = inline.index('type')
        return xcol, ycol, zcol, molcol, typecol

    def getmass(self, datfilename):
        # returns a dictionary of the mass of each atom type
        atommass = {}
        foundmass = False
        readingmasses = True
        atomnum = 1
        datfile = open(datfilename)
        for i in range(0, 4):
            datfile.readline()

        while foundmass == False:
            line = datfile.readline()
            line = line.split()

            if len(line) > 0:
                if line[0] == 'Masses':
                    foundmass = True
                    datfile.readline()

        while readingmasses == True:
            line = datfile.readline()
            line = line.split()
            if len(line) > 0:
                if int(line[0]) == atomnum:
                    atommass[int(line[0])] = float(line[1])
                    atomnum += 1

                else:
                    readingmasses = False

            else:
                readingmasses = False

        datfile.close()

        return atommass

    def readdata(self, trjfile, n, line, x, y, z, mol, atype, xcol, ycol, zcol,
                 molcol, typecol, i):
        # reads data from trjectory file into precreated arrays
        for j in range(0, 9):
            trjfile.readline()
        for a in range(0, n):
            inline = trjfile.readline()
            inline = inline.split()
            x[a] = inline[xcol]
            y[a] = inline[ycol]
            z[a] = inline[zcol]
            mol[a] = inline[molcol]
            atype[a] = inline[typecol]

        line[i] += n + 9
        return x, y, z, mol, atype, line

    def comprep(self, mol, n, atype, atommass, num_timesteps):
        # creates arrays to prepare for center of mass calculations
        nummol = int(max(mol))
        comx = [[0 for x in range(nummol)] for x in range(num_timesteps)]
        comy = [[0 for x in range(nummol)] for x in range(num_timesteps)]
        comz = [[0 for x in range(nummol)] for x in range(num_timesteps)]

        molmass = np.zeros(nummol)
        for atom in range(0, n):
            molmass[mol[atom] - 1] += atommass[atype[atom]]

        return nummol, comx, comy, comz, molmass

    def calccom(self, comx, comy, comz, x, y, z, mol, atype, atommass, molmass,
                Lx, Ly, Lz, Lx2, Ly2, Lz2, n, count, nummol):
        # calculates the center of mass for each molecule
        amass = np.zeros(n)
        for i in range(0, n):
            amass[i] = atommass[atype[i]]

        (comxt, comyt, comzt) = calccomf(n, nummol, x, y, z, mol,
                                         amass, molmass, Lx, Ly, Lz,
                                         Lx2, Ly2, Lz2)
        comx[count] += comxt
        comy[count] += comyt
        comz[count] += comzt
        count += 1

        return comx, comy, comz, count


class Conductivity(object):

    def calcConductivity(self, molcharges, trjfilename, logfilename,
                         datfilename, T, output):
        lrun = LammpsRun(datfilename, trjfilename, logfilename)
        dt = lrun.timestep
        tsjump = lrun.jump(trjfilename)
        (num_lines, n, num_timesteps, count, line) = self.getnum(trjfilename)
        atommass = self.getmass(datfilename)
        V = self.getdimensions(trjfilename[0])
        (vx, vy, vz, mol, atype, j, J) = self.createarrays(n, num_timesteps)
        (vxcol, vycol, vzcol, idcol, molcol, typecol) = self.getcolumns(
            trjfilename[0])
        for i in range(0, len(trjfilename)):
            trjfile = open(trjfilename[i])
            while line[i] < num_lines[i]:
                (vx, vy, vz, line, mol, atype) = self.readdata(trjfile, n,
                                                               line, vx, vy,
                                                               vz, vxcol,
                                                               vycol, vzcol,
                                                               idcol, i, mol,
                                                               molcol, atype,
                                                               typecol)
                if count == 0:
                    (comvx, comvy, comvz, nummol, molmass) = self.COMprep(mol,
                                                                          atype,
                                                                          atommass,
                                                                          n)
                (comvx, comvy, comvz) = self.calcCOMv(comvx, comvy, comvz, vx,
                                                      vy, vz, mol, atype,
                                                      atommass, molmass, n,
                                                      nummol)
                (j, count) = self.calcj(molcharges, comvx, comvy, comvz, j,
                                        count)
        J = self.clacJ(j, J)
        integral = self.integrateJ(J, tsjump * dt)
        time = []
        for i in range(0, len(integral)):
            time.append(i * tsjump * dt)
        popt = self.fitcurve(time, integral)
        cond = self.greenkubo(popt, T, V)
        output['Conductivity']['Green_Kubo'] = cond
        file = open('integral.dat', 'w')
        for i in range(0, len(integral)):
            file.write(str(integral[i]) + '\n')
        return output

    def getdimensions(self, trjfilename):
        # uses trjectory file to get the length of box sides
        trjfile = open(trjfilename)
        for i in range(0, 5):
            trjfile.readline()
        xbounds = trjfile.readline()
        xbounds = xbounds.split()
        ybounds = trjfile.readline()
        ybounds = ybounds.split()
        zbounds = trjfile.readline()
        zbounds = zbounds.split()
        Lx = float(xbounds[1]) - float(xbounds[0])
        Ly = float(ybounds[1]) - float(ybounds[0])
        Lz = float(zbounds[1]) - float(zbounds[0])
        V = Lx * Ly * Lz / 10 ** 30
        trjfile.close()
        return V

    def getnum(self, trjfilename):
        # uses the trjectory file and returns the number of lines and the number of atoms
        trjfile = open(trjfilename[0])
        for j in range(0, 3):
            trjfile.readline()
        n = int(trjfile.readline())
        trjfile.close()
        num_timesteps = 1
        num_lines = []
        for i in range(0, len(trjfilename)):
            num_lines.append(int(sum(1 for line in open(trjfilename[i]))))
            num_timesteps += int(num_lines[i] / (n + 9)) - 1
        line = [10 for x in trjfilename]
        for j in range(1, len(trjfilename)):
            line[j] += n + 9
        count = 0
        return num_lines, n, num_timesteps, count, line

    def createarrays(self, n, num_timesteps):
        vx = np.zeros(n)
        vy = np.zeros(n)
        vz = np.zeros(n)
        mol = np.zeros(n)
        atype = np.zeros(n)
        j = np.zeros((num_timesteps, 3))
        J = np.zeros(num_timesteps / 2)
        return vx, vy, vz, mol, atype, j, J

    def getcolumns(self, trjfilename):
        # defines the columns each data type is in in the trjectory file
        trjfile = open(trjfilename)
        for j in range(0, 8):
            trjfile.readline()
        inline = trjfile.readline()
        inline = inline.split()
        inline.remove('ITEM:')
        inline.remove('ATOMS')
        vxcol = inline.index('vx')
        vycol = inline.index('vy')
        vzcol = inline.index('vz')
        idcol = inline.index('id')
        molcol = inline.index('mol')
        typecol = inline.index('type')
        trjfile.close()
        return vxcol, vycol, vzcol, idcol, molcol, typecol

    def readdata(self, trjfile, n, line, vx, vy, vz, vxcol, vycol, vzcol,
                 idcol, i, mol, molcol, atype, typecol):
        # reads data from trjectory file into precreated arrays
        for j in range(0, 9):
            trjfile.readline()
        for a in range(0, n):
            inline = trjfile.readline()
            inline = inline.split()
            vx[int(inline[idcol]) - 1] = float(inline[vxcol])
            vy[int(inline[idcol]) - 1] = float(inline[vycol])
            vz[int(inline[idcol]) - 1] = float(inline[vzcol])
            mol[int(inline[idcol]) - 1] = int(inline[molcol])
            atype[int(inline[idcol]) - 1] = int(inline[typecol])

        line[i] += n + 9
        return vx, vy, vz, line, mol, atype

    def calcj(self, molcharges, comvx, comvy, comvz, j, count):
        j[count][0] = np.dot(comvx, molcharges)
        j[count][1] = np.dot(comvy, molcharges)
        j[count][2] = np.dot(comvz, molcharges)
        count += 1
        return j, count

    def clacJ(self, j, J):
        for i in range(0, int(len(j)/2)):
            for k in range(i, i + int((len(j)/2))):
                J[k - i] += np.dot(j[i], j[k])
        J *= float(2) / float(len(j))
        return J

    def integrateJ(self, J, dt):
        integral = cumtrapz(J, dx=dt)
        return integral

    def exp(self, x, a, b, c):
        return a * np.exp(-b * x) + c

    def fitcurve(self, time, integral):
        sigma = 1 / np.sqrt(time[1:])
        popt, pcov = curve_fit(self.exp, time[1:], integral[1:], sigma=sigma)
        return popt

    def greenkubo(self, popt, T, V):
        k = 1.38e-23
        el = 1.60217e-19
        cond = popt[2] / 3 / k / T / V * el ** 2 / 10 ** 5
        return cond

    def COMprep(self, mol, atype, atommass, n):
        nummol = max(mol)
        comvx = np.zeros(nummol)
        comvy = np.zeros(nummol)
        comvz = np.zeros(nummol)

        molmass = np.zeros(nummol)
        for atom in range(0, n):
            molmass[mol[atom] - 1] += atommass[atype[atom]]
        return comvx, comvy, comvz, nummol, molmass

    def getmass(self, datfilename):
        # returns a dictionary of the mass of each atom type
        atommass = {}
        foundmass = False
        readingmasses = True
        atomnum = 1
        datfile = open(datfilename)
        for i in range(0, 4):
            datfile.readline()

        while foundmass == False:
            line = datfile.readline()
            line = line.split()

            if len(line) > 0:
                if line[0] == 'Masses':
                    foundmass = True
                    datfile.readline()

        while readingmasses == True:
            line = datfile.readline()
            line = line.split()
            if len(line) > 0:
                if int(line[0]) == atomnum:
                    atommass[int(line[0])] = float(line[1])
                    atomnum += 1

                else:
                    readingmasses = False

            else:
                readingmasses = False
        datfile.close()
        return atommass

    def calcCOMv(self, comvx, comvy, comvz, vx, vy, vz, mol, atype, atommass,
                 molmass, n, nummol):
        amass = np.zeros(n)
        for i in range(0, n):
            amass[i] = atommass[atype[i]]
        (comvxt, comvyt, comvzt) = calccomf(n, nummol, vx, vy, vz, mol,
                                                    amass, molmass, 0, 0, 0,
                                                    100000, 100000, 100000)
        comvx = copy.deepcopy(comvxt)
        comvy = copy.deepcopy(comvyt)
        comvz = copy.deepcopy(comvzt)

        return comvx, comvy, comvz


class NernstEinsteinConductivity(object):
    def calcNEconductivity(self, output, molcharge, Lx, Ly, Lz, nummoltype,
                           moltypel, T):
        V = Lx * Ly * Lz * 10 ** -30
        e = 1.60217657e-19
        k = 1.3806488e-23
        NEcond = 0
        for i in range(0, len(moltypel)):
            q = float(molcharge[moltypel[i]])
            if q != 0:
                try:
                    D = float(output['Diffusivity'][moltypel[i]])
                except ValueError:
                    output[
                        'Nernst Einstien Conductivity in S/m'] = 'runtime not long enough'
                    return output
                N = int(nummoltype[i])
                NEcond += N * q ** 2 * D
        NEcond *= e ** 2 / k / T / V

        output['Conductivity']['Nernst_Einstein'] = NEcond

        return output


class MeanSquareDisplacement(object):
    """
    Calculates the MSD and diffusivity for all molecules in the system
    given a list of center of mass coordinates.

    Outputs are stored in a dictionary called output to later be stored
    in JSON format
    """

    def runMSD(self, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, moltype,
               moltypel, dt, tsjump, output):
        (comx, comy, comz) = self.unwrap(comx, comy, comz, Lx, Ly, Lz, Lx2,
                                         Ly2, Lz2)
        num_timesteps = len(comx)
        (MSDt, MSD, diffusivity) = self.gettimesteps(num_timesteps, moltypel)
        (molcheck, nummol) = self.setmolarray(moltype, moltypel)
        for i in range(0, MSDt):
            for j in range(i, i + MSDt):
                r2 = self.calcr2(comx, comy, comz, i, j)
                MSD = self.MSDadd(r2, MSD, molcheck, i, j)
        MSD = self.MSDnorm(MSD, MSDt, nummol)
        Time = self.createtime(dt, tsjump, MSDt)
        (lnMSD, lntime) = self.takelnMSD(MSD, Time)
        for molecule in range(0, len(moltypel)):
            firststep = self.findlinearregion(lnMSD, lntime, dt, molecule)
            self.getdiffusivity(Time, MSD, firststep, molecule, diffusivity)
        self.append_dict(MSD, moltypel, diffusivity, output, Time)
        return output

    def unwrap(self, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2):
        for i in range(1, len(comx)):
            for j in range(0, len(comx[i])):
                if (comx[i][j] - comx[i - 1][j]) > Lx2:
                    while (comx[i][j] - comx[i - 1][j]) > Lx2:
                        comx[i][j] -= Lx
                elif (comx[i][j] - comx[i - 1][j]) < (-Lx2):
                    while (comx[i][j] - comx[i - 1][j]) < (-Lx2):
                        comx[i][j] += Lx

                if (comy[i][j] - comy[i - 1][j]) > Ly2:
                    while (comy[i][j] - comy[i - 1][j]) > Ly2:
                        comy[i][j] -= Ly
                elif (comy[i][j] - comy[i - 1][j]) < (-Ly2):
                    while (comy[i][j] - comy[i - 1][j]) < (-Ly2):
                        comy[i][j] += Ly

                if (comz[i][j] - comz[i - 1][j]) > Lz2:
                    while (comz[i][j] - comz[i - 1][j]) > Lz2:
                        comz[i][j] -= Lz
                elif (comz[i][j] - comz[i - 1][j]) < (-Lz2):
                    while (comz[i][j] - comz[i - 1][j]) < (-Lz2):
                        comz[i][j] += Lz
        return comx, comy, comz

    def gettimesteps(self, num_timesteps, moltypel):
        MSDt = int(np.floor(num_timesteps / 2))
        MSD = np.zeros((len(moltypel), MSDt))
        diffusivity = []
        return MSDt, MSD, diffusivity

    def setmolarray(self, moltype, moltypel):
        molcheck = np.zeros((len(moltypel), len(moltype)))
        for i in range(0, len(moltype)):
            molcheck[moltype[i]][i] = 1
        nummol = np.zeros(len(moltypel))
        for i in range(0, len(nummol)):
            nummol[i] = np.sum(molcheck[i])
        return molcheck, nummol

    def calcr2(self, comx, comy, comz, i, j):
        r2 = (comx[j] - comx[i]) ** 2 + (comy[j] - comy[i]) ** 2 + (comz[j] -
                                                                    comz[
                                                                        i]) ** 2
        return r2

    def MSDadd(self, r2, MSD, molcheck, i, j):
        for k in range(0, len(molcheck)):
            sr2 = np.dot(r2, molcheck[k])
            MSD[k][j - i] += sr2
        return MSD

    def MSDnorm(self, MSD, MSDt, nummol):
        for i in range(0, len(nummol)):
            MSD[i] /= MSDt * nummol[i]
        return MSD

    def createtime(self, dt, tsjump, MSDt):
        Time = np.arange(0, MSDt, dtype=np.float64)
        Time *= dt * tsjump
        return Time

    def takelnMSD(self, MSD, Time):
        lnMSD = np.log(MSD[:, 1:])
        lntime = np.log(Time[1:])
        return lnMSD, lntime

    def findlinearregion(self, lnMSD, lntime, dt, molecule):
        timestepskip = int(500 / dt)
        linearregion = True
        maxtime = len(lnMSD[0])
        numskip = 1
        tolerance = 0.05
        while linearregion == True:
            if numskip * timestepskip + 1 > maxtime:
                firststep = maxtime - 1 - (numskip - 1) * timestepskip
                return firststep
                linearregion = False
            else:
                t1 = maxtime - 1 - (numskip - 1) * timestepskip
                t2 = maxtime - 1 - numskip * timestepskip
                slope = (lnMSD[molecule][t1] - lnMSD[molecule][t2]) / (
                    lntime[t1] - lntime[t2])
                if abs(slope - 1.) < tolerance:
                    numskip += 1
                else:
                    firststep = t1
                    return firststep
                    linearregion = False

    def getdiffusivity(self, Time, MSD, firststep, molecule, diffusivity):
        calctime = []
        calcMSD = []
        for i in range(firststep, len(Time)):
            calctime.append(Time[i])
            calcMSD.append(MSD[molecule][i])
        if len(calctime) == 1:
            diffusivity.append('runtime not long enough')
        else:
            line = stats.linregress(calctime, calcMSD)
            slope = line[0]
            diffusivity.append(slope / 600000)

    def append_dict(self, MSD, moltypel, diffusivity, output, Time):
        output['MSD'] = {}
        output['MSD']['units'] = 'Angstroms^2, fs'
        output['Diffusivity'] = {}
        output['Diffusivity']['units'] = 'm^2/s'
        for i in range(0, len(moltypel)):
            output['MSD'][moltypel[i]] = copy.deepcopy(MSD[i].tolist())
            output['Diffusivity'][moltypel[i]] = copy.deepcopy(diffusivity[i])
        output['MSD']['time'] = Time.tolist()


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
        g1 = siteradial(x, y, z, aindex1, aindex2, mol, Lx, Ly, Lz,
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