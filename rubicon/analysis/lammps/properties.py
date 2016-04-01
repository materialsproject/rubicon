# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines classes that uses the processed lammps data to compute
properties such as viscosity, center of mass, conductivities
"""

import copy
from multiprocessing import Pool

import numpy as np

from scipy.integrate import cumtrapz
from scipy.optimize import curve_fit
import scipy.integrate

from six.moves import range
from six.moves import zip

from rubicon.io.lammps.lammpsio import LammpsLog
from rubicon.analysis.lammps.process import TimeData
from rubicon.analysis.lammps._md_analyzer import calccom as calccomf

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
    print("Install matplotlib")


__author__ = "Michael Humbert, A Sharma, Kiran Mathew"


class Viscosity(object):
    """
    All the properties are evaluated based on the *properties input.
    Right now, pymatgen only supports viscosity (argument: viscosity) evaluation.
    Example: 'md_properties log.lammps viscosity' will return the viscosity
        of the system.
    """

    def __init__(self, log_filename, properties):
        self.lammps_log = LammpsLog(log_filename, *properties)

    def compute(self):
        self.lammps_log.parselog()
        self.wline = self.lammps_log.wline
        print('Done reading the log file. Starting Calculations...')
        NCORES = 8
        p = Pool(NCORES)
        if 'viscosity' in self.lammps_log.properties:
            a1 = self.lammps_log.LOG['pxy']
            a2 = self.lammps_log.LOG['pxz']
            a3 = self.lammps_log.LOG['pyz']
            a4 = self.lammps_log.LOG['pxx'] - self.lammps_log.LOG['pyy']
            a5 = self.lammps_log.LOG['pyy'] - self.lammps_log.LOG['pzz']
            a6 = self.lammps_log.LOG['pxx'] - self.lammps_log.LOG['pzz']
            array_array = [a1, a2, a3, a4, a5, a6]
            pv = p.map(self.autocorrelate, array_array)
            pcorr = (pv[0] + pv[1] + pv[2]) / 6 + (pv[3] + pv[4] + pv[5]) / 24
            visco = ( scipy.integrate.cumtrapz(pcorr,
                                               self.lammps_log.LOG['step'][:len(
                                pcorr)])) \
                    * self.lammps_log.timestep * 10 ** -15 * 1000 * 101325. ** 2 \
                    * self.lammps_log.LOG['vol'][-1] \
                    * 10 ** -30 / (1.38 * 10 ** -23 * self.lammps_log.temp)
            #plt.plot(np.array(self.lammps_log.LOG['step'][
            #                  :len(pcorr) - 1]) * self.lammps_log.timestep,
            #         visco)
            #plt.xlabel('Time (femtoseconds)')
            #plt.ylabel('Viscosity (cp)')
            #plt.savefig('viscosity_parallel.png')
            output = open('viscosity_parallel.txt', 'w')
            output.write(
                '#Time (fs), Average Pressure Correlation (atm^2), Viscosity (cp)\n')
            for line in zip(np.array(
                    self.lammps_log.LOG['step'][:len(
                        pcorr) - 1]) * self.lammps_log.timestep - self.lammps_log.cutoff,
                            pcorr,
                            visco):
                output.write(' '.join(str(x) for x in line) + '\n')
            output.close()
            print('Viscosity Calculation Complete!')

    def autocorrelate(self, a):
        b = np.concatenate((a, np.zeros(len(a))), axis=1)
        c = np.fft.ifft(np.fft.fft(b) * np.conjugate(np.fft.fft(b))).real
        d = c[:len(c) / 2]
        d = d / (np.array(list(range(len(a)))) + 1)[::-1]
        return d[:self.wline]


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
        (num_lines, n, num_timesteps, count, line) = self.getnum(trjfilename)
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

    def getnum(self, trjfilename):
        # uses the trjectory file and returns the number of lines and the number of atoms
        trjfile = open(trjfilename[0])
        for i in range(0, 3):
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
    g = TimeData()

    def calcConductivity(self, molcharges, trjfilename, logfilename,
                         datfilename, T, output):
        dt = self.g.dt(logfilename)
        tsjump = self.g.jump(trjfilename[0])
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