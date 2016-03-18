# -*- coding: utf-8 -*-
"""
Created on Fri May  8 10:17:55 2015

@author: mhumbert
"""
import copy

import calccomf

import numpy as np

from scipy.integrate import cumtrapz
from scipy.optimize import curve_fit

from gettimedata import gettimedata


class calcCond:
    g = gettimedata()

    def calcConductivity(self, molcharges, trjfilename, logfilename,
                         datfilename, T, output):
        dt = self.g.getdt(logfilename)
        tsjump = self.g.getjump(trjfilename[0])
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
        return (num_lines, n, num_timesteps, count, line)

    def createarrays(self, n, num_timesteps):
        vx = np.zeros(n)
        vy = np.zeros(n)
        vz = np.zeros(n)
        mol = np.zeros(n)
        atype = np.zeros(n)
        j = np.zeros((num_timesteps, 3))
        J = np.zeros(num_timesteps / 2)
        return (vx, vy, vz, mol, atype, j, J)

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
        return (vxcol, vycol, vzcol, idcol, molcol, typecol)

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
        return (vx, vy, vz, line, mol, atype)

    def calcj(self, molcharges, comvx, comvy, comvz, j, count):
        j[count][0] = np.dot(comvx, molcharges)
        j[count][1] = np.dot(comvy, molcharges)
        j[count][2] = np.dot(comvz, molcharges)
        count += 1
        return (j, count)

    def clacJ(self, j, J):
        for i in range(0, len(j) / 2):
            for k in range(i, i + (len(j) / 2)):
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
        return (comvx, comvy, comvz, nummol, molmass)

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
        (comvxt, comvyt, comvzt) = calccomf.calccom(n, nummol, vx, vy, vz, mol,
                                                    amass, molmass, 0, 0, 0,
                                                    100000, 100000, 100000)
        comvx = copy.deepcopy(comvxt)
        comvy = copy.deepcopy(comvyt)
        comvz = copy.deepcopy(comvzt)

        return (comvx, comvy, comvz)
