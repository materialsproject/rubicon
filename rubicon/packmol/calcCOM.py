# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:07:34 2015

@author: mhumbert
"""

import numpy as np
import linecache
import os
import time
import calccomf


class calcCOM:
    
    def calcCOM(self, trjfilename, datfilename):
        (num_lines, n, num_timesteps, count, line)=self.getnum(trjfilename)
        (Lx, Lx2, Ly, Ly2, Lz, Lz2) = self.getdimensions(trjfilename[0])  
        (x,y,z,mol,atype) = self.createarrays(n)
        (xcol, ycol, zcol, molcol, typecol) = self.getcolumns(trjfilename[0])
        atommass = self.getmass(datfilename)
        for i in range(0,len(trjfilename)):
            print i
            while line[i] < num_lines[i]:
                (x,y,z,mol,atype,line) = self.readdata(trjfilename[i], n, line, x, y, z, mol, atype, xcol, ycol, zcol, molcol, typecol,i)
                if count == 0:
                    (nummol, comx, comy, comz, molmass) = self.comprep(mol, n, atype, atommass, num_timesteps)
                (comx, comy, comz, count) = self.calccom(comx, comy, comz, x, y, z, mol, atype, atommass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2, n, count, nummol)
            linecache.clearcache()
        self.saveCOM(comx, comy, comz)
        return (comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2)
        
    def getnum(self,trjfilename):
        # uses the trjectory file and returns the number of lines and the number of atoms
        n = int(linecache.getline(trjfilename[0],4))
        num_timesteps=1
        num_lines=[]
        for i in range(0,len(trjfilename)):
            num_lines.append(int(sum(1 for line in open(trjfilename[i]))))
            num_timesteps += int(num_lines[i] / (n+9))-1
        line = [10 for x in trjfilename]
        for j in range(1,len(trjfilename)):
            line[j] += n+9
        count = 0
        return (num_lines, n, num_timesteps, count, line)
        
    def getdimensions(self,trjfilename):
        # uses trjectory file to get the length of box sides
        xbounds = linecache.getline(trjfilename, 6)
        xbounds = xbounds.split()
        ybounds = linecache.getline(trjfilename, 7)
        ybounds = ybounds.split()
        zbounds = linecache.getline(trjfilename, 8)
        zbounds = zbounds.split()
        Lx = float(xbounds[1])-float(xbounds[0])
        Lx2 = Lx/2
        Ly = float(ybounds[1])-float(ybounds[0])
        Ly2 = Ly/2
        Lz = float(zbounds[1])-float(zbounds[0])
        Lz2 = Lz/2 
        return (Lx, Lx2, Ly, Ly2, Lz, Lz2)
        
    def createarrays(self,n):
        #creates numpy arrays for data reading
        x = np.zeros(n)
        y = np.zeros(n)
        z = np.zeros(n)
        mol = np.zeros(n)
        atype = np.zeros(n)
        return (x,y,z,mol,atype)

    def getcolumns(self,trjfilename):
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
        return (xcol, ycol, zcol, molcol, typecol)
        
    def getmass(self, datfilename):
        # returns a dictionary of the mass of each atom type
        linenum=5
        atommass = {}
        foundmass= False
        readingmasses = True
        atomnum = 1
        
        while foundmass == False:
            line = linecache.getline(datfilename, linenum)
            line = line.split()
            
            if len(line) > 0:
                if line[0] == 'Masses':
                    foundmass = True
                    linenum += 2
                else:
                    linenum += 1
            else:
                linenum += 1
                
        while readingmasses == True:
            line = linecache.getline(datfilename, linenum)
            line = line.split()
            if len(line) > 0:
                if int(line[0]) == atomnum:
                    atommass[int(line[0])] = float(line[1])
                    atomnum += 1
                    linenum += 1
                    
                else:
                    readingmasses = False
                
            else:
                readingmasses = False
                    
        return atommass
        
    def readdata(self, trjfilename, n, line, x, y, z, mol, atype, xcol, ycol, zcol, molcol, typecol, i):
        # reads data from trjectory file into precreated arrays
        for a in range(0,n):
            inline = linecache.getline(trjfilename, line[i] + a)
            inline = inline.split()
            x[a]= inline[xcol]
            y[a]= inline[ycol]
            z[a]= inline[zcol]
            mol[a]= inline[molcol]
            atype[a]= inline[typecol]
            
        line[i] += n+9
        return (x,y,z,mol,atype,line)
        
    def comprep(self, mol, n, atype, atommass, num_timesteps):
        #creates arrays to prepare for center of mass calculations
        nummol = int(max(mol))
        comx = [[0 for x in range(nummol)]for x in range(num_timesteps)]
        comy = [[0 for x in range(nummol)]for x in range(num_timesteps)]
        comz = [[0 for x in range(nummol)]for x in range(num_timesteps)]

        molmass = np.zeros(nummol)
        for atom in range(0,n):
            molmass[mol[atom]-1] += atommass[atype[atom]]
            
        return (nummol, comx, comy, comz, molmass)
        
    def calccom(self, comx, comy, comz, x, y, z, mol, atype, atommass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2, n, count, nummol):
        #calculates the center of mass for each molecule
        amass = np.zeros(n)
        for i in range(0,n):
            amass[i] = atommass[atype[i]]
            
        (comxt, comyt, comzt)= calccomf.calccom(n, nummol, x, y, z, mol, amass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2)
        comx[count] += comxt
        comy[count] += comyt
        comz[count] += comzt
        count += 1
        print('timestep ' + str(count) + ' finished')
            
        return (comx, comy, comz, count)
        
    def saveCOM(self, comx, comy, comz):
        try:
            os.remove("COMx.dat")
        except OSError:
            pass        
        
        COMxfile = open("COMx.dat", "a")
        for timepoint in range(0,len(comx)):
            for molecule in range(0,len(comx[0])):
                COMxfile.write('      ' + str(comx[timepoint][molecule]))
            COMxfile.write('\n')
        COMxfile.close()
        
        try:
            os.remove("COMy.dat")
        except OSError:
            pass        
        
        COMyfile = open("COMy.dat", "a")
        for timepoint in range(0,len(comy)):
            for molecule in range(0,len(comy[0])):
                COMyfile.write('      ' + str(comy[timepoint][molecule]))
            COMyfile.write('\n')
        COMyfile.close()
        
        try:
            os.remove("COMz.dat")
        except OSError:
            pass        
        
        COMzfile = open("COMz.dat", "a")
        for timepoint in range(0,len(comz)):
            for molecule in range(0,len(comz[0])):
                COMzfile.write('      ' + str(comz[timepoint][molecule]))
            COMzfile.write('\n')
        COMzfile.close()