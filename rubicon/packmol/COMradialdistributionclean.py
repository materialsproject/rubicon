# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:07:34 2015

@author: mhumbert
"""

import numpy as np
import linecache
import matplotlib.pyplot as plt
import os
import comradial

class COMradialdistribution:
    
    def runradial(self, trjfilename, datfilename, mol1, mol2, firststep=1):
        (num_lines, n, num_timesteps, count)=self.getnum(trjfilename)
        (Lx, Lx2, Ly, Ly2, Lz, Lz2) = self.getdimensions(trjfilename)
        (maxr, binsize, numbins, line) = self.setgparam(Lx2,Ly2,Lz2,n,firststep)    
        (g,x,y,z,mol,atype) = self.createarrays(numbins,n)
        (xcol, ycol, zcol, molcol, typecol) = self.getcolumns(trjfilename)
        (linenum, nummoltype, moltypel, moltype) = self.getmoltype(datfilename)
        atommass = self.getmass(datfilename, linenum)
        (nummol1, nummol2, mol1, mol2) = self.getnummol(moltypel, nummoltype, mol1, mol2)
        while line < num_lines:
            (x,y,z,mol,atype,line) = self.readdata(trjfilename, n, line, x, y, z, mol, atype, xcol, ycol, zcol, molcol, typecol)
            if count == 0:
                (nummol, comx, comy, comz, xt, yt, zt, molmass) = self.comprep(mol, n, atype, atommass)
            (comx, comy, comz) = self.calccom(comx, comy, comz, xt, yt, zt, x, y, z, mol, atype, atommass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2, n)
            (g, count) = self.radialdistribution(g, mol1, mol2, nummol, moltype, comx, comy, comz, Lx, Ly, Lz, binsize, numbins, maxr, count)
            print('timestep ' + str(count+firststep-1) + ' of ' + str(num_timesteps-1) + ' finished')
        (radiuslist, g) = self.radialnormalization(numbins,binsize,Lx,Ly,Lz,nummol1,nummol2,count,g)
        self.plot(radiuslist, g)
        self.writetofile(numbins, radiuslist, g)
        
    def getnum(self,trjfilename):
        # uses the trjectory file and returns the number of lines and the number of atoms
        num_lines = int(sum(1 for line in open(trjfilename)))
        n = int(linecache.getline(trjfilename,4))
        num_timesteps = int(num_lines / (n+9))
        count = 0
        return (num_lines, n, num_timesteps, count)
        
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
        
    def setgparam(self,Lx2,Ly2,Lz2,n,firststep):
        # uses side lengths to set the maximum radius for box and number of bins
        #also sets the first line using data on firststep and numbre of atoms
        maxr = min(Lx2,Ly2,Lz2)
        binsize = 0.1
        numbins = int(np.ceil(maxr/binsize))
        line = 10 + int(firststep) * (n+9)
        return (maxr, binsize, numbins, line)
        
    def createarrays(self,numbins,n):
        #creates numpy arrays for data reading
        g = np.zeros(numbins)
        x = np.zeros(n)
        y = np.zeros(n)
        z = np.zeros(n)
        mol = np.zeros(n)
        atype = np.zeros(n)
        return (g,x,y,z,mol,atype)

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
        
    def getmoltype(self, datfilename):
        # determines molecule types and number of each molecule type
        #also creates a list of molecule type of each molecule
        linenum=3
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
                
        for i in range(0,len(moltypel)):
            for j in range(0,nummoltype[i]):
                moltype.append(int(i))
                
        return (linenum, nummoltype, moltypel, moltype)
        
    def getmass(self, datfilename, linenum):
        # returns a dictionary of the mass of each atom type
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
        
    def getnummol(self, moltypel, nummoltype, mol1, mol2):
        # returns number of each mpolecule type and converts the molecule type to an integer
        nummol1 = nummoltype[moltypel.index(mol1)]
        nummol2 = nummoltype[moltypel.index(mol2)]
        mol1 = int(moltypel.index(mol1))
        mol2 = int(moltypel.index(mol2))
        return (nummol1, nummol2, mol1, mol2)
        
    def readdata(self, trjfilename, n, line, x, y, z, mol, atype, xcol, ycol, zcol, molcol, typecol):
        # reads data from trjectory file into precreated arrays
        for a in range(0,n):
            inline = linecache.getline(trjfilename, line + a)
            inline = inline.split()
            x[a]= inline[xcol]
            y[a]= inline[ycol]
            z[a]= inline[zcol]
            mol[a]= inline[molcol]
            atype[a]= inline[typecol]
            
        line += n+9
            
        return (x,y,z,mol,atype,line)
        
    def comprep(self, mol, n, atype, atommass):
        #creates arrays to prepare for center of mass calculations
        nummol = int(max(mol))
        comx = np.zeros(nummol)
        comy = np.zeros(nummol)
        comz = np.zeros(nummol)
        xt = np.zeros(nummol)
        yt = np.zeros(nummol)
        zt = np.zeros(nummol)
        molmass = np.zeros(nummol)
        for atom in range(0,n):
            molmass[mol[atom]-1] += atommass[atype[atom]]
            
        return (nummol, comx, comy, comz, xt, yt, zt, molmass)
        
    def calccom(self, comx, comy, comz, xt, yt, zt, x, y, z, mol, atype, atommass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2, n):
        #calculates the center of mass for each molecule
        comx *= 0
        comy *= 0 
        comz *=0
        xt *= 0
        yt *= 0
        zt *= 0
        for atom in range(0,n):
            if xt[mol[atom]-1] == 0:
                xt[mol[atom]-1] = x[atom]
                yt[mol[atom]-1] = y[atom]
                zt[mol[atom]-1] = z[atom]
                
            else:
                if x[atom]-xt[mol[atom]-1] > Lx2:
                    x[atom] -= Lx
                    
                elif xt[mol[atom]-1]-x[atom] > Lx2:
                    x[atom] += Lx
                    
                if y[atom]-yt[mol[atom]-1] > Ly2:
                    y[atom] -= Ly
                    
                elif yt[mol[atom]-1]-y[atom] > Ly2:
                    y[atom] += Ly
                    
                if z[atom]-zt[mol[atom]-1] > Lz2:
                    z[atom] -= Lz
                    
                elif zt[mol[atom]-1]-z[atom] > Lz2:
                    z[atom] += Lz
                
            comx[mol[atom]-1] += x[atom]*atommass[atype[atom]]
            comy[mol[atom]-1] += y[atom]*atommass[atype[atom]]
            comz[mol[atom]-1] += z[atom]*atommass[atype[atom]]
            
        comx /= molmass
        comy /= molmass
        comz /= molmass
            
        return (comx, comy, comz)
        
    def radialdistribution(self, g, mol1, mol2, nummol, moltype, comx, comy, comz, Lx, Ly, Lz, binsize, numbins, maxr, count):
        #implements a FORTRAN code to calculate the number of molecules within each shell
        g1 = comradial.comradial(mol1, mol2, nummol, moltype, comx, comy, comz, Lx, Ly, Lz, binsize, numbins, maxr)
        g += g1
        print(g1[10])
        count += 1
        return (g, count)
        
    def radialnormalization(self,numbins,binsize,Lx,Ly,Lz,nummol1,nummol2,count,g):
        # normalizes g to box density
        radiuslist = (np.arange(numbins)+1)*binsize
        g *= Lx*Ly*Lz/nummol1/nummol2/4/np.pi/(radiuslist)**2/binsize/count
        
        return (radiuslist, g)
        
    def plot(self, radiuslist, g):
        # plots radial distribution function
        plt.figure()
        plt.plot(radiuslist,g)
        plt.xlabel('radius')
        plt.ylabel('g(r)')
        plt.title('Pair Distribution Function')
        #plt.savefig('Pairdist2.png')
        #plt.show()
        
    def writetofile(self, numbins, radiuslist, g):
        try:
            os.remove("gr.dat")
        except OSError:
            pass        
        
        gfile = open("gr.dat", "a")
        for radius in range(0,numbins):
            gfile.write(str(radiuslist[radius]) + "     " + str(g[radius]) + "\n")
        gfile.close()