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
    
    def runradial(self, datfilename, mol1, mol2, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, output, firststep=1):
        (maxr, binsize, numbins, count, g) = self.setgparam(Lx2,Ly2,Lz2,firststep)    
        (linenum, nummoltype, moltypel, moltype) = self.getmoltype(datfilename)
        (nummol1, nummol2, mol1, mol2) = self.getnummol(moltypel, nummoltype, mol1, mol2)
        while count < len(comx):
            (count) = self.radialdistribution(g, mol1, mol2, len(comx[1]), moltype, comx, comy, comz, Lx, Ly, Lz, binsize, numbins, maxr, count)
            print('timestep ' + str(count-firststep) + ' of ' + str(len(comx)-firststep) + ' finished')
        (radiuslist) = self.radialnormalization(numbins,binsize,Lx,Ly,Lz,nummol1,nummol2,count,g,firststep)
        #self.plot(radiuslist, g)
        #self.writetofile(numbins, radiuslist, g)
        self.append_dict(radiuslist, g, output, mol1, mol2, moltypel)
        return output
        
    def setgparam(self,Lx2,Ly2,Lz2,firststep):
        # uses side lengths to set the maximum radius for box and number of bins
        #also sets the first line using data on firststep and number of atoms
        maxr = min(Lx2,Ly2,Lz2)
        binsize = 0.1
        numbins = int(np.ceil(maxr/binsize))
        count = firststep
        g = np.zeros(numbins)
        return (maxr, binsize, numbins, count, g)
        
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
        
    def getnummol(self, moltypel, nummoltype, mol1, mol2):
        # returns number of each mpolecule type and converts the molecule type to an integer
        nummol1 = nummoltype[moltypel.index(mol1)]
        nummol2 = nummoltype[moltypel.index(mol2)]
        mol1 = int(moltypel.index(mol1))
        mol2 = int(moltypel.index(mol2))
        return (nummol1, nummol2, mol1, mol2)
        
    def radialdistribution(self, g, mol1, mol2, nummol, moltype, comx, comy, comz, Lx, Ly, Lz, binsize, numbins, maxr, count):
        #implements a FORTRAN code to calculate the number of molecules within each shell
        
        g1 = comradial.comradial(mol1, mol2, nummol, moltype, comx[count], comy[count], comz[count], Lx, Ly, Lz, binsize, numbins, maxr)
        g += g1

        count += 1
        return (count)
        
    def radialnormalization(self,numbins,binsize,Lx,Ly,Lz,nummol1,nummol2,count,g, firststep):
        # normalizes g to box density
        print('start normalization')
        radiuslist = (np.arange(numbins)+1)*binsize
        g *= Lx*Ly*Lz/nummol1/nummol2/4/np.pi/(radiuslist)**2/binsize/(count-firststep)
        print('end normalization')
        return (radiuslist)
        
    def plot(self, radiuslist, g):
        # plots radial distribution function
        print('begin plot')
        plt.plot(radiuslist,g)
        plt.xlabel('radius')
        plt.ylabel('g(r)')
        plt.title('Pair Distribution Function')
        plt.savefig('Pairdist2.png')
        #plt.show()
        print('end plot')        
        
    def writetofile(self, numbins, radiuslist, g):
        try:
            os.remove("gr2.dat")
        except OSError:
            pass        
        
        gfile = open("gr2.dat", "a")
        for radius in range(0,numbins):
            gfile.write(str(radiuslist[radius]) + "     " + str(g[radius]) + "\n")
        gfile.close()
        
    def append_dict(self, radiuslist, g, output, mol1, mol2, moltypel):
        output['g(r) for {0} and {1}'.format(moltypel[mol1], moltypel[mol2])] = g
        output['r list for g(r) in angstroms'] = radiuslist