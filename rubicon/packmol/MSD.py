# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:08:27 2015

@author: mhumbert
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats
import time
import copy

class MSD:
    
    def runMSD(self, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, moltype, moltypel, dt, tsjump, output):
        (comx, comy, comz) = self.unwrap(comx,comy,comz,Lx,Ly,Lz,Lx2,Ly2,Lz2)
        num_timesteps = len(comx)
        (MSDt, MSD, diffusivity) = self.gettimesteps(num_timesteps, moltypel)
        (molcheck,nummol) = self.setmolarray(moltype,moltypel)
        for i in range(0,MSDt):
            for j in range(i,i+MSDt):
                r2 = self.calcr2(comx, comy, comz, i, j)
                MSD = self.MSDadd(r2, MSD, molcheck, i, j)
        MSD = self.MSDnorm(MSD, MSDt, nummol)
        Time = self.createtime(dt, tsjump, MSDt)
        #self.writetofile(MSD, Time)
        (MSD, Time, lnMSD, lntime) = self.plotMSD(MSD, Time)
        for molecule in range(0,len(moltypel)):
            firststep = self.findlinearregion(lnMSD, lntime, dt, molecule)
            self.getdiffusivity(Time, MSD, firststep, molecule, diffusivity)
        self.append_dict(MSD, moltypel, diffusivity, output, Time)
        return output
    
    def unwrap(self, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2):
        for i in range(1,len(comx)):
            for j in range(0,len(comx[i])):
                if (comx[i][j]-comx[i-1][j])>Lx2:
                    while (comx[i][j]-comx[i-1][j])>Lx2:
                        comx[i][j] -= Lx
                elif (comx[i][j]-comx[i-1][j]) < (-Lx2):
                    while (comx[i][j]-comx[i-1][j]) < (-Lx2):
                        comx[i][j] += Lx
                    
                if (comy[i][j]-comy[i-1][j])>Ly2:
                    while (comy[i][j]-comy[i-1][j])>Ly2:
                        comy[i][j] -= Ly
                elif (comy[i][j]-comy[i-1][j]) < (-Ly2):
                    while (comy[i][j]-comy[i-1][j]) < (-Ly2):
                        comy[i][j] += Ly
                    
                if (comz[i][j]-comz[i-1][j])>Lz2:
                    while (comz[i][j]-comz[i-1][j])>Lz2:
                        comz[i][j] -= Lz
                elif (comz[i][j]-comz[i-1][j])< (-Lz2):
                    while (comz[i][j]-comz[i-1][j])< (-Lz2):
                        comz[i][j] += Lz
        return (comx, comy, comz)
        
    def gettimesteps(self, num_timesteps,moltypel):
        MSDt = int(np.floor(num_timesteps/2))
        MSD = np.zeros((len(moltypel), MSDt))
        diffusivity = []
        return (MSDt, MSD, diffusivity)
        
    def setmolarray(self, moltype, moltypel):
        molcheck = np.zeros((len(moltypel), len(moltype)))
        for i in range(0,len(moltype)):
            molcheck[moltype[i]][i]=1
        nummol = np.zeros(len(moltypel))
        for i in range(0,len(nummol)):
            nummol[i] = np.sum(molcheck[i])        
        return (molcheck,nummol)
        
    def calcr2(self, comx, comy, comz, i, j):
        r2 = (comx[j]-comx[i])**2 + (comy[j]-comy[i])**2 + (comz[j]-comz[i])**2
        
        return r2
        
    def MSDadd(self, r2, MSD, molcheck, i, j):
        for k in range(0,len(molcheck)):        
            sr2 = np.dot(r2, molcheck[k])
            MSD[k][j-i] += sr2
        return MSD
        
    def MSDnorm(self, MSD,MSDt,nummol):
        for i in range(0,len(nummol)):        
            MSD[i] /= MSDt*nummol[i]
        
        return MSD
        
    def createtime(self, dt, tsjump, MSDt):
        Time = np.arange(0,MSDt)
        Time *= dt*tsjump
        return Time
    
    def writetofile(self, MSD, Time):
        try:
            os.remove("MSD.dat")
        except OSError:
            pass        
        
        MSDfile = open("MSD.dat", "a")
        for timepoint in range(0,len(MSD[0])):
            MSDfile.write(str(Time[timepoint]))
            for molecule in range(0,len(MSD)):
                MSDfile.write('      ' + str(MSD[molecule][timepoint]))
            MSDfile.write('\n')
        MSDfile.close()
    
    def plotMSD(self, MSD, Time):
        
        lnMSD = np.log(MSD)
        lntime = np.log(Time)
        
        for i in range(0,len(MSD)):
            plt.figure()
            plt.plot(Time, MSD[i])
            #plt.savefig('MSD.png')
            #plt.show()
            
            plt.figure()
            plt.plot(lntime,lnMSD[i])
            #plt.show()
            
        return (MSD, Time, lnMSD, lntime)
        
    def findlinearregion(self, lnMSD, lntime, dt, molecule):
        timestepskip=int(500/dt)
        linearregion=True
        maxtime = len(lnMSD[0])
        numskip=1
        tolerance=0.05
        while linearregion == True:
            if numskip*timestepskip+1 > maxtime:
                firststep = maxtime-1-(numskip-1)*timestepskip
                return firststep
                linearregion= False
            else:
                t1=maxtime-1-(numskip-1)*timestepskip
                t2=maxtime-1-numskip*timestepskip
                slope = (lnMSD[molecule][t1]-lnMSD[molecule][t2])/(lntime[t1]-lntime[t2])
                if abs(slope-1.) < tolerance:
                    numskip += 1
                else:
                    firststep=t1
                    return firststep
                    linearregion = False
                
                
    def getdiffusivity(self, Time, MSD, firststep, molecule, diffusivity):
        calctime = []        
        calcMSD = []
        for i in range(firststep, len(Time)):
            calctime.append(Time[i])
            calcMSD.append(MSD[molecule][i])
        if len(calctime)==1:
            diffusivity.append('runtime not long enough')
        else:
            line = stats.linregress(calctime,calcMSD)
            slope = line[0]
            diffusivity.append(slope/600000)
            
    def append_dict(self, MSD, moltypel, diffusivity, output, Time):
        for i in range(0,len(moltypel)):
            output['MSD of {0} in angstroms^2'.format(moltypel[i])] = copy.deepcopy(MSD[i])
            output['Diffusivity of {0} in m^2/sec'.format(moltypel[i])] = copy.deepcopy(diffusivity[i])
            
        output['Time in fs'] = Time
            
        
