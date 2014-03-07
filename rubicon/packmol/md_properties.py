# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 01:16:09 2013

@author: asharma6
"""

import lammpsio
import sys
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from multiprocessing import Pool

def autocorrelate (a):
    b=np.concatenate((a,np.zeros(len(a))),axis=1)
    c= np.fft.ifft(np.fft.fft(b)*np.conjugate(np.fft.fft(b))).real
    d=c[:len(c)/2]
    d=d/(np.array(range(len(a)))+1)[::-1]
    return d[:wline]

if __name__=='__main__':
    
    """
    All the properties are evaluated based on the *properties input. Right now, pymatgen only supports viscosity (argument: viscosity) evaluation.
                Example: 'md_properties log.lammps viscosity' will return the viscosity of the system.
    """
    
    logfilename=sys.argv[1]
    properties=sys.argv[2:]
    l=lammpsio.LammpsLog(logfilename, *properties)
    l.parselog()
    wline=l.wline
    print 'Done reading the log file. Starting Calculations...'
    NCORES=8
    p=Pool(NCORES)
    
  
    if 'viscosity' in l.properties:
        a1=l.LOG['pxy']
        a2=l.LOG['pxz']
        a3=l.LOG['pyz']
        a4=l.LOG['pxx']-l.LOG['pyy']
        a5=l.LOG['pyy']-l.LOG['pzz']
        a6=l.LOG['pxx']-l.LOG['pzz']
        array_array=[a1,a2,a3,a4,a5,a6]
        pv=p.map(autocorrelate,array_array)
        pcorr = (pv[0]+pv[1]+pv[2])/6+(pv[3]+pv[4]+pv[5])/24
        
        
        visco = (scipy.integrate.cumtrapz(pcorr,l.LOG['step'][:len(pcorr)]))*l.timestep*10**-15*1000*101325.**2*l.LOG['vol'][-1]*10**-30/(1.38*10**-23*l.temp)  
        plt.plot(np.array(l.LOG['step'][:len(pcorr)-1])*l.timestep,visco)
        plt.xlabel('Time (femtoseconds)')
        plt.ylabel('Viscosity (cp)')
        plt.savefig('viscosity_parallel.png')
    
        output=open('viscosity_parallel.txt','w')
        output.write('#Time (fs), Average Pressure Correlation (atm^2), Viscosity (cp)\n')
        for line in zip(np.array(l.LOG['step'][:len(pcorr)-1])*l.timestep-l.cutoff,pcorr,visco):
          output.write(' '.join(str(x) for x in line)+'\n')
        output.close()
        print 'Viscosity Calculation Comlete!'