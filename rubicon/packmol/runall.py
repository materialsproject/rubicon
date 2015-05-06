# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:14:34 2015

@author: mhumbert
"""

from MSD import MSD
from calcCOM import calcCOM
from gettimedata import gettimedata
from getmoldata import getmoldata
from COMradial import COMradialdistribution
import time

c = calcCOM()
m = MSD()
gt = gettimedata()
gm = getmoldata()
crd = COMradialdistribution()

trjfile=['sample_files/NaSCN.lammpstrj']
datfile='sample_files/data.water_1NaSCN'
logfile='sample_files/mol.log'
output = {}

start = time.time()
tic = time.time()
tsjump = gt.getjump(trjfile[0])
(nummoltype, moltypel, moltype) = gm.getmoltype(datfile)
dt = gt.getdt(logfile)
(comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2) = c.calcCOM(trjfile,datfile)
toc = time.time()
print('calc com ' + str((toc-tic)/60)+ ' minutes')


tic = time.time()
output = m.runMSD(comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, moltype, moltypel, dt, tsjump, output)
toc = time.time()
print('MSD '+ str((toc-tic)/60)+ ' minutes')
tic = time.time()
for i in range(0,len(moltypel)):
    for j in range(i,len(moltypel)):
        output = crd.runradial(datfile, moltypel[i], moltypel[j], comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, output, nummoltype, moltypel, moltype, firststep=1)
toc = time.time()
print('radial ' + str((toc-tic)/60)+ ' minutes')
stop = time.time()
print('total ' + str((stop-start)/60)+ ' minutes')
print output.keys()
