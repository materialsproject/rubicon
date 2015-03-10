# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:07:34 2015

@author: mhumbert
"""

import time
import math
from COMradialdistributionclean import COMradialdistribution

tic = time.time()
crd = COMradialdistribution()

crd.runradial('sample_files/NaSCN.lammpstrj','sample_files/data.water_1NaSCN','H2O','H2O', firststep=0)

toc = time.time()
hours = math.floor((toc - tic)/3600)
minutes = math.floor((toc-tic-3600*hours)/60)
seconds = (toc-tic-3600*hours - 60*minutes)

print("length of run: " + str(int(hours)) + " hours " + str(int(minutes)) + " minutes " + str(int(round(seconds))) + " seconds")