# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:07:34 2015

@author: mhumbert
"""

import math
import time

from rubicon.analysis.lammps.atomradialdistributionclean import siteradialdistribution

tic = time.time()
ard = siteradialdistribution()

ard.runatomradial('sample_files/NaSCN.lammpstrj',
                  'sample_files/data.water_1NaSCN', 'H2O', 1, 1, 'H2O', 1, 1,
                  firststep=0)

toc = time.time()
hours = math.floor((toc - tic) / 3600)
minutes = math.floor((toc - tic - 3600 * hours) / 60)
seconds = (toc - tic - 3600 * hours - 60 * minutes)

print("length of run: " + str(int(hours)) + " hours " + str(
    int(minutes)) + " minutes " + str(int(round(seconds))) + " seconds")
