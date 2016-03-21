# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 10:41:46 2015

@author: mhumbert
"""


class gettimedata:
    '''
            Uses a lammps trajectory file and log file to determine the 
            timestep size and the trajectory print frequency
            
    '''

    def getdt(self, logfilename):
        logfile = open(logfilename)
        foundtimestep = False
        while foundtimestep == False:
            inline = logfile.readline()
            inline = inline.split()
            if len(inline) > 0:
                if inline[0] == 'timestep':
                    dt = float(inline[1])
                    foundtimestep = True
        logfile.close()

        return dt

    def getjump(self, trjfilename):
        trjfile = open(trjfilename)
        trjfile.readline()
        t1 = trjfile.readline()
        t1 = int(t1)
        trjfile.readline()
        n = int(trjfile.readline())
        for i in range(0, n + 6):
            trjfile.readline()
        t2 = int(trjfile.readline())
        tsjump = t2 - t1
        return tsjump
