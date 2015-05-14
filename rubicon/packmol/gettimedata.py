# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 10:41:46 2015

@author: mhumbert
"""

import linecache

class gettimedata:
    
    def getdt(self, logfilename):
        line = 1
        foundtimestep=False
        while foundtimestep==False:
            inline = linecache.getline(logfilename, line)
            inline = inline.split()
            
            if len(inline) == 0:
                line +=1
            
            elif inline[0] == 'timestep':
                dt = float(inline[1])
                foundtimestep=True
                
            else:
                line += 1 
        
        return dt
        
    def getjump(self, trjfilename):
        n = int(linecache.getline(trjfilename,4))
        tsjump = int(linecache.getline(trjfilename,n+11))-int(linecache.getline(trjfilename,2))
        return tsjump