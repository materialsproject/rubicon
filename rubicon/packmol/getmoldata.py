# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 10:36:25 2015

@author: mhumbert
"""

import linecache

class getmoldata:
    
    '''
    
             Determines molecule types and number of each molecule type and 
             creates a list of molecule type of each molecule
             
             Requires the following comments in the lammps data file starting 
            at the third line
            
            # "number" "molecule" molecules
            
            where "number" is the number of that molecule type and
            "molecule" is a name for that molecule
            
            Do not include blank lines in between the molecule types
            
    '''
    
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
                
        return (nummoltype, moltypel, moltype)
        
