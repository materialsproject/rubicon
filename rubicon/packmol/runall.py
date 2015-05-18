# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:14:34 2015

@author: mhumbert
"""
'''

        Sample driver file to calculate the MSD and diffusivity for all 
        molecules in a system as well as the center of mass radial distribution
        function for all pairs of molecules in the system
        
        Requires the following comments in the lammps data file starting 
        at the third line
            
        # "number" "molecule" molecules
            
        where "number" is the number of that molecule type and
        "molecule" is a name for that molecule
            
        Do not include blank lines in between the molecule types
        
        Outputs are stored in a dictionary called output to later be stored
        in JSON format

'''

if __name__ ==  '__main__':

    from MSD import MSD
    from calcCOM import calcCOM
    from gettimedata import gettimedata
    from getmoldata import getmoldata
    from COMradial import COMradialdistribution
    
    c = calcCOM()
    m = MSD()
    gt = gettimedata()
    gm = getmoldata()
    crd = COMradialdistribution()
    
    trjfile=['sample_files/NaSCN.lammpstrj']
    datfile='sample_files/data.water_1NaSCN'
    logfile='sample_files/mol.log'
    output = {}
    output['RDF'] = {}
    output['RDF']['units'] = 'unitless and angstroms'
    
    tsjump = gt.getjump(trjfile[0])
    (nummoltype, moltypel, moltype) = gm.getmoltype(datfile)
    dt = gt.getdt(logfile)
    (comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2) = c.calcCOM(trjfile,datfile)
    
    
    output = m.runMSD(comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, moltype, moltypel, dt, tsjump, output)
    for i in range(0,len(moltypel)):
        for j in range(i,len(moltypel)):
            output = crd.runradial(datfile, moltypel[i], moltypel[j], comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, output, nummoltype, moltypel, moltype, firststep=1)
