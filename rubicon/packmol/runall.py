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

if __name__ == '__main__':
    from MSD import MSD
    from calcCOM import calcCOM
    from gettimedata import gettimedata
    from getmoldata import getmoldata
    from COMradialnofort import COMradialdistribution
    from getatomcharges import getatomcharges
    from calcNEconductivity import calcNEconductivity
    from getcoordinationnumber import getcoordinationnumber
    
    c = calcCOM()
    m = MSD()
    gt = gettimedata()
    gm = getmoldata()
    crd = COMradialdistribution()
    gc = getatomcharges()
    ne = calcNEconductivity()
    cn = getcoordinationnumber()

    trjfile = 'sample_files/NaSCN.lammpstrj'
    datfile = 'sample_files/data.water_1NaSCN'
    logfile = 'sample_files/mol.log'
    output = {}
    output['RDF'] = {}
    output['RDF']['units'] = 'unitless and angstroms'
    output['Conductivity'] = {}
    output['Conductivity']['units'] = 'S/m'
    T = 298  # get from lammpsio

    tsjump = gt.getjump(trjfile)
    (nummoltype, moltypel, moltype) = gm.getmoltype(datfile)
    dt = gt.getdt(logfile)
    n = gc.findnumatoms(datfile)
    (molcharges, atomcharges, n) = gc.getmolcharges(datfile, n)
    molcharge = gc.molchargedict(molcharges, moltypel, moltype)
    (comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2) = c.calcCOM([trjfile],
                                                              datfile)

    output = m.runMSD(comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, moltype,
                      moltypel, dt, tsjump, output)

    output = ne.calcNEconductivity(output, molcharge, Lx, Ly, Lz, nummoltype,
                                   moltypel, T)
    # print(output)
    # print('Conductivity Finshed')
    output = crd.runradial(datfile, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2,
                           Lz2, output, nummoltype, moltypel, moltype,
                           firststep=1)


    # output = cn.calccoordinationnumber(output,nummoltype,moltypel,Lx*Ly*Lz)
    # outputfile=open('test.json', 'w')
    # json.dump(output,outputfile,indent=4)
    # outputfile.close()
