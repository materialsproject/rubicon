# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import


"""
Run tests

TODO: convert ot proper unittests
"""

import math
import time

from rubicon.analysis.lammps.utils import MeanSquareDisplacement
from rubicon.analysis.lammps.properties import CenterOfMass, \
    NernstEinsteinConductivity, Conductivity
from rubicon.analysis.lammps.radial_distribution import \
    RadialDistributionPure, SiteRadialDistribution
from rubicon.analysis.lammps.process import AtomicCharges, \
    CoordinationNumber, MoleculeData, TimeData
from rubicon.analysis.lammps.ion_pair import IonPair


__author__ = "mhumbert"


def test_site_radial_distribution():
    tic = time.time()
    ard = SiteRadialDistribution()

    ard.runatomradial('sample_files/NaSCN.lammpstrj',
                      'sample_files/data.water_1NaSCN', 'H2O', 1, 1, 'H2O', 1,
                      1,
                      firststep=0)

    toc = time.time()
    hours = math.floor((toc - tic) / 3600)
    minutes = math.floor((toc - tic - 3600 * hours) / 60)
    seconds = (toc - tic - 3600 * hours - 60 * minutes)

    print(("length of run: " + str(int(hours)) + " hours " + str(
        int(minutes)) + " minutes " + str(int(round(seconds))) + " seconds"))


def test_conductivity():
    g = AtomicCharges()
    gm = MoleculeData()
    cc = Conductivity()

    T = 350  # from lammpsio

    trjfilename = ['sample_files/NaSCN.lammpstrj']
    datfilename = 'sample_files/data.water_1NaSCN'
    logfilename = 'sample_files/mol.log'
    output = {}
    output['Conductivity'] = {}
    output['Conductivity']['units'] = 'S/m'

    (nummoltype, moltypel, moltype) = gm.getmoltype(datfilename)
    n = g.findnumatoms(datfilename)
    (molcharges, atomcharges, n) = g.getmolcharges(datfilename, n)
    output = cc.calcConductivity(molcharges, trjfilename, logfilename,
                                 datfilename, T, output)
    print((output['Conductivity']['Green_Kubo']))

def test_other():
    """
    Sample driver file to calculate the MSD and diffusivity for all
    molecules in a system as well as the center of mass radial distribution
    function for all pairs of molecules in the system. Will also integrate
    the rdf to get the coordination number and calculate the ion pair lifetime
    for the system

    Requires the following comments in the lammps data file starting
    at the third line

    # 'number' 'molecule' molecules

    where 'number' is the number of that molecule type and
    'molecule' is a name for that molecule

    Do not include blank lines in between the molecule types

    Outputs are stored in a dictionary called output to later be stored
    in JSON format
    """
    c = CenterOfMass()
    m = MeanSquareDisplacement()
    gt = TimeData()
    gm = MoleculeData()
    crd = RadialDistributionPure()
    gc = AtomicCharges()
    ne = NernstEinsteinConductivity()
    cn = CoordinationNumber()
    ip = IonPair()

    trjfile = 'rubicon/analysis/lammps/sample_files/NaSCN.lammpstrj'
    datfile = 'rubicon/analysis/lammps/sample_files/data.water_1NaSCN'
    logfile = 'rubicon/analysis/lammps/sample_files/mol.log'
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
    ip.runionpair(comx, comy, comz, Lx, Ly, Lz, moltypel, moltype, tsjump,
                  dt, output, skipframes=0)
    output = crd.runradial(datfile, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2,
                           Lz2, output, nummoltype, moltypel, moltype,
                           firststep=1)
    output = cn.calccoordinationnumber(output, nummoltype, moltypel,
                                       Lx * Ly * Lz)


if __name__ == '__main__':
    test_site_radial_distribution()
    test_conductivity()
    test_other()