# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:10:12 2015

@author: mhumbert
"""


class getatomcharges:
    def findnumatoms(self, datfilename):
        datfile = open(datfilename)
        foundnumatoms = False
        for j in range(0, 3):
            datfile.readline()
        while foundnumatoms == False:
            line = datfile.readline()
            line = line.split()
            if len(line) >= 2:
                if line[1] == 'atoms':
                    n = int(line[0])
                    foundnumatoms = True
        datfile.close()
        return n

    def getmolcharges(self, datfilename, n):
        datfile = open(datfilename)
        for j in range(0, 4):
            datfile.readline()
        atomcharges = [0 for x in range(0, n)]
        mol = [0 for x in range(0, n)]
        foundatoms = False
        readingcharges = True

        while foundatoms == False:
            line = datfile.readline()
            line = line.split()

            if len(line) > 0:
                if line[0] == 'Atoms':
                    foundatoms = True
                    datfile.readline()

        while readingcharges == True:
            line = datfile.readline()
            line = line.split()
            if len(line) >= 6:
                atomcharges[int(line[0]) - 1] = float(line[3])
                mol[int(line[0]) - 1] = int(line[1])

            else:
                readingcharges = False

        nummol = max(mol)
        molcharges = [0 for x in range(0, nummol)]
        for atom in range(0, n):
            molcharges[mol[atom] - 1] += atomcharges[atom]

        datfile.close()
        return (molcharges, atomcharges, n)

    def molchargedict(self, molcharges, moltypel, moltype):
        molcharge = {}
        for molecules in range(0, len(moltypel)):
            molcharge[moltypel[molecules]] = molcharges[
                moltype.index(molecules)]
        return molcharge
