
"""
This module implements input and output processing from Lampps.
"""

import copy
import re
import numpy as np
from string import Template
from pymatgen.core.structure import Molecule
from pymatgen.core.units import Energy
from pymatgen.serializers.json_coders import MSONable
from pymatgen.util.coord_utils import get_angle


__author__ = 'navnidhirajput'

class LMPInput():
    """
    write input file for lammps
    """
    def __init__(self):
        pass

    def write_data_file(self,basename=None):
        """Write LAMMPS data file using the specified file basename. Return
        complete filename."""

        with open("mol.data",'w') as f:
            f.write()


    def get_data_file(self):
        """Organize, format, and return a LAMMPS data (forcefield) file.
        """

        # create LAMMPS data file
        pass

    def set_atom(self):
        pass

    #def set_atom_types(self,atom_types):
    #    self.atoms=len(atom_types)



    def set_coeff(self,gff,top):

        lines=[]
        if gff is not None:
            lines.append('LAMMPS Data File')
            lines.append(' \n')
            lines.append("{} {}".format(len(gff.bonds),"atom type"))
            lines.append("{} {}".format(len(gff.bonds),"bond type"))
            lines.append("{} {}".format(len(gff.angles),"angle type"))
            lines.append("{} {}".format(len(gff.dihedrals),"dihedral type"))
            lines.append("{} {}".format(len(gff.imdihedrals),"improper dihedral type"))
            lines.append('\n')

            lines.append("{} {}".format(len(top.bonds),"atoms"))
            lines.append("{} {}".format(len(top.angles),"bonds"))
            lines.append("{} {}".format(len(top.dihedrals),"angles"))
            lines.append("{} {}".format(len(top.imdihedrals),"dihedrals"))
            lines.append('\n')

            if gff.masses is not None:
                lines.append('Masses')
                lines.append('\n')
                for i, v in enumerate(gff.masses.values()):
                    lines.append('{} {}'.format(i+1, v))
                lines.append('\n')

            if gff.vdws:
                lines.append('Pair Coeffs')
                lines.append('\n')
                for i, v in enumerate(gff.vdws.values()):
                    lines.append('{} {} {}'.format(i+1, v[0],v[1]))
                lines.append('\n')


            if gff.bonds is not None:
                lines.append('Bond Coeffs')
                lines.append('\n')
                for i, v in enumerate(gff.bonds.values()):
                    lines.append('{} {} {}'.format(i+1, v[0],v[1]))
                lines.append('\n')


            if gff.angles is not None:
                lines.append('Angle Coeffs')
                lines.append('\n')
                for i, v in enumerate(gff.angles.values()):
                    lines.append('{} {} {}'.format(i+1, v[0],v[1]))
                lines.append('\n')


            if gff.dihedrals is not None:
                lines.append('Dihedral Coeffs')
                lines.append('\n')
                for i, v in enumerate(gff.dihedrals.values()):
                    lines.append('{} {} {}'.format(i+1, v[0],v[1]))
                lines.append('\n')

            if gff.imdihedrals is not None:
                lines.append('Imp Dihedral Coeffs')
                lines.append('\n')
                for i, v in enumerate(gff.imdihedrals.values()):
                    lines.append('{} {} {}'.format(i+1, v[0],v[1]))
                lines.append('\n')

        return '\n'.join(lines)


    def set_atom(self):

        lines=[]
        lines.append('Atoms')

        return '\n'.join(lines)






    def set_dimension(self):
        pass
