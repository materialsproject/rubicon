
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

    def set_bond(self,bond):
        self.bonds=len(bond)

    def set_bond_types(self,bond_types):
        self.bonds=len(bond_types)

    def set_angle(self):
        pass

    def set_angle_types(self,angle_types):
        self.angles=len(angle_types)

    def set_dihedral(self):
        pass

    def set_dihedral_types(self,dihedral_types):
        self.dihedrals=len(dihedral_types)

    def set_improper(self):
        pass

    def set_improper_types(self,imdihedral_types):
        self.imdihedrals=len(imdihedral_types)

    def set_masses(self, atom_gaff ,gff_masses):

        lines=[]
        if gff_masses is not None:
            lines.append('Masses')
            for i, v in enumerate(gff_masses.values()):
                lines.append('{} {}'.format(i+1, v))
            print '\n'.join(lines)


    def set_pair_coeff(self, atom_gaff ,gff_vdw):

        lines=[]
        if gff_vdw is not None:
            lines.append('Pair Coeffs')
            for i, v in enumerate(gff_vdw.values()):
                lines.append('{} {} {}'.format(i+1, v[0],v[1]))
            print '\n'.join(lines)

    def set_pair_coeff(self, atom_gaff ,gff_vdw):

        lines=[]
        if gff_vdw is not None:
            lines.append('Pair Coeffs')
            for i, v in enumerate(gff_vdw.values()):
                lines.append('{} {} {}'.format(i+1, v[0],v[1]))
            print '\n'.join(lines)





    def set_dimension(self):
        pass
