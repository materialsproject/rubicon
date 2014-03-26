
"""
This module implements input and output processing from Lampps.
"""
from collections import defaultdict

import copy
import os
import re
import numpy as np
from string import Template
from pymatgen import read_mol
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
        #self.ant_atom_index={}
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



    def set_coeff(self,gff,top,pmr,my_ant):

        lines=[]
        num_mol=pmr.param_list[0]['number']+pmr.param_list[1]['number']
        num_dih=0
        num_total_dih=0
        for k,v in gff.dihedrals.iteritems():
                num_dih+=len(v)
        for k,v in my_ant.topdihedralFF.iteritems():
            num_total_dih+=len(v[1])

        if gff is not None:
            lines.append('LAMMPS Data File')
            lines.append(' \n')
            lines.append("{} {}".format(len(gff.masses),"atom type"))
            lines.append("{} {}".format(len(gff.bonds),"bond type"))
            lines.append("{} {}".format(len(gff.angles),"angle type"))

            lines.append("{} {}".format((num_dih),"dihedral type"))
            lines.append("{} {} {}".format(len(gff.imdihedrals),"improper type",'\n'))


            lines.append("{} {}".format((len(top.atoms)*num_mol),"atoms"))
            lines.append("{} {}".format(len((top.bonds)*num_mol),"bonds"))
            lines.append("{} {}".format((len(top.angles)*num_mol),"angles"))
            lines.append("{} {}".format((num_total_dih*num_mol),"dihedrals"))
            lines.append("{} {} {}".format((len(top.imdihedrals)*num_mol),"impropers",'\n'))


            lines.append("{} {} {}".format(pmr.param_list[0]['inside box'][0],pmr.param_list[0]['inside box'][3],"xlo  xhi"))
            lines.append("{} {} {}".format(pmr.param_list[0]['inside box'][1],pmr.param_list[0]['inside box'][4],"ylo  yhi"))
            lines.append("{} {} {} {}".format(pmr.param_list[0]['inside box'][2],pmr.param_list[0]['inside box'][5],"zlo  zhi",'\n'))



            if gff.masses is not None:
                lines.append('{} {}'.format('Masses','\n'))
                for i, v in enumerate(gff.masses.values()):
                    lines.append('{} {} {} {}'.format(i+1, v,'#',gff.masses.keys()[i]))
                lines.append('\n')

            if gff.vdws:
                lines.append('{} {}'.format('Pair Coeffs','\n'))
                for i, v in enumerate(gff.vdws.values()):
                    lines.append('{} {} {} {} {}'.format(i+1, v[0],v[1],'#',gff.vdws.keys()[i]))
                lines.append('\n')


            if gff.bonds is not None:
                lines.append('{} {}'.format('Bond Coeffs','\n'))
                for i, v in enumerate(gff.bonds.values()):
                    lines.append('{} {} {} {} {} {}'.format(i+1, v[0],v[1],'#',gff.bonds.keys()[i][0],gff.bonds.keys()[i][1]))
                lines.append('\n')


            if gff.angles is not None:
                lines.append('{} {}'.format('Angle Coeffs','\n'))
                for i, v in enumerate(gff.angles.values()):
                    lines.append('{} {} {} {} {} {} {}'.format(i+1, v[0],v[1],'#',gff.angles.keys()[i][0],gff.angles.keys()[i][1],gff.angles.keys()[i][2]))
                lines.append('\n')

            j=0
            if gff.dihedrals is not None:
                lines.append('{} {}'.format('Dihedral Coeffs','\n'))
                for i, v in enumerate(gff.dihedrals.values()):
                    for func_form,d in v.iteritems():
                        j+=1
                        lines.append('{} {} {} {} {} {} {} {} {}'.format(j,d[0],func_form,d[1],'#',gff.dihedrals.keys()[i][0],gff.dihedrals.keys()[i][1],gff.dihedrals.keys()[i][2],gff.dihedrals.keys()[i][3]))
                lines.append('\n')

            if gff.imdihedrals is not None:
                lines.append('{} {}'.format('Imp Dihedral Coeffs','\n'))
                for i, v in enumerate(gff.imdihedrals.values()):
                    lines.append('{} {} {} {} {}'.format(i+1, v[0],v[1],'#',gff.imdihedrals.keys()[i][0],gff.imdihedrals.keys()[i][1],gff.imdihedrals.keys()[i][2],gff.imdihedrals.keys()[i][3]))
                lines.append('\n')

        return '\n'.join(lines)


    def set_atom(self,filename,pmr,gff,ac):
        """
        TODO
        """
        lines=[]
        atom_type_index={}
        lines.append('{}{}{}'.format('\n','Atoms','\n'))
        s=read_mol(filename)
        self.box_mol_index=[]
        #self.box_atom_index=[]
        for i, v in enumerate(gff.masses):
            atom_type_index[gff.masses.keys()[i]]=i+1
        i = 0
        mol_index = 0
        #iterate over types of mol
        for mol, parm in zip(pmr.mols,pmr.param_list):
            num_atoms=len(mol)
            num_this_mol=parm['number']
            #iterate every molecule of molecule type
            for imol in range(num_this_mol):
                mol_coords = s.cart_coords[i:i+num_atoms]
                mol_index += 1
                #iterate over atoms in every molecule
                d={}
                for k, v in enumerate(mol_coords):
                    lines.append('{}  {}  {}  {}  {}  {} {} {}  {} {}'.format(k+i+1,mol_index,atom_type_index[ac.atom_index_gaff[k+1]],v[0],v[1],v[2],'#',mol_index,ac.atom_index_gaff[k+1],ac.atom_index[k+1]))
                    d[ac.atom_index[k+1]]=k+i+1
                i += num_atoms
                self.box_mol_index.append(d)


        return '\n'.join(lines)


    def set_bonds(self,pmr,gff,ac,top):
        """
        TODO
        """

        lines=[]
        bond_type_index={}
        lines.append('{}{}{}'.format('\n','Bonds','\n'))
        for i, v in enumerate(gff.bonds):
            bond_type_index[gff.bonds.keys()[i]]=i+1

        i = 0
        mol_index = 0
        #iterate over types of mol
        for mol, parm in zip(pmr.mols,pmr.param_list):
            num_atoms=len(mol)
            num_this_mol=parm['number']
            #iterate over first molecule
            for imol in range(num_this_mol):
                mol_bonds = top.bonds[i:i+num_atoms]
                mol_index += 1
                #iterate over bonds in first molecule
                for k, v in enumerate(top.bonds):
                    a=ac.atom_gaff[top.bonds[k][0]]
                    b = ac.atom_gaff[top.bonds[k][1]]
                    bond_label=tuple(sorted([a,b]))
                    lines.append('{}  {}  {}  {}  {}  {}  {}  {}'.format(i+k+1 ,bond_type_index[bond_label],self.box_mol_index[0][v[0]],self.box_mol_index[0][v[1]],'#',mol_index,top.bonds[k][0],top.bonds[k][1]))
                i += len(top.bonds)
        return '\n'.join(lines)

    def set_angles(self,pmr,gff,ac,top):
        """
        TODO
        """
        lines=[]
        angle_type_index={}
        lines.append('{}{}{}'.format('\n','Angles','\n'))
        for i, v in enumerate(gff.angles):
            angle_type_index[gff.angles.keys()[i]]=i+1

        i = 0
        mol_index = 0
        #iterate over types of mol
        for mol, parm in zip(pmr.mols,pmr.param_list):
            num_atoms=len(mol)
            num_this_mol=parm['number']
            #iterate over first molecule
            for imol in range(num_this_mol):
                mol_angles = top.angles[i:i+num_atoms]
                mol_index += 1
                #iterate over bonds in first molecule
                for k, v in enumerate(top.angles):
                    a = ac.atom_gaff[top.angles[k][0]]
                    b = ac.atom_gaff[top.angles[k][1]]
                    c = ac.atom_gaff[top.angles[k][2]]
                    angle_label=tuple(sorted([a,b,c]))
                    lines.append('{}  {}  {}  {}  {}  {}  {}  {}  {}  {}'
                    .format(i+k+1 ,angle_type_index[angle_label],self.box_mol_index[0][v[0]],self.box_mol_index[0][v[1]],self.box_mol_index[0][v[2]],'#',mol_index,top.angles[k][0],top.angles[k][1],top.angles[k][2]))

                i += len(top.angles)
        return '\n'.join(lines)

    def set_dihedrals(self,pmr,gff,ac,top,my_ant):
        """
        TODO
        """
        lines=[]
        dihedral_type_index={}
        lines.append('{}{}{}'.format('\n','Dihedrals','\n'))
        for i, v in enumerate(gff.dihedrals):
            dihedral_type_index[gff.dihedrals.keys()[i]]=i+1

        i = 0
        j=0
        l=0
        mol_index = 0
        #iterate over types of mol
        for mol, parm in zip(pmr.mols,pmr.param_list):
            num_atoms=len(mol)
            num_this_mol=parm['number']
            #iterate over first molecule
            for imol in range(num_this_mol):
                mol_bonds = top.dihedrals[i:i+num_atoms]
                mol_index += 1
                l+=1
                #iterate over bonds in first molecule
                for k, v in my_ant.topdihedralFF.iteritems():
                    A=k.split()[0]
                    B=k.split()[1]
                    C=k.split()[2]
                    D=k.split()[3]
                    a = ac.atom_gaff[k.split()[0]]
                    b = ac.atom_gaff[k.split()[1]]
                    c = ac.atom_gaff[k.split()[2]]
                    d = ac.atom_gaff[k.split()[3]]
                    dihedral_label= (a,b,c,d)
                    if dihedral_label[0]>dihedral_label[3]:
                        dihedral_label=tuple(reversed(list(dihedral_label)))
                    for func_form,d in v[1].iteritems():
                        j+=1
                        lines.append('{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}'
                    .format(j ,dihedral_type_index[dihedral_label],self.box_mol_index[0][A],self.box_mol_index[0][B],self.box_mol_index[0][C],self.box_mol_index[0][D],'#',mol_index,top.dihedrals[l][0],top.dihedrals[l][1],top.dihedrals[l][2],top.dihedrals[l][3]))
                i += len(top.dihedrals)
        return '\n'.join(lines)


    def set_imdihedrals(self,pmr,gff,ac,top):
        """
        TODO
        """

        lines=[]
        imdihedral_type_index={}
        lines.append('{}{}{}'.format('\n','Impropers','\n'))
        for i, v in enumerate(gff.imdihedrals):
            imdihedral_type_index[gff.imdihedrals.keys()[i]]=i+1
        i = 0
        j=0
        mol_index = 0
        #iterate over types of mol
        for mol, parm in zip(pmr.mols,pmr.param_list):
            num_atoms=len(mol)
            num_this_mol=parm['number']
            #iterate over first molecule
            for imol in range(num_this_mol):
                mol_bonds = top.imdihedrals[i:i+num_atoms]
                mol_index += 1
                #iterate over bonds in first molecule
                for k, v in enumerate(top.imdihedrals):
                    j+=1
                    a=ac.atom_gaff[top.imdihedrals[k][0]]
                    b = ac.atom_gaff[top.imdihedrals[k][1]]
                    c = ac.atom_gaff[top.imdihedrals[k][2]]
                    d = ac.atom_gaff[top.imdihedrals[k][3]]
                    imdihedral_label=tuple([a,b,c,d])
                    lines.append('{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}'
                    .format(j,imdihedral_type_index[imdihedral_label],self.box_mol_index[0][v[0]],self.box_mol_index[0][v[1]],self.box_mol_index[0][v[2]],self.box_mol_index[0][v[3]],'#',mol_index,top.imdihedrals[k][0],top.imdihedrals[k][1],top.imdihedrals[k][2],top.imdihedrals[k][3]))
                i += num_atoms
        return '\n'.join(lines)
