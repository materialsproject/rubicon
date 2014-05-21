"""
This module implements input and output processing from Lampps.
"""
import math
from pybel import Molecule
from pymatgen.packmol.packmol import PackmolRunner
from rubicon.gff.antechamberio import AntechamberRunner
from rubicon.gff.gff import Gff
from rubicon.gff.topology import TopMol

__author__ = 'navnidhirajput'


class LmpInput():
    """
    write lammps data input file
    """

    def __init__(self):
        self.lines = []

    def _set_gff_types(self,ffmol,mols_in_box):

        """
        set force field information about number of atom types, bond types etc.
        """
        lines = []
        num_dih = 0
        num_atoms_types = 0
        num_bonds_types = 0
        num_angles_types = 0
        num_dihedrals_types = 0
        num_impropers_types = 0
        lines.append("{} {}".format('LAMMPS Data File',' \n'))
        for k,v in enumerate(mols_in_box.param_list):
                lines.append("{} {} {} {} {}".format('#',v['number'],"mol",k+1,"molecule"))
        lines.append('\n')
        for gff,top in zip( ffmol.gff_list,ffmol.top_list):
            num_atoms_types += len(gff.masses)
            num_bonds_types += len(gff.bonds)
            num_angles_types += len(gff.angles)
            num_impropers_types += (len(gff.imdihedrals))
            for k, v in gff.dihedrals.iteritems():
                num_dih += len(v)
        num_dihedrals_types += num_dih

        lines.append("{} {}".format(num_atoms_types, "atom type"))
        lines.append("{} {}".format(num_bonds_types, "bond type"))
        lines.append("{} {}".format(num_angles_types, "angle type"))
        lines.append("{} {}".format((num_dihedrals_types), "dihedral type"))
        lines.append(
            "{} {} {}".format(num_impropers_types, "improper type", '\n'))

        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_top_types(self,ffmol,mols_in_box):

        """
        set force field information about number of atom types, bond types etc.
        """
        lines = []
        num_dih = 0
        num_atoms = 0
        num_bonds = 0
        num_angles = 0
        num_dihedrals = 0
        num_impropers = 0
        num_total_dih = 0
        num_mol = 0

        for k,v in enumerate(mols_in_box.param_list):
            num_mol += v['number']

        for gff,top in zip(ffmol.gff_list,ffmol.top_list):
            top._get_ff_dihedrals(gff.dihedrals,top.dihedrals,gff.atom_gaff)
            num_atoms += (len(top.atoms) * num_mol)
            num_bonds += len((top.bonds) * num_mol)
            num_angles += (len(top.angles) * num_mol)
            num_impropers += (len(top.imdihedrals) * num_mol)
            for k, v in gff.dihedrals.iteritems():
                num_dih += len(v)
            for k, v in top.topdihedralff.iteritems():
                num_total_dih += len(v[1])
        num_dihedrals += (num_total_dih * num_mol)
        lines.append("{} {}".format(num_atoms, "atoms"))
        lines.append("{} {}".format(num_bonds, "bonds"))
        lines.append("{} {}".format(num_angles, "angles"))
        lines.append("{} {}".format(num_dihedrals, "dihedrals"))
        lines.append(
            "{} {} {}".format(num_impropers, "impropers",
                              '\n'))

        self.lines.extend(lines)
        return '\n'.join(lines)


    def _set_box_dimensions(self,mols_in_box):

        """
        set force field information about number of atom types, bond types etc.
        """
        lines = []
        lines.append("{} {} {}".format(mols_in_box.param_list[0]['inside box'][0],
                                           mols_in_box.param_list[0]['inside box'][3],
                                           "xlo  xhi"))
        lines.append("{} {} {}".format(mols_in_box.param_list[0]['inside box'][1],
                                           mols_in_box.param_list[0]['inside box'][4],
                                           "ylo  yhi"))
        lines.append(
                "{} {} {} {}".format(mols_in_box.param_list[0]['inside box'][2],
                                     mols_in_box.param_list[0]['inside box'][5],
                                     "zlo  zhi", '\n'))

        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_masses(self,ffmol):

        lines = []
        num_atoms = 0
        lines.append('{} {}'.format('Masses', '\n'))
        for gff,top in zip(ffmol.gff_list,ffmol.top_list):
            if gff.masses is not None:
                    for i, v in enumerate(gff.masses.values()):
                        lines.append('{} {} {} {}'.format(num_atoms+1, v, '#',
                                                          gff.masses.keys()[i]))
                        num_atoms = num_atoms +1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)


    def _set_pair_coeffs(self,ffmol):

        lines = []
        num_atoms = 0
        i = 0

        lines.append('{} {}'.format('Pair Coeffs', '\n'))

        for gff,top in zip(ffmol.gff_list,ffmol.top_list):
            if gff.vdws:
                for i, v in enumerate(gff.vdws.values()):
                    lines.append(
                        '{} {} {} {} {}'.format(num_atoms+1, v[1],v[0] * 1.7818, '#',
                                                gff.vdws.keys()[i]))
                    num_atoms = num_atoms +1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_bond_coeffs(self,ffmol):

        lines = []
        num_atoms = 0
        lines.append('{} {}'.format('Bond Coeffs', '\n'))

        for gff,top in zip(ffmol.gff_list,ffmol.top_list):
                for i, v in enumerate(gff.bonds.values()):
                    lines.append(
                        '{} {} {} {} {} {}'.format(num_atoms+1, v[0], v[1], '#',
                                                   gff.bonds.keys()[i][0],
                                                   gff.bonds.keys()[i][1]))
                    num_atoms = num_atoms +1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_angle_coeffs(self,ffmol):

        lines = []
        num_atoms = 0

        lines.append('{} {}'.format('Angle Coeffs', '\n'))
        for gff,top in zip(ffmol.gff_list,ffmol.top_list):

                for i, v in enumerate(gff.angles.values()):
                    lines.append(
                        '{} {} {} {} {} {} {}'.format(num_atoms+1, v[0], v[1], '#',
                                                      gff.angles.keys()[i][0],
                                                      gff.angles.keys()[i][1],
                                                      gff.angles.keys()[i][2]))
                    num_atoms = num_atoms +1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_dihedral_coeffs(self,ffmol,):

        lines = []
        num_atoms = 0
        j=0
        lines.append('{} {}'.format('Dihedral Coeffs', '\n'))
        for gff,top in zip(ffmol.gff_list,ffmol.top_list):
            if gff.dihedrals is not None:
                for i, v in enumerate(gff.dihedrals.values()):
                    for func_form, d in v.iteritems():
                        j += 1
                        lines.append(
                            '{}  {}  {}  {} {} {} {} {} {}'.format(j, d[0],
                            func_form, d[1],'#',gff.dihedrals.keys()[i][0],
                            gff.dihedrals.keys()[i][1],gff.dihedrals.keys()[i]
                            [2],gff.dihedrals.keys()[i][3]))

        lines.append('\n')

        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_improper_coeffs(self,ffmol):
        lines = []
        num_atoms = 0
        lines.append('{} {}'.format('Imp Dihedral Coeffs', '\n'))

        for gff,top in zip(ffmol.gff_list,ffmol.top_list):
                for i, v in enumerate(gff.imdihedrals.values()):
                    lines.append('{} {}  {}  {} {} {} {} {} {}'.format
                        (num_atoms + 1, v[0], round(math.cos(math.degrees(v[1])),0),v[2], '#',
                                gff.imdihedrals.keys()[i][0],
                                gff.imdihedrals.keys()[i][1],
                                gff.imdihedrals.keys()[i][2],
                                gff.imdihedrals.keys()[i][3]))
                num_atoms = num_atoms +1
        self.lines.extend(lines)
        return '\n'.join(lines)



    # def set_coeff(self, gff, top, pmr):
    #     """
    #     give the force field, topology and box size information
    #     """
    #
    #     lines = []
    #     num_mol = 0
    #     for k,v in enumerate(pmr.param_list):
    #         num_mol += v['number']
    #     num_dih = 0
    #     num_total_dih = 0
    #     num_imdih = 0
    #     num_total_imdih = 0
    #
    #     top._get_ff_dihedrals(gff.dihedrals,top.dihedrals,gff.atom_gaff)
    #
    #     for k, v in gff.dihedrals.iteritems():
    #         num_dih += len(v)
    #     for k, v in top.topdihedralff.iteritems():
    #         num_total_dih += len(v[1])
    #
    #
    #     for k, v in gff.imdihedrals.iteritems():
    #         num_imdih += len(v)
    #     for k, v in top.topimdihedralff.iteritems():
    #         num_total_imdih += len(v[1])
    #
    #     if gff is not None:
    #         lines.append("{} {}".format('LAMMPS Data File',' \n'))
    #         for k,v in enumerate(pmr.param_list):
    #             lines.append("{} {} {} {} {} {}".format('#',v['number'],"mol",k+1,"molecule",' \n'))
    #         lines.append("{} {}".format(len(gff.masses), "atom type"))
    #         lines.append("{} {}".format(len(gff.bonds), "bond type"))
    #         lines.append("{} {}".format(len(gff.angles), "angle type"))
    #         lines.append("{} {}".format((num_dih), "dihedral type"))
    #         lines.append(
    #             "{} {} {}".format(len(gff.imdihedrals), "improper type", '\n'))
    #
    #         lines.append("{} {}".format((len(top.atoms) * num_mol), "atoms"))
    #         lines.append("{} {}".format(len((top.bonds) * num_mol), "bonds"))
    #         lines.append("{} {}".format((len(top.angles) * num_mol), "angles"))
    #         lines.append("{} {}".format((num_total_dih * num_mol), "dihedrals"))
    #         lines.append(
    #             "{} {} {}".format((len(top.imdihedrals) * num_mol), "impropers",
    #                               '\n'))
    #
    #         lines.append("{} {} {}".format(pmr.param_list[0]['inside box'][0],
    #                                        pmr.param_list[0]['inside box'][3],
    #                                        "xlo  xhi"))
    #         lines.append("{} {} {}".format(pmr.param_list[0]['inside box'][1],
    #                                        pmr.param_list[0]['inside box'][4],
    #                                        "ylo  yhi"))
    #         lines.append(
    #             "{} {} {} {}".format(pmr.param_list[0]['inside box'][2],
    #                                  pmr.param_list[0]['inside box'][5],
    #                                  "zlo  zhi", '\n'))
    #
    #         if gff.masses is not None:
    #             lines.append('{} {}'.format('Masses', '\n'))
    #             for i, v in enumerate(gff.masses.values()):
    #                 lines.append('{} {} {} {}'.format(i + 1, v, '#',
    #                                                   gff.masses.keys()[i]))
    #             lines.append('\n')
    #
    #         if gff.vdws:
    #             lines.append('{} {}'.format('Pair Coeffs', '\n'))
    #             for i, v in enumerate(gff.vdws.values()):
    #                 lines.append(
    #                     '{} {} {} {} {}'.format(i + 1, v[1],v[0] * 1.7818, '#',
    #                                             gff.vdws.keys()[i]))
    #             lines.append('\n')
    #
    #         if gff.bonds is not None:
    #             lines.append('{} {}'.format('Bond Coeffs', '\n'))
    #             for i, v in enumerate(gff.bonds.values()):
    #                 lines.append(
    #                     '{} {} {} {} {} {}'.format(i + 1, v[0], v[1], '#',
    #                                                gff.bonds.keys()[i][0],
    #                                                gff.bonds.keys()[i][1]))
    #             lines.append('\n')
    #
    #         if gff.angles is not None:
    #             lines.append('{} {}'.format('Angle Coeffs', '\n'))
    #             for i, v in enumerate(gff.angles.values()):
    #                 lines.append(
    #                     '{} {} {} {} {} {} {}'.format(i + 1, v[0], v[1], '#',
    #                                                   gff.angles.keys()[i][0],
    #                                                   gff.angles.keys()[i][1],
    #                                                   gff.angles.keys()[i][2]))
    #             lines.append('\n')
    #
    #         j = 0
    #         if gff.dihedrals is not None:
    #             lines.append('{} {}'.format('Dihedral Coeffs', '\n'))
    #             for i, v in enumerate(gff.dihedrals.values()):
    #                 for func_form, d in v.iteritems():
    #                     j += 1
    #                     lines.append(
    #                         '{}  {}  {}  {} {} {} {} {} {}'.format(j, d[0],
    #                         func_form, d[1],'#',gff.dihedrals.keys()[i][0],
    #                         gff.dihedrals.keys()[i][1],gff.dihedrals.keys()[i]
    #                         [2],gff.dihedrals.keys()[i][3]))
    #             lines.append('\n')
    #
    #         if gff.imdihedrals is not None:
    #             lines.append('{} {}'.format('Imp Dihedral Coeffs', '\n'))
    #             for i, v in enumerate(gff.imdihedrals.values()):
    #                 lines.append('{} {}  {}  {} {} {} {} {} {}'.format(i + 1, v[0], round(math.cos(math.degrees(v[1])),0),v[2], '#',
    #                             gff.imdihedrals.keys()[i][0],
    #                             gff.imdihedrals.keys()[i][1],
    #                             gff.imdihedrals.keys()[i][2],
    #                             gff.imdihedrals.keys()[i][3]))
    #             lines.append('\n')
    #     self.lines.extend(lines)
    #     return '\n'.join(lines)


    def set_atom(self, pmr, gff):
        """
        set the Atoms section in lammps data file
        """
        lines = []
        atom_type_index = {}
        lines.append('{}{}{}'.format('\n', 'Atoms', '\n'))
        mol_pack = pmr.run()
        self.box_mol_index = []
        for i, v in enumerate(gff.masses):
            atom_type_index[gff.masses.keys()[i]] = i + 1
        i = 0
        mol_index = 0
        #iterate over types of mol

        for mol, parm in zip(pmr.mols, pmr.param_list):
            num_atoms = len(mol)
            num_this_mol = parm['number']

            #iterate every molecule of molecule type
            for imol in range(num_this_mol):
                mol_coords = mol_pack.cart_coords[i:i + num_atoms]
                mol_index += 1
                #iterate over atoms in every molecule
                d = {}
                for k, v in enumerate(mol_coords):
                    lines.append(
                        '{}  {}  {}  {}  {}  {} {} {}  {} {}'.format(k + i + 1,
                         mol_index,atom_type_index[gff.atom_index_gaff
                            [k + 1]],v[0], v[1],v[2], '#',mol_index,
                            gff.atom_index_gaff[k + 1],
                            gff.atom_index[k + 1]))

                    d[gff.atom_index[k + 1]] = k + i + 1
                i += num_atoms

        self.box_mol_index.append(d)
        self.lines.extend(lines)
        return '\n'.join(lines)


    def set_bonds(self, pmr, gff, top):
        """
        set the Bonds section in lammps data file
        """

        lines = []
        bond_type_index = {}
        lines.append('{}{}{}'.format('\n', 'Bonds', '\n'))
        for i, v in enumerate(gff.bonds):
            bond_type_index[gff.bonds.keys()[i]] = i + 1

        i = 0
        mol_index = 0
        #iterate over types of mol
        for mol, parm in zip(pmr.mols, pmr.param_list):
            num_atoms = len(mol)
            num_this_mol = parm['number']
            #iterate every molecule of molecule type
            for imol in range(num_this_mol):
                mol_bonds = top.bonds[i:i + num_atoms]

                mol_index += 1
                #iterate over bonds in first molecule
                for k, v in enumerate(mol_bonds):

                    a = gff.atom_gaff[top.bonds[k][0]]
                    b = gff.atom_gaff[top.bonds[k][1]]

                    bond_label = tuple(sorted([a, b]))

                    lines.append(
                        '{}  {}  {}  {}  {}  {}  {}  {}'.format(i + k + 1,
                        bond_type_index[bond_label],self.box_mol_index[imol][v[0]],
                        self.box_mol_index[imol][v[1]],'#', mol_index,
                        top.bonds[k][0],top.bonds[k][1]))
                i += len(top.bonds)
        self.lines.extend(lines)
        return '\n'.join(lines)

    def set_angles(self, pmr, gff, top):
        """
        set the Angles section in lammps data file
        """
        lines = []
        angle_type_index = {}
        lines.append('{}{}{}'.format('\n', 'Angles', '\n'))
        for i, v in enumerate(gff.angles):
            angle_type_index[gff.angles.keys()[i]] = i + 1

        i = 0
        mol_index = 0
        #iterate over types of mol
        for mol, parm in zip(pmr.mols, pmr.param_list):
            num_atoms = len(mol)
            num_this_mol = parm['number']
            #iterate over first molecule
            for imol in range(num_this_mol):
                mol_angles = top.angles[i:i + num_atoms]
                mol_index += 1
                #iterate over bonds in first molecule
                for k, v in enumerate(top.angles):
                    a = gff.atom_gaff[top.angles[k][0]]
                    b = gff.atom_gaff[top.angles[k][1]]
                    c = gff.atom_gaff[top.angles[k][2]]
                    angle_label = tuple(sorted([a, b, c]))
                    lines.append('{}  {}  {}  {}  {}  {}  {}  {}  {}  {}'
                    .format(i + k + 1, angle_type_index[angle_label],
                            self.box_mol_index[0][v[0]],
                            self.box_mol_index[0][v[1]],
                            self.box_mol_index[0][v[2]], '#', mol_index,
                            top.angles[k][0], top.angles[k][1],
                            top.angles[k][2]))

                i += len(top.angles)
        self.lines.extend(lines)
        return '\n'.join(lines)

    def set_dihedrals(self, pmr, gff, top):
        """
        set the Dihedrals section in lammps data file
        """
        lines = []
        dihedral_type_index = {}
        lines.append('{}{}{}'.format('\n', 'Dihedrals', '\n'))
        for i, v in enumerate(gff.dihedrals):
            dihedral_type_index[gff.dihedrals.keys()[i]] = i + 1

        i = 0
        j = 0
        l = 0
        mol_index = 0

        top._get_ff_dihedrals(gff.dihedrals,top.dihedrals,gff.atom_gaff)
        #iterate over types of mol
        for mol, parm in zip(pmr.mols, pmr.param_list):
            num_atoms = len(mol)
            num_this_mol = parm['number']
            #iterate over first molecule
            for imol in range(num_this_mol):
                mol_index += 1
                l += 1
                #iterate over bonds in first molecule

                for k, v in top.topdihedralff.iteritems():
                #for k, v in enumerate(top.dihedrals):

                    A = k.split()[0]
                    B = k.split()[1]
                    C = k.split()[2]
                    D = k.split()[3]
                    a = gff.atom_gaff[k.split()[0]]
                    b = gff.atom_gaff[k.split()[1]]
                    c = gff.atom_gaff[k.split()[2]]
                    d = gff.atom_gaff[k.split()[3]]
                    dihedral_label = (a, b, c, d)
                    if dihedral_label[0] > dihedral_label[3]:
                        dihedral_label = tuple(reversed(list(dihedral_label)))
                    for func_form, d in v[1].iteritems():
                        j += 1
                        lines.append(
                            '{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}'
                            .format(j, dihedral_type_index[dihedral_label],
                                    self.box_mol_index[0][A],
                                    self.box_mol_index[0][B],
                                    self.box_mol_index[0][C],
                                    self.box_mol_index[0][D], '#', mol_index,
                                    top.dihedrals[l][0], top.dihedrals[l][1],
                                    top.dihedrals[l][2], top.dihedrals[l][3]))
                i += len(top.dihedrals)
        self.lines.extend(lines)
        return '\n'.join(lines)


    def set_imdihedrals(self, pmr, gff, top):
        """
        set the Improper Dihedral section in lammps data file
        """

        lines = []
        imdihedral_type_index = {}
        lines.append('{}{}{}'.format('\n', 'Impropers', '\n'))
        for i, v in enumerate(gff.imdihedrals):
            imdihedral_type_index[gff.imdihedrals.keys()[i]] = i + 1
        i = 0
        j = 0
        mol_index = 0
        #iterate over types of mol
        for mol, parm in zip(pmr.mols, pmr.param_list):
            num_atoms = len(mol)
            num_this_mol = parm['number']
            #iterate over first molecule
            for imol in range(num_this_mol):
                mol_index += 1
                #iterate over improper dihedrals in first molecule

                for k, v in enumerate(top.imdihedrals):

                    j += 1
                    a = gff.atom_gaff[top.imdihedrals[k][0]]
                    b = gff.atom_gaff[top.imdihedrals[k][1]]
                    c = gff.atom_gaff[top.imdihedrals[k][2]]
                    d = gff.atom_gaff[top.imdihedrals[k][3]]
                    imdihedral_label = tuple([a, b, c, d])
                    lines.append(
                        '{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}'
                        .format(j, imdihedral_type_index[imdihedral_label],
                                self.box_mol_index[0][v[0]],
                                self.box_mol_index[0][v[1]],
                                self.box_mol_index[0][v[2]],
                                self.box_mol_index[0][v[3]], '#', mol_index,
                                top.imdihedrals[k][0], top.imdihedrals[k][1],
                                top.imdihedrals[k][2], top.imdihedrals[k][3]))
                i += num_atoms
        self.lines.extend(lines)
        return '\n'.join(lines)

    def run(self, mols, mols_in_box, ffmol):

       my_lammps_list = []
       my_lampps = LmpInput()
       my_lammps_list.append(my_lampps._set_gff_types(ffmol, mols_in_box))
       my_lammps_list.append(my_lampps._set_top_types(ffmol, mols_in_box))
       my_lammps_list.append(my_lampps._set_box_dimensions(mols_in_box))
       my_lammps_list.append(my_lampps._set_masses(ffmol))
       my_lammps_list.append(my_lampps._set_pair_coeffs(ffmol))
       my_lammps_list.append(my_lampps._set_bond_coeffs(ffmol))
       my_lammps_list.append(my_lampps._set_angle_coeffs(ffmol))
       my_lammps_list.append(my_lampps._set_dihedral_coeffs(ffmol))
       my_lammps_list.append(my_lampps._set_improper_coeffs(ffmol))
       #my_lammps_list.append(my_lampps.set_atom(mols_in_box,ffmol))
       '''
       for mol, gff,top in zip(mols, ffmol.gff_list,ffmol.top_list):
            my_lampps = LmpInput()
            my_lammps_list.append(my_lampps.set_coeff(gff, top, mols_in_box))
            my_lammps_list.append(my_lampps.set_atom(mols_in_box, gff))
            my_lammps_list.append(my_lampps.set_bonds(mols_in_box, gff, top))
            my_lammps_list.append(my_lampps.set_angles(mols_in_box, gff, top))
            my_lammps_list.append(
                my_lampps.set_dihedrals(mols_in_box, gff, top))
            my_lammps_list.append(
                my_lampps.set_imdihedrals(mols_in_box, gff, top))
                '''

       return '\n'.join(my_lammps_list)



