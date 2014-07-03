"""
This module implements input and output processing from Lampps.
"""
import math
__author__ = 'navnidhirajput'


class LmpInput():
    """
    write lammps data input file
    """

    def __init__(self, ffmol, mols_in_box,):
        self.lines = []
        self.mols_in_box = mols_in_box
        self.ffmol = ffmol

    def _set_gff_types(self, ffmol, mols_in_box):

        """
        set force field information about number of atom types, bond types etc.
        """
        lines = []
        num_dih = 0
        num_atoms_types = 0
        num_atoms_types_sorted = 0
        num_bonds_types = 0
        num_angles_types = 0
        num_dihedrals_types = 0
        num_impropers_types = 0
        atom_type_list = []
        bond_type_list = []
        angle_type_list = []
        dihedral_type_list = []
        improper_type_list = []
        lines.append('LAMMPS Data File\n')
        for k, v in enumerate(mols_in_box.param_list):
            lines.append("{} {} {} {} {}".format('#', v['number'], "mol", k + 1,
                                                 "molecule"))
        lines.append('\n')
        for gff, top in zip(ffmol.gff_list, ffmol.top_list):
            atom_type_list.extend(gff.masses.keys())
            bond_type_list.extend(gff.bonds.keys())
            angle_type_list.extend(gff.angles.keys())
            dihedral_type_list.extend(gff.dihedrals.keys())
            improper_type_list.extend(gff.imdihedrals.keys())
            num_atoms_types = len(set(atom_type_list))
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
            "{} {}{}".format(num_impropers_types, "improper type",'\n'))
        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_top_types(self, ffmol, mols_in_box):
        """
        set force field information about number of atom types, bond types etc.
        """
        lines = []
        num_atoms = 0
        num_bonds = 0
        num_angles = 0
        num_dihedrals = 0
        num_impropers = 0
        num_total_dih = 0
        num_mol = 0
        for gff, top, v, in zip(ffmol.gff_list, ffmol.top_list,
                                mols_in_box.param_list):
            top._get_ff_dihedrals(gff.dihedrals, top.dihedrals, gff.atom_gaff)
            num_mol = v['number']
            num_atoms += (len(top.atoms)) * num_mol
            num_bonds += len((top.bonds) * num_mol)
            num_angles += (len(top.angles) * num_mol)
            num_impropers += (len(top.imdihedrals) * num_mol)

            for k, n in top.topdihedralff.iteritems():
                num_total_dih += len(n[1])
        num_dihedrals += (num_total_dih * num_mol)
        lines.append("{} {}".format(num_atoms, "atoms"))
        lines.append("{} {}".format(num_bonds, "bonds"))
        lines.append("{} {}".format(num_angles, "angles"))
        lines.append("{} {}".format(num_dihedrals, "dihedrals"))
        lines.append(
            "{} {}{}".format(num_impropers, "impropers",
                              '\n'))
        self.lines.extend(lines)
        return '\n'.join(lines)


    def _set_box_dimensions(self, mols_in_box):
        """
        set force field information about number of atom types, bond types etc.
        """
        lines = []
        lines.append(
            "{} {} {}".format(mols_in_box.param_list[0]['inside box'][0],
                              mols_in_box.param_list[0]['inside box'][3],
                              "xlo  xhi"))
        lines.append(
            "{} {} {}".format(mols_in_box.param_list[0]['inside box'][1],
                              mols_in_box.param_list[0]['inside box'][4],
                              "ylo  yhi"))
        lines.append(
            "{} {} {}{}".format(mols_in_box.param_list[0]['inside box'][2],
                                 mols_in_box.param_list[0]['inside box'][5],
                                 "zlo  zhi", '\n'))

        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_masses(self, ffmol, mols_in_box):

        lines = []
        num_atoms = 0
        mol_index = 1
        element_list = []
        #lines.append('{} {}'.format('Masses', '\n'))
        lines.append('Masses\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.masses is not None:
                for i, v in enumerate(gff.masses.values()):
                    if gff.masses.keys()[i] in element_list:
                        continue
                    element_list.append(gff.masses.keys()[i])
                    lines.append('{} {} {} {} {} {}'.format(num_atoms + 1, v, '#',
                                                         gff.masses.keys()[i],
                                                         mol_index,
                                            mol.site_properties["mol_name"][0]))
                    num_atoms = num_atoms + 1
            mol_index += 1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)


    def _set_pair_coeffs(self, ffmol, mols_in_box):

        lines = []
        num_atoms = 0
        mol_index = 1
        element_list = []
        lines.append('Pair Coeffs\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.vdws:
                for i, v in enumerate(gff.vdws.values()):
                    if gff.vdws.keys()[i] in element_list:
                        continue
                    element_list.append(gff.vdws.keys()[i])
                    lines.append(
                        '{} {} {} {} {} {} {}'.format(num_atoms + 1, v[1],
                                                   v[0] * 1.7818, '#',
                                                   gff.vdws.keys()[i],
                                                   mol_index,
                                                   mol.site_properties["mol_name"][0]))
                    num_atoms = num_atoms + 1
            mol_index += 1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_bond_coeffs(self, ffmol, mols_in_box):

        lines = []
        num_atoms = 0
        element_list = []
        mol_index = 1
        lines.append('Bond Coeffs\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.bonds:
                for i, v in enumerate(gff.bonds.values()):
                    lines.append(
                        '{} {}  {} {} {} {} {} {}'.format(num_atoms + 1, v[0],
                                                       v[1], '#',
                                                       gff.bonds.keys()[i][0],
                                                       gff.bonds.keys()[i][1],
                                                       mol_index,
                                                       mol.site_properties["mol_name"][0]))
                    num_atoms = num_atoms + 1
            mol_index += 1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_angle_coeffs(self, ffmol, mols_in_box):

        lines = []
        num_atoms = 0
        element_list = []
        mol_index = 1
        lines.append('Angle Coeffs\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.angles:
                for i, v in enumerate(gff.angles.values()):
                    lines.append(
                        '{} {} {} {} {} {} {} {} {}'.format(num_atoms + 1, v[0],
                                                         v[1], '#',
                                                         gff.angles.keys()[i][
                                                             0],
                                                         gff.angles.keys()[i][
                                                             1],
                                                         gff.angles.keys()[i][
                                                             2], mol_index,
                                                         mol.site_properties["mol_name"][0]))
                    num_atoms = num_atoms + 1
            mol_index += 1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_dihedral_coeffs(self, ffmol, mols_in_box):

        lines = []
        num_atoms = 0
        j = 0
        element_list = []
        mol_index = 1
        lines.append('Dihedral Coeffs\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.dihedrals is not None:
                for i, v in enumerate(gff.dihedrals.values()):
                    for func_form, d in v.iteritems():
                        j += 1
                        lines.append(
                            '{}  {}  {}  {} {} {} {} {} {} {} {}'.format(j, d[0],
                              func_form,
                              d[1], '#',
                              gff.dihedrals.keys()[
                                  i][0],
                              gff.dihedrals.keys()[
                                  i][1],
                              gff.dihedrals.keys()[
                                  i]
                              [2],
                              gff.dihedrals.keys()[
                                  i][3],
                              mol_index,
                              mol.site_properties["mol_name"][0]))
            mol_index += 1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_improper_coeffs(self, ffmol, mols_in_box):
        lines = []
        num_atoms = 0
        element_list = []
        mol_index = 1
        lines.append('Imp Dihedral Coeffs\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.imdihedrals:
                for i, v in enumerate(gff.imdihedrals.values()):
                    lines.append('{} {}  {}  {} {} {} {} {} {} {} {}'.format
                        (num_atoms + 1, v[0],
                         round(math.cos(math.degrees(v[1])), 0), v[2], '#',
                         gff.imdihedrals.keys()[i][0],
                         gff.imdihedrals.keys()[i][1],
                         gff.imdihedrals.keys()[i][2],
                         gff.imdihedrals.keys()[i][3], mol_index,
                         mol.site_properties["mol_name"][0]))
                    num_atoms = num_atoms + 1
            mol_index += 1
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)


    def _set_atom(self, ffmol, mols_in_box):
        """
        set the Atoms section in lammps data file
        """
        lines = []
        element_list = []
        num_atoms = 0
        masses_index = 0
        atom_type_index = {}
        lines.append('Atoms\n')
        mol_pack = mols_in_box.run()
        self.box_mol_index = []
        i = 0
        mol_index = 0
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.masses is not None:
                for m, v in enumerate(gff.masses.values()):
                    if gff.masses.keys()[m] in element_list:
                        continue
                    element_list.append(gff.masses.keys()[m])
                    masses_index = masses_index + 1
                    atom_type_index[gff.masses.keys()[m]] = masses_index
            num_atoms = len(mol)
            num_this_mol = parm['number']


            #iterate every molecule of molecule type
            for imol in range(num_this_mol):

                mol_coords = mol_pack.cart_coords[i:i + num_atoms]
                mol_index += 1

                d = {}
                for k, v in enumerate(mol_coords):
                    lines.append(
                        '{}  {}  {}  {}  {}  {} {} {}  {} {} {} {}'.format(k + i + 1,
                         mol_index,
                         atom_type_index[
                             gff.atom_index_gaff
                             [
                                 k + 1]],
                         gff.charges[mol.site_properties["mol_name"][0]][gff.atom_index[
                             k + 1]],
                         v[0], v[1],
                         v[2], '#',
                         mol_index,
                         gff.atom_index_gaff[
                             k + 1],
                         gff.atom_index[
                             k + 1],
                         mol.site_properties["mol_name"][0]))
                    d[gff.atom_index[k + 1]] = k + i + 1
                self.box_mol_index.append(d)
                i += num_atoms
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)


    def _set_bonds(self, ffmol, mols_in_box):
        """
        set the Bonds section in lammps data file
        """
        lines = []
        bond_index = 0
        bond_type_index = {}
        i = 0
        mol_index = 0
        lines.append('Bonds\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.bonds:
                for m, v in enumerate(gff.bonds.values()):
                    bond_index = bond_index + 1
                    bond_type_index[gff.bonds.keys()[m]] = bond_index
            num_this_mol = parm['number']
            for imol in range(num_this_mol):
                mol_bonds = top.bonds
                mol_index += 1
                for k, v in enumerate(mol_bonds):
                    a = gff.atom_gaff[top.bonds[k][0]]
                    b = gff.atom_gaff[top.bonds[k][1]]
                    bond_label = tuple(sorted([a, b]))
                    lines.append('{} {} {} {} {} {} {} {} {}'.format((i + k + 1),
                                  bond_type_index[
                                      bond_label],
                                  self.box_mol_index[
                                      mol_index - 1][
                                      v[0]],
                                  self.box_mol_index[
                                      mol_index - 1][
                                      v[1]],
                                  '#',
                                  mol_index,
                                  top.bonds[k][
                                      0],
                                  top.bonds[k][
                                      1],
                                  mol.site_properties["mol_name"][0]))
                i += len(top.bonds)
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)


    def _set_angles(self, ffmol, mols_in_box):
        """
        set the Angles section in lammps data file
        """
        lines = []
        angle_type_index = {}
        angle_index = 0
        i = 0
        mol_index = 0
        lines.append('Angles\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.angles:
                for m, v in enumerate(gff.angles.values()):
                    angle_index = angle_index + 1
                    angle_type_index[gff.angles.keys()[m]] = angle_index
                #num_atoms = len(mol)
            num_this_mol = parm['number']
            #iterate over first molecule
            for imol in range(num_this_mol):
                mol_angles = top.angles
                mol_index += 1
                #iterate over bonds in first molecule
                for k, v in enumerate(top.angles):
                    a = gff.atom_gaff[top.angles[k][0]]
                    b = gff.atom_gaff[top.angles[k][1]]
                    c = gff.atom_gaff[top.angles[k][2]]
                    angle_label = (a, b, c)
                    if angle_label[0] > angle_label[2]:
                        angle_label = tuple(reversed(list(angle_label)))
                    lines.append('{}  {}  {}  {}  {}  {}  {}  {}  {}  {} {}'
                    .format(i + k + 1, angle_type_index[angle_label],
                            self.box_mol_index[mol_index - 1][v[0]],
                            self.box_mol_index[mol_index - 1][v[1]],
                            self.box_mol_index[mol_index - 1][v[2]], '#',
                            mol_index,
                            top.angles[k][0], top.angles[k][1],
                            top.angles[k][2],
                            mol.site_properties["mol_name"][0]))

                i += len(top.angles)
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)

    def _set_dihedrals(self, ffmol, mols_in_box):
        """
        set the Dihedrals section in lammps data file
        """
        lines = []
        dihedral_index = 0
        dihedral_type_index = {}
        i = 0
        j = 0
        l = 0
        mol_index = 0
        lines.append('Dihedrals\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.dihedrals is not None:
                for m, v in enumerate(gff.dihedrals.values()):
                    dihedral_index = dihedral_index + 1
                    dihedral_type_index[
                        gff.dihedrals.keys()[m]] = dihedral_index

                top._get_ff_dihedrals(gff.dihedrals, top.dihedrals,
                                      gff.atom_gaff)
                num_this_mol = parm['number']
                #iterate over first molecule
                for imol in range(num_this_mol):
                    mol_dihedrals = top.topdihedralff
                    mol_index += 1
                    l += 1
                    #iterate over bonds in first molecule
                    for k, v in top.topdihedralff.iteritems():
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
                            dihedral_label = tuple(
                                reversed(list(dihedral_label)))
                        for func_form, d in v[1].iteritems():
                            j += 1
                            lines.append(
                                '{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {} {}'
                                .format(j, dihedral_type_index[dihedral_label],
                                        self.box_mol_index[mol_index - 1][A],
                                        self.box_mol_index[mol_index - 1][B],
                                        self.box_mol_index[mol_index - 1][C],
                                        self.box_mol_index[mol_index - 1][D],
                                        '#', mol_index,
                                        top.dihedrals[l][0],
                                        top.dihedrals[l][1],
                                        top.dihedrals[l][2],
                                        top.dihedrals[l][3],
                                        mol.site_properties["mol_name"][0]))
                    i += len(top.dihedrals)
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)


    def _set_imdihedrals(self, ffmol, mols_in_box):
        """
        set the Improper Dihedral section in lammps data file
        """
        lines = []
        imdihedral_type_index = {}
        imdihedrals_index = 0
        i = 0
        j = 0
        mol_index = 0
        lines.append('Impropers\n')
        for gff, top, mol, parm in zip(ffmol.gff_list, ffmol.top_list,
                                       mols_in_box.mols,
                                       mols_in_box.param_list):
            if gff.imdihedrals is not None:
                for m, v in enumerate(gff.imdihedrals.values()):
                    imdihedrals_index = imdihedrals_index + 1
                    imdihedral_type_index[
                        gff.imdihedrals.keys()[m]] = imdihedrals_index

                #iterate over types of mol
                num_this_mol = parm['number']
                #iterate over first molecule
                for imol in range(num_this_mol):
                    mol_impropers = top.imdihedrals
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
                            '{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {} {}'
                            .format(j, imdihedral_type_index[imdihedral_label],
                                    self.box_mol_index[mol_index - 1][v[0]],
                                    self.box_mol_index[mol_index - 1][v[1]],
                                    self.box_mol_index[mol_index - 1][v[2]],
                                    self.box_mol_index[mol_index - 1][v[3]],
                                    '#', mol_index,
                                    top.imdihedrals[k][0],
                                    top.imdihedrals[k][1],
                                    top.imdihedrals[k][2],
                                    top.imdihedrals[k][3],
                                    mol.site_properties["mol_name"][0]))
                    i += len(top.imdihedrals)
        lines.append('\n')
        self.lines.extend(lines)
        return '\n'.join(lines)

    def get_lammps_string(self):
        """
        returns a string of lammps data input file
        """
        my_lammps_list = []
        my_lammps_list.append(self._set_gff_types(self.ffmol, self.mols_in_box))
        my_lammps_list.append(self._set_top_types(self.ffmol, self.mols_in_box))
        my_lammps_list.append(self._set_box_dimensions(self.mols_in_box))
        my_lammps_list.append(self._set_masses(self.ffmol, self.mols_in_box))
        my_lammps_list.append(
            self._set_pair_coeffs(self.ffmol, self.mols_in_box))
        my_lammps_list.append(
            self._set_bond_coeffs(self.ffmol, self.mols_in_box))
        my_lammps_list.append(
            self._set_angle_coeffs(self.ffmol, self.mols_in_box))
        my_lammps_list.append(
            self._set_dihedral_coeffs(self.ffmol, self.mols_in_box))
        my_lammps_list.append(
            self._set_improper_coeffs(self.ffmol, self.mols_in_box))
        my_lammps_list.append(self._set_atom(self.ffmol, self.mols_in_box, ))
        my_lammps_list.append(self._set_bonds(self.ffmol, self.mols_in_box, ))
        my_lammps_list.append(self._set_angles(self.ffmol, self.mols_in_box, ))
        my_lammps_list.append(
            self._set_dihedrals(self.ffmol, self.mols_in_box, ))
        my_lammps_list.append(
            self._set_imdihedrals(self.ffmol, self.mols_in_box, ))

        return '\n'.join(my_lammps_list)

    def write_lammps_data(self, filename=None):
        """
        write lammps data input file
        """

        with open(filename, 'w') as f:
            f.write(self.get_lammps_string())

