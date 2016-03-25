# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import six
from six.moves import range

__author__ = 'navnidhirajput'


class TopMol(object):
    """
    gff_atom: GAFF atom name in .rtf file
    atom: List to map atom names with gaff atom name, [['C1', 'c'],...]
    bonds: List of bonds in .rtf files(atom name)
    angles: List of angles in .rtf files (atom name)
    dihedrals: List of dihedrals in .rtf files (atom name)
    improper dihedrals: List of improper dihedrlas in .rtf files (atom name)
    """

    def __init__(self, atoms, charges, bonds, angles, dihedrals, imdihedrals,
                 num_bonds, num_gaff_atoms, gaff_atoms):

        self.atoms = atoms
        self.bonds = bonds
        self.charges = dict() if charges is None else  charges
        self.angles = angles
        self.dihedrals = dihedrals
        self.imdihedrals = imdihedrals
        self.num_bonds = num_bonds
        self.topbondff = dict()
        self.topangleff = dict()
        self.topdihedralff = dict()
        self.topimdihedralff = dict()
        self.num_gaff_atoms = num_gaff_atoms
        self.gaff_atoms = gaff_atoms

    @classmethod
    def from_file(cls, filename, continue_on_corrupt_file=False):

        """
        read the .rtf file created by antechamber and stores
        the atoms, charges, bonds, angles, dihedrals and improper dihedrals of
        a molecule
        """
        atoms = []
        charges = []
        bonds = []
        angles = []
        dihedrals = []
        imdihedrals = []
        number_gaff_atoms = []
        gaff_atoms = []
        with open(filename) as f:
            lines = f.readlines()
            for line in lines:
                if 'need revision' in line and not continue_on_corrupt_file:
                    raise TopCorruptionException
                if len(line.strip()) == 0:
                    continue
                token = line.split()
                if token[0] == 'MASS':
                    gaff_atoms.append(token[1:3])
                    number_gaff_atoms.append(token[1])
                if token[0] == 'ATOM':
                    atoms.append(token[1:3])
                    charges.append(token[3:4])
                if token[0] == 'BOND':
                    bonds.append(token[1:3])
                elif token[0] == 'ANGL':
                    angles.append(token[1:4])
                elif token[0] == 'DIHE':
                    dihedrals.append(token[1:5])
                elif token[0] == 'IMPH':
                    imdihedrals.append(token[1:5])

            topology = TopMol(atoms, charges, bonds, angles, dihedrals,
                              imdihedrals, len(bonds), number_gaff_atoms,
                              gaff_atoms)
            return topology

    def _get_ff_dihedrals(self, gff_dihedrals, top_dihedral, atom_gaff):

        self.gaff_info = []
        for keys, values in six.iteritems(gff_dihedrals):
            self.gaff_info = [keys, values]

        for item in top_dihedral:
            d1 = item[0] + ' ' + item[1] + ' ' + item[2] + ' ' + item[3]
            a1, a2, a3, a4 = atom_gaff[item[0]], atom_gaff[item[1]], atom_gaff[
                item[2]], atom_gaff[item[3]]
            dihedral_label = (a1, a2, a3, a4)
            if dihedral_label[0] > dihedral_label[3]:
                dihedral_label = tuple(reversed(list(dihedral_label)))
            if dihedral_label in gff_dihedrals:
                self.topdihedralff[d1] = (
                    dihedral_label, gff_dihedrals[dihedral_label])
            self.num_dih_types = len(set(self.topdihedralff.keys()))

    def _get_ff_bonds(self, gff_bonds, top_bond, atom_gaff):

        self.gaff_info = []
        for keys, values in six.iteritems(gff_bonds):
            self.gaff_info = [keys, values]
        for item in top_bond:
            d1 = item[0] + ' ' + item[1]
            a1, a2 = atom_gaff[item[0]], atom_gaff[item[1]]
            if (str(a1), str(a2)) in gff_bonds:
                self.topbondff[d1] = (
                    (str(a1), str(a2)), gff_bonds[(str(a1), str(a2))])
            elif (str(a2), str(a1)) in gff_bonds:
                bond_type = tuple(sorted((str(a1), str(a2))))
                self.topbondff[d1] = ((str(a1), str(a2)), gff_bonds[bond_type])
                self.num_bond_types = len(set(self.topbondff.keys()))

    def _get_ff_angles(self, gff_angles, top_angle, atom_gaff):

        self.gaff_info = []
        for keys, values in six.iteritems(gff_angles):
            self.gaff_info = [keys, values]
        for item in top_angle:
            d1 = item[0] + ' ' + item[1] + ' ' + item[2]
            a1, a2, a3 = atom_gaff[item[0]], atom_gaff[item[1]], atom_gaff[
                item[2]]
            if (str(a1), str(a2), str(a3)) in gff_angles:
                self.topangleff[d1] = ((str(a1), str(a2), str(a3)),
                                       gff_angles[(str(a1), str(a2), str(a3))])
            elif (str(a3), str(a2), str(a1)) in gff_angles:
                angle_type = tuple(sorted((str(a1), str(a2), str(a3))))
                self.topangleff[d1] = (
                    (str(a1), str(a2), str(a3)), gff_angles[angle_type])
            self.num_ang_types = len(set(self.topangleff.keys()))

    def _get_ff_imdihedrals(self, gff_imdihedrals, top_imdihedral, atom_gaff):

        self.gaff_info = []
        for keys, values in six.iteritems(gff_imdihedrals):
            self.gaff_info = [keys, values]

        for item in top_imdihedral:
            d1 = item[0] + ' ' + item[1] + ' ' + item[2] + ' ' + item[3]
            a1, a2, a3, a4 = atom_gaff[item[0]], atom_gaff[item[1]], atom_gaff[
                item[2]], atom_gaff[item[3]]
            imdihedral_label = (a1, a2, a3, a4)
            if imdihedral_label[0] > imdihedral_label[3]:
                imdihedral_label = tuple(reversed(list(imdihedral_label)))
            if imdihedral_label in gff_imdihedrals:
                self.topimdihedralff[d1] = (
                    imdihedral_label, gff_imdihedrals[imdihedral_label])
            self.num_imdih_types = len(set(self.topimdihedralff.keys()))

    @classmethod
    def gaff_fromfile(self, filename=None):
        """
        read the gaff.txt file to store
        atoms, bonds, angles, dihedrals and improper dihedrals of
        a molecule
        """
        atom_section = False
        atoms = {}
        bonds = []
        angles = []
        dihedrals = []
        imdihedrals = []
        num_gaff_atoms = []
        gaff_atoms = []

        with open(filename) as f:
            for line in f.readlines():
                if line.startswith('atom_name'):
                    atom_section = True
                    continue
                if atom_section:
                    if len(line.strip()) == 0:
                        atom_section = False
                        continue
                    atom_name = line[0:2]
                    atom_mass = line[3:9]
                    atoms[atom_name] = atom_mass
                if len(line.strip()) == 0:
                    continue
                token = line.split()
                if token[0] == 'ATOM':
                    atoms.append(token[1:3])
                if token[0] == 'BOND':
                    bonds.append(token[1:3])
                elif token[0] == 'ANGL':
                    angles.append(token[1:4])
                elif token[0] == 'DIHE':
                    dihedrals.append(token[1:5])
                elif token[0] == 'IMPH':
                    imdihedrals.append(token[1:5])
            topology = TopMol(atoms, bonds, angles, dihedrals, imdihedrals,
                              len(bonds), num_gaff_atoms, gaff_atoms)
            return topology


class TopCorruptionException(Exception):
    pass


def correct_corrupted_top_files(corrupted_file=None, gaff_file=None):
    rtf_lines = []
    rtf_lines.append('{}{}'.format('* Topology File.\n', '*\n'))

    top = TopMol.from_file(corrupted_file, True)
    top_gaff = TopMol.gaff_fromfile(gaff_file)
    atom_index = 1

    for x in range(0, len(top.atoms)):
        if top.atoms[x][1].lower() in list(top_gaff.atoms.keys()):
            rtf_lines.append('{}  {}  {}   {}'.format('MASS', atom_index,
                                                      top.atoms[x][1].lower(),
                                                      top_gaff.atoms[
                                                          top.atoms[x][
                                                              1].lower()]))
            rtf_lines.append('\n')
        atom_index += 1
    rtf_lines.append('{}{}'.format('\nRESI MOL  0.000\n', 'GROUP\n'))

    for atoms in top.atoms:
        rtf_lines.append(
            '{}  {}  {}'.format('ATOM', atoms[0], atoms[1].lower()))
        rtf_lines.append('\n')
    rtf_lines.append('\n')

    for bonds in top.bonds:
        rtf_lines.append('{}  {}  {}'.format('BOND', bonds[0], bonds[1]))
        rtf_lines.append('\n')

    rtf_lines.append('\n')

    for angles in top.angles:
        rtf_lines.append(
            '{}  {}  {}  {}'.format('ANGL', angles[0], angles[1], angles[2]))
        rtf_lines.append('\n')

    rtf_lines.append('\n')

    for dihedrals in top.dihedrals:
        rtf_lines.append(
            '{}  {}  {}  {}  {}'.format('DIHE', dihedrals[0], dihedrals[1],
                                        dihedrals[2], dihedrals[3]))
        rtf_lines.append('\n')

    with open('mol.rtf', 'w') as f:
        f.writelines(rtf_lines)
