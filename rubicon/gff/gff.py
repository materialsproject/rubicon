from collections import defaultdict
import json
import os

__author__ = 'navnidhirajput'

from pymatgen.serializers.json_coders import MSONable


class Gff(MSONable):
    """
    A force field library. Reads the output file from
    AntechamberRunner and populate the FF library

    Args:
        bonds: store the bond distance (A) and spring constant (Kcal/molA2)
         in a dict.
        angles: store the bond distance (A) and spring constant
        (Kcal/mol*radian2) in a dict.
        dihedral: store the magnitude of torsion (Kcal/mol), phase offset in
        degree and the periodicity of torsion in a dict.
        imdihedrals: store improper dihedral information which include
        the magnitude of torsion (Kcal/mol), phase
        offset in degree and the periodicity of torsion in a dict
        vdws: store the van der waal radius (A) and van der wall depth for a
        given atom (Kcal/mol) in a dict.

        All the dictionaries are then converted into json serializable object
            Example:
            {'bonds': {'c-o': (648.0, 1.214)
            'angles': {'os-c-os': (76.45, 111.38)
            'dihedrals': {'X-c-os-X': (2.7, 180.0)
            'imdihedrals': {'o-os-c-os': (1.1, 180.0)}
            'vdws': {'c': (1.908, 0.086)}
    """


    def __init__(self, atoms, bonds, angles, dihedrals, imdihedrals, vdws, masses, charges):

        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.imdihedrals = imdihedrals
        self.vdws = vdws
        self.masses = masses
        self.charges = dict() if charges is None else charges
        self.atom_index=dict()
        self.atom_index_gaff=dict()
        self.atom_gaff=dict()



    @classmethod
    def from_forcefield_para(cls, filename=None, continue_on_corrupt_file = False):
        """
        reads ANTECHAMBER.FRCMOD and stores the force field parameters for
        bonds, angles, dihedrals, improper dihedrals, van der waals and Masses
        """
        atoms = {}
        bonds = {}
        angles = {}
        dihedrals = defaultdict(dict)
        imdihedrals = {}
        vdws = {}



        with open(filename) as f:
            bond_section = False
            angle_section = False
            dihedral_section = False
            imdihedral_section = False
            vdw_section = False
            mass_section = False

            for line in f.readlines():
                if 'need revision' in line and not continue_on_corrupt_file:
                    raise FFCorruptionException

                if line.startswith('MASS'):
                    mass_section = True
                    continue
                if mass_section:
                    if len(line.strip()) == 0:
                        mass_section = False
                        continue

                    atom_type = line[0:2].strip()
                    mass = float(line[3:9])
                    atoms[atom_type.lower()] = (mass)

                if line.startswith('BOND'):
                    bond_section = True
                    continue
                if bond_section:
                    if len(line.strip()) == 0:
                        bond_section = False
                        continue
                    bond_type = (line[0:2].strip(), line[3:5].strip())
                    bond_type = tuple(sorted(bond_type))
                    bond_constant = float(line[7:13])
                    bond_distance = float(line[16:21])
                    bonds[bond_type] = (bond_constant, bond_distance)

                if line.startswith('ANGLE'):
                    angle_section = True
                    continue
                if angle_section:
                    if len(line.strip()) == 0:
                        angle_section = False
                        continue
                    angle_type = line[0:2].strip(), line[3:5].strip(), \
                                 line[6:8].strip()

                    if line[0:2].strip()> line[6:8].strip():
                        angle_type=tuple(reversed(angle_type))

                    angle_constant = float(line[11:17])
                    angle_distance = float(line[22:29])
                    angles[angle_type] = (angle_constant, angle_distance)


                if line.startswith('DIHE'):
                    dihedral_section = True
                    continue
                if dihedral_section:
                    if len(line.strip()) == 0:
                        dihedral_section = False
                        continue
                    dihedral_type = line[0:2].strip(), line[3:5].strip(), \
                                    line[6:8].strip(), line[9:11].strip()
                    if dihedral_type[0] > dihedral_type[3]:
                        dihedral_type = tuple(reversed(list(dihedral_type)))
                    dihedral_func_type = (line[49:50])
                    dihedral_constant = float(line[19:24])
                    dihedral_angle = float(line[31:38])
                    dihedrals[(dihedral_type)][dihedral_func_type] = (
                    dihedral_constant, dihedral_angle)


                if line.startswith('IMPROPER'):
                    imdihedral_section = True
                    continue
                if imdihedral_section:
                    if len(line.strip()) == 0:
                        imdihedral_section = False
                        continue
                    imdihedral_type = line[0:2].strip(), line[3:5].strip(), \
                                      line[6:8].strip(), line[9:11].strip()
                    imdihedral_distance = float(line[19:24])
                    imdihedral_angle = float(line[31:38])
                    imdihedral_function= float(line[47:50])
                    imdihedrals[imdihedral_type] = (
                    imdihedral_distance, imdihedral_angle, imdihedral_function)


                if line.startswith('NONBON'):
                    vdw_section = True
                    continue
                if vdw_section:
                    if len(line.strip()) == 0:
                        vdw_section = False
                        continue
                    token = line.split()
                    vdw_type = line[2:4].strip()
                    vdw_type = line[2:4].strip()
                    sigma = float(line[14:20])
                    epsilon = abs(float(line[22:28]))
                    vdws[vdw_type] = (sigma, epsilon)


            return atoms, bonds,angles,dihedrals,imdihedrals,vdws




    @classmethod
    def from_gaff_para(cls, filename=None):
        """
        reads gaff_nidhi.dat and stores the force field parameters for
        bonds, angles, dihedrals, improper dihedrals, van der waals and Masses
        """
        atom_section = False
        atoms={}
        bonds = {}
        angles = {}
        general_dihedrals = {}
        specific_dihedrals = {}
        dihedrals = {}
        vdws = {}


        with open('gaff_nidhi.txt') as f:
            atom_section = False
            bond_section = False
            angle_section = False
            dihedral_section = False
            vdw_section = False

            for line in f.readlines():

                if line.startswith('atom_name'):
                    atom_section = True
                    continue
                if atom_section:
                    if len(line.strip())==0:
                        atom_section = False
                        continue
                    atom_name = line[0:2]
                    atom_mass = line[3:9]
                    atoms[atom_name] = atom_mass

                if line.startswith('bond_type'):
                    bond_section = True
                    continue
                if bond_section:
                    if len(line.strip())==0:
                        bond_section = False
                        continue
                    bond_type = (line[0:2].strip(), line[3:5].strip())
                    bond_type = tuple(sorted(bond_type))
                    bond_constant = float(line[7:13])
                    bond_distance = float(line[16:21])
                    bonds[bond_type] = (bond_constant, bond_distance)


                if line.startswith('angle_type'):
                        angle_section = True
                        continue
                if angle_section:
                    if len(line.strip()) == 0:
                        angle_section = False
                        continue
                    angle_type = line[0:2].strip(), line[3:5].strip(),line[6:8].strip()

                    if line[0:2].strip()> line[6:8].strip():
                        angle_type=tuple(reversed(angle_type))

                    angle_constant = float(line[11:17])
                    angle_distance = float(line[22:29])
                    angles[angle_type] = (angle_constant, angle_distance)


                if line.startswith('dihedral_type'):
                    dihedral_section = True
                    continue
                if dihedral_section:
                    if len(line.strip()) == 0:
                        dihedral_section = False
                        continue
                    dihedral_type = line[0:2].strip(), line[3:5].strip(), \
                                    line[6:8].strip(), line[9:11].strip()
                    if 'X' in dihedral_type:

                                    general_dihedral_type = line[3:5].strip(),line[6:8].strip()
                                    general_dihedral_func_type = (line[14:15])
                                    general_dihedral_constant = float(line[18:24])
                                    general_dihedral_angle = float(line[31:38])
                                    general_dihedrals[(general_dihedral_type)] = (general_dihedral_func_type,
                                    general_dihedral_constant, general_dihedral_angle)


                    else:
                        specific_dihedral_type = line[0:2].strip(), line[3:5].strip(), \
                                line[6:8].strip(), line[9:11].strip()
                        specific_dihedral_func_type = (line[14:15])
                        specific_dihedral_constant = float(line[18:24])
                        specific_dihedral_angle = float(line[31:38])
                        specific_dihedrals[(specific_dihedral_type)] = (specific_dihedral_func_type,
                        specific_dihedral_constant, specific_dihedral_angle)



                if line.startswith('non_bonded'):
                    vdw_section = True
                    continue
                if vdw_section:
                    if len(line.strip()) == 0:
                        vdw_section = False
                        continue
                    token = line.split()
                    vdw_type = token[0]
                    sigma = token[1]
                    epsilon = token[2]
                    vdws[vdw_type] = (sigma, epsilon)

            return atoms,bonds,angles,specific_dihedrals,general_dihedrals,vdws



    def read_atom_index(self,mol,filename=None):

        """
        read ANTECHAMBER_AC.AC to store the antechamber atom name
        and GAFF atom name and index of atoms in a dict
        """

        with open(filename) as f:
            for line in f.readlines():
                token = line.split()
                if token[0]=='ATOM':
                    index=int(token[1])
                    atom_name=token[2]
                    gaff_name=token[-1]
                    self.atom_gaff[atom_name]=gaff_name
                    self.atom_index[index]=atom_name
                    self.atom_index_gaff[index]=gaff_name
            self.atom_gaff.update(self.atom_gaff)
        self.num_types = len(set(self.atom_gaff.values()))




    def read_charges(self):
        """
        reads charges.json to read charges from GAFF atom
        """
        filename = os.path.join(os.path.dirname(__file__),'charges.json')
        jsonfile = open(filename)
        self.charges = json.load(jsonfile, encoding="utf-8")
        return self.charges



    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "bonds": self.bonds,
                "angles": self.angles,
                "dihedrals": self.dihedrals,
                "imdihedrals": self.imdihedrals,
                "vdws": self.vdws}


    @classmethod
    def from_dict(cls, d):
        return Gff(bonds=d["bonds"],
                   angles=d["angles"],
                   dihedrals=d["dihedrals"],
                   imdihedrals=d["imdihedrals"],
                   vdws=d["vdws"]
        )

class FFCorruptionException(Exception):
    pass


def correct_corrupted_frcmod_files(corrupted_file = None, gaff_file = None):

    frc_lines = []
    frc_lines.append('{}{}'.format('remark goes here\n', 'MASS\n'))


    gff = Gff.from_forcefield_para(corrupted_file, True)
    gff_para = Gff.from_gaff_para(gaff_file)

    gaff_atom = gff[0].keys()
    for ant_atom_name in gaff_atom:
        print gaff_atom[0].lower()
        if ant_atom_name in gff_para[0].keys():
             frc_lines.append('{}  {}'.format(ant_atom_name,gff_para[0][ant_atom_name]))
             frc_lines.append('\n')

    frc_lines.append('\nBOND\n')
    for bond in gff[1]:
        if bond in gff_para[1].keys():

            frc_lines.append('{}{}{}    {}     {}'.format(bond[0], '-', bond[1],gff_para[1][bond][0],gff_para[1][bond][1]))
            frc_lines.append('\n')

    frc_lines.append('\nANGLE\n')
    for angle in gff[2]:
        if angle in gff_para[2].keys():

            frc_lines.append('{}{}{}{}{}   {}    {}'.format(angle[0], '-', angle[1], '-', angle[2],gff_para[2][angle][0],gff_para[2][angle][1]))
            frc_lines.append('\n')

    frc_lines.append('\nDIHE\n')
    for dihedral in gff[3]:
        if dihedral in gff_para[3].keys():
            frc_lines.append('{}{}{}{}{}{}{}'.format(dihedral[0], '-', dihedral[1], '-', dihedral[2], '-', dihedral[3]))
            frc_lines.append('\n')
        elif dihedral[1:3]  in gff_para[4].keys():

            frc_lines.append('{}{}{}{}{}{}{}'.format(dihedral[0], '-', dihedral[1], '-', dihedral[2], '-', dihedral[3]))
            frc_lines.append('\n')

        elif dihedral[::-1][1:3] in gff_para[4].keys():

            frc_lines.append('{}{}{}{}{}{}{}   {}       {}      {}      {}'.format(dihedral[0], '-', dihedral[1], '-', dihedral[2], '-', dihedral[3],
                                                                   '1',gff_para[4][dihedral[::-1][1:3]][1],gff_para[4][dihedral[::-1][1:3]][2],gff_para[4][dihedral[::-1][1:3]][0]))
            frc_lines.append('\n')


    frc_lines.append('\nIMPROPER\n')
    for improper in gff[4]:
        if improper in gff_para[3].keys():
            frc_lines.append('{}{}{}{}{}{}{}'.format(improper[0], '-', improper[1], '-', improper[2], '-', improper[3]))
            frc_lines.append('\n')
        elif improper[1:3]  in gff_para[4].keys():

            frc_lines.append('{}{}{}{}{}{}{}'.format(improper[0], '-', improper[1], '-', improper[2], '-', improper[3]))
            frc_lines.append('\n')

        elif improper[::-1][1:3] in gff_para[4].keys():

            frc_lines.append('{}{}{}{}{}{}{}   {}       {}      {}      {}'.format(improper[0], '-', improper[1], '-', improper[2], '-', improper[3],
                                                                   '1',gff_para[4][improper[::-1][1:3]][1],gff_para[4][improper[::-1][1:3]][2],gff_para[4][improper[::-1][1:3]][0]))
            frc_lines.append('\n')


    frc_lines.append('\nNONBON\n')
    for vdw in gff[5]:
        if vdw in gff_para[5].keys():
            frc_lines.append('{}  {}    {}'.format(vdw,gff_para[5][vdw][0],gff_para[5][vdw][1]))


        frc_lines.append('\n')

    with open('mol.frcmod', 'w') as f:
        f.writelines(frc_lines)








