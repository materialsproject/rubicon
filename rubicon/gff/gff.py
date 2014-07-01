from collections import defaultdict
import json
import os

__author__ = 'navnidhirajput'

from pymatgen.serializers.json_coders import MSONable


class Gff(MSONable):
    """
    A force field library. Right now reads the output file from
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


    def __init__(self, bonds, angles, dihedrals, imdihedrals, vdws, masses,charges):

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
    def from_forcefield_para(cls, filename=None):

        bonds = dict()
        angles = dict()
        dihedrals = defaultdict(dict)
        imdihedrals = dict()
        vdws = dict()
        masses = dict()


        with open(filename) as f:
            bond_section = False
            angle_section = False
            dihedral_section = False
            imdihedral_section = False
            vdw_section = False
            mass_section = False

            for line in f.readlines():

                if line.startswith('MASS'):
                    mass_section = True
                    continue
                if mass_section:
                    if len(line.strip()) == 0:
                        mass_section = False
                        continue

                    atom_type = line[0:2].strip()
                    mass = float(line[3:9])
                    masses[atom_type] = (mass)

                if line.startswith('BOND'):
                    bond_section = True
                    continue
                if bond_section:
                    if len(line.strip()) == 0:
                        bond_section = False
                        continue
                    bond_type = (line[0:2].strip(), line[3:5].strip())
                    bond_type = tuple(sorted(bond_type))
                    bond_k_distance = float(line[7:13])
                    bond_distance = float(line[16:21])
                    bonds[bond_type] = (bond_k_distance, bond_distance)

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

                    angle_k_distance = float(line[11:17])
                    angle_distance = float(line[22:29])
                    angles[angle_type] = (angle_k_distance, angle_distance)


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
                    dihedral_k_distance = float(line[19:24])
                    dihedral_angle = float(line[31:38])
                    dihedrals[(dihedral_type)][dihedral_func_type] = (
                    dihedral_k_distance, dihedral_angle)


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
                    vdw_type = line[2:4].strip()
                    sigma = float(line[14:20])
                    epsilon = abs(float(line[22:28]))
                    vdws[vdw_type] = (sigma, epsilon)

            return Gff(bonds,angles,dihedrals,imdihedrals,vdws,masses,None)


    def read_atom_index(self,mol,filename=None):

        with open(filename) as f:

            for line in f.readlines():
                token = line.split()
                if token[0]=='ATOM':
                    index=int(token[1])
                    atom_name=token[2]
                    atom_gaff=token[9]
                    atom_name=token[2]
                    gaff_name=token[-1]
                    self.atom_gaff[atom_name]=gaff_name
                    self.atom_index[index]=atom_name
                    self.atom_index_gaff[index]=atom_gaff
                    # print (self.atom_index.values())
                    # mol.add_site_property("atomname",self.atom_index.values())
            self.atom_gaff.update(self.atom_gaff)
            # for atom_name in (self.atom_index.values()):
            #     atom_name.append(atom_name)
            # mol.add_site_property("atomname",atom_name)
            # print  mol.site_properties["atomname"]
        self.num_types = len(set(self.atom_gaff.values()))




    def read_atomType(self,filename=None):

        with open(filename) as f:

            for line in f.readlines():
                token = line.split()
                if token[0]=='ATOM':
                    atom_name=token[2]
                    gaff_name=token[-1]
                    self.atom_gaff[atom_name]=gaff_name
            self.atom_gaff.update(self.atom_gaff)
        self.num_types = len(set(self.atom_gaff.values()))

    def read_charges(self):
        filename = os.path.join(os.path.dirname(__file__),'charges.json')
        jsonfile = open(filename)
        self.charges = json.load(jsonfile, encoding="utf-8")
        #print self.charges["TFN"]['O']






    @property
    def to_dict(self):
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









