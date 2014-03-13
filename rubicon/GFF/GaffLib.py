__author__ = 'navnidhirajput'

from pymatgen.serializers.json_coders import MSONable

class GaffLib(MSONable):

    def __init__(self,bonds={}, angles={}, dihedrals={},
                 imdihedrals={}, vdws={},masses={}):

        self.bonds = bonds
        self.angles=angles
        self.dihedrals=dihedrals
        self.imdihedrals=imdihedrals
        self.vdws=vdws
        self.masses=masses




    def read_forcefield_para(self,filename=None):

        bonds = dict()
        angles = dict()
        dihedrals = dict()
        imdihedrals = dict()
        vdws = dict()
        masses=dict()


        with open(filename) as f:
            bond_section = False
            angle_section = False
            dihedral_section = False
            imdihedral_section = False
            vdw_section = False
            mass_section=False

            for line in f.readlines():
                if line.startswith('BOND'):
                    bond_section = True
                    continue
                if bond_section:
                    if len(line.strip())==0:
                        bond_section = False
                        continue
                    token = line.split()

                    bond_type=token[0]
                    bond_k_distance=float(token[1])
                    bond_distance=float(token[2])
                    bonds[bond_type]=(bond_k_distance,bond_distance)


                if line.startswith('ANGLE'):
                    angle_section = True
                    continue
                if angle_section:
                    if len(line.strip())==0:
                        angle_section = False
                        continue
                    token = line.split()
                    angle_type=token[0]
                    angle_k_distance=float(token[1])
                    angle_distance=float(token[2])
                    angles[angle_type]=(angle_k_distance,angle_distance)


                if line.startswith('DIHEDRAL'):
                    dihedral_section = True
                    continue
                if dihedral_section:
                    if len(line.strip())==0:
                        dihedral_section = False
                        continue
                    token = line.split()
                    dihedral_type=token[0]
                    dihedral_k_distance=float(token[2])
                    dihedral_func_type=float(token[1])
                    dihedral_angle=float(token[3])
                    dihedrals[dihedral_type]=(dihedral_k_distance,dihedral_angle)


                if line.startswith('IMPHI'):
                    imdihedral_section = True
                    continue
                if imdihedral_section:
                    if len(line.strip())==0:
                        imdihedral_section = False
                        continue
                    token = line.split()
                    imdihedral_type=token[0]
                    imdihedral_distance=float(token[1])
                    imdihedral_func_type=float(token[3])
                    imdihedral_angle=float(token[2])
                    imdihedrals[imdihedral_type]=(imdihedral_distance,imdihedral_angle)

                '''
                if line.startswith('NONBONDED'):
                    vdw_section = True
                    continue
                if vdw_section:
                    if len(line.strip())==0:
                        vdw_section = False
                        continue
                    if line.startswith("CUTNB"):
                        continue
                    if line.startswith("!"):
                        continue
                    token = line.split()
                    vdw_type=token[0]
                    epsilon=abs(float(token[2]))
                    sigma=float(token[3])
                    vdws[vdw_type]=(sigma,epsilon)
            self.bonds.update(bonds)
            self.angles.update(angles)
            self.dihedrals.update(dihedrals)
            self.imdihedrals.update(imdihedrals)
            '''

            self.bonds.update(bonds)
            self.angles.update(angles)
            self.dihedrals.update(dihedrals)
            self.imdihedrals.update(imdihedrals)
            print bonds
            print angles
            print dihedrals
            print imdihedrals




    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                 "bonds": self.bonds,
                 "angles": self.angles,
                 "dihedrals": self.dihedrals,
                 "imdihedrals": self.imdihedrals,
                 "vdws":self.vdws}


    @classmethod
    def from_dict(cls, d):
        return GaffLib(bonds=d["bonds"],
                    angles=d["angles"],
                    dihedrals=d["dihedrals"],
                    imdihedrals=d["imdihedrals"],
                    vdws=d["vdws"]
                    )


my_gaff_lib=GaffLib()
my_gaff_lib.read_forcefield_para("gaff.dat")
my_gaff_lib.to_dict