__author__ = 'navnidhirajput'

from pymatgen.serializers.json_coders import MSONable


class TopMol(MSONable):

    """
    Reads the topology of a molecule from Antechamber output(.rtf) files.
    """

    def __init__(self):


        self.bond_type = []
        self.angle_type=[]
        self.dihedral_type=[]
        self.imdihedral_type=[]


    def read_topology(self,filename=None):


        with open(filename) as f:
            bond_section = False

            for line in f.readlines():
                if len(line.strip())==0:
                    bond_section = False
                    continue

                token = line.split()
                if token[0]=='BOND':
                    self.bond_type.append(token[1]+'-'+token[2])

                if token[0]=='ANGL':
                    self.angle_type.append(token[1]+'-'+token[2]+'-'+token[3])


                if token[0]=='DIHE':
                    self.dihedral_type.append(token[1]+'-'+token[2]+'-'+token[3]+'-'+token[4])


                if token[0]=='IMPH':
                    self.imdihedral_type.append(token[1]+'-'+token[2]+'-'+token[3]+'-'+token[4])


    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                 "bond_type": self.bond_type,
                 "angle_type": self.angle_type,
                 "dihedral_type": self.dihedral_type,
                 "imdihedral_type": self.imdihedral_type}


    @classmethod
    def from_dict(cls, d):
        return TopMol(bond_type=d["bond_type"],
                    angle_type=d["angle_type"],
                    dihedral_type=d["dihedral_type"],
                    imdihedral_type=d["imdihedral_type"])


my_top=TopMol()
my_top.read_topology("mol.rtf")
d1=my_top.to_dict
print d1

