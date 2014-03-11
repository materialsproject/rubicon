__author__ = 'navnidhirajput'




class AC():

    """
    load topology data from antechamber(.rtf) file
    """


    def __init__(self):

        self.atom_index=dict()
        self.atom_gaff=dict()

    def read_atomIndex(self,filename=None):

        with open(filename) as f:

            for line in f.readlines():
                token = line.split()
                if token[0]=='ATOM':
                    self.index=int(token[1])
                    self.atom_name=token[2]
                    self.atom_index[self.atom_name]=self.index
            self.atom_index.update(self.atom_index)


    def read_atomType(self,filename=None):

        with open(filename) as f:

            for line in f.readlines():
                token = line.split()
                if token[0]=='ATOM':
                    self.atom_name=token[2]
                    self.gaff_name=token[-1]
                    self.atom_gaff[self.atom_name]=self.gaff_name
            self.atom_gaff.update(self.atom_gaff)
        self.num_types = len(set(self.atom_gaff.values()))


class TopMol():

    def __init__(self):
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.imdihedrals = []


    def get_bonds(self,filename):
        with open(filename) as f:
            bond_section = False

            for line in f.readlines():
                if len(line.strip())==0:
                    bond_section = False
                    continue

                token = line.split()
                if token[0]=='BOND':
                    self.atom1=token[1]
                    self.atom2=token[2]
                    #self.bonds.add(tuple([self.atom1, self.atom2]))
                    self.bonds.append([self.atom1, self.atom2])


    def get_angles(self,filename):
        with open(filename) as f:
            angle_section = False

            for line in f.readlines():
                if len(line.strip())==0:
                    angle_section = False
                    continue

                token = line.split()
                if token[0]=='ANGL':
                    self.atom1=token[1]
                    self.atom2=token[2]
                    self.atom3=token[3]
                    self.angles.append([self.atom1,self.atom2,self.atom3])


    def get_dihedrals(self,filename):
        with open(filename) as f:
            dihedral_section = False

            for line in f.readlines():
                if len(line.strip())==0:
                    dihedral_section = False
                    continue

                token = line.split()
                if token[0]=='DIHE':
                    self.atom1=token[1]
                    self.atom2=token[2]
                    self.atom3=token[3]
                    self.atom4=token[4]
                    self.dihedrals.append([self.atom1,self.atom2,self.atom3,self.atom4])


    def get_imdihedrals(self,filename):
        with open(filename) as f:
            imdihedral_section = False

            for line in f.readlines():
                if len(line.strip())==0:
                    imdihedral_section = False
                    continue

                token = line.split()
                if token[0]=='IMPH':
                    self.atom1=token[1]
                    self.atom2=token[2]
                    self.atom3=token[3]
                    self.atom4=token[4]
                    self.imdihedrals.append([self.atom1,self.atom2,self.atom3,
                                             self.atom4])


















