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

    def __init__(self,bonds,angles,dihedrals,imdihedrals):
        #self.bonds = []
        #self.angles = []
        #self.dihedrals = []
        #self.imdihedrals = []
        self.bonds=bonds
        self.angles=angles
        self.dihedrals=dihedrals
        self.imdihedrals= imdihedrals

    @classmethod
    def get_bonds(cls,lines):
        bond_section = False
        bonds=[]

        for line in lines:
            if len(line.strip())==0:
                bond_section = False
                continue

            token = line.split()
            if token[0]=='BOND':
                atom1=token[1]
                atom2=token[2]
                bonds.append([atom1, atom2])
        return bonds

    @classmethod
    def get_angles(self,lines):
        angle_section = False
        angles=[]

        for line in lines:
            if len(line.strip())==0:
                angle_section = False
                continue

            token = line.split()
            if token[0]=='ANGL':
                atom1=token[1]
                atom2=token[2]
                atom3=token[3]
                angles.append([atom1,atom2,atom3])
        return angles

    @classmethod
    def get_dihedrals(self,lines):

        dihedral_section = False
        dihedrals=[]

        for line in lines:
            if len(line.strip())==0:
                dihedral_section = False
                continue

            token = line.split()
            if token[0]=='DIHE':
                atom1=token[1]
                atom2=token[2]
                atom3=token[3]
                atom4=token[4]
                dihedrals.append([atom1,atom2,atom3,atom4])
        return dihedrals

    @classmethod
    def get_imdihedrals(self,lines):
        imdihedral_section = False
        imdihedrals=[]
        for line in lines:
            if len(line.strip())==0:
                imdihedral_section = False
                continue

            token = line.split()
            if token[0]=='IMPH':
                atom1=token[1]
                atom2=token[2]
                atom3=token[3]
                atom4=token[4]
                imdihedrals.append([atom1,atom2,atom3,atom4])
        return imdihedrals

    @classmethod
    def from_file(cls,filename):
        with open(filename) as f :
            lines=f.readlines()
            #print lines
            bonds=cls.get_bonds(lines)
            angles=cls.get_angles(lines)
            dihedrals=cls.get_dihedrals(lines)
            imdihedrals=cls.get_imdihedrals(lines)
            topology=TopMol(bonds,angles,dihedrals,imdihedrals)
        return topology


















