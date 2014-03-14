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
                    index=int(token[1])
                    atom_name=token[2]
                    self.atom_index[atom_name]=index
            self.atom_index.update(self.atom_index)


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



class TopMol():

    def __init__(self,bonds,angles,dihedrals,imdihedrals):
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
                bonds.append(token[1:3])
        print bonds
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
                angles.append(token[1:4])
        print angles
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
                dihedrals.append(token[1:5])
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
                imdihedrals.append(token[1:5])
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


















