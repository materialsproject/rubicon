from pymatgen import CovalentBond

__author__ = 'navnidhirajput'




class AC():

    """
    load topology data from antechamber(.rtf) file
    """


    def __init__(self):

        self.atom_index=dict()
        self.atom_name=0

    def read_atomType(self,filename=None):

        with open(filename) as f:

            for line in f.readlines():
                token = line.split()
                if token[0]=='ATOM':
                    self.index=int(token[1])
                    self.atom_name=token[2]
                    self.atom_index[self.atom_name]=self.index
            self.atom_index.update(self.atom_index)



class TopBond():

    def __init__(self):
        self.bonds = []

    def get_bonds(self,filename, atom_index, mol):
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
                    #print self.atom1
                    id1 = atom_index[self.atom1] - 1
                    id2 = atom_index[self.atom2] - 1
                    self.bonds.append([mol.sites[id1], mol.sites[id2]])



class TopAngle():

    def __init__(self):
        self.angles = []

    def get_angles(self,filename, atom_index, mol):
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
                    #id1 = atom_index[self.atom1] - 1
                    #id2 = atom_index[self.atom2] - 1
                    #id3 = atom_index[self.atom3] - 1
                    #self.angles.append([mol.sites[id1], mol.sites[id2],
                    #                   mol.sites[id3]])


class TopDihedral():

    def __init__(self):
        self.dihedrals = []

    def get_dihedrals(self,filename, atom_index, mol):
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
                    self.angles.append([self.atom1,self.atom2,self.atom3,self.atom4])
                    #id1 = atom_index[self.atom1] - 1
                    #id2 = atom_index[self.atom2] - 1
                    #id3 = atom_index[self.atom3] - 1
                    #id4 = atom_index[self.atom4] - 1
                    #self.dihedrals.append([mol.sites[id1], mol.sites[id2],
                    #                    mol.sites[id3], mol.sites[id4]])



class TopImDihedral():

    def __init__(self):
        self.imdihedrals = []

    def get_imdihedrals(self,filename, atom_index, mol):
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
                    self.angles.append([self.atom1,self.atom2,self.atom3,self.atom4])
                    #id1 = atom_index[self.atom1] - 1
                    #id2 = atom_index[self.atom2] - 1
                    #id3 = atom_index[self.atom3] - 1
                    #id4 = atom_index[self.atom4] - 1
                    #self.imdihedrals.append([mol.sites[id1], mol.sites[id2],
                    #                    mol.sites[id3], mol.sites[id4]])




