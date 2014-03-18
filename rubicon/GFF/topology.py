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

    def __init__(self,bonds,angles,dihedrals,imdihedrals,num_bonds):
        self.bonds=bonds
        self.angles=angles
        self.dihedrals=dihedrals
        self.imdihedrals= imdihedrals
        self.num_bonds=num_bonds
        


    @classmethod
    def from_file(cls,filename):
        bonds=[]
        angles=[]
        dihedrals=[]
        imdihedrals=[]
        with open(filename) as f :
            lines=f.readlines()
            for line in lines:
                if len(line.strip())==0:
                    continue
                token = line.split()
                if token[0]=='BOND':
                    bonds.append(token[1:3])


                elif token[0]=='ANGL':
                    angles.append(token[1:4])
                elif token[0]=='DIHE':
                    dihedrals.append(token[1:5])
                elif token[0]=='IMPH':
                    imdihedrals.append(token[1:5])


            topology=TopMol(bonds,angles,dihedrals,imdihedrals,len(bonds))
            return topology


















