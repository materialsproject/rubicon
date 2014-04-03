from pymatgen import Molecule

__author__ = 'navnidhirajput'


class TopMol(Molecule):

    def __init__(self,atoms,bonds,angles,dihedrals,imdihedrals,num_bonds):
        self.atoms=atoms
        self.bonds=bonds
        self.angles=angles
        self.dihedrals=dihedrals
        self.imdihedrals= imdihedrals
        self.num_bonds=num_bonds
        self.atom_index=dict()
        self.atom_index_gaff=dict()
        self.atom_gaff=dict()


    @classmethod
    def from_file(cls,filename):
        atoms=[]
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
                if token[0]=='ATOM':
                    atoms.append(token[1:3])
                if token[0]=='BOND':
                    bonds.append(token[1:3])
                elif token[0]=='ANGL':
                    angles.append(token[1:4])
                elif token[0]=='DIHE':
                    dihedrals.append(token[1:5])
                elif token[0]=='IMPH':
                    imdihedrals.append(token[1:5])
            topology=TopMol(atoms,bonds,angles,dihedrals,imdihedrals,len(bonds))
            return topology




















