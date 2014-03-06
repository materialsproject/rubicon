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




class TopBond():

    def __init__(self):
        #self.bonds = set()
        self.bonds = []


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




class TopAngle():

    def __init__(self):
        self.angles = []

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



class TopDihedral():

    def __init__(self):
        self.dihedrals = []

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




class TopImDihedral():

    def __init__(self):
        self.imdihedrals = []

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

class TopBondFF():

    """
    takes topology and FF parameters of bonds

    """
    def __init__(self):
        self.topbondFF=dict()


    def get_FF_bonds(self,bonds,top_bond,atom_gaff):

        self.gaff_info=[]
        #print bonds['ca-ca'][0]

        for keys, values in bonds.iteritems():
                self.gaff_info=[keys,values]

        for item in top_bond:
            d1=item[0]+'-'+ item[1]
            a1=atom_gaff[item[0]]
            a2=atom_gaff[item[1]]
            if a1+'-'+a2 in bonds:
                self.topbondFF[d1]=(a1+'-'+a2,bonds[a1+'-'+a2])
            else:
                self.topbondFF[d1]=(a2+'-'+a1,bonds[a2+'-'+a1])
       



class TopAngleFF():

    """
    takes topology and FF parameters of Angles

    """
    def __init__(self):
        self.topangleFF=dict()


    def get_FF_angles(self,angles,top_angle,atom_gaff):

        self.gaff_info=[]
        #print bonds['ca-ca'][0]

        for keys, values in angles.iteritems():
                self.gaff_info=[keys,values]
        #print self.gaff_info

        for item in top_angle:
            d1=item[0]+'-'+ item[1]+'-'+item[2]
            a1=atom_gaff[item[0]]
            a2=atom_gaff[item[1]]
            a3=atom_gaff[item[2]]
            if a1+'-'+a2+'-'+a3 in angles:
                self.topangleFF[d1]=(a1+'-'+a2+'-'+a3,angles[a1+'-'+a2+'-'+a3])

            elif a1+'-'+a2+'-'+a3:
                self.topangleFF[d1]=(a1+'-'+a2+'-'+a3,angles[a1+'-'+a2+'-'+a3])
















