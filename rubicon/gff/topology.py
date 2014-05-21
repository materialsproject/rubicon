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
        self.topbondff = dict()
        self.topangleff = dict()
        self.topdihedralff = dict()
        self.topimdihedralff = dict()



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



    def _get_ff_dihedrals(self, gff_dihedrals, top_dihedral, atom_gaff):

        self.gaff_info = []
        for keys, values in gff_dihedrals.iteritems():
            self.gaff_info = [keys, values]

        for item in top_dihedral:
            d1 = item[0] + ' ' + item[1] + ' ' + item[2] + ' ' + item[3]
            a1, a2, a3, a4 = atom_gaff[item[0]], atom_gaff[item[1]], atom_gaff[
                item[2]], atom_gaff[item[3]]
            dihedral_label = (a1, a2, a3, a4)
            if dihedral_label[0] > dihedral_label[3]:
                dihedral_label = tuple(reversed(list(dihedral_label)))
            if dihedral_label in gff_dihedrals:
                self.topdihedralff[d1] = (
                dihedral_label, gff_dihedrals[dihedral_label])
            self.num_dih_types = len(set(self.topdihedralff.keys()))




    def _get_ff_bonds(self, gff_bonds, top_bond, atom_gaff):

        self.gaff_info = []
        for keys, values in gff_bonds.iteritems():
            self.gaff_info = [keys, values]
        for item in top_bond:
            d1 = item[0] + ' ' + item[1]
            a1, a2 = atom_gaff[item[0]], atom_gaff[item[1]]
            if (str(a1), str(a2)) in gff_bonds:
                self.topbondff[d1] = (
                (str(a1), str(a2)), gff_bonds[(str(a1), str(a2))])
            elif (str(a2), str(a1)) in gff_bonds:
                bond_type = tuple(sorted((str(a1), str(a2))))
                self.topbondff[d1] = ((str(a1), str(a2)), gff_bonds[bond_type])
                self.num_bond_types = len(set(self.topbondff.keys()))
                #print "TOPBONDFF",self.topbondff


    def _get_ff_angles(self, gff_angles, top_angle, atom_gaff):

        self.gaff_info = []
        for keys, values in gff_angles.iteritems():
            self.gaff_info = [keys, values]
        for item in top_angle:
            d1 = item[0] + ' ' + item[1] + ' ' + item[2]
            a1, a2, a3 = atom_gaff[item[0]], atom_gaff[item[1]], atom_gaff[
                item[2]]
            if (str(a1), str(a2), str(a3)) in gff_angles:
                self.topangleff[d1] = (((str(a1), str(a2), str(a3))),
                                       gff_angles[(str(a1), str(a2), str(a3))])
            elif (str(a3), str(a2), str(a1)) in gff_angles:
                angle_type = tuple(sorted((str(a1), str(a2), str(a3))))
                self.topangleff[d1] = (
                ((str(a1), str(a2), str(a3))), gff_angles[angle_type])
            self.num_ang_types = len(set(self.topangleff.keys()))

    def _get_ff_imdihedrals(self, gff_imdihedrals, top_imdihedral, atom_gaff):

        self.gaff_info = []
        for keys, values in gff_imdihedrals.iteritems():
            self.gaff_info = [keys, values]

        for item in top_imdihedral:
            d1 = item[0] + ' ' + item[1] + ' ' + item[2] + ' ' + item[3]
            a1, a2, a3, a4 = atom_gaff[item[0]], atom_gaff[item[1]], atom_gaff[
                item[2]], atom_gaff[item[3]]
            imdihedral_label = (a1, a2, a3, a4)
            if imdihedral_label[0] > imdihedral_label[3]:
                imdihedral_label = tuple(reversed(list(imdihedral_label)))
            if imdihedral_label in gff_imdihedrals:
                self.topimdihedralff[d1] = (
                imdihedral_label, gff_imdihedrals[imdihedral_label])
            self.num_imdih_types = len(set(self.topimdihedralff.keys()))























