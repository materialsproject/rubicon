#from rubicon.gff.antechamberio import AntechamberRunner


class AtomTyper():

    def get_ffmol(self, mols):
        pass

"""
class AntechamberAtomTyper(AtomTyper):
    def __init__(self, filename):
        self.atom_types
        pass

    def get_ffmol(self, mols):
        ant = AntechamberRunner(mols)
        gff_list, top_list = ant._run_antechamber('mol.pdb', mols)
        #run antechamber
        #save atom types to self.atom_types

    def get_atom_types(self):
        return self.atom_types

    def read_atom_index(self,filename=None):

        with open(filename) as f:

            for line in f.readlines():
                token = line.split()
                if token[0]=='ATOM':
                    index=int(token[1])
                    atom_name=token[2]
                    atom_gaff=token[9]
                    self.atom_index[index]=atom_name

                    self.atom_index_gaff[index]=atom_gaff


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


    def _convert_to_pdb(self, molecule, filename=None):

        generate pdb file for a given molecule

        write_mol(molecule, filename)

    def _run_parmchk(self, filename=None):

        run parmchk using ANTECHAMBER_AC.AC file

        Args:
            filename = pdb file of the molecule

        command_parmchk = (
        'parmchk -i ' + filename + ' -f ac -o mol.frcmod -a Y')
        return_cmd = subprocess.call(shlex.split(command_parmchk))
        return return_cmd


    def _run_antechamber(self, filename=None, mols=[]):

        generate and run antechamber command for specified pdb file

        Args:
            filename = pdb file of the molecule

        scratch = tempfile.gettempdir()
        return_cmd = None

        with ScratchDir(scratch,
                        copy_to_current_on_exit=ANTECHAMBER_DEBUG) as d:
            gff_list = []
            top_list=[]

            for mol in mols:

                self._convert_to_pdb(mol, 'mol.pdb')
                command = (
                'antechamber -i ' + filename + ' -fi pdb -o ' + filename[:-4] +
                " -fo charmm")

                return_cmd = subprocess.call(shlex.split(command))
                self.molname = filename.split('.')[0]
                self._run_parmchk('ANTECHAMBER_AC.AC')

                #gff = self._parse_output()
                self.read_atom_index('ANTECHAMBER_AC.AC')
                self.read_atomType('ANTECHAMBER_AC.AC')
                top = TopMol.from_file('mol.rtf')
                #my_gff.read_forcefield_para('mol.frcmod')


                my_gff = Gff.from_forcefield_para('mol.frcmod')

                self.read_atomType('ANTECHAMBER_AC.AC')
                #my_gff = self._parse_output()
                #self._get_ff_bonds(my_gff.bonds, top.bonds, self.atom_gaff)
                #self._get_ff_angles(my_gff.angles, top.angles,
                #                    self.atom_gaff)
                #self._get_ff_dihedrals(my_gff.dihedrals, top.dihedrals,
                #                       self.atom_gaff)
                #self._get_ff_imdihedrals(my_gff.imdihedrals, top.imdihedrals,
                #                         self.atom_gaff)
                print "MYGFFBONDS",my_gff.bonds


                #gff_list.append(copy.deepcopy(my_gff))
                print '*', type(my_gff)
                print '**', gff_list
                gff_list.append(my_gff)
                print '***', [x.bonds for x in gff_list]
                #print "GFFLIST0",gff_list[0].bonds,gff_list
                top_list.append(top)

            #print "GFFLIST0",gff_list[0].bonds
            #print "TOPLIST0",top_list[2].bonds

            return  gff_list,top_list
"""