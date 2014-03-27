

__author__ = 'navnidhirajput'


""" A wrapper for Antechamber which generates force field files
    for a specified molecule using pdb file as input
"""


import subprocess
from pymatgen import write_mol
import shlex
from gff import GFF
from monty.io import ScratchDir
import tempfile
from pymatgen.packmol.packmol import PackmolRunner
from rubicon.gff.lamppsio import LMPInput
from rubicon.gff.topology import AC, TopMol

ANTECHAMBER_DEBUG = False



class Antechamber():
    """
    A wrapper for Antechamber software

    """

    def __init__(self,molecule=None):
        """
        Args:
        molecule: The input molecule. Default is None
        """
        self.molecule=molecule
        self.topbondFF=dict()
        self.topangleFF=dict()
        self.topdihedralFF=dict()
        self.topimdihedralFF=dict()


    def convert_to_pdb(self,molecule,filename):
        """
        generate pdb file for a given molecule
        """
        write_mol(molecule,filename)
        self.filename=filename

    def run_antechamber(self,filename=None,mols=[]):
        """
        generate and run antechamber command for specified pdb file

        Args:
            filename = pdb file of the molecule
        """
        scratch = tempfile.gettempdir()
        return_cmd=None
        my_gff=None
        my_ant=None
        with ScratchDir(scratch,copy_from_current_on_enter = ANTECHAMBER_DEBUG) as d:

            gff_list = []
            my_lammps_list=[]
            for mol in mols:
                #my_ant = Antechamber(mol)
                self.convert_to_pdb(mol, 'mol.pdb')
                command=('antechamber -i ' +filename +' -fi pdb -o ' +filename[:-4]+
                 " -fo charmm")
                return_cmd=subprocess.call(shlex.split(command))
                self.molname = filename.split('.')[0]
                self.run_parmchk('ANTECHAMBER_AC.AC')
                gff = self.parse_output()
                ac = AC()
                ac.read_atomIndex('ANTECHAMBER_AC.AC')
                ac.read_atomType('ANTECHAMBER_AC.AC')
                pmr = PackmolRunner([mol, mol], [{"number":1,"inside box":[0.,0.,0.,40.,40.,40.]}, {"number":2}])

                top = TopMol.from_file('mol.rtf')

                my_gff = GFF()
                my_gff.read_forcefield_para('mol.frcmod')

                atom_gaff = AC()
                atom_gaff.read_atomType('ANTECHAMBER_AC.AC')

                self.get_FF_bonds(my_gff.bonds, top.bonds, atom_gaff.atom_gaff)
                self.get_FF_angles(my_gff.angles, top.angles, atom_gaff.atom_gaff)
                self.get_FF_dihedrals(my_gff.dihedrals, top.dihedrals, atom_gaff.atom_gaff)
                self.get_FF_imdihedrals(my_gff.imdihedrals, top.imdihedrals, atom_gaff.atom_gaff)

                gff_list.append(gff)

                my_lampps=LMPInput()

                my_lampps.set_coeff(my_gff,top,pmr,self)
                my_lampps.set_atom(pmr,my_gff,ac)
                my_lampps.set_bonds(pmr,my_gff,ac,top)
                my_lampps.set_angles(pmr,my_gff,ac,top)
                my_lampps.set_dihedrals(pmr,my_gff,ac,top,self)
                my_lampps.set_imdihedrals(pmr,my_gff,ac,top)
                my_lammps_list.append(my_lampps)


        return return_cmd,my_gff,my_ant,my_lammps_list,gff_list

    def run_parmchk(self,filename=None):
        """
        run parmchk using ANTECHAMBER_AC.AC file

        Args:
            filename = pdb file of the molecule
        """
        command_parmchk=('parmchk -i '+ filename +' -f ac -o mol.frcmod -a Y')
        return_cmd=subprocess.call(shlex.split(command_parmchk))
        return return_cmd


    def parse_output(self):
        """
        create an object of forcefield_type_molecule
        and call the function
        """
        ffm=GFF()
        ffm.read_forcefield_para(self.molname+".frcmod")
        return ffm


    def get_FF_bonds(self,gff_bonds,top_bond,atom_gaff):
        self.gaff_info=[]
        for keys, values in gff_bonds.iteritems():
                self.gaff_info=[keys,values]
        for item in top_bond:
            d1=item[0]+' '+ item[1]
            a1,a2 = atom_gaff[item[0]],atom_gaff[item[1]]
            if (str(a1),str(a2)) in gff_bonds:
                self.topbondFF[d1]=((str(a1),str(a2)),gff_bonds[(str(a1),str(a2))])
            elif (str(a2),str(a1)) in gff_bonds:
                bond_type=tuple(sorted((str(a1),str(a2))))
                self.topbondFF[d1]=((str(a1),str(a2)),gff_bonds[bond_type])
                self.num_bond_types = len(set(self.topbondFF.keys()))


    def get_FF_angles(self,gff_angles,top_angle,atom_gaff):
        self.gaff_info=[]
        for keys, values in gff_angles.iteritems():
                self.gaff_info=[keys,values]
        for item in top_angle:
            d1=item[0]+' '+ item[1]+' '+item[2]
            a1,a2,a3 = atom_gaff[item[0]],atom_gaff[item[1]],atom_gaff[item[2]]
            if (str(a1),str(a2),str(a3)) in gff_angles:
                self.topangleFF[d1]=(((str(a1),str(a2),str(a3))),gff_angles[(str(a1),str(a2),str(a3))])
            elif (str(a3),str(a2),str(a1)) in gff_angles:
                angle_type = tuple(sorted((str(a1),str(a2),str(a3))))
                self.topangleFF[d1]=(((str(a1),str(a2),str(a3))),gff_angles[angle_type])
            self.num_ang_types = len(set(self.topangleFF.keys()))



    def get_FF_dihedrals(self,gff_dihedrals,top_dihedral,atom_gaff):
        self.gaff_info=[]
        for keys, values in gff_dihedrals.iteritems():
                self.gaff_info=[keys,values]
        for item in top_dihedral:
            d1=item[0]+' '+ item[1]+' '+item[2]+ ' '+item[3]
            a1,a2,a3,a4 = atom_gaff[item[0]],atom_gaff[item[1]],atom_gaff[item[2]],atom_gaff[item[3]]
            dihedral_label=(a1,a2,a3,a4)
            if dihedral_label[0]>dihedral_label[3]:
                        dihedral_label=tuple(reversed(list(dihedral_label)))
            if dihedral_label in gff_dihedrals:
                self.topdihedralFF[d1]=(dihedral_label,gff_dihedrals[dihedral_label])
            self.num_dih_types = len(set(self.topdihedralFF.keys()))



    def get_FF_imdihedrals(self,gff_imdihedrals,top_imdihedral,atom_gaff):
        self.gaff_info=[]
        for keys, values in gff_imdihedrals.iteritems():
                self.gaff_info=[keys,values]
        for item in top_imdihedral:
            d1=item[0]+' '+ item[1]+' '+item[2]+ ' '+item[3]
            a1,a2,a3,a4 = atom_gaff[item[0]],atom_gaff[item[1]],atom_gaff[item[2]],atom_gaff[item[3]]
            imdihedral_label=(a1,a2,a3,a4)
            if imdihedral_label[0]> imdihedral_label[3]:
                        imdihedral_label=tuple(reversed(list(imdihedral_label)))
            if imdihedral_label in gff_imdihedrals:
                self.topimdihedralFF[d1]=(imdihedral_label,gff_imdihedrals[imdihedral_label])
            self.num_imdih_types = len(set(self.topimdihedralFF.keys()))

    @staticmethod
    def run(mols):

        gff_list = []
        for mol in mols:
            my_ant = Antechamber(mol)
            my_ant.convert_to_pdb(mol, 'mol.pdb')
            my_ant.run_antechamber('mol.pdb')
            my_ant.run_parmchk('ANTECHAMBER_AC.AC')
            gff = my_ant.parse_output()
            ac = AC()
            ac.read_atomIndex('ANTECHAMBER_AC.AC')
            ac.read_atomType('ANTECHAMBER_AC.AC')
            pmr = PackmolRunner([mol, mol], [{"number":1,"inside box":[0.,0.,0.,40.,40.,40.]}, {"number":2}])

            top = TopMol.from_file('mol.rtf')

            my_gff = GFF()
            my_gff.read_forcefield_para('mol.frcmod')

            atom_gaff = AC()
            atom_gaff.read_atomType('ANTECHAMBER_AC.AC')

            my_ant.get_FF_bonds(my_gff.bonds, top.bonds, atom_gaff.atom_gaff)
            my_ant.get_FF_angles(my_gff.angles, top.angles, atom_gaff.atom_gaff)
            my_ant.get_FF_dihedrals(my_gff.dihedrals, top.dihedrals, atom_gaff.atom_gaff)
            my_ant.get_FF_imdihedrals(my_gff.imdihedrals, top.imdihedrals, atom_gaff.atom_gaff)


            gff_list.append(gff)

            my_lampps=LMPInput()

            my_lampps.set_coeff(my_gff,top,pmr,my_ant)
            my_lampps.set_atom(pmr,my_gff,ac)
            my_lampps.set_bonds(pmr,my_gff,ac,top)
            my_lampps.set_angles(pmr,my_gff,ac,top)
            my_lampps.set_dihedrals(pmr,my_gff,ac,top,my_ant)
            my_lampps.set_imdihedrals(pmr,my_gff,ac,top)
            return my_gff,my_ant,my_lampps








