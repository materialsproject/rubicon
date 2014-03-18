from pymatgen.core import bonds

__author__ = 'navnidhirajput'



from pymatgen.serializers.json_coders import MSONable
import glob

import os



class GFF(MSONable):

    """
    A force field library. Right now reads the output file from Antechamber
    and populate the FF library

    Args:
        bonds: store the bond distance (A) and spring constant (Kcal/molA2)
         in a dict.
        angles: store the bond distance (A) and spring constant
        (Kcal/mol*radian2) in a dict.
        dihedral: store the magnitude of torsion (Kcal/mol), phase offset in
        degree and the periodicity of torsion in a dict.
        imdihedrals: store improper dihedral information which include
        the magnitude of torsion (Kcal/mol), phase
        offset in degree and the periodicity of torsion in a dict
        vdws: store the van der waal radius (A) and van der wall depth for a
        given atom (Kcal/mol) in a dict.

        All the dictionaries are then converted into json serializable object
            Example:
            {'bonds': {'c-o': (648.0, 1.214)
            'angles': {'os-c-os': (76.45, 111.38)
            'dihedrals': {'X-c-os-X': (2.7, 180.0)
            'imdihedrals': {'o-os-c-os': (1.1, 180.0)}
            'vdws': {'c': (1.908, 0.086)}
    """



    def __init__(self,bonds={}, angles={}, dihedrals={},
                 imdihedrals={}, vdws={},masses={}):

        self.bonds = bonds
        self.angles=angles
        self.dihedrals=dihedrals
        self.imdihedrals=imdihedrals
        self.vdws=vdws
        self.masses=masses




    def read_forcefield_para(self,filename=None):

        bonds = dict()
        angles = dict()
        dihedrals = dict()
        imdihedrals = dict()
        vdws = dict()


        with open(filename) as f:
            bond_section = False
            angle_section = False
            dihedral_section = False
            imdihedral_section = False
            vdw_section = False

            for line in f.readlines():
                if line.startswith('BOND'):
                    bond_section = True
                    continue
                if bond_section:
                    if len(line.strip())==0:
                        bond_section = False
                        continue
                    token = line.split()
                    bond_type=token[0]+'-'+token[1]
                    bond_k_distance=float(token[2])
                    bond_distance=float(token[3])
                    bonds[bond_type]=(bond_k_distance,bond_distance)

                if line.startswith('ANGLE'):
                    angle_section = True
                    continue
                if angle_section:
                    if len(line.strip())==0:
                        angle_section = False
                        continue
                    token = line.split()
                    angle_type=token[0]+'-'+token[1]+'-'+token[2]
                    angle_k_distance=float(token[3])
                    angle_distance=float(token[4])
                    angles[angle_type]=(angle_k_distance,angle_distance)


                if line.startswith('DIHEDRAL'):
                    dihedral_section = True
                    continue
                if dihedral_section:
                    if len(line.strip())==0:
                        dihedral_section = False
                        continue
                    token = line.split()
                    dihedral_type=token[0]+'-'+token[1]+'-'+token[2]+'-'+token[3]
                    dihedral_k_distance=float(token[4])
                    dihedral_func_type=float(token[5])
                    dihedral_angle=float(token[6])
                    dihedrals[dihedral_type]=(dihedral_k_distance,dihedral_angle)


                if line.startswith('IMPROPER'):
                    imdihedral_section = True
                    continue
                if imdihedral_section:
                    if len(line.strip())==0:
                        imdihedral_section = False
                        continue


                    imdihedral_type=line[0:11]
                    imdihedral_distance=float(line[18:28])
                    imdihedral_angle=float(line[29:41])
                    imdihedrals[imdihedral_type]=(imdihedral_distance,imdihedral_angle)




                if line.startswith('NONBONDED'):
                    vdw_section = True
                    continue
                if vdw_section:
                    if len(line.strip())==0:
                        vdw_section = False
                        continue
                    if line.startswith("CUTNB"):
                        continue
                    if line.startswith("!"):
                        continue
                    token = line.split()
                    vdw_type=token[0]
                    epsilon=abs(float(token[2]))
                    sigma=float(token[3])
                    vdws[vdw_type]=(sigma,epsilon)
            self.bonds.update(bonds)
            self.angles.update(angles)
            self.dihedrals.update(dihedrals)
            self.imdihedrals.update(imdihedrals)
            self.vdws.update(vdws)



    def read_mass(self,filename=None):
        masses=dict()
        with open(filename) as f:
            mass_section = True
            for line in f.readlines():
                if line.startswith("*"):
                        continue
                if line.startswith(" "):
                        continue
                if mass_section:
                    if len(line.strip())==0:
                        mass_section = False
                        continue
                    token = line.split()
                    atom_type=token[2]
                    mass=float(token[3])
                    masses[atom_type]=(mass)
            self.masses.update(masses)




    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                 "bonds": self.bonds,
                 "angles": self.angles,
                 "dihedrals": self.dihedrals,
                 "imdihedrals": self.imdihedrals,
                 "vdws":self.vdws}


    @classmethod
    def from_dict(cls, d):
        return GFF(bonds=d["bonds"],
                    angles=d["angles"],
                    dihedrals=d["dihedrals"],
                    imdihedrals=d["imdihedrals"],
                    vdws=d["vdws"]
                    )



class GFF_library(GFF):

    def append_gff(self):
      """
        this will append the FF library after reading the parameters from
        different molecules.
        """
      gff_sec=GFF()

      self.bonds.update(gff_sec.bonds)
      self.angles.update(gff_sec.angles)
      self.dihedrals.update(gff_sec.dihedrals)
      self.imdihedrals.update(gff_sec.imdihedrals)
      self.vdws.update(gff_sec.vdws)
      self.masses.update(gff_sec.masses)

      #print gff_sec.bonds
      #print gff_sec.angles
      #print gff_sec.dihedrals
      #print gff_sec.masses

      return gff_sec


#class TopFF()


































