import math
import itertools
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from pymatgen.io.babelio import BabelMolAdaptor
from rubicon.utils.ion_arranger.energy_evaluator import EnergyEvaluator
import openbabel as ob
import numpy as np

__author__ = 'xiaohuiqu'

class AtomicRadiusUtils(object):
    metals = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V',
              'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr',
              'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
              'In', 'Sn', 'Sb', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
              'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
              'Hf', 'Ta', 'W', 'Re',  'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
              'Pb']
    angstrom2au = 1.0/0.52917721092
    def __init__(self, covalent_radius_scale=2.0, metal_radius_scale=0.5):
        self.covalent_radius_scale = covalent_radius_scale
        self.metal_radius_scale = metal_radius_scale

    def get_radius(self, mol):
        radius = []
        ref_radius = CovalentRadius.radius
        num_atoms = mol.NumAtoms()
        element_table = ob.OBElementTable()
        for i in range(1, num_atoms + 1):
            a = mol.GetAtom(i)
            atomic_num = a.GetAtomicNum()
            symbol = element_table.GetSymbol(atomic_num)
            if symbol in self.metals:
                scale = self.metal_radius_scale
            else:
                scale = self.covalent_radius_scale
            rad = ref_radius[symbol] * self.angstrom2au * scale
            radius.append(rad)
        return radius

class CoulumbEnergyEvaluator(EnergyEvaluator):

    def __init__(self, mol_coords, mol_charges, cation_charges, anion_charges,
                 mol_radius, cation_radius, anion_radius):
        super(CoulumbEnergyEvaluator, self).__init__(mol_coords)
        self.mol_charges = mol_charges
        self.cation_charges = cation_charges
        self.anion_charges = anion_charges
        self.mol_radius = mol_radius
        self.cation_radius = cation_radius
        self.anion_radius = anion_radius

    def calc_energy(self, cation_coords, anion_coords):
        energy = 0.0
        for frag_coords in cation_coords:
            energy += self._pair_energy(
                self.mol_coords, self.mol_charges, self.mol_radius,
                frag_coords, self.cation_charges, self.cation_radius)
        for frag_coords in anion_coords:
            energy += self._pair_energy(
                self.mol_coords, self.mol_charges, self.mol_radius,
                frag_coords, self.anion_charges, self.anion_radius)
        for cc in cation_coords:
            for ac in anion_coords:
                energy += self._pair_energy(
                    cc, self.cation_charges, self.cation_radius,
                    ac, self.anion_charges, self.anion_radius)
        return energy

    @staticmethod
    def _pair_energy(coords1, charges1, radius1, coords2, charges2,
                    radius2):
        energy = 0.0
        for coord1, charge1, rad1 in zip(coords1, charges1, radius1):
            for coord2, charge2, rad2 in zip(coords2, charges2, radius2):
                distance = math.sqrt(sum([(x1-x2)**2 for x1, x2
                                          in zip(coord1, coord2)]))
                if distance <= rad1 + rad2:
                    continue
                else:
                    electrostatic_energy = charge1 * charge2 / distance
                    energy += electrostatic_energy
        return energy

class HardSphereEnergyEvaluator(EnergyEvaluator):
    overlap_energy = 1.0E4

    def __init__(self, mol_coords, mol_radius, cation_radius, anion_radius):
        super(HardSphereEnergyEvaluator, self).__init__(mol_coords)
        self.mol_radius = mol_radius
        self.cation_radius = cation_radius
        self.anion_radius = anion_radius

    def calc_energy(self, cation_coords, anion_coords):
        energies = []
        for frag_coords in cation_coords:
            energies.append(self._pair_energy(self.mol_coords, self.mol_radius,
                frag_coords, self.cation_radius))
        for frag_coords in anion_coords:
            energies.append(self._pair_energy(self.mol_coords, self.mol_radius,
                frag_coords, self.anion_radius))
        for cc in cation_coords:
            for ac in anion_coords:
                energies.append(self._pair_energy(cc, self.cation_radius,
                    ac, self.anion_radius))
        energy = min(energies)
        return energy

    @classmethod
    def _pair_energy(cls, coords1, radius1, coords2, radius2):
        energy = 0.0
        for coord1, rad1 in zip(coords1, radius1):
            for coord2, rad2 in zip(coords2, radius2):
                distance = math.sqrt(sum([(x1-x2)**2 for x1, x2
                                          in zip(coord1, coord2)]))
                if distance <= rad1 + rad2:
                    energy += cls.overlap_energy
        return energy

class ContactDetector(object):
    def __init__(self, mol_coords, mol_radius, cation_radius, anion_radius, cap=0.0):
        self.mol_coords = mol_coords
        self.mol_radius = mol_radius
        self.cation_radius = cation_radius
        self.anion_radius = anion_radius

    def is_contact(self, cation_coords, anion_coords):
        pass


    def _get_contact_matrix(self, cation_coords, anion_coords):
        fragments = [(self.mol_coords, self.mol_radius)]
        fragments.extend(zip(cation_coords, [self.cation_radius] * len(cation_coords)))
        fragments.extend(zip(anion_coords, [self.anion_radius] * len(anion_coords)))
        num_frag = len(fragments)
        fragments = [(c, r, i) for i, (c, r) in enumerate(fragments)]
        contact_matrix = np.identity(num_frag, dtype=int)
        frag_pair = itertools.combinations(fragments, r=2)
        for p in frag_pair:
            ((c1s, r1s, i1), (c2s, r2s, i2)) = p
            contact = 0
            for c1, r1, c2, r2 in zip(c1s, r1s, c2s, r2s):
                distance = math.sqrt(sum([(x1-x2)**2 for x1, x2
                                          in zip(c1, c2)]))
                if distance <= r1 + r2:
                    contact = 1
                    break
            contact_matrix[i1, i2] = contact
            contact_matrix[i2, i1] = contact
        return contact_matrix



class GravitationalEnergyEvaluator(EnergyEvaluator):

    gravity = 0.001

    def __init__(self, mol_coords, mol_radius, cation_radius, anion_radius):
        super(GravitationalEnergyEvaluator, self).__init__(mol_coords)
        self.mol_radius = mol_radius
        self.cation_radius = cation_radius
        self.anion_radius = anion_radius

    def calc_energy(self, cation_coords, anion_coords):
        energy = 0.0
        for frag_coords in cation_coords:
            energy += self._gravitational_energy(frag_coords, self.cation_radius)
        for frag_coords in anion_coords:
            energy += self._gravitational_energy(frag_coords, self.anion_radius)
        return energy

    def _gravitational_energy(self, ion_coords, ion_radius):
        grav_distance = 10000.0
        for mol_coord, mol_rad in zip(self.mol_coords, self.mol_radius):
            for ion_coord, ion_rad in zip(ion_coords, ion_radius):
                distance = math.sqrt(sum([(x1-x2)**2 for x1, x2
                                          in zip(mol_coord, ion_coord)]))
                if distance > mol_rad + ion_rad:
                    gd = distance - (mol_rad + ion_rad)
                    if gd < grav_distance:
                        grav_distance = gd
                else:
                    return 0.0
        if grav_distance > 9999.0:
            return 0.0
        else:
            return self.gravity * grav_distance

class HardSphereElectrostaticEnergyEvaluator(EnergyEvaluator):
    def __init__(self, mol_coords, mol_radius, cation_radius, anion_radius,
                 mol_charges, cation_charges, anion_charges):
        super(HardSphereElectrostaticEnergyEvaluator, self).__init__(mol_coords)
        self.hardsphere = HardSphereEnergyEvaluator(mol_coords, mol_radius, cation_radius, anion_radius)
        self.coulumb = CoulumbEnergyEvaluator(mol_coords, mol_charges, cation_charges, anion_charges,
                                              mol_radius, cation_radius, anion_radius)
        self.gravitation = GravitationalEnergyEvaluator(mol_coords, mol_radius, cation_radius, anion_radius)
        self.calculators = [self.hardsphere, self.coulumb, self.gravitation]

    def calc_energy(self, cation_coords, anion_coords):
        energy = 0.0
        for c in self.calculators:
            energy += c.calc_energy(cation_coords, anion_coords)
            if energy > HardSphereEnergyEvaluator.overlap_energy * 0.9:
                return energy
        return energy

    @classmethod
    def from_qchem_output(cls, mol_qcout, cation_qcout, anion_qcout,
                          covalent_radius_scale=2.0, metal_radius_scale=0.5):
        from rubicon.utils.ion_arranger.ion_arranger import IonPlacer
        mol_charges = mol_qcout.data[0]["charges"]["chelpg"]
        cation_charges = cation_qcout.data[0]["charges"]["chelpg"]
        anion_charges = anion_qcout.data[0]["charges"]["chelpg"]
        pymatgen_mol_molecule = mol_qcout.data[0]["molecules"][-1]
        pymatgen_mol_cation = cation_qcout.data[0]["molecules"][-1]
        pymatgen_mol_anion = anion_qcout.data[0]["molecules"][-1]
        # noinspection PyProtectedMember
        obmol_mol = BabelMolAdaptor(pymatgen_mol_molecule)._obmol
        # noinspection PyProtectedMember
        obmol_cation = BabelMolAdaptor(pymatgen_mol_cation)._obmol
        # noinspection PyProtectedMember
        obmol_anion = BabelMolAdaptor(pymatgen_mol_anion)._obmol
        mol_coords = IonPlacer.normalize_molecule(obmol_mol)
        rad_util = AtomicRadiusUtils(covalent_radius_scale, metal_radius_scale)
        mol_radius = rad_util.get_radius(obmol_mol)
        cation_radius = rad_util.get_radius(obmol_cation)
        anion_radius = rad_util.get_radius(obmol_anion)
        evaluator = HardSphereElectrostaticEnergyEvaluator(
            mol_coords, mol_radius, cation_radius, anion_radius,
            mol_charges, cation_charges, anion_charges)
        return evaluator
