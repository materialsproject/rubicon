import os
from random import Random
from time import time
import math
import itertools

import inspyred
import openbabel as ob
from pymatgen.core.structure import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.qchemio import QcOutput
from pymatgen.io.smartio import write_mol
import numpy as np

from rubicon.utils.ion_arranger.hard_sphere_energy_evaluators import HardSphereElectrostaticEnergyEvaluator, \
    AtomicRadiusUtils
from rubicon.utils.ion_arranger.semi_emprical_qm_energy_evaluator import SemiEmpricalQuatumMechanicalEnergyEvaluator


class IonPlacer():





    def __init__(self, molecule, fragments, nums_fragments, energy_evaluator,
                 prng=None, seed=None):
        self.prng = prng if prng else Random()
        self.seed = seed if seed else time()
        self.prng.seed(self.seed)
        self.final_pop = None
        self.best = None
        self.ea = inspyred.swarm.PSO(self.prng)
        self.ea.terminator = inspyred.ec.terminators.evaluation_termination
        self.ea.topology = inspyred.swarm.topologies.ring_topology
        self.ea.observer = inspyred.ec.observers.best_observer
        self.molecule = molecule
        self.fragments = fragments
        self.nums_fragments = nums_fragments
        self.mol_coords = self.normalize_molecule(self.molecule)
        for frag in self.fragments:
            self.normalize_molecule(frag)
        self.bounder = self.get_bounder(self.mol_coords, self.fragments, self.nums_fragments)
        self.best_pymatgen_mol = None
        self.playing_time = None
        self.energy_evaluator = energy_evaluator

    @classmethod
    def get_bounder(cls, mol_coords, fragments, nums_fragments):
        lower_bound = []
        upper_bound = []
        frag_dimensions = []
        mol_radius = max([math.sqrt(sum([x**2 for x in c]))
                          for c in mol_coords])
        frag_dimensions.append(tuple([mol_radius, 1]))
        frag_atom_counts = []
        for frag, num_frag in zip(fragments, nums_fragments):
            num_atoms = frag.NumAtoms()
            frag_coords = cls.get_mol_coords(frag)
            frag_radius = max([math.sqrt(sum([x**2 for x in c]))
                               for c in frag_coords])
            frag_dimensions.append(tuple([frag_radius, num_frag]))
            frag_atom_counts.append(num_atoms)
        max_length = 2 * sum([r*n for r, n in frag_dimensions])
        x_min = max_length + 5.0 * sum(nums_fragments)
        for num_frag, num_atoms in zip(nums_fragments, frag_atom_counts):
            single_frag_lower_bound = [-x_min] * 3
            single_frag_upper_bound = [x_min] * 3
            if num_atoms > 1:
                single_frag_lower_bound.extend([0.0, -math.pi])
                single_frag_upper_bound.extend([math.pi, math.pi])
            lower_bound.extend(single_frag_lower_bound * num_frag)
            upper_bound.extend(single_frag_upper_bound * num_frag)
        bounder = inspyred.ec.Bounder(lower_bound, upper_bound)
        return bounder

    @staticmethod
    def normalize_molecule(mol):
        mol.Center()
        coords = IonPlacer.get_mol_coords(mol)
        return coords

    @staticmethod
    def get_mol_coords(mol):
        coords = []
        num_atoms = mol.NumAtoms()
        for i in range(1, num_atoms + 1):
            a = mol.GetAtom(i)
            coords.append([a.GetX() * AtomicRadiusUtils.angstrom2au,
                           a.GetY() * AtomicRadiusUtils.angstrom2au,
                           a.GetZ() * AtomicRadiusUtils.angstrom2au])
        return coords

    @staticmethod
    def get_mol_species(mol):
        species = []
        num_atoms = mol.NumAtoms()
        element_table = ob.OBElementTable()
        for i in range(1, num_atoms + 1):
            a = mol.GetAtom(i)
            atomic_num = a.GetAtomicNum()
            symbol = element_table.GetSymbol(atomic_num)
            species.append(symbol)
        return species

    @staticmethod
    def rotate(mol, theta, phi):
        m_theta = np.array([[1.0, 0.0, 0.0],
                            [0.0, math.cos(theta), -math.sin(theta)],
                            [0.0, math.sin(theta), math.cos(theta)]])
        m_phi = np.array([[math.cos(phi), 0.0, math.sin(phi)],
                          [0.0, 1.0, 0.0],
                          [-math.sin(phi), 0.0, math.cos(phi)]])
        m = np.dot(m_phi, m_theta)
        m_list = list(itertools.chain(*m))
        mol.Rotate(ob.double_array(m_list))
        return mol

    def decode_solution(self, x):
        fragments_coords = []
        fragments_atom_counts = [frag.NumAtoms() for frag in self.fragments]
        xi = 0
        for frag, num_frag, num_atoms in zip(self.fragments, self.nums_fragments, fragments_atom_counts):
            for i in range(num_frag):
                cc = ob.OBMol(frag)
                tx, ty, tz = [t/AtomicRadiusUtils.angstrom2au for t in x[xi: xi+3]]
                xi += 3
                if num_atoms > 1:
                    theta = x[xi]
                    xi += 1
                    phi = x[xi]
                    xi += 1
                    self.rotate(cc, theta, phi)
                cc.Translate(ob.vector3(tx, ty, tz))
                fragments_coords.append(self.get_mol_coords(cc))
        return fragments_coords

    @staticmethod
    def generate_conformers(random, args):
        # generator
        lower_bound = args['_ec'].bounder.lower_bound
        upper_bound = args['_ec'].bounder.upper_bound
        return [random.uniform(l, u) for l, u in zip(lower_bound, upper_bound)]

    # noinspection PyUnusedLocal
    def evaluate_conformers(self, candidates, args):
        # evaluator
        fitness = []
        for c in candidates:
            fragments_coords = self.decode_solution(c)
            energy = self.energy_evaluator.calc_energy(fragments_coords)
            fitness.append(energy)
        return fitness

    def place(self, max_evaluations=30000, pop_size=100, neighborhood_size=5):
        t1 = time()
        self.final_pop = self.ea.evolve(generator=self.generate_conformers,
                                        evaluator=self.evaluate_conformers,

                                        pop_size=pop_size,
                                        bounder=self.bounder,
                                        maximize=False,
                                        max_evaluations=max_evaluations,
                                        neighborhood_size=neighborhood_size)
        self.best = max(self.final_pop)
        # max means best, not necessarily smallest
        best_fragments_coords = self.decode_solution(
            self.best.candidate)
        mol_elements = self.get_mol_species(self.molecule)
        fragments_elements = []
        for frag, num_frag in zip(self.fragments, self.nums_fragments):
            fragments_elements.extend([self.get_mol_species(frag)]*num_frag)
        species = []
        coords_au = []
        species.extend(mol_elements)
        coords_au.extend(self.mol_coords)
        for elements, c in zip(fragments_elements, best_fragments_coords):
            species.extend(elements)
            coords_au.extend(c)
        coords_ang = [[x/AtomicRadiusUtils.angstrom2au for x in c] for c in coords_au]
        self.best_pymatgen_mol = Molecule(species, coords_ang)
        t2 = time()
        self.playing_time = t2 - t1
        return self.best_pymatgen_mol


def main():
    def gcd(a, b):
        if b == 0:
                return a
        else:
                return gcd(b, a % b)

    def lcm(a, b):
            return a * b / gcd(a, b)
    import argparse
    parser = argparse.ArgumentParser(
        description="Place salt around a molecule")
    parser.add_argument("-m", "--molecule", dest="molecule", type=str,
                        required=True,
                        help="the file name of molecule")
    parser.add_argument("-l", "--ligand", dest="fragments", type=str, nargs='+',
                        required=True,
                        help="the list of fragment file names to to be placed around the molecule")
    parser.add_argument("-n", "--nums_fragments", dest="nums_fragments", type=int, nargs='+',
                        required=True,
                        help="the number of each fragment, the order must be the same with FRAGMENTS")
    parser.add_argument("-c", "--charge", dest="charge", type=int,
                        required=True,
                        help="total charge of the system")
    parser.add_argument("-o", "--outputfile", dest="outputfile", type=str,
                        required=True,
                        help="the file name of the aligned conformer")
    parser.add_argument("-i", "--iterations", dest="iterations", type=int,
                        default=600,
                        help="maximum number of evaluations")
    parser.add_argument("-s", "--size", dest="size", type=int, default=15,
                        help="population size")
    parser.add_argument("-k", "--num_neighbours", dest="num_neighbours", type=int, default=2,
                        help="number of neighbours")
    parser.add_argument("-e", "--evaluator", dest="evaluator", type=str, default="hardsphere",
                        choices=["hardsphere", "sqm"], help="Energy Evaluator")
    options = parser.parse_args()
    if options.evaluator == 'hardsphere':
        qcout_molecule = QcOutput(options.molecule)
        qcout_cation = QcOutput(options.cation)
        qcout_anion = QcOutput(options.anion)
        total_charge_cation = qcout_cation.data[0]["molecules"][-1].charge
        total_charge_anion = qcout_anion.data[0]["molecules"][-1].charge
        total_charge_mol = qcout_molecule.data[0]["molecules"][-1].charge
        num_lcm = lcm(total_charge_cation, -total_charge_anion)
        num_cation = num_lcm/total_charge_cation
        num_anion = num_lcm/-total_charge_anion
        pymatgen_mol_molecule = qcout_molecule.data[0]["molecules"][-1]
        pymatgen_mol_cation = qcout_cation.data[0]["molecules"][-1]
        pymatgen_mol_anion = qcout_anion.data[0]["molecules"][-1]
        # noinspection PyProtectedMember
        molecule = BabelMolAdaptor(pymatgen_mol_molecule)._obmol
        # noinspection PyProtectedMember
        obmol_cation = BabelMolAdaptor(pymatgen_mol_cation)._obmol
        # noinspection PyProtectedMember
        obmol_anion = BabelMolAdaptor(pymatgen_mol_anion)._obmol
        energy_evaluator = HardSphereElectrostaticEnergyEvaluator.from_qchem_output(
            qcout_molecule, qcout_cation, qcout_anion)
        fragments = [obmol_cation, obmol_anion]
    else:
        # noinspection PyProtectedMember
        molecule = BabelMolAdaptor.from_file(options.molecule,
                                             os.path.splitext(options.molecule)[1][1:])._obmol
        fragments = []
        for frag_file in options.fragments:
            file_format = os.path.splitext(frag_file)[1][1:]
            # noinspection PyProtectedMember
            fragments.append(BabelMolAdaptor.from_file(frag_file, file_format)._obmol)
        energy_evaluator = SemiEmpricalQuatumMechanicalEnergyEvaluator()
    placer = IonPlacer(molecule=molecule, fragments=fragments, nums_fragments=options.nums_fragments,
                       energy_evaluator=energy_evaluator)
    placer.place(max_evaluations=options.iterations,
                 pop_size=options.size,
                 neighborhood_size=options.num_neighbours)
    print 'It took {:.1f} seconds to place the salt'.format(placer
                                                            .playing_time)
    write_mol(placer.best_pymatgen_mol, options.outputfile)

if __name__ == '__main__':
    main()