import copy
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
import numpy as np
import simplerandom.random as srr

from rubicon.utils.ion_arranger.hard_sphere_energy_evaluators import HardSphereElectrostaticEnergyEvaluator, \
    AtomicRadiusUtils
from rubicon.utils.ion_arranger.semi_emprical_qm_energy_evaluator import SemiEmpricalQuatumMechanicalEnergyEvaluator


class IonPlacer():





    def __init__(self, molecule, fragments, nums_fragments, energy_evaluator,
                 prng="kiss", random_seed=None, taboo_tolerance_ang=1.0, taboo_tolerance_particle_ratio=0.5,
                 topology="ring", initial_guess="breadth", bound_setter="chain", always_write_best=False,
                 max_generations_each_conformer=100):
        if prng == "kiss":
            self.prng = srr.KISS()
        elif prng == "python":
            self.prng = Random()
        else:
            raise Exception("Random number generator is not supported")
        self.seed = random_seed if random_seed else int(time())
        self.prng.seed(self.seed)
        print "Random Seed:", self.seed
        self.final_pop = None
        self.best = None
        self.ea = inspyred.swarm.PSO(self.prng)
        self.ea.terminator = inspyred.ec.terminators.evaluation_termination
        if topology == "ring":
            self.ea.topology = inspyred.swarm.topologies.ring_topology
        elif topology == "star":
            self.ea.topology = inspyred.swarm.topologies.star_topology
        else:
            raise ValueError("only RING and STAR topology are supported")
        self.initial_guess = initial_guess
        self.ea.observer = self.fitness_observer
        self.molecule = molecule
        self.fragments = fragments
        self.nums_fragments = nums_fragments
        self.mol_coords = self.normalize_molecule(self.molecule)
        for frag in self.fragments:
            self.normalize_molecule(frag)
        self.bounder = self.get_bounder(self.mol_coords, self.fragments, self.nums_fragments, bound_setter=bound_setter)
        if self.initial_guess == "center":
            self.max_radius = self.get_max_radius(self.mol_coords, self.fragments, self.nums_fragments)
        elif self.initial_guess == "volume":
            self.max_radius = self.get_equivalent_radius(self.mol_coords, self.fragments, self.nums_fragments)
        else:
            self.max_radius = None
        self.best_pymatgen_mol = None
        self.playing_time = None
        self.energy_evaluator = energy_evaluator
        self.taboo_tolerance_au = taboo_tolerance_ang * AtomicRadiusUtils.angstrom2au
        self.taboo_tolerance_particle_ratio = taboo_tolerance_particle_ratio
        self.bound_setter = bound_setter
        self.always_write_best = always_write_best
        self.reevaluate_fitness = False
        self.max_generations_each_conformer = max_generations_each_conformer
        self.conformer_staring_generation = None

    def fitness_observer(self, population, num_generations, num_evaluations, args):
        current_fitness = [p.fitness for p in population]
        short_current_fitness = [round(x, 2) for x in current_fitness]
        print "BEST FITNESS", min(short_current_fitness)
        print "CURRENT FITNESS", short_current_fitness
        archieved_fitness = [p.fitness for p in self.ea.archive]
        short_archieved_fitness = [round(x, 2) for x in archieved_fitness]
        print "ARCHIEVED FITNESS", short_archieved_fitness
        print

    @classmethod
    def get_max_radius(cls, mol_coords, fragments, nums_fragments):
        mol_radius = max([math.sqrt(sum([x**2 for x in c]))
                          for c in mol_coords])
        max_radius = mol_radius
        for frag, num_frag in zip(fragments, nums_fragments):
            frag_coords = cls.get_mol_coords(frag)
            frag_radius = max([math.sqrt(sum([x**2 for x in c]))
                               for c in frag_coords])
            if frag_radius > max_radius:
                max_radius = frag_radius
        return max_radius

    @classmethod
    def get_equivalent_radius(cls, mol_coords, fragments, nums_fragments):
        extra_atom_radius = 1.5
        total_volume = 0.0
        mol_radius = max([math.sqrt(sum([x**2 for x in c]))
                          for c in mol_coords])
        safe_radius = mol_radius + extra_atom_radius
        mol_volume = (4.0/3.0) * math.pi * (safe_radius ** 3)
        total_volume += mol_volume
        for frag, num_frag in zip(fragments, nums_fragments):
            frag_coords = cls.get_mol_coords(frag)
            frag_radius = max([math.sqrt(sum([x**2 for x in c]))
                               for c in frag_coords])
            safe_radius = frag_radius + extra_atom_radius
            frag_volume = (4.0/3.0) * math.pi * (safe_radius ** 3)
            total_volume += frag_volume * num_frag
        equivalent_radius = (total_volume * (3.0/4.0) / math.pi) ** (1.0/3.0)
        return equivalent_radius


    @classmethod
    def get_bounder(cls, mol_coords, fragments, nums_fragments, bound_setter):
        lower_bound = []
        upper_bound = []
        if bound_setter == "chain":
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
        elif bound_setter == "volume":
            margin = 3.0
            super_mol_radius = cls.get_equivalent_radius(mol_coords, fragments, nums_fragments) + margin
            frag_atom_counts = []
            for frag in fragments:
                num_atoms = frag.NumAtoms()
                frag_atom_counts.append(num_atoms)
            for num_frag, num_atoms in zip(nums_fragments, frag_atom_counts):
                single_frag_lower_bound = [-super_mol_radius] * 3
                single_frag_upper_bound = [super_mol_radius] * 3
                if num_atoms > 1:
                    single_frag_lower_bound.extend([0.0, -math.pi])
                    single_frag_upper_bound.extend([math.pi, math.pi])
                lower_bound.extend(single_frag_lower_bound * num_frag)
                upper_bound.extend(single_frag_upper_bound * num_frag)
        else:
            raise ValueError("only chain and volume bound setter is supported")
        bounder = inspyred.ec.Bounder(lower_bound, upper_bound)
        return bounder

    @staticmethod
    def normalize_molecule(mol):
        mol.Center()
        mol.ToInertialFrame()
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

    def generate_conformers(self, random, args):
        # generator
        lower_bound = copy.deepcopy(args['_ec'].bounder.lower_bound)
        upper_bound = copy.deepcopy(args['_ec'].bounder.upper_bound)
        if self.initial_guess == "center":
            fragments_atom_counts = [frag.NumAtoms() for frag in self.fragments]
            xi = 0
            for num_frag, num_atoms in zip(self.nums_fragments, fragments_atom_counts):
                for i in range(num_frag):
                    lower_bound[xi: xi+3] = [-self.max_radius] * 3
                    upper_bound[xi: xi+3] = [self.max_radius] *3
                    xi += 3
                    if num_atoms > 1:
                        xi += 2
        return [random.uniform(l, u) for l, u in zip(lower_bound, upper_bound)]

    def clean_swarm_memory(self):
        self.ea._previous_population = []
        self.ea.archive = []
        self.reevaluate_fitness = True


    def taboo_current_solution(self, coords_fitness):
        best_fitness, best_index = self._get_best_index_and_fitness(coords_fitness)
        self.energy_evaluator.calc_energy(fragments_coords=coords_fitness[best_index][0])
        self.energy_evaluator.taboo_current_position()

        # clean swarm memory
        self.clean_swarm_memory()
        self.conformer_staring_generation = None

    def _get_best_index_and_fitness(self, coords_fitness):
        best_index = 0
        best_fitness = coords_fitness[best_index][1]
        for i, (c, fitness) in enumerate(coords_fitness):
            if fitness < best_fitness:
                best_index = i
                best_fitness = fitness
        return best_fitness, best_index

    def is_conformer_located(self, coords_fitness):
        best_fitness, best_index = self._get_best_index_and_fitness(coords_fitness)
        mopac_energy_threshold = -3.0
        if best_fitness > mopac_energy_threshold:
            # Even not optimized by MOPAC
            return False
        best_coords = list(itertools.chain(*coords_fitness[best_index][0]))
        distances_to_best = [max([math.sqrt(sum([(x1-x2)**2
                                                 for x1, x2 in zip(c1, c2)]))
                                  for c1, c2 in zip(list(itertools.chain(*particle)), best_coords)])
                             for particle, fitness in coords_fitness]
        num_particle_in_range = 0
        for d in distances_to_best:
            if d <= self.taboo_tolerance_au:
                num_particle_in_range += 1
        ratio_particle_in_range = float(num_particle_in_range)/len(distances_to_best)
        if ratio_particle_in_range >= self.taboo_tolerance_particle_ratio:
            return True
        else:
            return False

    # noinspection PyUnusedLocal
    def evaluate_conformers(self, candidates, args):
        # evaluator
        fitness = []
        all_coords = []
        for c in candidates:
            fragments_coords = self.decode_solution(c)
            energy = self.energy_evaluator.calc_energy(fragments_coords)
            fitness.append(energy)
            all_coords.append(fragments_coords)
        coords_fitness = zip(all_coords, fitness)
        if self.is_conformer_located(coords_fitness):
            self.taboo_current_solution(coords_fitness)
        if self.always_write_best:
            best_fitness, best_index = self._get_best_index_and_fitness(coords_fitness)
            self.write_structure(candidates[best_index], "current_best.xyz")
        if self.reevaluate_fitness:
            self.reevaluate_fitness = False
            fitness = self.evaluate_conformers(candidates, args)
        mopac_energy_threshold = -3.0
        if min(fitness) < mopac_energy_threshold and self.conformer_staring_generation is None:
            self.conformer_staring_generation = self.ea.num_evaluations
        if self.conformer_staring_generation is not None and \
                        self.ea.num_evaluations > self.conformer_staring_generation + \
                        self.max_generations_each_conformer:
            coords_fitness = zip(all_coords, fitness)
            self.taboo_current_solution(coords_fitness)
        return fitness

    def write_structure(self, candidate, filename="result.xyz"):
        fragments_coords = self.decode_solution(candidate)
        mol_elements = self.get_mol_species(self.molecule)
        fragments_elements = []
        for frag, num_frag in zip(self.fragments, self.nums_fragments):
            fragments_elements.extend([self.get_mol_species(frag)] * num_frag)
        species = []
        coords_au = []
        species.extend(mol_elements)
        coords_au.extend(self.mol_coords)
        for elements, c in zip(fragments_elements,fragments_coords):
            species.extend(elements)
            coords_au.extend(c)
        coords_ang = [[x / AtomicRadiusUtils.angstrom2au for x in c] for c in coords_au]
        pmg_mol = Molecule(species, coords_ang)
        file_format = os.path.splitext(filename)[1][1:]
        pmg_mol.to(file_format, filename)

    def place(self, max_evaluations=30000, pop_size=100, neighborhood_size=5):
        t1 = time()
        self.final_pop = self.ea.evolve(generator=self.generate_conformers,
                                        evaluator=self.evaluate_conformers,

                                        pop_size=pop_size,
                                        bounder=self.bounder,
                                        maximize=False,
                                        max_evaluations=max_evaluations,
                                        neighborhood_size=neighborhood_size,

                                        inertia=0.721,
                                        cognitive_rate=1.193,
                                        social_rate=1.193)
        self.best = max(self.final_pop)
        # max means best, not necessarily smallest
        self.write_structure(self.best.candidate, filename="result.xyz")
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
    parser.add_argument("-t", "--taboo_tolerance", dest="taboo_tolerance", type=float,
                        default=1.0,
                        help="The radius to taboo a solution (in Angstrom)")
    parser.add_argument("-r", "--ratio_taboo_particles", dest="ratio_taboo_particles", type=float,
                        default=0.5,
                        help="ratio of particle within the tolerance to consider taboo current solution")
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
    parser.add_argument("--force_ordered_fragment", dest="force_ordered_fragment", action="store_true",
                        help="set this option to keep the fragment of the same in the order of input along the X-axis")
    parser.add_argument("--topology", dest="topology", choices=["ring", "star"], type=str, default="ring",
                        help="the topology of the PSO information network")
    parser.add_argument("--initial_guess", dest="initial_guess", choices=["breadth", "center", "volume"],
                        default="breadth",
                        help="where should particles should be initially put")
    parser.add_argument("--bound_setter", dest="bound_setter", choices=["chain", "volume"], default="chain",
                        help="method to set the bound conditions of PSO")
    parser.add_argument("--always_write_best", dest="always_write_best", action="store_true",
                        help="enable this option to output the best structure at every iteration")
    parser.add_argument("--random_seed", dest="random_seed", default=None, type=int,
                        help="random seed for PSO, an integer is expected")
    parser.add_argument("--max_generations_each_conformer", dest="max_generations_each_conformer", default=100, type=int,
                        help="maximum generations for each conformer")
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
        energy_evaluator = SemiEmpricalQuatumMechanicalEnergyEvaluator(
            molecule, fragments, options.nums_fragments, total_charge=options.charge,
            taboo_tolerance_ang=options.taboo_tolerance, force_order_fragment=options.force_ordered_fragment,
            bound_setter=options.bound_setter)
    if len(fragments) != len(options.nums_fragments):
        raise ValueError("you must specify the duplicated count for every fragment")
    placer = IonPlacer(molecule=molecule, fragments=fragments, nums_fragments=options.nums_fragments,
                       energy_evaluator=energy_evaluator, taboo_tolerance_ang=options.taboo_tolerance,
                       taboo_tolerance_particle_ratio=options.ratio_taboo_particles, topology=options.topology,
                       initial_guess=options.initial_guess, bound_setter=options.bound_setter,
                       always_write_best=options.always_write_best, random_seed=options.random_seed,
                       max_generations_each_conformer=options.max_generations_each_conformer)
    energy_evaluator.arranger = placer
    placer.place(max_evaluations=options.iterations,
                 pop_size=options.size,
                 neighborhood_size=options.num_neighbours)
    print 'It took {:.1f} seconds to place the salt'.format(placer
                                                            .playing_time)


if __name__ == '__main__':
    main()