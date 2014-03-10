from random import Random
from time import time
import inspyred
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
import openbabel as ob


class HardSphereIonPlacer():
    def __init__(self, molecule, cation, anion, prng=None, radius_scale=2.0):
        lower_bound = [-100.0, -50.0]
        upper_bound = [100.0, 50.0]
        self.bounder = inspyred.ec.Bounder(lower_bound, upper_bound)
        self.prng = prng if prng else Random()
        self.radius_scale = radius_scale
        self.seed = time()
        self.prng.seed(self.seed)
        self.final_pop = None
        self.best = None
        self.ea = inspyred.swarm.PSO(self.prng)
        self.ea.terminator = inspyred.ec.terminators.evaluation_termination
        self.ea.topology = inspyred.swarm.topologies.ring_topology
        self.molecule = molecule
        self.cation = cation
        self.anion = anion


    @staticmethod
    def normalize_molecule(mol, radius_scale):
        mol.Center()
        coords = []
        radius = []
        ref_radius = CovalentRadius.radius
        num_atoms = mol.NumAtoms()
        element_table = ob.OBElementTable()
        for i in range(1, num_atoms + 1):
            a = mol.GetAtom(i)
            atomic_num = a.GetAtomicNum()
            symbol = element_table.GetSymbol(atomic_num)
            rad = ref_radius[symbol] * radius_scale
            radius.append(rad)
            coords.append([a.GetX(), a.GetY(), a.GetZ()])
        return coords, radius



    @staticmethod
    def generate_conformers(random, args):
        # generator
        lower_bound = args['_ec'].bounder.lower_bound
        upper_bound = args['_ec'].bounder.upper_bound
        return [random.uniform(l, u) for l, u in zip(lower_bound, upper_bound)]

    @staticmethod
    def evaluate_conformers(candidates, args):
        # evaluator
        fitness = []
        for c in candidates:
            fit = (c[0] - 1.5)**2 + (c[1] + 35)**2
            fitness.append(fit)
        return fitness

    def place(self):
        self.final_pop = self.ea.evolve(generator=self.generate_conformers,
                                        evaluator=self.evaluate_conformers,
                                        pop_size=100,
                                        bounder=self.bounder,
                                        maximize=False,
                                        max_evaluations=30000,
                                        neighborhood_size=5)
        self.best = max(self.final_pop)
        # max means best, not necessarily smallest


if __name__ == '__main__':
    placer = HardSphereIonPlacer()
    placer.place()
    print('Best Solution: \n{0}'.format(str(placer.best)))