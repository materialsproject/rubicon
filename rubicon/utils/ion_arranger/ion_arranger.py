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





    def __init__(self, molecule, cation, anion, num_cation, num_anion, energy_evaluator,
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
        self.cation = cation
        self.anion = anion
        self.num_cation = num_cation
        self.num_anion = num_anion
        self.mol_coords = self.normalize_molecule(self.molecule)
        self.normalize_molecule(self.cation)
        self.normalize_molecule(self.anion)
        self.bounder = self.get_bounder(self.mol_coords, cation, anion, num_cation, num_anion)
        self.best_pymatgen_mol = None
        self.playing_time = None
        self.energy_evaluator = energy_evaluator

    @classmethod
    def get_bounder(cls, mol_coords, ob_cation, ob_anion, num_cation, num_anion):
        lower_bound = []
        upper_bound = []
        cation_num_atoms = ob_cation.NumAtoms()
        anion_num_atoms = ob_anion.NumAtoms()
        mol_radius = max([math.sqrt(sum([x**2 for x in c]))
                          for c in mol_coords])
        cation_coords = cls.get_mol_coords(ob_cation)
        cation_radius = max([math.sqrt(sum([x**2 for x in c]))
                             for c in cation_coords])
        anion_coords = cls.get_mol_coords(ob_anion)
        anion_radius = max([math.sqrt(sum([x**2 for x in c]))
                            for c in anion_coords])
        max_length = 2 * (mol_radius + cation_radius * num_cation +
                          anion_radius * num_anion)
        x_min = max_length + 5.0 * (num_cation + num_anion)
        for i in range(num_cation):
            lower_bound.extend([-x_min] * 3)
            upper_bound.extend([x_min] * 3)
            if cation_num_atoms > 1:
                lower_bound.extend([0.0, -math.pi])
                upper_bound.extend([math.pi, math.pi])
        for i in range(num_anion):
            lower_bound.extend([-x_min] * 3)
            upper_bound.extend([x_min] * 3)
            if anion_num_atoms > 1:
                lower_bound.extend([0.0, -math.pi])
                upper_bound.extend([math.pi, math.pi])
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
        cation_coords = []
        anion_coords = []
        cation_num_atoms = self.cation.NumAtoms()
        anion_num_atoms = self.anion.NumAtoms()
        xi = 0
        for i in range(self.num_cation):
            cc = ob.OBMol(self.cation)
            tx, ty, tz = [t/AtomicRadiusUtils.angstrom2au for  t in x[xi: xi+3]]
            xi += 3
            if cation_num_atoms > 1:
                theta = x[xi]
                xi += 1
                phi = x[xi]
                xi += 1
                self.rotate(cc, theta, phi)
            cc.Translate(ob.vector3(tx, ty, tz))
            cation_coords.append(self.get_mol_coords(cc))
        for i in range(self.num_anion):
            ac = ob.OBMol(self.anion)
            tx, ty, tz = [t/AtomicRadiusUtils.angstrom2au for t in x[xi: xi+3]]
            xi += 3
            if anion_num_atoms > 1:
                theta = x[xi]
                xi += 1
                phi = x[xi]
                xi += 1
                IonPlacer.rotate(ac, theta, phi)
            ac.Translate(ob.vector3(tx, ty, tz))
            anion_coords.append(IonPlacer.get_mol_coords(ac))
        return cation_coords, anion_coords

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
            cation_coords, anion_coords = self.decode_solution(c)
            energy = self.energy_evaluator.calc_energy(cation_coords, anion_coords)
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
        best_cation_coords, best_anion_coords = self.decode_solution(
            self.best.candidate)
        mol_elements = self.get_mol_species(self.molecule)
        cation_elements = self.get_mol_species(self.cation)
        anion_elements = self.get_mol_species(self.anion)
        species = []
        coords_au = []
        species.extend(mol_elements)
        coords_au.extend(self.mol_coords)
        for c in best_cation_coords:
            species.extend(cation_elements)
            coords_au.extend(c)
        for c in best_anion_coords:
            species.extend(anion_elements)
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
                        help="the QChem output file name of molecule")
    parser.add_argument("-c", "--cation", dest="cation", type=str,
                        required=True,
                        help="the QChem output file name of cation")
    parser.add_argument("-a", "--anion", dest="anion", type=str,
                        required=True,
                        help="the QChem output file name of anion")
    parser.add_argument("-o", "--outputfile", dest="outputfile", type=str,
                        required=True,
                        help="the file name of the aligned conformer")
    parser.add_argument("-n", "--num_iter", dest="num_iter", type=int,
                        default=600,
                        help="maximum number of evaluations")
    parser.add_argument("-s", "--size", dest="size", type=int, default=15,
                        help="population size")
    parser.add_argument("-k", "--num_neighbours", dest="num_neighbours", type=int, default=2,
                        help="number of neighbours")
    parser.add_argument("-e", "--evaluator", dest="evaluator", type=str, default="hardsphere",
                        choices=["hardsphere", "sqm"], help="Energy Evaluator")
    options = parser.parse_args()
    qcout_molecule = QcOutput(options.molecule)
    qcout_cation = QcOutput(options.cation)
    qcout_anion = QcOutput(options.anion)
    total_charge_cation = qcout_cation.data[0]["molecules"][-1].charge
    total_charge_anion = qcout_anion.data[0]["molecules"][-1].charge
    num_lcm = lcm(total_charge_cation, -total_charge_anion)
    num_cation = num_lcm/total_charge_cation
    num_anion = num_lcm/-total_charge_anion
    pymatgen_mol_molecule = qcout_molecule.data[0]["molecules"][-1]
    pymatgen_mol_cation = qcout_cation.data[0]["molecules"][-1]
    pymatgen_mol_anion = qcout_anion.data[0]["molecules"][-1]
    # noinspection PyProtectedMember
    obmol_molecule = BabelMolAdaptor(pymatgen_mol_molecule)._obmol
    # noinspection PyProtectedMember
    obmol_cation = BabelMolAdaptor(pymatgen_mol_cation)._obmol
    # noinspection PyProtectedMember
    obmol_anion = BabelMolAdaptor(pymatgen_mol_anion)._obmol
    hardsphere_evaluator = lambda: HardSphereElectrostaticEnergyEvaluator.from_qchem_output(
        qcout_molecule, qcout_cation, qcout_anion)
    sqm_evaluator = lambda: SemiEmpricalQuatumMechanicalEnergyEvaluator(
        obmol_molecule, obmol_cation, obmol_anion, total_charge=0, num_cation=num_cation, num_anion=num_anion)
    evaluator_creators = {"hardsphere": hardsphere_evaluator, "sqm": sqm_evaluator}
    creator = evaluator_creators[options.evaluator]
    energy_evaluator = creator()
    placer = IonPlacer(
        obmol_molecule, obmol_cation, obmol_anion, num_cation, num_anion, energy_evaluator)
    placer.place(max_evaluations=options.num_iter,
                 pop_size=options.size,
                 neighborhood_size=options.num_neighbours)
    print 'It took {:.1f} seconds to place the salt'.format(placer
                                                            .playing_time)
    write_mol(placer.best_pymatgen_mol, options.outputfile)

if __name__ == '__main__':
    main()