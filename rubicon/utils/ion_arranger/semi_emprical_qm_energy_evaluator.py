import copy
import os
import shutil

from custodian import Custodian
import itertools
import math
from monty.tempfile import ScratchDir
from pymatgen.core.structure import Molecule
from pymatgen.core.units import Energy
from pymatgen.io.babelio import BabelMolAdaptor

from rubicon.io.mopacio.custodian.mopac_error_handlers import MopacErrorHandler
from rubicon.io.mopacio.custodian.mopacjob import MopacJob
from rubicon.io.mopacio.mopacio import MopOutput, MopTask
from rubicon.utils.ion_arranger.energy_evaluator import EnergyEvaluator
from rubicon.utils.ion_arranger.hard_sphere_energy_evaluators import HardSphereEnergyEvaluator, AtomicRadiusUtils, \
    ContactDetector, LargestContactGapEnergyEvaluator, OrderredLayoutEnergyEvaluator, UmbrellarForceEnergyEvaluator, \
    ContactGapRMSDEnergyEvaluator


__author__ = 'xiaohuiqu'


class SemiEmpricalQuatumMechanicalEnergyEvaluator(EnergyEvaluator):

    def __init__(self, ob_mol, ob_fragments, nums_fragments, total_charge,
                 lower_covalent_radius_scale=2.0, lower_metal_radius_scale=0.8,
                 upper_covalent_radius_scale=3.0, upper_metal_radius_scale=1.5,
                 taboo_tolerance_ang=1.0, force_order_fragment=False, bound_setter="chain"):
        from rubicon.utils.ion_arranger.ion_arranger import IonPlacer
        self.mol_coords = IonPlacer.normalize_molecule(ob_mol)
        super(SemiEmpricalQuatumMechanicalEnergyEvaluator, self).__init__(self.mol_coords)
        self.total_charge = total_charge
        self.force_ordered_fragment = force_order_fragment
        self.lower_sphere = self._construct_hardsphere_energy_evaluator(
            lower_covalent_radius_scale, lower_metal_radius_scale,
            self.mol_coords, ob_mol, ob_fragments, nums_fragments)
        self.contact_detector = self._construct_contact_detector(
            upper_covalent_radius_scale, upper_metal_radius_scale,
            self.mol_coords, ob_mol, ob_fragments, nums_fragments)
        self.layout_order = OrderredLayoutEnergyEvaluator(self.mol_coords, nums_fragments)
        grav = "rmsd_gap"
        if grav == "rmsd_gap":
            self.gravitation = self._construct_contact_cap_rmsd_energy_evaluator(
                upper_covalent_radius_scale, upper_metal_radius_scale,
                self.mol_coords, ob_mol, ob_fragments, nums_fragments)
        elif grav == "largest_gap":
            self.gravitation = self._construct_largest_cap_energy_evaluator(
                upper_covalent_radius_scale, upper_metal_radius_scale,
                self.mol_coords, ob_mol, ob_fragments, nums_fragments, bound_setter=bound_setter)
        self.mol_species = IonPlacer.get_mol_species(ob_mol)
        self.fragments_species = [IonPlacer.get_mol_species(frag) for frag in ob_fragments]
        self.ob_fragments = ob_fragments
        self.nums_fragments = nums_fragments
        self.run_number = 1
        self.global_best_energy = 0.0
        self.memory_position_tolerance_au = 0.3 * AtomicRadiusUtils.angstrom2au
        self.memory_positions = []
        self.memory_size = 1000
        self.current_raw_position = None
        self.current_optimized_position = None
        self.taboo_tolerance_au = taboo_tolerance_ang * AtomicRadiusUtils.angstrom2au
        self.tabooed_raw_positions = []
        self.tabooed_optimizated_positions = []
        self.current_best_raw_position = None
        self.current_best_optimized_position = None
        self.best_umbrella_ratio = 3.0
        self.best_energy = None
        self.best_mol = None
        self.best_run_number = None
        self.sqm_umbrella_centers = []
        self.sqm_umbrella = UmbrellarForceEnergyEvaluator(self.mol_coords, self.sqm_umbrella_centers, self.taboo_tolerance_au,
                                                          self.best_umbrella_ratio)
        self.gap_umbrella_centers = []
        self.gap_umbrella = None
        self.smallest_gap = 9999999999.9
        self.arranger = None


    def save_current_best_position(self):
        conformer_dir = "conformers"
        if not os.path.exists(conformer_dir):
            os.mkdir(conformer_dir)
        mol_file_name = os.path.join(conformer_dir,
                                     "step_{}_energy_{}.xyz".format(self.best_run_number, self.best_energy))
        BabelMolAdaptor(self.best_mol).write_file(mol_file_name, file_format="xyz")

    def taboo_current_position(self, raw_position_only=False):
        self.tabooed_raw_positions.append(list(itertools.chain(*self.current_raw_position)))
        if not raw_position_only:
            if self.current_optimized_position is None:
                mol = self._get_super_molecule(self.current_raw_position)
                energy, final_mol = self.run_mopac(mol)
                from rubicon.utils.ion_arranger.ion_arranger import IonPlacer
                final_coords = IonPlacer.normalize_molecule(final_mol)
                self.current_optimized_position = final_coords
            self.tabooed_optimizated_positions.append(self.current_optimized_position)

            # taboo is called from ion arranger
            self.save_current_best_position()
            self.memory_positions = []
            self.best_energy = None
            self.current_best_optimized_position = None
            self.current_best_raw_position = None
            self.best_mol = None
            self.best_run_number = None
            del self.sqm_umbrella_centers[:]

    def is_current_position_tabooed(self, position_type):
        if position_type not in ["raw", "optimized"]:
            raise ValueError("Position type must be either raw or optimized")
        current_pos = self.current_raw_position if position_type == "raw" \
            else self.current_optimized_position
        if position_type == "raw":
            current_pos = list(itertools.chain(*current_pos))
        tabooed_positions = self.tabooed_raw_positions if position_type == "raw" \
            else self.tabooed_optimizated_positions
        for p in tabooed_positions:
            distance = max([math.sqrt(sum([(x1-x2)**2 for x1, x2 in zip(c1, c2)]))
                            for c1, c2 in
                            zip(p, current_pos)])
            if distance < self.taboo_tolerance_au:
                return True
        return False

    def query_memory_positions(self, fragments_coords):
        current_pos = list(itertools.chain(*fragments_coords))
        for p, energy in self.memory_positions:
            distance = max([math.sqrt(sum([(x1-x2)**2 for x1, x2 in zip(c1, c2)]))
                            for c1, c2 in
                            zip(p, current_pos)])
            if distance < self.memory_position_tolerance_au:
                return energy
        return None

    def append_position_to_memory(self, fragments_coords, energy):
        self.memory_positions.insert(0, tuple([list(itertools.chain(*fragments_coords)), energy]))
        if len(self.memory_positions) > self.memory_size:
            self.memory_positions.pop()

    def is_raw_position_inside_best_umbrella(self, fragments_coords):
        # replaced by Umbrella force
        if self.current_best_raw_position is None:
            return False
        current_pos = list(itertools.chain(*fragments_coords))
        distance = max([math.sqrt(sum([(x1-x2)**2 for x1, x2 in zip(c1, c2)]))
                        for c1, c2 in
                        zip(self.current_best_raw_position, current_pos)])
        if distance < self.taboo_tolerance_au * self.best_umbrella_ratio:
            return True
        else:
            return False

    def is_optimized_position_inside_the_best_promixity(self, fragments_coords):
        if self.current_best_optimized_position is None:
            return True
        distance = max([math.sqrt(sum([(x1-x2)**2 for x1, x2 in zip(c1, c2)]))
                        for c1, c2 in
                        zip(self.current_best_optimized_position, fragments_coords)])
        if distance < self.taboo_tolerance_au:
            return True
        else:
            return False

    def build_or_extend_gap_umbrella(self, current_gap, fragments_coords):
        smallest_gap_umbrella_radius = self.taboo_tolerance_au * self.best_umbrella_ratio
        umbrella_radius = min([smallest_gap_umbrella_radius, current_gap])
        if self.gap_umbrella is None:
            self.gap_umbrella = UmbrellarForceEnergyEvaluator(self.mol_coords, self.gap_umbrella_centers,
                                                              self.taboo_tolerance_au,
                                                              self.best_umbrella_ratio)
            self.arranger.clean_swarm_memory()
            self.smallest_gap = 9999999999.9
        if current_gap < self.smallest_gap:
            umbrella_coords = list(itertools.chain(*fragments_coords))
            self.sqm_umbrella_centers.append(umbrella_coords)
            self.smallest_gap = current_gap
            self.gap_umbrella.umbrella_radius = umbrella_radius


    def calc_energy(self, fragments_coords):
        # Energy orders:
        # Tabooed(1.0E5) > Umbrella(7.0E4) > Layout (6.0E4) > HardSphere (5.0E4)  >
        # gravity (4.0E4) > PM7 (-1.0E1)

        tabooed_energy = 1.0E5

        self.current_raw_position = fragments_coords
        self.current_optimized_position = None

        if self.is_current_position_tabooed(position_type="raw"):
            return tabooed_energy

        if self.gap_umbrella is not None:
            gap_umbrella_energy = self.gap_umbrella.calc_energy(fragments_coords)
            if gap_umbrella_energy > self.gap_umbrella.base_energy - 100.0:
                return gap_umbrella_energy

        sqm_umbrella_energy = self.sqm_umbrella.calc_energy(fragments_coords)
        if sqm_umbrella_energy > self.sqm_umbrella.base_energy - 100.0:
            return sqm_umbrella_energy

        if self.force_ordered_fragment:
            energy = self.layout_order.calc_energy(fragments_coords)
            if energy > self.layout_order.base_energy - 100.0:
                return energy

        energy = self.lower_sphere.calc_energy(fragments_coords)
        if energy > self.lower_sphere.base_energy - 100.0:
            return energy

        if not self.contact_detector.is_contact(fragments_coords):
            energy = self.gravitation.calc_energy(fragments_coords)
            if energy > self.gravitation.base_energy - 100.0:
                self.build_or_extend_gap_umbrella(current_gap=energy-self.gravitation.base_energy,
                                                  fragments_coords=fragments_coords)
                return energy

        memorized_energy = self.query_memory_positions(fragments_coords)
        if memorized_energy is not None:
            energy = memorized_energy
        else:
            mol = self._get_super_molecule(fragments_coords)
            energy, final_mol = self.run_mopac(mol)
            from rubicon.utils.ion_arranger.ion_arranger import IonPlacer
            final_coords = IonPlacer.normalize_molecule(final_mol)
            self.current_optimized_position = final_coords
            if self.is_current_position_tabooed(position_type="optimized"):
                self.taboo_current_position(raw_position_only=True)
                return tabooed_energy
            energy = round(energy, 3)
            if self.is_optimized_position_inside_the_best_promixity(final_coords):
                if energy < self.best_energy or self.current_best_optimized_position is None:
                    self.best_energy = energy
                    self.current_best_optimized_position = final_coords
                    self.current_best_raw_position = list(itertools.chain(*fragments_coords))
                    if len(self.sqm_umbrella_centers) == 0 and self.arranger is not None:
                        self.arranger.clean_swarm_memory()
                        self.gap_umbrella = None
                        self.gap_umbrella_centers = []
                    self.sqm_umbrella_centers.append(self.current_best_raw_position)
                    self.best_mol = final_mol
                    self.best_run_number = self.run_number
            else:
                energy = tabooed_energy
            self.append_position_to_memory(fragments_coords, energy)
        # coarse grained energy,
        # make potential energy surface simpler
        return energy

    def run_mopac(self, mol):
        cur_dir = os.getcwd()
        all_errors = set()
        energy = 0.0
        with ScratchDir(rootpath=cur_dir, copy_from_current_on_enter=False, copy_to_current_on_exit=True):
            order_text = ["st", "nd", "th"]
            title = "Salt Alignment {}{} Calculation".format(
                self.run_number, order_text[self.run_number-1 if self.run_number < 3 else 2])
            self.run_number += 1
            mop = MopTask(mol, self.total_charge, "opt", title, "PM7", {"CYCLES": 2000})
            mop.write_file("mol.mop")
            job = MopacJob()
            handler = MopacErrorHandler()
            c = Custodian(handlers=[handler], jobs=[job], max_errors=50)
            custodian_out = c.run()
            mopout = MopOutput("mol.out")
            final_pmg_mol = mopout.data["molecules"][-1]
            final_ob_mol = BabelMolAdaptor(final_pmg_mol)._obmol
            energy = Energy(mopout.data["energies"][-1][-1], "eV").to("Ha")
            if energy < self.global_best_energy:
                self.global_best_energy = energy
                shutil.copy("mol.out", os.path.join(cur_dir, "best_mol.out"))
            for run in custodian_out:
                for correction in run['corrections']:
                    all_errors.update(correction['errors'])
            if len(all_errors) > 0:
                out_error_logs_file = os.path.join(cur_dir, "out_error_logs.txt")
                bak_cmd = "cat {} >> {}".format("mol.out", out_error_logs_file)
                os.system(bak_cmd)
        return energy, final_ob_mol

    def _get_super_molecule(self, fragments_coords):
        super_mol_species = []
        super_mol_coords_au = []
        super_mol_species.extend(copy.deepcopy(self.mol_species))
        super_mol_coords_au.extend(copy.deepcopy(self.mol_coords))
        duplicated_fragments_species = []
        for frag_sp, num_frag in zip(self.fragments_species, self.nums_fragments):
            duplicated_fragments_species.extend([frag_sp]*num_frag)
        for coords, sp in zip(fragments_coords, duplicated_fragments_species):
            super_mol_coords_au.extend(coords)
            super_mol_species.extend(sp)
        super_mol_coords_ang = [[x/AtomicRadiusUtils.angstrom2au for x in coord] for coord in super_mol_coords_au]
        return Molecule(super_mol_species, super_mol_coords_ang)

    @staticmethod
    def _construct_hardsphere_energy_evaluator(covalent_radius_scale, metal_radius_scale,
                                               mol_coords, ob_mol, ob_fragments, nums_fragments):
        rad_util = AtomicRadiusUtils(covalent_radius_scale, metal_radius_scale)
        mol_radius = rad_util.get_radius(ob_mol)
        fragments_atom_radius = [rad_util.get_radius(frag) for frag in ob_fragments]
        sphere = HardSphereEnergyEvaluator(
            mol_coords, mol_radius, fragments_atom_radius, nums_fragments)
        return sphere

    @staticmethod
    def _construct_contact_detector(covalent_radius_scale, metal_radius_scale,
                                    mol_coords, ob_mol, ob_fragments, nums_fragments):
        rad_util = AtomicRadiusUtils(covalent_radius_scale, metal_radius_scale)
        mol_radius = rad_util.get_radius(ob_mol)
        fragments_atom_radius = [rad_util.get_radius(frag) for frag in ob_fragments]
        detector = ContactDetector(mol_coords, mol_radius, fragments_atom_radius, nums_fragments)
        return detector

    @staticmethod
    def _construct_largest_cap_energy_evaluator(covalent_radius_scale, metal_radius_scale,
                                                mol_coords, ob_mol, ob_fragments, nums_fragments,
                                                bound_setter="chain"):
        rad_util = AtomicRadiusUtils(covalent_radius_scale, metal_radius_scale)
        mol_radius = rad_util.get_radius(ob_mol)
        fragments_atom_radius = [rad_util.get_radius(frag) for frag in ob_fragments]
        from rubicon.utils.ion_arranger.ion_arranger import IonPlacer
        bounder = IonPlacer.get_bounder(mol_coords, ob_fragments, nums_fragments, bound_setter)
        max_cap = max(bounder.upper_bound) * 2.0 / AtomicRadiusUtils.angstrom2au
        evaluator = LargestContactGapEnergyEvaluator(
            mol_coords, mol_radius, fragments_atom_radius, nums_fragments, max_cap, threshold=0.01)
        return evaluator

    @staticmethod
    def _construct_contact_cap_rmsd_energy_evaluator(covalent_radius_scale, metal_radius_scale,
                                                mol_coords, ob_mol, ob_fragments, nums_fragments):
        rad_util = AtomicRadiusUtils(covalent_radius_scale, metal_radius_scale)
        mol_radius = rad_util.get_radius(ob_mol)
        fragments_atom_radius = [rad_util.get_radius(frag) for frag in ob_fragments]
        evaluator = ContactGapRMSDEnergyEvaluator(mol_coords, mol_radius, fragments_atom_radius, nums_fragments)
        return evaluator