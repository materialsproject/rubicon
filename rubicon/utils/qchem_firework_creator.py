# coding=utf-8
import copy
import json
import os
import math

from fireworks import Firework
from fireworks.core.firework import Tracker
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.qchemio import QcTask, QcInput

from rubicon.dupefinders.dupefinder_eg import DupeFinderEG
from rubicon.firetasks.qchem_task import QChemTask


__author__ = 'xiaohuiqu'


class QChemFireWorkCreator():
    def __init__(self, mol, molname, mission, additional_user_tags=None,
                 dupefinder=None, priority=1, update_spec=None, large=False):
        self.molname = molname
        self.mol = mol
        self.large = large
        initial_inchi = self.get_inchi(mol)
        user_tags = {'mission': mission,
                     "molname": molname}
        if additional_user_tags:
            user_tags.update(additional_user_tags)
        spec = dict()
        spec['user_tags'] = user_tags
        spec['_priority'] = priority
        spec['_dupefinder'] = dupefinder.to_dict() if dupefinder \
            else DupeFinderEG().to_dict()
        tracker_out = Tracker("mol.qcout", nlines=20)
        tracker_std = Tracker("mol.qclog", nlines=10)
        tracker_joberr = Tracker("FW_job.error", nlines=20)
        tracker_jobout = Tracker("FW_job.out", nlines=20)
        spec["_trackers"] = [tracker_out, tracker_std, tracker_joberr,
                             tracker_jobout]
        spec['run_tags'] = dict()
        spec['implicit_solvent'] = {}
        spec['inchi'] = initial_inchi
        spec['num_atoms'] = len(mol)
        if update_spec:
            spec.update(update_spec)
        self.base_spec = lambda: copy.deepcopy(spec)

    @staticmethod
    def get_inchi(mol):
        bb = BabelMolAdaptor(mol)
        pbmol = bb.pybel_mol
        return pbmol.write("inchi").strip()

    @staticmethod
    def get_state_name(charge, spin_multiplicity):
        charge_state = {-2: "anion_2", -1: "anion", 0: "neutral", 1: "cation", 2: "cation_2"}
        spin_state = {1: "singlet", 2: "doublet", 3: "triplet"}
        return spin_state[spin_multiplicity] + " " + charge_state[charge]

    @staticmethod
    def get_exchange_correlation_basis_auxbasis_remparams(method):
        if method is None:
            method = "B3LYP/6-31+G*"
        aux_basis = None
        correlation = None
        rem_params = None
        theoretical_level, basis_set = method.split('/')
        basis_set = basis_set.lower()
        if theoretical_level.lower() == "b3lyp-xdm":
            exchange = 'b3lyp'
            rem_params = {"dftvdw_jobnumber": 1,
                          "dftvdw_method": 1,
                          "dftvdw_print": 1,
                          "dftvdw_kai": 800,
                          "dftvdw_use_ele_drv": 1}
        elif theoretical_level.lower() in ["xygjos", "lxygjos"]:
            exchange = theoretical_level.lower()
            if basis_set == "6-31+g*":
                aux_basis = "rimp2-aug-cc-pvdz"
            else:
                aux_basis = "rimp2-aug-cc-pvtz"
            if exchange == "lxygjos":
                rem_params = {"omega": 200}
        elif theoretical_level.lower() == "pbe-d3":
            exchange = 'pbe'
            correlation = 'pbe'
            rem_params = {"dft_d": "empirical_grimme3",
                          "dft_d3_s6": 1000,
                          "dft_d3_rs6": 1217,
                          "dft_d3_s8": 722,
                          "dft_d3_3body": True}
        elif theoretical_level.lower() == "blyp-d3":
            exchange = 'b'
            correlation = 'lyp'
            rem_params = {"dft_d": "empirical_grimme3",
                          "dft_d3_s6": 1000,
                          "dft_d3_rs6": 1094,
                          "dft_d3_s8": 1682,
                          "dft_d3_3body": True}
        elif theoretical_level.lower() == "b3lyp-d3":
            exchange = 'b3lyp'
            rem_params = {"dft_d": "empirical_grimme3",
                          "dft_d3_s6": 1000,
                          "dft_d3_rs6": 1261,
                          "dft_d3_s8": 1703,
                          "dft_d3_3body": True}
        else:
            exchange = theoretical_level.lower()
        method_token = [t for t in [basis_set, exchange, aux_basis, correlation, rem_params]
                        if t]
        return exchange, correlation, basis_set,  aux_basis, rem_params, method_token


    def geom_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db,
                priority=None, method=None, task_type_prefix=None):
        task_type = "geometry optimization"
        if task_type_prefix:
            task_type = task_type_prefix + " " + task_type
        state_name = self.get_state_name(charge, spin_multiplicity)
        if not method:
            if self.large:
                method = "PBE-D3/6-31+G*"
            else:
                method = "B3lYP/6-31+G*"
        title = self.molname + " " + state_name + " " + method + " " + task_type
        exchange, correlation, basis_set,  aux_basis, rem_params, method_token = self.\
            get_exchange_correlation_basis_auxbasis_remparams(method)
        qctask = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                        jobtype="opt", title=title, exchange=exchange, correlation=correlation,
                        basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        if self.large:
            qctask.set_geom_max_iterations(200)
            qctask.set_scf_algorithm_and_iterations(iterations=100)
            qctask.scale_geom_opt_threshold(gradient=1.0,
                                            displacement=10.0,
                                            energy=10.0)
        qcinp = QcInput([qctask])
        spec = self.base_spec()
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        spec['run_tags']['methods'] = method_token
        spec["qm_method"] = method
        if priority:
            spec['_priority'] = priority
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemGeomOptDBInsertionTask
        fw_geom_cal = Firework([QChemTask()],
                               spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_geom_db = Firework([QChemGeomOptDBInsertionTask()],
                              spec=spec_db, name=task_name_db, fw_id=fw_id_db)

        return fw_geom_cal, fw_geom_db

    def freq_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db,
                priority=None, method=None):
        if not method:
            if self.large:
                method = "PBE-D3/6-31+G*"
            else:
                method = "B3lYP/6-31+G*"
        task_type = "vibrational frequency"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + method + " " + task_type
        exchange, correlation, basis_set,  aux_basis, rem_params, method_token = self. \
            get_exchange_correlation_basis_auxbasis_remparams(method)
        if exchange.lower() in ["xygjos"]:
            rem_params["IDERIV"] = 1
        qctask = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                        jobtype="freq", title=title, exchange=exchange, correlation=correlation,
                        basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        if self.large:
            qctask.set_scf_algorithm_and_iterations(iterations=100)
        qcinp = QcInput([qctask])
        spec = self.base_spec()
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        spec['run_tags']['methods'] = method_token
        spec["qm_method"] = method
        if priority:
            spec['_priority'] = priority
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemFrequencyDBInsertionTask
        fw_freq_cal = Firework([QChemTask()],
                               spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_freq_db = Firework([QChemFrequencyDBInsertionTask()],
                              spec=spec_db, name=task_name_db, fw_id=fw_id_db)
        return fw_freq_cal, fw_freq_db

    @staticmethod
    def get_dielectric_constant(solvent, use_vdW_surface):
        dielec_const_file = os.path.join(os.path.dirname(__file__),
                                         "../utils/data", "pcm_params.json")
        with open(dielec_const_file) as f:
            dielec_consts = json.load(f)
        if isinstance(solvent, dict):
            # e. g. {"components": {"ethylenecarbonate": 3.0,
            #                       "ethylmethylcarbonate": 7.0}
            #        "metrics": "molarity"}
            # supported metrics: molarity and volume
            # probe_radius mix by volume
            # dielectric constants mix according to
            # Peiming Wang, Andrzej Anderko; Fluid Phase Equilibria 186 (2001) 103–122
            components = sorted([(comp.lower(), float(amount))
                                 for comp, amount in solvent["components"].items()])
            compounds = [comp for comp, amount in components]
            for comp in compounds:
                if comp not in dielec_consts:
                    raise Exception("Don't know the dielectric constants for "
                                    "solvent '{}' in mixtures".format(comp))
            rads = [dielec_consts[comp]["solvent_probe_radius"] for comp in compounds]
            dielecs = [dielec_consts[comp]["dielectric"] for comp in compounds]
            if solvent["metrics"] == "volume":
                total_volume = sum([amount for comp, amount in components])
                vol_ratio = [amount/total_volume for comp, amount in components]
            elif solvent["metrics"] == "molarity":
                volumes = [amount * (r**3) for (comp, amount), r in zip(components, rads)]
                total_volume = sum(volumes)
                vol_ratio = [v/total_volume for v in volumes]
            else:
                raise Exception("Only molartiy and volume ratio is supported")
            moles = [v / (r**3) for v, r in zip(vol_ratio, rads)]
            total_mole = sum(moles)
            volume_per_molecule = 1.0/total_mole
            probe_radius = volume_per_molecule ** (1.0/3)
            polarizations = [((epsilon - 1) * (2 * epsilon + 1))/(9 * epsilon) for epsilon in dielecs]
            total_polarization = sum([p * v for p, v in zip(polarizations, vol_ratio)])
            pol_per_vol = total_polarization / 1.0
            a, b, c = 2, - (1 + 9 * pol_per_vol), -1
            dielectric_constant = (-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a)
            solvent_name = ', '.join(['{:.0%} {:s}'.format(vol, name)
                                      for vol, name in zip(vol_ratio, compounds)]) \
                           + ' in volume'
        elif isinstance(solvent, str) or isinstance(solvent, unicode):
            solvent_pure = solvent.lower()
            if solvent_pure not in dielec_consts:
                raise Exception("Don't know the dielectric constants for "
                                "solvent '{}'".format(solvent_pure))
            dielectric_constant = dielec_consts[solvent_pure]["dielectric"]
            probe_radius = dielec_consts[solvent_pure]["solvent_probe_radius"]
            solvent_name = solvent
        else:
            raise Exception("Please use a string for pure solvent, and use"
                            "dict for mixtures")
        if use_vdW_surface:
            probe_radius = 0.0
        return dielectric_constant, probe_radius, solvent_name

    @staticmethod
    def get_smx_solvent(implicit_solvent, solvent, solvent_method):
        smx_data_file = os.path.join(os.path.dirname(__file__),
                                     "../utils/data", "smx_data.json")
        with open(smx_data_file) as f:
            smx_data = json.load(f)
        if solvent not in smx_data["builtin_solvent"]:
            smx_solvent = "other"
            if solvent not in smx_data["custom_solvent"]:
                raise Exception("Don't know the SMx parameters for "
                                "solvent '{}'".format(solvent))
            implicit_solvent['solvent_data'] = \
                smx_data["custom_solvent"][solvent]
        else:
            smx_solvent = solvent
        rem_options = dict(solvent_method=solvent_method.lower(),
                           smx_solvent=smx_solvent.lower())
        return rem_options, smx_solvent

    def set_solvent_method(self, qctask_sol, solvent, solvent_method, use_vdW_surface=False):
        implicit_solvent = dict()
        implicit_solvent['solvent_name'] = solvent
        if solvent_method.lower() in ['cpcm', 'ief-pcm']:
            if solvent_method.lower() == 'ief-pcm':
                solvent_theory = 'ssvpe'
            else:
                solvent_theory = 'cpcm'
            dielectric_constant, probe_radius, solvent_name = self.get_dielectric_constant(solvent, use_vdW_surface)
            qctask_sol.use_pcm(solvent_params={"Dielectric": dielectric_constant},
                               pcm_params={'Theory': solvent_theory,
                                           "SASrad": probe_radius})
            implicit_solvent['model'] = "{}_at_surface{:.2f}".format(solvent_method.lower(), probe_radius)
            implicit_solvent['dielectric_constant'] = dielectric_constant
            implicit_solvent['solvent_probe_radius'] = probe_radius
            implicit_solvent['radii'] = 'uff'
            implicit_solvent['vdwscale'] = 1.1
            implicit_solvent['solvent_name'] = solvent_name
        elif solvent_method.lower() == 'cosmo':
            dielectric_constant, dummy, solvent_name = self.get_dielectric_constant(solvent, use_vdW_surface)
            qctask_sol.use_cosmo(dielectric_constant)
            implicit_solvent['model'] = 'cosmo'
            implicit_solvent['dielectric_constant'] = dielectric_constant
            implicit_solvent['solvent_name'] = solvent_name
        elif solvent_method.lower() in ['sm12mk', 'sm12chelpg', 'sm12',
                                        'sm8']:
            implicit_solvent['model'] = solvent_method.lower()
            rem_options, smx_solvent = self.get_smx_solvent(implicit_solvent, solvent, solvent_method)
            qctask_sol.params['rem'].update(rem_options)
            implicit_solvent['smx_solvent'] = smx_solvent
        else:
            raise Exception("Don't know how to setup solvent model '{}'".
                            format(solvent_method))
        if not self.large:
            qctask_sol.set_dft_grid(128, 302)
            qctask_sol.set_integral_threshold(12)
            qctask_sol.set_scf_convergence_threshold(8)
        else:
            qctask_sol.set_scf_algorithm_and_iterations(iterations=100)
            if solvent_method.lower() in ['cpcm', 'ief-pcm']:
                qctask_sol.params['pcm'].update({"hpoints": 194,
                                                 "heavypoints": 194})
        return implicit_solvent

    def sp_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db,
              solvent_method="ief-pcm", use_vdW_surface=False, solvent="water", priority=None, qm_method=None,
              population_method=None, task_type_name=None):
        if not qm_method:
            qm_method = "B3LYP/6-31+G*"
        spec = self.base_spec()
        if priority:
            spec['_priority'] = priority
        task_type = task_type_name if task_type_name else "single point energy"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + qm_method + " " + task_type
        title += "\n Gas Phase"
        exchange, correlation, basis_set,  aux_basis, rem_params, method_token = self. \
            get_exchange_correlation_basis_auxbasis_remparams(qm_method)
        if population_method:
            if not rem_params:
                rem_params = dict()
            if population_method.lower() == "nbo":
                rem_params["nbo"] = 1
            elif population_method.lower() == "chelpg":
                rem_params["chelpg"] = True
        qctask_vac = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title, exchange=exchange, correlation=correlation,
                            basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        if not self.large:
            qctask_vac.set_dft_grid(128, 302)
            qctask_vac.set_integral_threshold(12)
            qctask_vac.set_scf_convergence_threshold(8)
        else:
            qctask_vac.set_scf_algorithm_and_iterations(iterations=100)

        title = " Solution Phase, {}".format(solvent)
        qctask_sol = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title, exchange=exchange, correlation=correlation,
                            basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        qctask_sol.set_scf_initial_guess(guess="read")
        implicit_solvent = self.set_solvent_method(qctask_sol, solvent, solvent_method, use_vdW_surface)

        qcinp = QcInput([qctask_vac, qctask_sol])
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        spec['run_tags']['methods'] = method_token
        spec["qm_method"] = qm_method
        spec['implicit_solvent'] = implicit_solvent
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemSinglePointEnergyDBInsertionTask
        fw_sp_cal = Firework([QChemTask()],
                             spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_sp_db = Firework([QChemSinglePointEnergyDBInsertionTask()],
                            spec=spec_db, name=task_name_db, fw_id=fw_id_db)
        return fw_sp_cal, fw_sp_db

    def vacuum_only_sp_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db,
                          priority=None, qm_method=None, population_method=None,
                          mixed_basis_generator=None, mixed_aux_basis_generator=None,
                          super_mol_snlgroup_id=None, super_mol_egsnl=None,
                          super_mol_inchi_root=None, ghost_atoms=None, bs_overlap=False):
        if not qm_method:
            qm_method = "B3LYP/6-31+G*"
        spec = self.base_spec()
        if priority:
            spec['_priority'] = priority
        if super_mol_snlgroup_id:
            from rubicon.workflows.bsse_wf import BSSEFragments
            task_type = "bsse {} fragment".format(BSSEFragments.OVERLAPPED if bs_overlap else BSSEFragments.ISOLATED)
        else:
            task_type = "vacuum only single point energy"
        if mixed_basis_generator or mixed_aux_basis_generator:
            population_method = population_method if population_method else "nbo"
            task_type = "atomic charge"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + qm_method + " " + task_type
        title += "\n Gas Phase"
        exchange, correlation, basis_set,  aux_basis, rem_params, method_token = self. \
            get_exchange_correlation_basis_auxbasis_remparams(qm_method)
        if population_method:
            if not rem_params:
                rem_params = dict()
            if population_method.lower() == "nbo":
                rem_params["nbo"] = 1
            elif population_method.lower() == "chelpg":
                rem_params["chelpg"] = True
            elif population_method.lower() == "hirshfeld":
                rem_params["hirshfeld"] = True
        ga = ghost_atoms if bs_overlap else None
        qctask_vac = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title, exchange=exchange, correlation=correlation,
                            basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params,
                            ghost_atoms=ga)
        if (not self.large) and (mixed_basis_generator is None and mixed_aux_basis_generator is None):
            qctask_vac.set_dft_grid(128, 302)
            qctask_vac.set_integral_threshold(12)
            qctask_vac.set_scf_convergence_threshold(8)
        else:
            qctask_vac.set_scf_algorithm_and_iterations(iterations=100)

        qcinp = QcInput([qctask_vac])
        spec["qcinp"] = qcinp.as_dict()
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        spec['run_tags']['methods'] = method_token
        spec["qm_method"] = qm_method
        if super_mol_snlgroup_id:
            spec["run_tags"]["super_mol_snlgroup_id"] = super_mol_snlgroup_id
            spec["snlgroup_id"] = super_mol_snlgroup_id
            spec["egsnl"] = super_mol_egsnl
            spec["inchi_root"] = super_mol_inchi_root
        if ghost_atoms:
            spec["run_tags"]["ghost_atoms"] = sorted(set(ghost_atoms))
            from rubicon.workflows.bsse_wf import BSSEFragments
            spec["run_tags"]["bsse_fragment_type"] = BSSEFragments.OVERLAPPED if bs_overlap else BSSEFragments.ISOLATED
        if mixed_basis_generator:
            spec["_mixed_basis_set_generator"] = mixed_basis_generator
        if mixed_aux_basis_generator:
            spec["_mixed_aux_basis_set_generator"] = mixed_aux_basis_generator
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemSinglePointEnergyDBInsertionTask
        fw_sp_cal = Firework([QChemTask()],
                             spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_sp_db = Firework([QChemSinglePointEnergyDBInsertionTask()],
                            spec=spec_db, name=task_name_db, fw_id=fw_id_db)
        return fw_sp_cal, fw_sp_db

    def aimd_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db,
                num_steps, time_step, temperature, priority=None, qm_method=None):
        if not qm_method:
            qm_method = "B3LYP/6-31+G*"
        spec = self.base_spec()
        if priority:
            spec['_priority'] = priority
        task_type = "ab initio molecule dynamics"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + qm_method + " " + task_type
        exchange, correlation, basis_set,  aux_basis, rem_params, method_token = self. \
            get_exchange_correlation_basis_auxbasis_remparams(qm_method)
        if rem_params is None:
            rem_params = dict()
        rem_params["aimd_method"] = "bomd"
        rem_params["time_step"] = time_step
        rem_params["aimd_steps"] = num_steps
        rem_params["aimd_init_veloc"] = "thermal"
        rem_params["aimd_temp"] = temperature
        rem_params["fock_extrap_order"] = 6
        rem_params["fock_extrap_points"] = 12
        qctask_vac = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                            jobtype="aimd", title=title, exchange=exchange, correlation=correlation,
                            basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        qctask_vac.set_scf_algorithm_and_iterations(iterations=100)

        qcinp = QcInput([qctask_vac])
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        spec['run_tags']['methods'] = method_token
        spec['run_tags']['rem_params'] = rem_params
        spec["qm_method"] = qm_method
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemAIMDDBInsertionTask
        fw_sp_cal = Firework([QChemTask()],
                             spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_sp_db = Firework([QChemAIMDDBInsertionTask()],
                            spec=spec_db, name=task_name_db, fw_id=fw_id_db)
        return fw_sp_cal, fw_sp_db