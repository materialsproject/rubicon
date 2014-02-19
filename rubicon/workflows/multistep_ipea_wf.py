import copy
import itertools
from fireworks.core.firework import Workflow, FireWork, Tracker
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.qchemio import QcInput, QcTask
from rubicon.dupefinders.dupefinder_eg import DupeFinderEG
from rubicon.firetasks.qchem_task import QChemTask

__author__ = 'xiaohuiqu'


class QChemFireWorkCreator():
    def __init__(self, mol, molname, mission, additional_user_tags=None,
                 dupefinder=None, priority=1, update_spec=None, large=False):
        self.bs = '6-31+G*'
        self.dft = 'B3LYP'
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
        spec['run_tags']['methods'] = [self.bs.lower(), self.dft.lower()]
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
        charge_state = {-1: "anion", 0: "neutral", 1: "cation"}
        spin_state = {1: "singlet", 2: "doublet", 3: "triplet"}
        return spin_state[spin_multiplicity] + " " + charge_state[charge]

    def geom_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db):
        task_type = "geometry optimization"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + self.dft + " " + \
            self.bs + " " + task_type
        if self.large:
            title = self.molname + " " + state_name + " PBE-D3/6-31+G* " + \
                    task_type
            qctask = QcTask(self.mol, charge=charge,
                            spin_multiplicity=spin_multiplicity,
                            jobtype="opt", title=title,
                            exchange='pbe', correlation='pbe',
                            basis_set=self.bs,
                            rem_params={"DFT_D": "EMPIRICAL_GRIMME3",
                                        "DFT_D3_S6": 1000,
                                        "DFT_D3_RS6": 1217,
                                        "DFT_D3_S8": 722,
                                        "DFT_D3_3BODY": False})
            qctask.set_geom_max_iterations(200)
            qctask.set_scf_algorithm_and_iterations(iterations=100)
            qctask.scale_geom_opt_threshold(gradient=1.0,
                                            displacement=10.0,
                                            energy=10.0)
        else:
            qctask = QcTask(self.mol, charge=charge,
                            spin_multiplicity=spin_multiplicity,
                            jobtype="opt", title=title,
                            exchange=self.dft, basis_set=self.bs)
        qctask.set_memory(total=1100)
        qcinp = QcInput([qctask])
        spec = self.base_spec()
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemGeomOptDBInsertionTask
        fw_geom_cal = FireWork([QChemTask()],
                               spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_geom_db = FireWork([QChemGeomOptDBInsertionTask()],
                              spec=spec_db, name=task_name_db, fw_id=fw_id_db)

        return fw_geom_cal, fw_geom_db

    def freq_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db):
        task_type = "vibrational frequency"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + self.dft + " " +\
            self.bs + " " + task_type
        if self.large:
            title = self.molname + " " + state_name + " PBE-D3/6-31+G* " + \
                    task_type
            qctask = QcTask(self.mol, charge=charge,
                            spin_multiplicity=spin_multiplicity,
                            jobtype="freq", title=title,
                            exchange='pbe', correlation='pbe',
                            basis_set=self.bs,
                            rem_params={"DFT_D": "EMPIRICAL_GRIMME3",
                                        "DFT_D3_S6": 1000,
                                        "DFT_D3_RS6": 1217,
                                        "DFT_D3_S8": 722,
                                        "DFT_D3_3BODY": False})
            qctask.set_scf_algorithm_and_iterations(iterations=100)
        else:
            qctask = QcTask(self.mol, charge=charge,
                            spin_multiplicity=spin_multiplicity,
                            jobtype="freq", title=title,
                            exchange=self.dft, basis_set=self.bs)
        qctask.set_memory(total=1100)
        qcinp = QcInput([qctask])
        spec = self.base_spec()
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemFrequencyDBInsertionTask
        fw_freq_cal = FireWork([QChemTask()],
                               spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_freq_db = FireWork([QChemFrequencyDBInsertionTask()],
                              spec=spec_db, name=task_name_db, fw_id=fw_id_db)
        return fw_freq_cal, fw_freq_db

    def sp_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db):
        task_type = "single point energy"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + self.dft + " " + \
            self.bs + " " + task_type
        title += "\n Gas Phase"
        qctask_vac = QcTask(self.mol, charge=charge,
                            spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title,
                            exchange=self.dft, basis_set=self.bs,
                            rem_params={"CHELPG": True})
        if self.large:
            qctask_vac.set_scf_algorithm_and_iterations(iterations=100)
        qctask_vac.set_memory(total=1100)
        qctask_vac.set_dft_grid(128, 302)
        qctask_vac.set_integral_threshold(12)
        qctask_vac.set_scf_convergence_threshold(8)
        title = " Solution Phase"
        qctask_sol = QcTask(self.mol, charge=charge,
                            spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title,
                            exchange=self.dft, basis_set=self.bs,
                            rem_params={"CHELPG": True})
        if self.large:
            qctask_sol.set_scf_algorithm_and_iterations(iterations=100)
        qctask_sol.use_pcm()
        qctask_sol.set_scf_initial_guess(guess="read")
        qctask_sol.set_memory(total=1100)
        qctask_sol.set_dft_grid(128, 302)
        qctask_sol.set_integral_threshold(12)
        qctask_sol.set_scf_convergence_threshold(8)
        qcinp = QcInput([qctask_vac, qctask_sol])
        spec = self.base_spec()
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        implicit_solvent = dict()
        implicit_solvent['model'] = 'ief-pcm'
        implicit_solvent['dielectric_constant'] = 78.3553
        implicit_solvent['radii'] = 'uff'
        implicit_solvent['vdwscale'] = 1.1
        spec['implicit_solvent'] = implicit_solvent
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemSinglePointEnergyDBInsertionTask
        fw_sp_cal = FireWork([QChemTask()],
                             spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_sp_db = FireWork([QChemSinglePointEnergyDBInsertionTask()],
                            spec=spec_db, name=task_name_db, fw_id=fw_id_db)
        return fw_sp_cal, fw_sp_db


def multistep_ipea_fws(mol, name, mission, dupefinder=None, priority=1,
                       parent_fwid=None):
    large = False
    if len(mol) > 50:
        large = True
    fw_creator = QChemFireWorkCreator(mol=mol, molname=name, mission=mission,
                                      dupefinder=dupefinder,
                                      priority=priority, large=large)
    fwid_base = 1
    if parent_fwid:
        if not (isinstance(parent_fwid, int) or isinstance(parent_fwid, list)):
            raise ValueError("Parent FireWork ID must be integer or list")
        parent_fwid = parent_fwid if isinstance(parent_fwid, list) \
            else [parent_fwid]
        fwid_base = max(parent_fwid) + 1
    fireworks = []
    links_dict = dict()
    # the task in the order of anion, neutral, cation
    cgi_cal, ngi_cal, agi_cal = (None, None, None)
    cfi_db, nfi_db, afi_db = (None, None, None)
    cgi_db, ngi_db, agi_db = (None, None, None)
    charge = (-1, 0, 1)
    spin_multiplicity = (2, 1, 2)
    if len(mol) > 1:
        fw_ids = zip(* [iter(range(fwid_base + 0, fwid_base + 6))] * 2)
        fws = (fw_creator.geom_fw(ch, spin, fwid_cal, fwid_db)
               for ch, spin, (fwid_cal, fwid_db)
               in zip(charge, spin_multiplicity, fw_ids))
        (cgi_cal, cgi_db), (ngi_cal, ngi_db), (agi_cal, agi_db) = fw_ids
        fireworks.extend(itertools.chain.from_iterable(fws))
        links_dict.update(dict(fw_ids))

        if not large:
            fw_ids = zip(* [iter(range(fwid_base + 6, fwid_base + 6 + 6))] * 2)
            fws = (fw_creator.freq_fw(ch, spin, fwid_cal, fwid_db)
                   for ch, spin, (fwid_cal, fwid_db)
                   in zip(charge, spin_multiplicity, fw_ids))
            (cfi_cal, cfi_db), (nfi_cal, nfi_db), (afi_cal, afi_db) = fw_ids
            fireworks.extend(itertools.chain.from_iterable(fws))
            links_dict.update(dict(fw_ids))
            links_dict.update({cgi_db: cfi_cal,
                               ngi_db: nfi_cal,
                               agi_db: afi_cal})

    fw_ids = zip(* [iter(range(fwid_base + 12, fwid_base + 12 + 6))] * 2)
    fws = (fw_creator.sp_fw(ch, spin, fwid_cal, fwid_db)
           for ch, spin, (fwid_cal, fwid_db)
           in zip(charge, spin_multiplicity, fw_ids))
    (cspi_cal, cspi_db), (nspi_cal, nspi_db), (aspi_cal, aspi_db) = fw_ids
    links_dict.update((dict(fw_ids)))
    fireworks.extend(itertools.chain.from_iterable(fws))
    if len(mol) > 1:
        if large:
            links_dict.update({cgi_db: cspi_cal, ngi_db: nspi_cal,
                               agi_db: aspi_cal})
        else:
            links_dict.update({cfi_db: cspi_cal, nfi_db: nspi_cal,
                               afi_db: aspi_cal})
        links_dict.update({nspi_db: [cgi_cal, agi_cal]})
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = ngi_cal
    else:
        links_dict.update({nspi_db: [cspi_cal, aspi_cal]})
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = nspi_cal
    return fireworks, links_dict


def mol_to_ipea_wf(mol, name, mission, dupefinder=None, priority=1,
                   parent_fwid=None):
    fireworks, links_dict = multistep_ipea_fws(mol, name, mission, dupefinder,
                                               priority, parent_fwid)
    return Workflow(fireworks, links_dict, name)