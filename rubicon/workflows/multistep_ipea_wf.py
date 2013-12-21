import copy
from fireworks.core.firework import Workflow, FireWork, Tracker
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.qchemio import QcInput, QcTask
from rubicon.dupefinders.dupefinder_eg import DupeFinderEG
from rubicon.firetasks.qchem_task import QChemTask

__author__ = 'xiaohuiqu'


class QChemFireWorkCreator():
    def __init__(self, mol, molname, mission, additional_user_tags=None,
                 dupefinder=DupeFinderEG, priority=1, update_spec=None):
        self.bs = '6-31+G*'
        self.dft = 'B3LYP'
        self.molname = molname
        self.mol = mol
        initial_inchi = self.get_inchi(mol)
        user_tags = {'mission': mission,
                     "molname": molname}
        if additional_user_tags:
            user_tags.update(additional_user_tags)
        spec = dict()
        spec['user_tags'] = user_tags
        spec['_priority'] = priority
        spec['_dupefinder'] = dupefinder().to_dict()
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

    def geom_fw(self, charge, spin_multiplicity, fw_id):
        task_type = "geometry optimization"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + self.dft + " " + \
            self.bs + " " + task_type
        qctask = QcTask(self.mol, charge=charge,
                        spin_multiplicity=spin_multiplicity,
                        jobtype="opt", title=title,
                        exchange=self.dft, basis_set=self.bs)
        qctask.set_memory(total=28000, static=3000)
        qcinp = QcInput([qctask])
        spec = self.base_spec()
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemGeomOptDBInsertionTask
        fw_geom = FireWork([QChemTask(),
                            QChemGeomOptDBInsertionTask()],
                           spec=spec, name=task_name, fw_id=fw_id)
        return fw_geom

    def freq_fw(self, charge, spin_multiplicity, fw_id):
        task_type = "vibrational frequency"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + self.dft + " " +\
            self.bs + " " + task_type
        qctask = QcTask(self.mol, charge=charge,
                        spin_multiplicity=spin_multiplicity,
                        jobtype="freq", title=title,
                        exchange=self.dft, basis_set=self.bs)
        qctask.set_memory(total=28000, static=3000)
        qcinp = QcInput([qctask])
        spec = self.base_spec()
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemFrequencyDBInsertionTask
        fw_freq = FireWork([QChemTask(),
                            QChemFrequencyDBInsertionTask()],
                           spec=spec, name=task_name, fw_id=fw_id)
        return fw_freq

    def sp_fw(self, charge, spin_multiplicity, fw_id):
        task_type = "single point energy"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + self.dft + " " + \
            self.bs + " " + task_type
        title += "\n Gas Phase"
        qctask_vac = QcTask(self.mol, charge=charge,
                            spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title,
                            exchange=self.dft, basis_set=self.bs)
        qctask_vac.set_memory(total=28000, static=3000)
        qctask_vac.set_dft_grid(128, 302)
        qctask_vac.set_integral_threshold(12)
        qctask_vac.set_scf_convergence_threshold(8)
        title = " Solution Phase"
        qctask_sol = QcTask(self.mol, charge=charge,
                            spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title,
                            exchange=self.dft, basis_set=self.bs)
        qctask_sol.use_pcm()
        qctask_sol.set_scf_initial_guess(guess="read")
        qctask_sol.set_memory(total=28000, static=3000)
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
        fw_freq = FireWork([QChemTask(),
                            QChemSinglePointEnergyDBInsertionTask()],
                           spec=spec, name=task_name, fw_id=fw_id)
        return fw_freq


def multistep_ipea_fws(mol, name, mission, dupefinder=DupeFinderEG, priority=1,
                       parent_fwid=None):
    fw_creator = QChemFireWorkCreator(mol, name, mission, None, dupefinder,
                                      priority)
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
    cg_fwid, ng_fwid, ag_fwid = (None, None, None)
    cf_fwid, nf_fwid, af_fwid = (None, None, None)
    charge = (-1, 0, 1)
    spin_multiplicity = (2, 1, 2)
    if len(mol) > 1:
        fw_ids = range(fwid_base + 0, fwid_base + 3)
        fws = (fw_creator.geom_fw(ch, spin, fwid)
               for ch, spin, fwid in zip(charge, spin_multiplicity, fw_ids))
        cg_fwid, ng_fwid, ag_fwid = fw_ids
        fireworks.extend(fws)

        fw_ids = range(fwid_base + 3, fwid_base + 3 + 3)
        fws = (fw_creator.freq_fw(ch, spin, fwid)
               for ch, spin, fwid in zip(charge, spin_multiplicity, fw_ids))
        cf_fwid, nf_fwid, af_fwid = fw_ids
        fireworks.extend(fws)
        links_dict.update({cg_fwid: cf_fwid, ng_fwid: nf_fwid,
                           ag_fwid: af_fwid})
    fw_ids = range(fwid_base + 6, fwid_base + 6 + 3)
    fws = (fw_creator.sp_fw(ch, spin, fwid)
           for ch, spin, fwid in zip(charge, spin_multiplicity, fw_ids))
    csp_fwid, nsp_fwid, asp_fwid = fw_ids
    fireworks.extend(fws)
    if len(mol) > 1:
        links_dict.update({cf_fwid: csp_fwid, nf_fwid: nsp_fwid,
                           af_fwid: asp_fwid})
        links_dict.update({nsp_fwid: [cg_fwid, ag_fwid]})
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = ng_fwid
    else:
        links_dict.update({nsp_fwid: [csp_fwid, asp_fwid]})
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = nsp_fwid
    return fireworks, links_dict


def mol_to_ipea_wf(mol, name, mission, dupefinder=DupeFinderEG, priority=1,
                   parent_fwid=None):
    fireworks, links_dict = multistep_ipea_fws(mol, name, mission, dupefinder,
                                               priority, parent_fwid)
    return Workflow(fireworks, links_dict, name)