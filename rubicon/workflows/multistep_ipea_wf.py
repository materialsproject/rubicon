import copy
from fireworks.core.firework import Workflow, FireWork, Tracker
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.qchemio import QcBatchInput, QcInput
from rubicon.dupefinders.dupefinder_eg import DupeFinderEG
from rubicon.firetasks.qchem_task import QChemTask

__author__ = 'xiaohuiqu'


class QChemFireWorkCreator():
    def __init__(self, mol, molname, mission, additional_user_tags=None,
                 dupefinder=DupeFinderEG, priority=1, update_spec=None):
        initial_inchi = self.get_inchi(mol)
        user_tags = {'mission': mission,
                     "initial_inchi": initial_inchi,
                     "molname": molname}
        if additional_user_tags:
            user_tags.update(additional_user_tags)
        self.bs = '6-31+G*'
        self.dft = 'B3LYP'
        self.molname = molname
        self.ut = lambda: copy.deepcopy(user_tags)
        self.mol = mol
        self.dupefinder = dupefinder
        self.priority = priority
        self.update_spec = update_spec

    @staticmethod
    def get_inchi(mol):
        bb = BabelMolAdaptor(mol)
        pbmol = bb.pybel_mol
        return pbmol.write("inchi").strip()

    def geom_fw(self, charge_shift, fw_id):
        charge = self.mol.charge + charge_shift
        charge_state_name = {0: "original", 1: "cation", -1: "anion"}
        title = self.molname + " " + self.dft + " " + self.bs + \
            charge_state_name[charge_shift] + " Geometry Optimization"
        qcinp = QcInput(self.mol, charge=charge, jobtype="opt", title=title,
                        exchange=self.dft, basis_set=self.bs)
        qcbat = QcBatchInput([qcinp])
        spec = dict()
        spec["qcinp"] = qcbat.to_dict
        spec['user_tags'] = self.ut()
        spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
        spec['user_tags']['charge_shift'] = charge_shift
        spec['_priority'] = self.priority
        if self.dupefinder:
            spec['_dupefinder'] = self.dupefinder().to_dict()
        tracker_out = Tracker("mol.qcout", nlines=20)
        tracker_std = Tracker("mol.qclog", nlines=10)
        spec["_trackers"] = [tracker_out.to_dict(), tracker_std.to_dict()]
        spec['user_tags']['methods'] = "B3LYP/" + self.bs
        spec['task_type'] = charge_state_name[charge_shift] + \
            ' Geometry Optimization'
        if self.update_spec:
            spec.update(self.update_spec)
        task_name = charge_state_name[charge_shift] + ' geom opt'
        from rubicon.firetasks.multistep_qchem_task \
            import QChemGeomOptDBInsertionTask
        fw_geom = FireWork([QChemTask(),
                            QChemGeomOptDBInsertionTask()],
                           spec=spec, name=task_name, fw_id=fw_id)
        return fw_geom

    def freq_fw(self, charge_shift, fw_id):
        charge = self.mol.charge + charge_shift
        charge_state_name = {0: "original", 1: "cation", -1: "anion"}
        title = self.molname + " " + self.dft + " " + self.bs + \
            charge_state_name[charge_shift] + " Vibrational Frequency Analysis"
        qcinp = QcInput(self.mol, charge=charge, jobtype="freq", title=title,
                        exchange=self.dft, basis_set=self.bs)
        qcbat = QcBatchInput([qcinp])
        spec = dict()
        spec["qcinp"] = qcbat.to_dict
        spec['user_tags'] = self.ut()
        spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
        spec['user_tags']['charge_shift'] = charge_shift
        spec['_priority'] = self.priority
        if self.dupefinder:
            spec['_dupefinder'] = self.dupefinder().to_dict()
        tracker_out = Tracker("mol.qcout", nlines=20)
        tracker_std = Tracker("mol.qclog", nlines=10)
        spec["_trackers"] = [tracker_out.to_dict(), tracker_std.to_dict()]
        spec['user_tags']['methods'] = "B3LYP/" + self.bs
        spec['task_type'] = charge_state_name[charge_shift] + \
            ' Vibrational Frequency'
        if self.update_spec:
            spec.update(self.update_spec)
        task_name = charge_state_name[charge_shift] + ' Vibrational Frequency'
        from rubicon.firetasks.multistep_qchem_task \
            import QChemFrequencyDBInsertionTask
        fw_freq = FireWork([QChemTask(),
                            QChemFrequencyDBInsertionTask()],
                           spec=spec, name=task_name, fw_id=fw_id)
        return fw_freq

    def sp_fw(self, charge_shift, fw_id):
        charge = self.mol.charge + charge_shift
        charge_state_name = {0: "original", 1: "cation", -1: "anion"}
        title = self.molname + " " + self.dft + " " + self.bs + \
            charge_state_name[charge_shift] + " Single Point Energy"
        title += "\n Gas Phase"
        qcinp_vac = QcInput(self.mol, charge=charge, jobtype="sp", title=title,
                            exchange=self.dft, basis_set=self.bs)
        title = " Solution Phase"
        qcinp_sol = QcInput(self.mol, charge=charge, jobtype="sp", title=title,
                            exchange=self.dft, basis_set=self.bs)
        qcinp_sol.use_pcm()
        qcinp_sol.set_scf_initial_guess(guess="read")
        qcbat = QcBatchInput([qcinp_vac, qcinp_sol])
        spec = dict()
        spec["qcinp"] = qcbat.to_dict
        spec['user_tags'] = self.ut()
        spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
        spec['user_tags']['charge_shift'] = charge_shift
        spec['_priority'] = self.priority
        if self.dupefinder:
            spec['_dupefinder'] = self.dupefinder().to_dict()
        tracker_out = Tracker("mol.qcout", nlines=20)
        tracker_std = Tracker("mol.qclog", nlines=10)
        spec["_trackers"] = [tracker_out.to_dict(), tracker_std.to_dict()]
        spec['user_tags']['methods'] = "B3LYP/" + self.bs
        spec['task_type'] = charge_state_name[charge_shift] + \
            ' Single Point Energy'
        if self.update_spec:
            spec.update(self.update_spec)
        task_name = charge_state_name[charge_shift] + ' Single Point Energy'
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
    if len(mol) > 1:
        charge_shifts = (-1, 0, 1)
        fw_ids = range(fwid_base + 0, fwid_base + 3)
        fws = (fw_creator.geom_fw(cs, fwid)
               for cs, fwid in zip(charge_shifts, fw_ids))
        cg_fwid, ng_fwid, ag_fwid = fw_ids
        fireworks.extend(fws)

        fw_ids = range(fwid_base + 3, fwid_base + 3 + 3)
        fws = (fw_creator.freq_fw(cs, fwid)
               for cs, fwid in zip(charge_shifts, fw_ids))
        cf_fwid, nf_fwid, af_fwid = fw_ids
        fireworks.extend(fws)
        links_dict.update({cg_fwid: cf_fwid, ng_fwid: nf_fwid,
                           ag_fwid: af_fwid})
    charge_shifts = (-1, 0, 1)
    fw_ids = range(fwid_base + 6, fwid_base + 6 + 3)
    fws = (fw_creator.sp_fw(cs, fwid)
           for cs, fwid in zip(charge_shifts, fw_ids))
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