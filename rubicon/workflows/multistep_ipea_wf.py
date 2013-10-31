import copy
from fireworks.core.firework import Workflow, FireWork
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.nwchemio import NwTask, NwInput
from pymatgen.symmetry.pointgroup import PointGroupAnalyzer
from rubicon.dupefinder.dupefinder_eg import DupeFinderEG
from rubicon.firetasks.nwchem_task import NWChemTask

__author__ = 'xiaohuiqu'


class NWChemFireWorkCreator():
    def __init__(self, mol, molname, mission, additional_user_tags=None, dupefinder=DupeFinderEG,
                 priority=1):
        theory_directive = {"iterations": 300, "vectors": "atomic"}
        symmetry_options = None
        pga = PointGroupAnalyzer(mol)
        if pga.sch_symbol == 'D*h' or pga.sch_symbol == 'C*v':
            # linear molecule, turn off symmetry
            symmetry_options = ['c1']
        initial_inchi = self.get_inchi(mol)
        user_tags = {'mission': mission,
                     "initial_inchi": initial_inchi,
                     "molname": molname}
        if additional_user_tags:
            user_tags.update(additional_user_tags)

        self.td = lambda: copy.deepcopy(theory_directive)
        self.sym = lambda: copy.deepcopy(symmetry_options)
        self.bs = '6-31+G*'
        self.ut = lambda: copy.deepcopy(user_tags)
        self.mol = mol
        self.dupefinder = dupefinder
        self.priority = priority


    def get_inchi(self, mol):
        bb = BabelMolAdaptor(mol)
        pbmol = bb.pybel_mol
        return pbmol.write("inchi").strip()


    def geom_fw(self, charge_shift, fw_id):
        charge = self.mol.charge + charge_shift
        tasks_geom = [NwTask.dft_task(self.mol, charge=charge,
                                      operation="optimize",
                                      xc="b3lyp", basis_set=self.bs,
                                      theory_directives=self.td(),
                                      alternate_directives={"driver":
                                                            {"maxiter": 300}})]
        nwi = NwInput(self.mol, tasks_geom, symmetry_options=self.sym())
        spec = nwi.to_dict
        spec['user_tags'] = self.ut()
        charge_state_name = {0: "original", 1: "cation", -1: "anion"}
        spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
        spec['user_tags']['charge_shift'] = charge_shift
        spec['_priority'] = self.priority
        spec['_dupefinder'] = self.dupefinder().to_dict()
        spec['user_tags']['methods'] = "B3LYP/" + self.bs
        spec['task_type'] = 'Geometry Optimization'
        task_name = charge_state_name[charge_shift] + ' geom opt'
        from rubicon.firetasks.multistep_nwchem_task \
            import NWChemGeomOptDBInsertionTask
        fw_geom = FireWork([NWChemTask(),
                            NWChemGeomOptDBInsertionTask()],
                           spec=spec, name=task_name, fw_id=fw_id)
        return fw_geom


    def freq_fw(self, charge_shift, fw_id):
        charge = self.mol.charge + charge_shift
        tasks_geom = [NwTask.dft_task(self.mol, charge=charge,
                                      operation="freq",
                                      xc="b3lyp", basis_set=self.bs,
                                      theory_directives=self.td())]
        nwi = NwInput(self.mol, tasks_geom, symmetry_options=self.sym())
        spec = nwi.to_dict
        spec['user_tags'] = self.ut()
        charge_state_name = {0: "original", 1: "cation", -1: "anion"}
        spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
        spec['user_tags']['charge_shift'] = charge_shift
        spec['_priority'] = self.priority
        spec['_dupefinder'] = self.dupefinder().to_dict()
        spec['user_tags']['methods'] = "B3LYP/" + self.bs
        spec['task_type'] = 'Vibrational Frequency'
        task_name = charge_state_name[charge_shift] + ' freq'
        from rubicon.firetasks.multistep_nwchem_task \
            import NWChemFrequencyDBInsertionTask
        fw_freq = FireWork([NWChemTask(),
                            NWChemFrequencyDBInsertionTask()],
                           spec=spec, name=task_name, fw_id=fw_id)
        return fw_freq


    def sp_fw(self, charge_shift, fw_id):
        charge = self.mol.charge + charge_shift
        tasks_geom = [NwTask.dft_task(self.mol, charge=charge,
                                      operation="energy",
                                      xc="b3lyp", basis_set=self.bs,
                                      theory_directives=self.td()),
                      NwTask.dft_task(self.mol, charge=charge,
                                      operation="energy",
                                      xc="b3lyp", basis_set=self.bs,
                                      theory_directives=self.td(),
                                      alternate_directives={'cosmo':
                                                            {"dielec": 78.0}})]
        nwi = NwInput(self.mol, tasks_geom, symmetry_options=self.sym())
        spec = nwi.to_dict
        spec['user_tags'] = self.ut()
        charge_state_name = {0: "original", 1: "cation", -1: "anion"}
        spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
        spec['user_tags']['charge_shift'] = charge_shift
        spec['_priority'] = self.priority
        spec['_dupefinder'] = self.dupefinder().to_dict()
        spec['user_tags']['methods'] = "B3LYP/" + self.bs
        spec['task_type'] = 'Single Point Energy'
        task_name = charge_state_name[charge_shift] + ' single point energy'
        from rubicon.firetasks.multistep_nwchem_task \
            import NWChemSinglePointEnergyDBInsertionTask
        fw_sp = FireWork([NWChemTask(),
                          NWChemSinglePointEnergyDBInsertionTask()],
                         spec=spec, name=task_name, fw_id=fw_id)
        return fw_sp


def multistep_ipea_fws(mol, name, mission, dupefinder=DupeFinderEG, priority=1, parent_fwid = None):
    fw_creator = NWChemFireWorkCreator(mol, name, mission, dupefinder, priority)
    fwid_base = 1
    if parent_fwid:
        if not (isinstance(parent_fwid, int) or isinstance(parent_fwid, list)):
            raise ValueError("Parent FireWork ID must be integer or list")
        parent_fwid = parent_fwid if isinstance(parent_fwid, list) else [parent_fwid]
        fwid_base = max(parent_fwid)
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


def mol_to_ipea_wf(mol, name, mission, dupefinder=DupeFinderEG, priority=1, parent_fwid = None):
    fireworks, links_dict = multistep_ipea_fws(mol, name, mission, dupefinder, priority,
                                               parent_fwid)
    return Workflow(fireworks, links_dict, name)