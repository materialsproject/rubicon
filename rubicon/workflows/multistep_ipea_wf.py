import copy
from fireworks.core.firework import FireWork, Workflow
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.nwchemio import NwTask, NwInput
from pymatgen.symmetry.pointgroup import PointGroupAnalyzer
from rubicon.firetasks.multistep_nwchem_task import NWChemGeomOptDBInsertionTask, NWChemFrequencyDBInsertionTask, NWChemSinglePointEnergyDBInsertionTask
from rubicon.firetasks.nwchem_task import NWChemTask

__author__ = 'xiaohuiqu'

def get_inchi(mol):
    bb = BabelMolAdaptor(mol)
    pbmol = bb.pybel_mol
    return pbmol.write("inchi")


def mol_to_geom_fw(mol, bs, charge_shift, fw_id, name, sym, td, user_tags):
    charge = mol.charge + charge_shift
    tasks_geom = [NwTask.dft_task(mol, charge=charge, operation="optimize",
                                  xc="b3lyp", basis_set=bs,
                                  theory_directives=td(),
                                  alternate_directives={"driver":
                                                        {"maxiter": 300}})]
    nwi = NwInput(mol, tasks_geom, symmetry_options=sym())
    spec = nwi.to_dict
    spec['user_tags'] = copy.deepcopy(user_tags)
    charge_state_name = {0: "original", 1: "cation", -1: "anion"}
    spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
    task_name = name + ' ' + charge_state_name[charge_shift] + ' geom opt'
    fw_geom = FireWork([NWChemTask(), NWChemGeomOptDBInsertionTask()],
                       spec=spec, name=task_name, fw_id=fw_id)
    return fw_geom


def mol_to_freq_fw(mol, bs, charge_shift, fw_id, name, sym, td, user_tags):
    charge = mol.charge + charge_shift
    tasks_geom = [NwTask.dft_task(mol, charge=charge, operation="freq",
                                  xc="b3lyp", basis_set=bs,
                                  theory_directives=td())]
    nwi = NwInput(mol, tasks_geom, symmetry_options=sym())
    spec = nwi.to_dict
    spec['user_tags'] = copy.deepcopy(user_tags)
    charge_state_name = {0: "original", 1: "cation", -1: "anion"}
    spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
    task_name = name + ' ' + charge_state_name[charge_shift] + ' freq'
    fw_freq = FireWork([NWChemTask(), NWChemFrequencyDBInsertionTask()],
                       spec=spec, name=task_name, fw_id=fw_id)
    return fw_freq


def mol_to_sp_fw(mol, bs, charge_shift, fw_id, name, sym, td, user_tags):
    charge = mol.charge + charge_shift
    tasks_geom = [NwTask.dft_task(mol, charge=charge, operation="energy",
                                  xc="b3lyp", basis_set=bs,
                                  theory_directives=td()),
                  NwTask.dft_task(mol, charge=charge, operation="energy",
                                  xc="b3lyp", basis_set=bs,
                                  theory_directives=td(),
                                  alternate_directives={'cosmo':
                                                        {"dielec": 78.0}})]
    nwi = NwInput(mol, tasks_geom, symmetry_options=sym())
    spec = nwi.to_dict
    spec['user_tags'] = copy.deepcopy(user_tags)
    charge_state_name = {0: "original", 1: "cation", -1: "anion"}
    spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
    task_name = name + ' ' + charge_state_name[charge_shift] + ' single point energy'
    fw_sp = FireWork([NWChemTask(), NWChemSinglePointEnergyDBInsertionTask()],
                     spec=spec, name=task_name, fw_id=fw_id)
    return fw_sp



def mol_to_ipea_wf(mol, name, mission):

    theory_directive = {"iterations": 300, "vectors": "atomic"}
    symmetry_options = None
    pga = PointGroupAnalyzer(mol)
    if pga.sch_symbol == 'D*h' or pga.sch_symbol == 'C*v':
        # linear molecule, turn off symmetry
        symmetry_options = ['c1']
    initial_inchi = get_inchi(mol)

    td = lambda: copy.deepcopy(theory_directive)
    sym = lambda: copy.deepcopy(symmetry_options)
    bs = '6-31+G*'
    user_tags = {'mission': mission, "initial_inchi": initial_inchi}

    fireworks = []
    links_dict = dict()

    # the task in the order of anion, neutral, cation
    charge_shifts = (-1, 0, 1)
    fw_ids = range(0, 3)
    fws = (mol_to_geom_fw(mol, '6-31+G*', cs, fwid, name, sym, td, user_tags)
           for cs, fwid in zip(charge_shifts, fw_ids))
    cg_fwid, ng_fwid, ag_fwid = fw_ids
    fireworks.extend(fws)

    cf_fwid, nf_fwid, af_fwid = (None, None, None)
    if len(mol) > 1:
        fw_ids = range(3, 3+3)
        fws = (mol_to_freq_fw(mol, '6-31+G*', cs, fwid, name, sym, td,
                              user_tags)
               for cs, fwid in zip(charge_shifts, fw_ids))
        cf_fwid, nf_fwid, af_fwid = fw_ids
        fireworks.extend(fws)
        links_dict.update({cg_fwid: cf_fwid, ng_fwid: nf_fwid,
                           ag_fwid: af_fwid})

    charge_shifts = (-1, 0, 1)
    fw_ids = range(6, 6+3)
    fws = (mol_to_sp_fw(mol, '6-31+G*', cs, fwid, name, sym, td, user_tags)
           for cs, fwid in zip(charge_shifts, fw_ids))
    csp_fwid, nsp_fwid, asp_fwid = fw_ids
    fireworks.extend(fws)
    if len(mol) > 1:
        links_dict.update({cf_fwid: csp_fwid, nf_fwid: nsp_fwid,
                           af_fwid: asp_fwid})
    else:
        links_dict.update({cg_fwid: csp_fwid, ng_fwid: nsp_fwid,
                           ag_fwid: asp_fwid})
    links_dict.update({nsp_fwid: [cg_fwid, ag_fwid]})

    return Workflow(fireworks, links_dict, name)