from fireworks.core.firework import Workflow
from rubicon.firetasks.multistep_nwchem_task import NWChemFireWorkCreator

__author__ = 'xiaohuiqu'


def mol_to_ipea_wf(mol, name, mission):

    fw_creator = NWChemFireWorkCreator(mol, name, mission)

    fireworks = []
    links_dict = dict()

    # the task in the order of anion, neutral, cation
    cg_fwid, ng_fwid, ag_fwid = (None, None, None)
    cf_fwid, nf_fwid, af_fwid = (None, None, None)
    if len(mol) > 1:
        charge_shifts = (-1, 0, 1)
        fw_ids = range(0, 3)
        fws = (fw_creator.geom_fw(cs, fwid)
               for cs, fwid in zip(charge_shifts, fw_ids))
        cg_fwid, ng_fwid, ag_fwid = fw_ids
        fireworks.extend(fws)

        fw_ids = range(3, 3+3)
        fws = (fw_creator.freq_fw(cs, fwid)
               for cs, fwid in zip(charge_shifts, fw_ids))
        cf_fwid, nf_fwid, af_fwid = fw_ids
        fireworks.extend(fws)
        links_dict.update({cg_fwid: cf_fwid, ng_fwid: nf_fwid,
                           ag_fwid: af_fwid})

    charge_shifts = (-1, 0, 1)
    fw_ids = range(6, 6+3)
    fws = (fw_creator.sp_fw(cs, fwid)
           for cs, fwid in zip(charge_shifts, fw_ids))
    csp_fwid, nsp_fwid, asp_fwid = fw_ids
    fireworks.extend(fws)
    if len(mol) > 1:
        links_dict.update({cf_fwid: csp_fwid, nf_fwid: nsp_fwid,
                           af_fwid: asp_fwid})
        links_dict.update({nsp_fwid: [cg_fwid, ag_fwid]})
    else:
        links_dict.update({nsp_fwid: [csp_fwid, asp_fwid]})


    return Workflow(fireworks, links_dict, name)