from fireworks import Workflow
from rubicon.workflows.multistep_ipea_wf import QChemFireWorkCreator

__author__ = 'xiaohuiqu'


def single_point_energy_fws(mol, name, mission, parameters, dupefinder=None, priority=1,
                            parent_fwid=None, additional_user_tags=None):
    large = False
    if len(mol) > 50:
        large = True
    solvent = parameters.get("solvent", "water")
    method = parameters.get("method", None)
    population_method = parameters.get("population_method", None)
    energy_method, geom_method = method.split("//")
    charge = mol.charge
    spin_multiplicity = mol.spin_multiplicity
    fw_creator = QChemFireWorkCreator(mol=mol, molname=name, mission=mission,
                                      dupefinder=dupefinder,
                                      priority=priority, large=large, additional_user_tags=additional_user_tags)
    fwid_base = 1
    if parent_fwid:
        if not (isinstance(parent_fwid, int) or isinstance(parent_fwid, list)):
            raise ValueError("Parent FireWork ID must be integer or list")
        parent_fwid = parent_fwid if isinstance(parent_fwid, list) \
            else [parent_fwid]
        fwid_base = max(parent_fwid) + 1
    fws = []
    links_dict = dict()
    geom_cal_fwid = None
    freq_db_fwid = None
    if len(mol) > 1:
        geom_cal_fwid, geom_db_fwid = fwid_base + 0, fwid_base + 1
        fw_geom = fw_creator.geom_fw(
            charge, spin_multiplicity, geom_cal_fwid, geom_db_fwid, priority, geom_method)
        fws.append(fw_geom)
        links_dict[geom_cal_fwid] = geom_db_fwid

        freq_cal_fwid, freq_db_fwid = fwid_base + 2, fwid_base + 3
        fw_freq = fw_creator.freq_fw(
            charge, spin_multiplicity, freq_cal_fwid, freq_db_fwid, priority, geom_method)
        fws.append(fw_freq)
        links_dict[geom_db_fwid] = freq_cal_fwid
        links_dict[freq_cal_fwid] = freq_db_fwid

    sp_cal_fwid, sp_db_fwid = fwid_base + 4, fwid_base + 5
    fw_sp = fw_creator.sp_fw(
        charge, spin_multiplicity, sp_cal_fwid, sp_db_fwid, solvent_method="sm12mk", solvent=solvent,
        priority=priority, method=energy_method, population_method=population_method
    )
    fws.append(fw_sp)
    links_dict[sp_cal_fwid] = sp_db_fwid
    if len(mol) > 1:
        links_dict[freq_db_fwid] = sp_cal_fwid
        links_dict[parent_fwid] = geom_cal_fwid
    else:
        links_dict[parent_fwid] = sp_cal_fwid

    return fws, links_dict


def single_point_energy_wf(mol, name, mission, parameters, dupefinder=None, priority=1,
                           parent_fwid=None, additional_user_tags=None):
    fws, links_dict = single_point_energy_fws(mol, name, mission, parameters, dupefinder, priority,
                                              parent_fwid, additional_user_tags)
    return Workflow(fws, links_dict, name)