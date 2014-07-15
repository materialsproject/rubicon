from fireworks import Workflow

from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator


__author__ = 'xiaohuiqu'


def single_point_energy_fws(mol, name, mission, solvent, solvent_method, qm_method, pop_method, dupefinder=None, priority=1,
                            parent_fwid=None, additional_user_tags=None):
    large = False
    if len(mol) > 50:
        large = True
    energy_method, sol_qm_method, geom_method = qm_method.split("//")
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
        fws.extend(fw_geom)
        links_dict[geom_cal_fwid] = geom_db_fwid

        freq_cal_fwid, freq_db_fwid = fwid_base + 2, fwid_base + 3
        fw_freq = fw_creator.freq_fw(
            charge, spin_multiplicity, freq_cal_fwid, freq_db_fwid, priority, geom_method)
        fws.extend(fw_freq)
        links_dict[geom_db_fwid] = freq_cal_fwid
        links_dict[freq_cal_fwid] = freq_db_fwid

    sol_cal_fwid, sol_db_fwid = fwid_base + 4, fwid_base + 5
    fw_sol = fw_creator.sp_fw(
        charge, spin_multiplicity, sol_cal_fwid, sol_db_fwid, solvent_method=solvent_method, solvent=solvent,
        priority=priority, qm_method=sol_qm_method, population_method=pop_method, task_type_name="solvation energy")
    fws.extend(fw_sol)
    links_dict[sol_cal_fwid] = sol_db_fwid

    vac_sp_cal_fwid, vac_sp_db_fwid = fwid_base + 6, fwid_base + 7
    fw_vac_sp = fw_creator.vacuum_only_sp_fw(
        charge, spin_multiplicity, vac_sp_cal_fwid, vac_sp_db_fwid, priority=priority, qm_method=energy_method)
    fws.extend(fw_vac_sp)
    links_dict[vac_sp_cal_fwid] = vac_sp_db_fwid
    links_dict[sol_db_fwid] = vac_sp_cal_fwid

    if len(mol) > 1:
        links_dict[freq_db_fwid] = sol_cal_fwid
        for pfw_id in parent_fwid:
            links_dict[pfw_id] = geom_cal_fwid
    else:
        for pfw_id in parent_fwid:
            links_dict[pfw_id] = sol_cal_fwid
    return fws, links_dict


def single_point_energy_wf(mol, name, **kwargs):
    fws, links_dict = single_point_energy_fws(mol, name, **kwargs)
    return Workflow(fws, links_dict, name)