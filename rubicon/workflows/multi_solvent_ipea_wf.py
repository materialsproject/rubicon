# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import itertools

from fireworks import Workflow
from six.moves import range
from six.moves import zip

from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator


def multi_solvent_ipea_fws(mol, name, mission, solvents, solvent_method,
                           use_vdW_surface, ref_charge,
                           spin_multiplicities=(2, 1, 2),
                           dupefinder=None, priority=1, parent_fwid=None,
                           additional_user_tags=None, qm_method=None):
    large = False
    if len(mol) > 50:
        large = True
    energy_method, geom_method = qm_method.split("//") if qm_method else (
        None, None)
    fw_creator = QChemFireWorkCreator(mol=mol, molname=name, mission=mission,
                                      dupefinder=dupefinder,
                                      priority=priority, large=large,
                                      additional_user_tags=additional_user_tags)
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
    charges = [ref_charge + i for i in (-1, 0, 1)]
    if len(mol) > 1:
        fw_ids = list(
            zip(*[iter(list(range(fwid_base + 0, fwid_base + 6)))] * 2))
        fws = (fw_creator.geom_fw(ch, spin, fwid_cal, fwid_db,
                                  priority=priority + 1, method=geom_method)
               for ch, spin, (fwid_cal, fwid_db)
               in zip(charges, spin_multiplicities, fw_ids))
        (cgi_cal, cgi_db), (ngi_cal, ngi_db), (agi_cal, agi_db) = fw_ids
        fireworks.extend(itertools.chain.from_iterable(fws))
        links_dict.update(dict(fw_ids))

        if not large:
            fw_ids = list(zip(
                *[iter(list(range(fwid_base + 6, fwid_base + 6 + 6)))] * 2))
            fws = (fw_creator.freq_fw(ch, spin, fwid_cal, fwid_db,
                                      priority=priority + 1,
                                      method=geom_method)
                   for ch, spin, (fwid_cal, fwid_db)
                   in zip(charges, spin_multiplicities, fw_ids))
            (cfi_cal, cfi_db), (nfi_cal, nfi_db), (afi_cal, afi_db) = fw_ids
            fireworks.extend(itertools.chain.from_iterable(fws))
            links_dict.update(dict(fw_ids))
            links_dict.update({cgi_db: cfi_cal,
                               ngi_db: nfi_cal,
                               agi_db: afi_cal})

    num_solvents = len(solvents)
    if num_solvents < 1:
        raise ValueError("You must provide at least one solvent")

    sp_fw_ids = []
    for sol_id, solvent in enumerate(solvents):
        fwid_start = fwid_base + 12 + (sol_id * 6)
        fwid_end = fwid_base + 12 + (sol_id * 6) + 6
        fw_ids = list(zip(*[iter(list(range(fwid_start, fwid_end)))] * 2))
        sp_fw_ids.append(fw_ids)
        fws = (fw_creator.sp_fw(ch, spin, fwid_cal, fwid_db,
                                solvent=solvent,
                                use_vdw_surface=use_vdW_surface,
                                solvent_method=solvent_method,
                                qm_method=energy_method)
               for ch, spin, (fwid_cal, fwid_db)
               in zip(charges, spin_multiplicities, fw_ids))
        links_dict.update(dict(fw_ids))
        fireworks.extend(itertools.chain.from_iterable(fws))
    cspi_cal_list = []
    nspi_cal_list = []
    aspi_cal_list = []
    for (cspi_cal, cspi_db), (nspi_cal, nspi_db), (aspi_cal, aspi_db) in \
            sp_fw_ids:
        cspi_cal_list.append(cspi_cal)
        nspi_cal_list.append(nspi_cal)
        aspi_cal_list.append(aspi_cal)
    if len(mol) > 1:
        if large:
            links_dict.update({cgi_db: cspi_cal_list, ngi_db: nspi_cal_list,
                               agi_db: aspi_cal_list})
            links_dict[ngi_db].extend([cgi_cal, agi_cal])
        else:
            links_dict.update({cfi_db: cspi_cal_list, nfi_db: nspi_cal_list,
                               afi_db: aspi_cal_list})
            links_dict[nfi_db].extend([cgi_cal, agi_cal])
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = ngi_cal
    else:
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = nspi_cal_list + aspi_cal_list + \
                                     cspi_cal_list
    return fireworks, links_dict


def mol_to_solvent_ipea_wf(mol, name, **kwargs):
    fireworks, links_dict = multi_solvent_ipea_fws(mol, name, **kwargs)
    return Workflow(fireworks, links_dict, name)
