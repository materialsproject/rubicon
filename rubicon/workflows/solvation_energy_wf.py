# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import itertools

from fireworks.core.firework import Workflow
from six.moves import range
from six.moves import zip

from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator

__author__ = 'xiaohuiqu'


def solvation_energy_fws(mol, name, mission, solvents, solvent_method,
                         use_vdW_surface, dupefinder=None, priority=1,
                         parent_fwid=None,
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
    geom_cal_id = None
    freq_db_id = None
    geom_db_id = None

    if len(mol) > 1:
        geom_cal_id, geom_db_id = fwid_base + 0, fwid_base + 1
        fws_geom = fw_creator.geom_fw(
            charge=0, spin_multiplicity=1, fw_id_cal=geom_cal_id,
            fw_id_db=geom_db_id, method=geom_method)
        fireworks.extend(fws_geom)
        links_dict[geom_cal_id] = geom_db_id

        if not large:
            freq_cal_id, freq_db_id = fwid_base + 2, fwid_base + 3
            fws_freq = fw_creator.freq_fw(
                charge=0, spin_multiplicity=1, fw_id_cal=freq_cal_id,
                fw_id_db=freq_db_id, method=geom_method)
            fireworks.extend(fws_freq)
            links_dict[freq_cal_id] = freq_db_id
            links_dict[geom_db_id] = freq_cal_id

    num_solvents = len(solvents)
    if num_solvents < 1:
        raise ValueError("You must provide at least one solvent")

    sp_fw_ids = list(zip(*[iter(list(range(fwid_base + 4,
                                           fwid_base + 4 + num_solvents * 2)))] * 2))
    sp_fws = (fw_creator.sp_fw(
        charge=0, spin_multiplicity=1, fw_id_cal=fwid_cal,
        fw_id_db=fwid_db, solvent_method=solvent_method,
        use_vdw_surface=use_vdW_surface, solvent=solvent,
        qm_method=energy_method)
              for (fwid_cal, fwid_db), solvent in zip(sp_fw_ids, solvents))
    sp_cal_ids, sp_db_ids = list(zip(*sp_fw_ids))
    links_dict.update((dict(sp_fw_ids)))
    fireworks.extend(itertools.chain.from_iterable(sp_fws))
    if len(mol) > 1:
        if large:
            links_dict[geom_db_id] = list(sp_cal_ids)
        else:
            links_dict[freq_db_id] = list(sp_cal_ids)
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = geom_cal_id
    else:
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = list(sp_cal_ids)
    return fireworks, links_dict


def mol_to_solvation_energy_wf(mol, name, **kwargs):
    fireworks, links_dict = solvation_energy_fws(
        mol, name, **kwargs)
    return Workflow(fireworks, links_dict, name)
