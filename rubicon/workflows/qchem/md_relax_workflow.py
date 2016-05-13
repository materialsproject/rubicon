# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import copy
import itertools

from fireworks import Workflow
from six.moves import range
from six.moves import zip

from rubicon.utils.atomic_charge_mixed_basis_set_generator import \
    AtomicChargeMixedBasisSetGenerator
from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator

__author__ = 'xiaohuiqu'

fs_to_au = 1.0E-15 / 2.418884326505E-17


def md_relax_fws(mol, name, mission, qm_method, high_temperature=323,
                 low_temperature=273,
                 md_steps=500, time_step=1.0, md_runs=3, normal_basis="6-31G*",
                 diffuse_basis="6-31+G*",
                 charge_threshold=-0.5, dupefinder=None, priority=1,
                 parent_fwid=None, additional_user_tags=None):
    """

    :param mol:
    :param name:
    :param mission:
    :param qm_method:
    :param high_temperature:
    :param low_temperature:
    :param md_steps: simulation step per MD episode.
    :param time_step (int): in fs.
    :param md_runs: times to run MD simulation at every temperature
    :param dupefinder:
    :param priority:
    :param parent_fwid:
    :param additional_user_tags:
    :return:
    """
    time_step_au = int(time_step * fs_to_au)
    mixed_basis_generator = AtomicChargeMixedBasisSetGenerator(
        charge_threshold, normal_basis, diffuse_basis)
    charge = mol.charge
    spin_multiplicity = mol.spin_multiplicity
    large = len(mol) > 50
    fw_creator = QChemFireWorkCreator(mol=mol, molname=name, mission=mission,
                                      dupefinder=dupefinder,
                                      priority=priority,
                                      additional_user_tags=additional_user_tags,
                                      large=large)
    fwid_base = 1
    if parent_fwid:
        if not (isinstance(parent_fwid, int) or isinstance(parent_fwid, list)):
            raise ValueError("Parent FireWork ID must be integer or list")
        parent_fwid = parent_fwid if isinstance(parent_fwid, list) \
            else [parent_fwid]
        fwid_base = max(parent_fwid) + 1
    fws = []
    links_dict = dict()

    sp1_cal_fwid, sp1_db_fwid = fwid_base + 0, fwid_base + 1
    fw_sp1 = fw_creator.vacuum_only_sp_fw(charge, spin_multiplicity,
                                          sp1_cal_fwid, sp1_db_fwid, priority,
                                          qm_method,
                                          population_method="nbo",
                                          mixed_basis_generator=mixed_basis_generator)
    fws.extend(fw_sp1)
    links_dict[sp1_cal_fwid] = sp1_db_fwid
    if parent_fwid:
        for pfw_id in parent_fwid:
            links_dict[pfw_id] = sp1_cal_fwid

    geom1_cal_fwid, geom1_db_fwid = fwid_base + 2, fwid_base + 3
    fw_geom1 = fw_creator.geom_fw(
        charge, spin_multiplicity, geom1_cal_fwid, geom1_db_fwid, priority,
        qm_method, task_type_prefix="pre-md")
    fws.extend(fw_geom1)
    links_dict[geom1_cal_fwid] = geom1_db_fwid
    links_dict[sp1_db_fwid] = geom1_cal_fwid

    sp2_cal_fwid, sp2_db_fwid = fwid_base + 4, fwid_base + 5
    fw_sp2 = fw_creator.vacuum_only_sp_fw(charge, spin_multiplicity,
                                          sp2_cal_fwid, sp2_db_fwid, priority,
                                          qm_method,
                                          population_method="nbo",
                                          mixed_basis_generator=copy.deepcopy(
                                              mixed_basis_generator))
    fws.extend(fw_sp2)
    links_dict[sp2_cal_fwid] = sp2_db_fwid
    links_dict[geom1_db_fwid] = sp2_cal_fwid

    temperatures = list(itertools.chain(
        [int(high_temperature), int(low_temperature)] * md_runs))
    md_fw_ids = list(zip(
        *[iter(
            list(range(fwid_base + 6, fwid_base + 6 + md_runs * 2 * 2)))] * 2))
    md_fws = []
    for (md_cal_fwid, md_db_fwid), temperature in zip(md_fw_ids, temperatures):
        fw_md = fw_creator.aimd_fw(charge, spin_multiplicity, md_cal_fwid,
                                   md_db_fwid, md_steps, time_step_au,
                                   temperature, priority, qm_method)

        md_fws.extend(fw_md)
        fws.extend(fw_md)
        links_dict[md_cal_fwid] = md_db_fwid
        links_dict[md_cal_fwid - 1] = md_cal_fwid
    links_dict[sp2_db_fwid] = md_fw_ids[0][0]
    last_md_fwid = md_fw_ids[-1][-1]

    geom2_cal_fwid, geom2_db_fwid = last_md_fwid + 1, last_md_fwid + 2
    fw_geom2 = fw_creator.geom_fw(
        charge, spin_multiplicity, geom2_cal_fwid, geom2_db_fwid, priority,
        qm_method, task_type_prefix="post-md")
    fws.extend(fw_geom2)
    links_dict[geom2_cal_fwid] = geom2_db_fwid
    links_dict[last_md_fwid] = geom2_cal_fwid

    return fws, links_dict


def md_relax_wf(mol, name, **kwargs):
    fws, links_dict = md_relax_fws(mol, name, **kwargs)
    return Workflow(fws, links_dict, name)
