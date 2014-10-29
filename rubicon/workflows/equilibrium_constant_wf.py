import json
import logging
import os
import sys

from pymatgen.matproj.snl import StructureNL
from pymongo import MongoClient

from rubicon.workflows.bsse_wf import counterpoise_correction_generation_fw, BSSEFragment
from rubicon.workflows.single_point_energy_wf import single_point_energy_fws


__author__ = 'xiaohuiqu'


def get_reactions_collection():
    db_dir = os.environ['DB_LOC']
    db_path = os.path.join(db_dir, 'tasks_db.json')
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('QChemDrone')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setLevel(getattr(logging, 'INFO'))
    logger.addHandler(sh)
    with open(db_path) as f:
        db_creds = json.load(f)
    conn = MongoClient(db_creds['host'], db_creds['port'],)
    db = conn[db_creds['db']]
    if db_creds['username']:
        db.authenticate(db_creds['username'], db_creds['password'])
    coll = db['reactions']
    return coll

def equilibrium_constant_fws(name, mission, solvent, solvent_method, use_vdW_surface, qm_method, reaction_id,
                             dupefinder=None, priority=1, parent_fwid=None, additional_user_tags=None,
                             depend_on_parent=False):
    energy_method, sol_qm_method, geom_method = qm_method.split("//")
    if '||' in energy_method:
        sp_qm_method, bsse_qm_method = energy_method.split("//")
        qm_method = "//".join([energy_method, sol_qm_method, geom_method])
    else:
        bsse_qm_method = energy_method
    coll = get_reactions_collection()
    reaction_doc = coll.find_one({"reaction_id": reaction_id})
    reactant_snls = [StructureNL.from_dict(s) for s in reaction_doc["reactant_snls"]]
    product_snls = [StructureNL.from_dict(s) for s in reaction_doc["product_snls"]]
    reactant_nicknames = reaction_doc['reactant_nicknames']
    product_nicknames = reaction_doc['product_nicknames']
    reactant_charges = reaction_doc['reactant_charges']
    product_charges = reaction_doc['product_charges']
    reactant_spin_multiplicities = reaction_doc['reactant_spin_multiplicities']
    product_spin_multiplicities = reaction_doc['product_spin_multiplicities']
    reactant_fragments = [BSSEFragment.from_dict(frag) for frag in reaction_doc['reactant_fragments']]
    product_fragments = [BSSEFragment.from_dict(frag) for frag in reaction_doc['product_fragments']]

    fwid_base = 1
    if parent_fwid:
        if not (isinstance(parent_fwid, int) or isinstance(parent_fwid, list)):
            raise ValueError("Parent FireWork ID must be integer or list")
        parent_fwid = parent_fwid if isinstance(parent_fwid, list) \
            else [parent_fwid]
        fwid_base = max(parent_fwid) + 1

    current_fwid = fwid_base
    fws = []
    links_dict = dict()

    for snl, nick_name, charge, spin, fragments in reactant_snls + product_snls, \
                                        reactant_nicknames + product_nicknames, \
                                        reactant_charges + product_charges, \
                                        reactant_spin_multiplicities+ product_spin_multiplicities, \
                                        reactant_fragments + product_fragments:
        mol = snl.structure
        mol.set_charge_and_spin(charge, spin)
        sp_fws, sp_links_dict = single_point_energy_fws(
            mol, name=nick_name, mission=mission, solvent=solvent, solvent_method=solvent_method,
            use_vdW_surface=use_vdW_surface, qm_method=qm_method, pop_method=None, dupefinder=dupefinder,
            priority=priority, parent_fwid=current_fwid, additional_user_tags=additional_user_tags,
            depend_on_parent_fw=False)
        fws.extend(sp_fws)
        for k, v2 in sp_links_dict:
            v1 = links_dict.get(k, [])
            links_dict[k] = v1 + v2
        current_fwid = max([fw.fw_id for fw in sp_fws]) + 1

        sp_fws, sp_links_dict = counterpoise_correction_generation_fw(
            molname=nick_name, charge=charge, spin_multiplicity=spin, qm_method=bsse_qm_method, fragments=fragments,
            mission=mission, priority=priority, parent_fwid=current_fwid, additional_user_tags=additional_user_tags)
        fws.extend(sp_fws)
        for k, v2 in sp_links_dict:
            v1 = links_dict.get(k, [])
            links_dict[k] = v1 + v2
        current_fwid = max([fw.fw_id for fw in sp_fws]) + 1

    if depend_on_parent:
        all_fwids = [fw.fw_id for fw in fws]
        for p_fwid in parent_fwid:
            links_dict[p_fwid] = all_fwids

    return fws, links_dict
