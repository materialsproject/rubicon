import json
import logging
import os
from pymongo import MongoClient
import sys


def get_molecules_collection(db_dir):
    db_path = os.path.join(db_dir, 'molecules_db.json')
    with open(db_path) as f:
        db_creds = json.load(f)
    conn = MongoClient(db_creds['host'], db_creds['port'])
    db = conn[db_creds['database']]
    if db_creds['admin_user']:
        db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
    ipea_coll = db['molecules']
    return ipea_coll

def transform_molecule_doc(mol1):
    mol2 = dict()
    mol2["inchi_root"] = mol1["inchi_root"]
    mol2["can"] = mol1['vacuum_properties']["can"]['neutral']
    mol2["smiles"] = mol1['vacuum_properties']["smiles"]['neutral']
    mol2["elements"] = mol1["elements"]
    mol2["nelements"] = mol1["nelements"]
    mol2["user_tags"] = mol1["user_tags"]
    mol2["run_tags"] = mol1["run_tags"]
    if "base_molecule" in mol1:
        mol2["base_molecule"] = mol1["base_molecule"]
    if "functional_groups" in mol1:
        mol2["functional_groups"] = mol1["functional_groups"]
    mol2["reduced_cell_formula_abc"] = mol1["reduced_cell_formula_abc"]
    if "water" in mol1["solvated_properties"]:
        mol2["implicit_solvent"] = mol1["solvated_properties"]["water"]["implicit_solvent"]
        mol2["electrode_potentials"] = mol1["solvated_properties"]["water"]["electrode_potentials"]
    else:
        mol2["implicit_solvent"] = {}
    mol2["pretty_formula"] = mol1["pretty_formula"]
    mol2["formula"] = mol1["formula"]
    mol2["pointgroup"] = mol1["pointgroup"]

    mol2["task_id"] = mol1["vacuum_properties"]["task_id"]["neutral"]
    mol2["task_id_deprecated"] = mol1["vacuum_properties"]["task_id_deprecated"]["neutral"]
    mol2["snlgroup_id_final"] = mol1["vacuum_properties"]["snlgroup_id_final"]["neutral"]
    mol2["charge"] = mol1["charge"]
    mol2["spin_multiplicity"] = mol1["vacuum_properties"]["spin_multiplicity"]["neutral"]
    mol2["snl_final"] = mol1["vacuum_properties"]["snl_final"]["neutral"]
    mol2["molecule"] = mol1["vacuum_properties"]["molecule"]["neutral"]
    mol2["xyz"] = mol1["vacuum_properties"]["xyz"]["neutral"]
    mol2["inchi"] = mol1["vacuum_properties"]["inchi"]["neutral"]

    solname = None
    for prefer_k in ["water_ief-pcm_at_surface0_00",
                     "water_ief-pcm_at_surface",
                     "water_ief-pcm",
                     "ief-pcm_at_surface0_00",
                     "ief-pcm_at_surface",
                     "ief-pcm"]:
        for k in mol1["solvated_properties"].keys():
            if prefer_k in k:
                solname = k
                break
        if solname:
            break
    else:
        for k in mol1["solvated_properties"].keys():
            solname = k
            break
    if solname:
        if "solvated_properties" in mol1:
            if "IP" in mol1["solvated_properties"][solname]:
                mol2["IE"] = mol1["solvated_properties"][solname]["IP"]
            if "EA" in mol1["solvated_properties"][solname]:
                mol2["EA"] = mol1["solvated_properties"][solname]["EA"]

    mol2["svg"] = mol1["svg"]
    return mol2

def copy_collections():
    source_cred_dir = os.environ['DB_LOC_SRC']
    dest_cred_dir = os.environ['DB_LOC_DEST']
    coll_src = get_molecules_collection(source_cred_dir)
    coll_dest = get_molecules_collection(dest_cred_dir)
    molecules = coll_src.find()
    for mol_db in molecules:
        mol_web = transform_molecule_doc(mol_db)
        if "molname" in mol_web["user_tags"]:
            molname = mol_web["user_tags"]["molname"]
        else:
            molname = mol_web["inchi"]
        mol_doc = coll_dest.find_one({"inchi_root": mol_web["inchi_root"],
                                      "charge": mol_web["charge"]})
        if mol_doc:
            logging.info("Updating molecule \"{}\"".format(molname))
            coll_dest.update({"inchi_root": mol_web["inchi_root"]}, mol_web,
                             upsert=True)
        else:
            logging.info("INSERT molecule \"{}\"".format(molname))
            coll_dest.insert(mol_web)

def main():
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('CopyMolecules')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setLevel(getattr(logging, 'INFO'))
    logger.addHandler(sh)
    copy_collections()

if __name__ == '__main__':
    main()
