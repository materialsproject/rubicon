import json
import os
import bson
from pymatgen import Composition
from pymongo import MongoClient

__author__ = 'xiaohuiqu'


def get_db_collection():
    db_dir = os.environ['DB_LOC']
    db_path = os.path.join(db_dir, 'molecules_db.json')
    with open(db_path) as f:
        db_creds = json.load(f)
    host = db_creds['host']
    port = db_creds['port']
    user = db_creds['readonly_user']
    password = db_creds['readonly_password']
    database_name = db_creds['database']
    collection_name = db_creds['collection']
    conn = MongoClient(host=host, port=port)
    db = conn[database_name]
    if user:
        db.authenticate(user, password)
    collection = db[collection_name]
    return collection, db

def rename_task_id(collection):
    collection.update({}, {'$rename': {'task_id': 'task_id_deprecated'}},
                      multi=True)

def add_new_task_id(db):
    for doc in db.molecules.find():
        doc["task_id"] = "mol-" + str(doc["task_id_deprecated"])
        db.molecules.save(doc)

def add_field_reduced_cell_formula_abc(db):
    for doc in db.molecules.find():
        doc["reduced_cell_formula_abc"] = Composition(doc["pretty_formula"])\
            .alphabetical_formula
        db.molecules.save(doc)

def mark_g3_success(collection):
    g3_missions = ['G2-97 Test Set Benchmark (Shyue Scheme)',
                   'G2-97 Test Set Benchmark (Larry Scheme)']
    for mission in g3_missions:
        collection.update({"user_tags.mission": mission},
                          {"$set": {"state": "successful"}},
                          multi=True)

def mark_imaginary_freq_error(collection):
    mol_names = '''anthrachinon_wfs_16_trichloromethyl
quinoxaline_wfs_7_cyano
quinoxaline_wfs_7_ethynyl
viologen_wfs_6_hydroxyl
anthrachinon_wfs_15_carboxyl
anthrachinon_wfs_15_fluoro
anthrachinon_wfs_15_methyl
anthrachinon_wfs_15_nitro
anthrachinon_wfs_16_cyano
anthrachinon_wfs_16_methyl
quinoxaline_wfs_13_benzene
quinoxaline_wfs_13_methyl
quinoxaline_wfs_13_trichloromethyl
quinoxaline_wfs_6_benzene
quinoxaline_wfs_7_benzene
quinoxaline_wfs_7_carboxyl
quinoxaline_wfs_7_methyl
quinoxaline_wfs_7_trichloromethyl
thiane_wfs_2_amine
thiane_wfs_5_cyano
thiane_wfs_5_ethynyl
thiane_wfs_5_fluoro
thiophene_wfs_5_methyl
thiophene_wfs_5_trichloromethyl
thiophene_wfs_6_benzene
thiophene_wfs_6_ethanamide
viologen_wfs_6_benzene
viologen_wfs_6_carboxyl
viologen_wfs_6_cyano
viologen_wfs_6_fluoro
viologen_wfs_6_methyl
viologen_wfs_6_trichloromethyl
viologen_wfs_6_vinyl
viologen_wfs_7_carboxyl
viologen_wfs_7_fluoro'''
    for mol in mol_names.split('\n'):
        collection.update({"user_tags.fw_name": mol},
                          {"$set": {"state": "error"}})

def mark_geom_failed_error(collection):
    mol_names = '''anthrachinon_wfs_15_benzene
anthrachinon_wfs_15_dimethylamine
anthrachinon_wfs_15_ethanamide
anthrachinon_wfs_15_ethyl
anthrachinon_wfs_15_methylamine
anthrachinon_wfs_16_amine
anthrachinon_wfs_16_benzene
anthrachinon_wfs_16_carboxyl
anthrachinon_wfs_16_dimethylamine
anthrachinon_wfs_16_ethanamide
anthrachinon_wfs_16_ethyl
anthrachinon_wfs_16_methoxyl
anthrachinon_wfs_16_methylamine
anthrachinon_wfs_16_vinyl
quinoxaline_wfs_13_dimethylamine
quinoxaline_wfs_13_ethyl
quinoxaline_wfs_13_methoxyl
quinoxaline_wfs_13_methylamine
quinoxaline_wfs_6_dimethylamine
quinoxaline_wfs_6_ethanamide
quinoxaline_wfs_6_ethyl
quinoxaline_wfs_6_methoxyl
quinoxaline_wfs_6_methylamine
quinoxaline_wfs_7_dimethylamine
quinoxaline_wfs_7_ethanamide
quinoxaline_wfs_7_ethyl
quinoxaline_wfs_7_methoxyl
quinoxaline_wfs_7_methylamine
thiane_wfs_10_benzene
thiane_wfs_10_carboxyl
thiane_wfs_10_cyano
thiane_wfs_10_dimethylamine
thiane_wfs_10_ethanamide
thiane_wfs_10_ethyl
thiane_wfs_10_ethynyl
thiane_wfs_10_fluoro
thiane_wfs_10_hydroxyl
thiane_wfs_10_methoxyl
thiane_wfs_10_nitro
thiane_wfs_10_trichloromethyl
thiane_wfs_10_vinyl
thiane_wfs_2_benzene
thiane_wfs_2_carboxyl
thiane_wfs_2_cyano
thiane_wfs_2_dimethylamine
thiane_wfs_2_ethanamide
thiane_wfs_2_ethynyl
thiane_wfs_2_fluoro
thiane_wfs_2_hydroxyl
thiane_wfs_2_methoxyl
thiane_wfs_2_methyl
thiane_wfs_2_methylamine
thiane_wfs_2_nitro
thiane_wfs_2_trichloromethyl
thiane_wfs_2_vinyl
thiane_wfs_5_amine
thiane_wfs_5_benzene
thiane_wfs_5_carboxyl
thiane_wfs_5_dimethylamine
thiane_wfs_5_ethanamide
thiane_wfs_5_ethyl
thiane_wfs_5_hydroxyl
thiane_wfs_5_methoxyl
thiane_wfs_5_methyl
thiane_wfs_5_methylamine
thiane_wfs_5_nitro
thiane_wfs_5_trichloromethyl
thiane_wfs_5_vinyl
thiophene_wfs_5_benzene
thiophene_wfs_5_dimethylamine
thiophene_wfs_5_ethanamide
thiophene_wfs_5_ethyl
thiophene_wfs_5_methoxyl
thiophene_wfs_6_dimethylamine
thiophene_wfs_6_ethyl
thiophene_wfs_6_methoxyl
thiophene_wfs_6_methylamine
viologen_wfs_6_amine
viologen_wfs_6_dimethylamine
viologen_wfs_6_ethanamide
viologen_wfs_6_ethyl
viologen_wfs_6_ethynyl
viologen_wfs_6_methoxyl
viologen_wfs_6_methylamine
viologen_wfs_6_nitro
viologen_wfs_7_amine
viologen_wfs_7_cyano
viologen_wfs_7_dimethylamine
viologen_wfs_7_ethanamide
viologen_wfs_7_ethyl
viologen_wfs_7_methoxyl
viologen_wfs_7_methyl
viologen_wfs_7_methylamine
viologen_wfs_7_nitro
viologen_wfs_7_trichloromethyl'''
    for mol in mol_names.split('\n'):
        collection.update({"user_tags.fw_name": mol},
                          {"$set": {"state": "error"}})

if __name__ == '__main__':
    collection, db = get_db_collection()
    #rename_task_id(collection)
    #add_new_task_id(db)
    #mark_g3_success(collection)
    #mark_imaginary_freq_error(collection)
    #mark_geom_failed_error(collection)
    add_field_reduced_cell_formula_abc(db)