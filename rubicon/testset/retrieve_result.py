from collections import defaultdict
import csv
import json
import os
from pymongo import MongoClient
import pymongo

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
    conn = MongoClient(host=host, port=port, connect=True)
    db = conn[database_name]
    if user:
        db.authenticate(user, password)
    collection = db[collection_name]
    return collection

def get_quinoxaline_result(db_collection, mission_tag):
    result_cursor = db_collection.find(filter={"user_tags.mission": mission_tag},
                                       projection=['user_tags', 'IE', 'EA',
                                                   'inchi'])
    result_cursor.sort("user_tags.fw_name", pymongo.ASCENDING)
    calc_result = list(result_cursor)
    return calc_result

def write_csv(result, filename):
    headings =["Molecule", "IP", "EA", "Inchi"]
    with open(filename, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=headings)
        writer.writeheader()
        for row in result:
            row_dict = defaultdict(str)
            row_dict["Molecule"] = row['user_tags']['fw_name']
            row_dict["Inchi"] = row['inchi']
            row_dict['IP'] = row.get('IE')
            row_dict['EA'] = row.get('EA')
            writer.writerow(row_dict)

def write_path_map(db_collection, mission_tag):
    result_cursor = db_collection.find(filter={"user_tags.mission": mission_tag},
                                       projection=['user_tags', 'path'])
    result_cursor.sort("user_tags.fw_name", pymongo.ASCENDING)
    calc_result = list(result_cursor)
    name2path = {row['user_tags']['fw_name']: row['path'] for row in
                 calc_result}
    with open('name2path.json', 'w') as f:
        json.dump(name2path, f, indent=4)

if __name__ == '__main__':
    collection = get_db_collection()
    mission = 'Quinoxaline Derivatives'
    calc_result = get_quinoxaline_result(collection, mission)
    write_csv(calc_result, "Quinoxaline.csv")
    write_path_map(collection, mission)