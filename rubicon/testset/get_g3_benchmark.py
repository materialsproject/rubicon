from collections import defaultdict
import copy
import os
import json
from pymongo.mongo_client import MongoClient
import csv

__author__ = 'xiaohuiqu'


KCAL_TO_EV = 0.0433634

def get_g3_bench_collection():
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
    return collection



def get_calcualtion_result(mission_tag, bench_key_name, bench_dict, no_expt_bench_dict, db_collection):
    result_cursor = db_collection.find({"user_tags.mission": mission_tag},
                             fields=['pretty_formula', 'IE', 'EA', 'charge',
                                     'user_tags.fw_name'])
    calc_result = list(result_cursor)

    with open("gauname2refname.json") as f:
        gau2web_name_map = json.load(f)

    for m in calc_result:
        fw_name = m['user_tags']['fw_name']
        if fw_name in gau2web_name_map:
            web_name = gau2web_name_map[fw_name]
            d = bench_dict[web_name]
        else:
            if not fw_name in no_expt_bench_dict:
                no_expt_bench_dict[fw_name] = {}
            d = no_expt_bench_dict[fw_name]

        if 'IE' in m:
            if 'IP' in d:
                d['IP'][bench_key_name] = m['IE']
            else:
                d['IP'] = {bench_key_name: m['IE']}
        if 'EA' in m:
            if 'EA' in d:
                d['EA'][bench_key_name] = m['EA']
            else:
                d['EA'] = {bench_key_name: m['EA']}


def write_dict(bench_dict, writer):
    mols = sorted(bench_dict.keys())
    if 'unit' in mols:
        mols.remove('unit')
    for m in mols:
        row_dict = defaultdict(str)
        row_dict["Molecule"] = m
        for property, source_names in bench_dict[m].items():
            for source in source_names:
                col_name = '{}-{}'.format(property, source)
                row_dict[col_name] = bench_dict[m][property][source]
        writer.writerow(row_dict)


def write_csv(bench_dict, no_expt_bench_dict):
    source_names = set()
    property_names = set()
    for m in bench_dict.items():
        if m[0] != 'unit':
            for i, v in m[1].items():
                property_names.add(i)
                source_names = source_names | set(v.keys())
    headings = ['Molecule'] + ['{}-{}'.format(i, j)
                               for i in sorted(list(property_names))
                               for j in sorted(list(source_names))]
    with open("G3_bench.csv", 'w') as f:
        writer = csv.DictWriter(f, fieldnames=headings)
        writer.writeheader()
        write_dict(bench_dict, writer)
        write_dict(no_expt_bench_dict, writer)


if __name__ == '__main__':

    with open('G3_ref.json') as f:
        ref_data = json.load(f)
    bench = copy.deepcopy(ref_data)
    for m in bench.items():
        if m[0] != 'unit':
            for i, v in m[1].items():
                v['Expt'] *= KCAL_TO_EV
                v['G3'] *= KCAL_TO_EV
    no_expt_bench = {}

    collection = get_g3_bench_collection()

    get_calcualtion_result("G2-97 Test Set Benchmark (Shyue Scheme)", "Shyue", bench, no_expt_bench, collection)
    get_calcualtion_result("G2-97 Test Set Benchmark (Larry Scheme)", "Larry", bench, no_expt_bench, collection)

    all_bench = {}
    all_bench.update(bench)
    all_bench.update(no_expt_bench)
    with open("G3_bench.json", 'w') as f:
        json.dump(all_bench, f, indent=4, sort_keys=True)

    write_csv(bench, no_expt_bench)
