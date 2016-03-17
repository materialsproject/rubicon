import copy
import csv
import json
import os
from collections import defaultdict

from pymongo.mongo_client import MongoClient

from pymatgen.core.units import Energy

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
    conn = MongoClient(host=host, port=port,
                       connect=False)
    db = conn[database_name]
    if user:
        db.authenticate(user, password)
    collection = db[collection_name]
    return collection


def get_1st_round_calcualtion_result(mission_tag, bench_key_name, bench_dict,
                                     no_expt_bench_dict,
                                     db_collection):
    result_cursor = db_collection.find(
        filter={"user_tags.mission": mission_tag},
        projection=['pretty_formula', 'IE', 'EA', 'charge',
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


def get_2nd_round_calcualtion_result(mission_tag, bench_key_name, bench_dict,
                                     no_expt_bench_dict,
                                     db_collection):
    result_cursor = db_collection.find(
        filter={"user_tags.mission": mission_tag},
        projection=['user_tags', 'IE', 'EA', 'calculations',
                    'inchi'])
    calc_result = list(result_cursor)
    for m in calc_result:
        ref_name = m['user_tags']['fw_name']
        if ref_name in bench_dict:
            d = bench_dict[ref_name]
        else:
            if not ref_name in no_expt_bench_dict:
                no_expt_bench_dict[ref_name] = {}
            d = no_expt_bench_dict[ref_name]

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
        if 'scf' in m['calculations']:
            if 'neutral_energy' in d:
                d['neutral_energy'][bench_key_name] = float(
                    Energy(m['calculations']['scf']['energies'][-1], 'eV').to(
                        'Ha'))
            else:
                d['neutral_energy'] = {bench_key_name:
                                           float(Energy(
                                               m['calculations']['scf'][
                                                   'energies'][-1],
                                               'eV').to('Ha'))}
        if 'scf_IE' in m['calculations']:
            if 'cation_energy' in d:
                d['cation_energy'][bench_key_name] = float(
                    Energy(m['calculations']['scf_IE']['energies'][-1],
                           'eV').to('Ha'))
            else:
                d['cation_energy'] = {bench_key_name:
                                          float(Energy(
                                              m['calculations']['scf_IE'][
                                                  'energies'][-1],
                                              'eV').to('Ha'))}
        if 'scf_EA' in m['calculations']:
            if 'anion_energy' in d:
                d['anion_energy'][bench_key_name] = float(
                    Energy(m['calculations']['scf_EA']['energies'][-1],
                           'eV').to('Ha'))
            else:
                d['anion_energy'] = {bench_key_name:
                                         float(Energy(
                                             m['calculations']['scf_EA'][
                                                 'energies'][-1],
                                             'eV').to('Ha'))}


def write_dict(bench_dict, writer):
    mols = sorted(bench_dict.keys())
    if 'unit' in mols:
        mols.remove('unit')
    for m in mols:
        row_dict = defaultdict(str)
        row_dict["Molecule"] = m
        for property, source_names in bench_dict[m].items():
            if property != 'inchi':
                for source in source_names:
                    col_name = '{}-{}'.format(property, source)
                    row_dict[col_name] = bench_dict[m][property][source]
        writer.writerow(row_dict)


def write_csv(bench_dict, no_expt_bench_dict):
    headings = set(['Molecule'])
    for m in bench_dict.items():
        if m[0] != 'unit':
            for property, v in m[1].items():
                if property != 'inchi':
                    for source in v.keys():
                        headings.add('{}-{}'.format(property, source))
    with open("G3_bench.csv", 'w') as f:
        writer = csv.DictWriter(f, fieldnames=sorted(list(headings)))
        writer.writeheader()
        write_dict(bench_dict, writer)
        write_dict(no_expt_bench_dict, writer)


def get_1st_round_benchmark():
    with open('G3_ref_with_inchi.json') as f:
        ref_data = json.load(f)
    bench = copy.deepcopy(ref_data)
    for m in bench.items():
        if m[0] != 'unit':
            for i, v in m[1].items():
                v['Expt'] *= KCAL_TO_EV
                v['G3'] *= KCAL_TO_EV
    no_expt_bench = {}
    collection = get_g3_bench_collection()
    get_1st_round_calcualtion_result("G2-97 Test Set Benchmark (Shyue Scheme)",
                                     "Shyue", bench,
                                     no_expt_bench,
                                     collection)
    get_1st_round_calcualtion_result("G2-97 Test Set Benchmark (Larry Scheme)",
                                     "Larry", bench,
                                     no_expt_bench,
                                     collection)
    all_bench = {}
    all_bench.update(bench)
    all_bench.update(no_expt_bench)
    with open("G3_bench.json", 'w') as f:
        json.dump(all_bench, f, indent=4, sort_keys=True)
    write_csv(bench, no_expt_bench)


def get_2nd_round_mp2_geom_benchmark():
    with open('G3_ref_with_inchi.json') as f:
        ref_data = json.load(f)
    bench = copy.deepcopy(ref_data)
    bench.pop('unit')
    for n, m in bench.items():
        for p in ['IP', 'EA']:
            if p in m:
                v = m[p]
                v['Expt'] *= KCAL_TO_EV
                v['G3'] *= KCAL_TO_EV
    no_expt_bench = {}
    collection = get_g3_bench_collection()
    get_2nd_round_calcualtion_result("G2-97 MP2 Geom Benchmark (Larry Scheme)",
                                     "Larry", bench,
                                     no_expt_bench,
                                     collection)
    get_2nd_round_calcualtion_result("G2-97 MP2 Geom Benchmark (Shyue Scheme)",
                                     "Shyue", bench,
                                     no_expt_bench,
                                     collection)

    all_bench = {}
    all_bench.update(bench)
    all_bench.update(no_expt_bench)
    with open("G3_bench.json", 'w') as f:
        json.dump(all_bench, f, indent=4, sort_keys=True)
    write_csv(bench, no_expt_bench)


def get_3th_round_dft_methods_benchmark():
    with open('G3_ref_with_inchi.json') as f:
        ref_data = json.load(f)
    bench = copy.deepcopy(ref_data)
    bench.pop('unit')
    for n, m in bench.items():
        for p in ['IP', 'EA']:
            if p in m:
                v = m[p]
                v['Expt'] *= KCAL_TO_EV
                v['G3'] *= KCAL_TO_EV
    no_expt_bench = {}
    collection = get_g3_bench_collection()

    for xc in ['xctpssh', 'pbe0', 'm06', 'm06-2x', 'pw6b95', 'becke97-d',
               'b3lyp']:
        mission = "DFT Test (" + xc + ")"
        get_2nd_round_calcualtion_result(mission, xc, bench, no_expt_bench,
                                         collection)

    all_bench = {}
    all_bench.update(bench)
    all_bench.update(no_expt_bench)
    with open("G3_bench.json", 'w') as f:
        json.dump(all_bench, f, indent=4, sort_keys=True)
    write_csv(bench, no_expt_bench)


if __name__ == '__main__':
    get_3th_round_dft_methods_benchmark()
