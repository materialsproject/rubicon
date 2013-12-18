import json
import os
from pymongo import MongoClient
from rubicon.submission.submission_mongo import SubmissionMongoAdapter

__author__ = 'xiaohuiqu'


def get_calculation_db():
    db_dir = os.environ['DB_LOC']
    db_path = os.path.join(db_dir, 'qchem_calc_db.json')
    with open(db_path) as f:
        db_creds = json.load(f)
    conn = MongoClient(db_creds['host'], db_creds['port'])
    db = conn[db_creds['database']]
    if db_creds['admin_user']:
        db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
    calc_coll = db[db_creds['collection']]
    ipea_coll = db[db_creds['ipea_collection']]
    return calc_coll, ipea_coll


def build_ipea_db():
    sma = SubmissionMongoAdapter.auto_load()
    jobs = sma.jobs
    calc_coll, ipea_coll = get_calculation_db()
    ipea_coll.remove()
    for j in jobs.find():
        neutral_tid = None
        cation_tid = None
        anion_tid = None
        neutral_vaccum_energy = None
        neutral_sol_energy = None
        cation_vaccum_energy = None
        cation_sol_energy = None
        anion_vaccum_energy = None
        anion_sol_energy = None
        ip_vaccum = None
        ip_sol = None
        ea_vaccum = None
        ea_sol = None
        print j['submission_id']
        if 'original Single Point Energy' in j['task_dict']:
            neutral_tid = j['task_dict']['original Single Point Energy']
        if 'cation Single Point Energy' in j['task_dict']:
            cation_tid = j['task_dict']['cation Single Point Energy']
        if 'anion Single Point Energy' in j['task_dict']:
            anion_tid = j['task_dict']['anion Single Point Energy']

        if neutral_tid:
            neutral_energy_dict = calc_coll.find(
                {'task_id': neutral_tid},
                fields={"_id": False,
                        "calculations.scf.energies": True,
                        "calculations.scf_pcm.energies": True})[0]
            neutral_vaccum_energy = neutral_energy_dict['calculations'][
                'scf']['energies'][-1][-1]
            neutral_sol_energy = neutral_energy_dict['calculations'][
                'scf_pcm']['energies'][-1][-1]

        if cation_tid:
            cation_energy_dict = calc_coll.find(
                {'task_id': cation_tid},
                fields={"_id": False,
                        "calculations.scf.energies": True,
                        "calculations.scf_pcm.energies": True})[0]
            cation_vaccum_energy = cation_energy_dict['calculations']['scf'][
                'energies'][-1][-1]
            cation_sol_energy = cation_energy_dict['calculations']['scf_pcm'][
                'energies'][-1][-1]

        if anion_tid:
            anion_energy_dict = calc_coll.find(
                {'task_id': anion_tid},
                fields={"_id": False,
                        "calculations.scf.energies": True,
                        "calculations.scf_pcm.energies": True})[0]
            anion_vaccum_energy = anion_energy_dict['calculations']['scf'][
                'energies'][-1][-1]
            anion_sol_energy = anion_energy_dict['calculations']['scf_pcm'][
                'energies'][-1][-1]
        if neutral_tid:
            if cation_tid:
                ip_vaccum = cation_vaccum_energy - neutral_vaccum_energy
                ip_sol = cation_sol_energy - neutral_sol_energy
            if anion_tid:
                ea_vaccum = neutral_vaccum_energy - anion_vaccum_energy
                ea_sol = neutral_sol_energy - anion_sol_energy

        d = dict()
        d['submission_id'] = j['submission_id']
        if neutral_tid:
            if cation_tid:
                d['ip'] = {'vaccum': ip_vaccum, 'sol': ip_sol}
            if anion_tid:
                d['ea'] = {'vaccum': ea_vaccum, 'sol': ea_sol}
        if 'ip' in d or 'ea' in d:
            ipea_coll.update({'submission_id': j['submission_id']},
                             {"$set": d},
                             upsert=True)


if __name__ == '__main__':
    build_ipea_db()