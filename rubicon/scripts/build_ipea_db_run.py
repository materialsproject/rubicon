import json
import os
from pymongo import MongoClient
from rubicon.submission.submission_mongo_eg import SubmissionMongoAdapterEG

__author__ = 'xiaohuiqu'


def get_calculation_db(credential_file):
    db_dir = os.environ['DB_LOC']
    db_path = os.path.join(db_dir, credential_file)
    with open(db_path) as f:
        db_creds = json.load(f)
    conn = MongoClient(db_creds['host'], db_creds['port'])
    db = conn[db_creds['database']]
    if db_creds['admin_user']:
        db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
    coll = db[db_creds['collection']]
    return coll


def build_ipea_db():
    sma = SubmissionMongoAdapterEG.auto_load()
    jobs = sma.jobs
    tasks_coll = get_calculation_db("tasks_db.json")
    ipea_coll = get_calculation_db("molecules_db.json")
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
        if 'singlet neutral single point energy' in j['task_dict']:
            neutral_tid = j['task_dict'][
                'singlet neutral single point energy']
        if 'doublet cation single point energy' in j['task_dict']:
            cation_tid = j['task_dict'][
                'doublet cation single point energy']
        if 'doublet anion single point energy' in j['task_dict']:
            anion_tid = j['task_dict']['doublet anion single point energy']

        if neutral_tid:
            neutral_energy_dict = tasks_coll.find(
                {'task_id': neutral_tid},
                fields={"_id": False,
                        "state": True,
                        "calculations.scf.energies": True,
                        "calculations.scf_pcm.energies": True})[0]
            if neutral_energy_dict["state"] == "successful":
                neutral_vaccum_energy = neutral_energy_dict['calculations'][
                    'scf']['energies'][-1][-1]
                neutral_sol_energy = neutral_energy_dict['calculations'][
                    'scf_pcm']['energies'][-1][-1]
            else:
                neutral_tid = None

        if cation_tid:
            cation_energy_dict = tasks_coll.find(
                {'task_id': cation_tid},
                fields={"_id": False,
                        "state": True,
                        "calculations.scf.energies": True,
                        "calculations.scf_pcm.energies": True})[0]
            if cation_energy_dict["state"] == "successful":
                cation_vaccum_energy = cation_energy_dict['calculations']['scf'][
                    'energies'][-1][-1]
                cation_sol_energy = cation_energy_dict['calculations']['scf_pcm'][
                    'energies'][-1][-1]
            else:
                cation_tid = None

        if anion_tid:
            anion_energy_dict = tasks_coll.find(
                {'task_id': anion_tid},
                fields={"_id": False,
                        "state": True,
                        "calculations.scf.energies": True,
                        "calculations.scf_pcm.energies": True})[0]
            if anion_energy_dict["state"] == "successful":
                anion_vaccum_energy = anion_energy_dict['calculations']['scf'][
                    'energies'][-1][-1]
                anion_sol_energy = anion_energy_dict['calculations']['scf_pcm'][
                    'energies'][-1][-1]
            else:
                anion_tid = None
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