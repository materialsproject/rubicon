import json
import os
from pymongo import MongoClient
from rubicon.submission.submission_mongo import SubmissionMongoAdapter

__author__ = 'xiaohuiqu'

def get_calculation_db():
    db_dir = os.environ['DB_LOC']
    db_path = os.path.join(db_dir, 'nwchem_calc_db.json')
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
        print j['submission_id']
        neutral_tid = j['task_dict']['original Single Point Energy']
        cation_tid = j['task_dict']['cation Single Point Energy']
        anion_tid = j['task_dict']['anion Single Point Energy']

        neutral_energy_dict = calc_coll.find({'task_id': neutral_tid},
                                                fields={"_id":False,
                                                        "calculations.scf.energies": True,
                                                        "calculations.sol_scf.energies": True})[0]
        neutral_vaccum_energy = neutral_energy_dict['calculations']['scf']['energies'][0]
        neutral_sol_energy = neutral_energy_dict['calculations']['sol_scf'] \
                                                ['energies'][0]['sol phase']

        cation_energy_dict = calc_coll.find({'task_id': cation_tid},
                                                fields={"_id":False,
                                                        "calculations.scf.energies": True,
                                                        "calculations.sol_scf.energies": True})[0]
        cation_vaccum_energy = cation_energy_dict['calculations']['scf']['energies'][0]
        cation_sol_energy = cation_energy_dict['calculations']['sol_scf'] \
                                              ['energies'][0]['sol phase']

        anion_energy_dict = calc_coll.find({'task_id': anion_tid},
                                                fields={"_id":False,
                                                        "calculations.scf.energies": True,
                                                        "calculations.sol_scf.energies": True})[0]
        anion_vaccum_energy = anion_energy_dict['calculations']['scf']['energies'][0]
        anion_sol_energy = anion_energy_dict['calculations']['sol_scf'] \
                                              ['energies'][0]['sol phase']

        ip_vaccum = cation_vaccum_energy - neutral_vaccum_energy
        ip_sol = cation_sol_energy - neutral_sol_energy

        ea_vaccum = neutral_vaccum_energy - anion_vaccum_energy
        ea_sol = neutral_sol_energy - anion_sol_energy

        d = {'ip': {'vaccum': ip_vaccum, 'sol': ip_sol},
             'ea': {'vaccum': ea_vaccum, 'sol': ea_sol},
             'submission_id':j['submission_id']}
        ipea_coll.update({'submission_id':j['submission_id']}, {"$set": d}, upsert=True)



if __name__ == '__main__':
    build_ipea_db()