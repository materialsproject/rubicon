
__author__ = 'navnidhirajput'


import shlex

import numpy
from monty import json
from pymongo import MongoClient

from rubicon.io.lammps.lammpsio import LammpsLog


class LammpsProperties():
    """
    Store lammps properties in output dictionary

    Right now RDF, MSD, Bulk diffusivity


    """

    _fw_name = "Lammps Properties Writer"

    def _insert_doc(self, fw_spec=None):
        db_dir = shlex.os.environ['DB_LOC']
        db_path = shlex.os.path.join(db_dir, 'tasks_db.json')
        with open(db_path) as f:
            db_creds = json.load(f)
        conn = MongoClient(db_creds['host'], db_creds['port'], )
        db = conn[db_creds['database']]
        if db_creds['admin_user']:
            db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
            coll = db['lammps_output']
        parse_lammps = LammpsLog.from_file('mol.log')
        docs = parse_lammps.llog
        docs = {k: list(v) if isinstance(v, numpy.ndarray) else v for k, v in
                docs.items()}
        coll.insert(docs)
        coll.update(docs)

    def run_task(self, fw_spec):
        self._insert_doc()
