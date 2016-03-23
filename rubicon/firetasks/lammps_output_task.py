# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import json
import shlex

import numpy
from fireworks import FireTaskBase, explicit_serialize
from pymongo import MongoClient

from rubicon.io.lammps.lammpsio import LammpsLog

__author__ = 'navnidhirajput'


@explicit_serialize
class WritelammpsOutputTask(FireTaskBase):
    """
    Writes LAMMPS Output files.

    Required params:


    Optional params:

    """

    _fw_name = "Lammps Output Writer"

    def _insert_doc(self, fw_spec=None, filename="mol.log"):
        db_dir = shlex.os.environ['DB_LOC']
        db_path = shlex.os.path.join(db_dir, 'tasks_db.json')
        with open(db_path) as f:
            db_creds = json.load(f)
        conn = MongoClient(db_creds['host'], db_creds['port'], )
        db = conn[db_creds['database']]
        print("+++++", db_creds)
        if db_creds['admin_user']:
            db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
            coll = db['lammps_output']

        parse_lammps = LammpsLog.from_file(filename)
        docs = parse_lammps.llog
        docs = {k: list(v) if isinstance(v, numpy.ndarray) else v for k, v in
                docs.items()}
        coll.insert(docs)

    def run_task(self, fw_spec):
        mol_log_file = fw_spec["prev_lammps_log"]
        self._insert_doc(filename=mol_log_file)
