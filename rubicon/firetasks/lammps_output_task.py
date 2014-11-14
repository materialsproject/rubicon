from pymatgen.packmol.lammpsio import LammpsLog
from pymongo import MongoClient

__author__ = 'navnidhirajput'

import shlex
import subprocess
from monty import logging, json
from pymatgen import Molecule
from pymatgen.packmol.packmol import PackmolRunner
from rubicon.gff.boxmol import BoxMol
from rubicon.gff.lammpsin import DictLammpsInputSet
from rubicon.gff.lamppsio import LmpInput
from rubicon.gff.antechamberio import AntechamberRunner


__author__ = 'navnidhirajput'


from fireworks import FireTaskBase, explicit_serialize, Firework, Workflow


@explicit_serialize
class WritelammpsOutputTask(FireTaskBase):
    """
    Writes LAMMPS Output files.

    Required params:


    Optional params:

    """

    _fw_name = "Lammps Output Writer"


    def _insert_doc(self, fw_spec = None, update_duplicates = True):
        db_dir = shlex.os.environ['DB_LOC']
        db_path = shlex.os.path.join(db_dir, 'tasks_db.json')
        with open(db_path) as f:
            db_creds = json.load(f)
        conn = MongoClient(db_creds['host'], db_creds['port'],)
        db = conn[db_creds['database']]
        if db_creds['admin_user']:
            db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
        coll = db[lammps_data]
        conn.close()


        # parsing code
        #qcout = QcOutput(zpath(path))

                    #dn update/insertion code
                    # coll.update({"super_mol_snlgroup_id": fw_spec["snlgroup_id"],
                    #      "fragments_def": fw_spec["fragments"]},
                    #     {"$set": d},
                    #     upsert=True)
        llog=LammpsLog()
        print llog.avgs



test_yong = LammpsLog.from_file('mol.log')
docs = test_yong.llog



# if db.counter.find({"_id": "mol_taskid"}).count() == 0:
#         db.counter.insert({"_id": "mol_taskid", "c": 1})
#             conn.close()

docs["task_id"] = 348