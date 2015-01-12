import json
from pymongo import MongoClient
import shlex
import numpy

try:
    # just a walkaround before the packmol is merged to master branch
    # after packmol is merged to master branch, the try...catch block
    # should be removed
    from pymatgen.packmol.packmol import PackmolRunner
    from pymatgen.packmol.lammpsio import LammpsLog
except:
    pass
from fireworks import FireTaskBase, explicit_serialize, Firework, Workflow

__author__ = 'navnidhirajput'



# try:
#     # just a walkaround before the packmol is merged to master branch
#     # after packmol is merged to master branch, the try...catch block
#     # should be removed
#     import pymatgen
#     if 'packmol' in pymatgen.__dict__:
#         from pymatgen.packmol.packmol import PackmolRunner
#         from pymatgen.packmol.lammpsio import LammpsLog
# except:
#     pass



__author__ = 'navnidhirajput'



@explicit_serialize
class WritelammpsOutputTask(FireTaskBase):
    """
    Writes LAMMPS Output files.

    Required params:


    Optional params:

    """

    # _fw_name = "Lammps Output Writer"


    # def _insert_doc(self, fw_spec = None):
    #     db_dir = shlex.os.environ['DB_LOC']
    #     db_path = shlex.os.path.join(db_dir, 'tasks_db.json')
    #     with open(db_path) as f:
    #         db_creds = json.load(f)
    #     conn = MongoClient(db_creds['host'], db_creds['port'],)
    #     db = conn[db_creds['database']]
    #     if db_creds['admin_user']:
    #         db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
    #         coll = db['lammps_output']
    #     parse_lammps = LammpsLog.from_file('mol.log')
    #     docs = parse_lammps.llog
    #     docs = {k: list(v) if isinstance(v, numpy.ndarray) else v for k, v in docs.items()}
    #     coll.insert(docs)
    #     #coll.update(docs)
    #
    # def run_task(self, fw_spec):
    #     self._insert_doc()






