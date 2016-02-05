import shlex
from fireworks import explicit_serialize, FireTaskBase, FWAction
from flask import json
import numpy
from pymatgen.io import gaussian
from pymongo import MongoClient

__author__ = 'navnidhirajput'

@explicit_serialize
class GaussianGeomOptDBInsertionTask(FireTaskBase):
    _fw_name = "QChem Geometry Optimization DB Insertion Task"

    def _insert_doc(self, fw_spec = None, filename = "mol_geo.out"):
        db_dir = shlex.os.environ['DB_LOC']
        db_path = shlex.os.path.join(db_dir, 'tasks_db.json')
        with open(db_path) as f:
            db_creds = json.load(f)
        conn = MongoClient(db_creds['host'], db_creds['port'],)
        db = conn[db_creds['database']]
        if db_creds['admin_user']:
            db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
            coll = db['Gaussian_Geo']
        gaus_geo_parser = gaussian.GaussianOutput(filename)
        if gaus_geo_parser.properly_terminated == True:
            docs = gaus_geo_parser.final_structure
            coll.insert(docs.as_dict())






    def run_task(self, fw_spec):
        mol_geo_file = fw_spec["prev_gaussian_geo"]
        self._insert_doc(filename=mol_geo_file)
        #prev_gaussian_geo = shlex.os.path.join(shlex.os.getcwd(), 'mol_geo.out')


        file_path = fw_spec["prev_gaussian_geo"]
        update_spec = {'prev_gaussian_geo': file_path}
        return FWAction(update_spec=update_spec)




#
# if __name__ == '__main__':
#      task_geo = GaussianGeomOptDBInsertionTask()
#      task_geo._insert_doc()