#!/usr/bin/env python

"""
TODO: Modify module doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "6/15/13"


import os
import logging
import datetime

from pymongo import MongoClient

from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.io.nwchemio import NwOutput
from pymatgen.util.io_utils import clean_json
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.xyzio import XYZ

from rubicon.testset.parse_mol import get_nih_names

import json

logger = logging.getLogger(__name__)


class DeltaSCFNwChemToDbTaskDrone(AbstractDrone):
    """
    Assimilates a delta scf nwchem run and inserts it into the database.
    """

    #Version of this db creator document.
    __version__ = "0.1.0"

    def __init__(self, host="127.0.0.1", port=27017, database="mg_core_dev",
                 user=None, password=None,  collection="mol_tasks",
                 simulate_mode=False, update_duplicates=True):
        """
        Args:
            host:
                Hostname of database machine. Defaults to 127.0.0.1 or
                localhost.
            port:
                Port for db access. Defaults to mongo's default of 27017.
            database:
                Actual database to access. Defaults to "mg_core_dev".
            user:
                User for db access. Requires write access. Defaults to None,
                which means no authentication.
            password:
                Password for db access. Requires write access. Defaults to
                None, which means no authentication.
            collection:
                Collection to query. Defaults to "mol_tasks".
            simulate_mode:
                Allows one to simulate db insertion without actually performing
                the insertion.
        """
        self.host = host
        self.database = database
        self.user = user
        self.password = password
        self.collection = collection
        self.port = port
        self.simulate = simulate_mode
        self.update_duplicates = update_duplicates
        if not simulate_mode:
            conn = MongoClient(self.host, self.port, j=True)
            db = conn[self.database]
            if self.user:
                db.authenticate(self.user, self.password)
            if db.counter.find({"_id": "mol_taskid"}).count() == 0:
                db.counter.insert({"_id": "mol_taskid", "c": 1})
            conn.close()

    def assimilate(self, path):
        """
        Parses nwchem runs. Then insert the result into the db. and return the
        mol_id or doc of the insertion.

        Returns:
            If in simulate_mode, the entire doc is returned for debugging
            purposes. Else, only the task_id of the inserted doc is returned.
        """
        try:
            d = self.get_task_doc(path)
            tid = self._insert_doc(d)
            return tid, d
        except Exception as ex:
            import traceback
            print traceback.format_exc(ex)
            logger.error(traceback.format_exc(ex))
            return False

    @classmethod
    def get_user_tags(cls, path):
        '''
        Parse the user_tags from the FW.json file.
        The user_tags can be set in the creation of FireWork
        '''
        fwjsonfile = os.path.join(os.path.dirname(path), 'FW.json')
        user_tags = {}
        with open(fwjsonfile) as f:
            d = json.load(f)
        if 'user_tags' in d['spec'].keys():
            user_tags = d['spec']['user_tags']
        if 'name' in d:
            user_tags['fw_name'] = d['name']
        if len(user_tags)>0:
            return user_tags
        else:
            return None

    @classmethod
    def get_task_doc(cls, path):
        """
        Get the entire task doc for a path, including any post-processing.
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        nwo = NwOutput(path)
        data = nwo.data
        mol = data[0]["molecules"][-1]
        bb = BabelMolAdaptor(mol)
        pbmol = bb.pybel_mol
        xyz = XYZ(mol)
        smiles = pbmol.write("smi").split()[0]
        can = pbmol.write("can").split()[0]
        inchi = pbmol.write("inchi")
        svg = pbmol.write("svg")
        comp = mol.composition
        initial_mol = data[0]["molecules"][0]
        charge = data[0]["charge"]
        spin_mult = data[0]["spin_multiplicity"]
        data_dict = {}
        from pymatgen.symmetry.pointgroup import PointGroupAnalyzer

        pga = PointGroupAnalyzer(mol)
        sch_symbol = pga.sch_symbol
        stationary_type = None
        for d in data:
            if d["job_type"] == "NWChem Geometry Optimization":
                data_dict["geom_opt"] = d
            elif d["job_type"] == "NWChem Nuclear Hessian and Frequency " \
                                  "Analysis":
                data_dict["freq"] = d
                if d['frequencies'][0][0] < -0.00:
                    # it is stupied that -0.00 is less than 0.00
                    stationary_type = "non-minimum"
                else:
                    stationary_type = "minimum"
            elif d["job_type"] == "NWChem DFT Module":
                if d["charge"] == charge:
                    data_dict["scf"] = d
                elif d["charge"] == charge + 1:
                    data_dict["scf_IE"] = d
                elif d["charge"] == charge - 1:
                    data_dict["scf_EA"] = d
            elif  d["job_type"] == "NWChem DFT Module COSMO":
                if d["charge"] == charge:
                    data_dict["sol_scf"] = d
                elif d["charge"] == charge + 1:
                    data_dict["sol_IE"] = d
                elif d["charge"] == charge - 1:
                    data_dict["sol_EA"] = d



        data = data_dict


        d = {"path": os.path.abspath(path),
             "folder": os.path.basename(os.path.dirname(os.path.abspath(
                 path))),
             "calculations": data,
             "initial_molecule": initial_mol.to_dict,
             "final_molecule": mol.to_dict,
             "pointgroup": sch_symbol,
             "pretty_formula": comp.reduced_formula,
             "reduced_cell_formula_abc": comp.alphabetical_formula,
             "formula": comp.formula,
             "charge": charge,
             "spin_mult": spin_mult,
             "composition": comp.to_dict,
             "elements": list(comp.to_dict.keys()),
             "nelements": len(comp),
             "smiles": smiles, "can": can, "inchi": inchi, "svg": svg,
             "xyz": str(xyz),
             "names": get_nih_names(smiles)}

        if "scf_EA" in data_dict and \
                (not data_dict["scf_EA"]["has_error"]) and \
                (not data_dict["scf"]["has_error"]):
            d["EA"] = (data["scf"]["energies"][-1]
                       - data["scf_EA"]["energies"][-1])
        if "scf_IE" in data_dict and \
                (not data_dict["scf_IE"]["has_error"]) and \
                (not data_dict["scf"]["has_error"]):
            d["IE"] = (data["scf_IE"]["energies"][-1]
                       - data["scf"]["energies"][-1])

        if stationary_type:
            d['stationary_type'] = stationary_type

        user_tags = cls.get_user_tags(path)
        if user_tags:
            d['user_tags'] = user_tags
            if "initial_inchi" in user_tags:
                initial_inchi = user_tags['initial_inchi']
                if initial_inchi != d['inchi']:
                    d['state'] = 'rejected'
                    d['reject_reason'] = 'structural change'
        if "state" not in d:
            for scf in data_dict.values():
                if scf['has_error']:
                    d["state"] = "error"
        if "state" not in d:
            d["state"] = "successful"

        return clean_json(d)

    def _insert_doc(self, d):
        if not self.simulate:
            # Perform actual insertion into db. Because db connections cannot
            # be pickled, every insertion needs to create a new connection
            # to the db.
            conn = MongoClient(self.host, self.port)
            db = conn[self.database]
            if self.user:
                db.authenticate(self.user, self.password)
            coll = db[self.collection]

            result = coll.find_one({"path": d["path"]},
                                   fields=["path", "task_id"])
            if result is None or self.update_duplicates:
                d["last_updated"] = datetime.datetime.today()
                if result is None:
                    if ("task_id" not in d) or (not d["task_id"]):
                        id_num = db.counter.find_and_modify(
                            query={"_id": "mol_taskid"},
                            update={"$inc": {"c": 1}}
                        )["c"]
                        d["task_id"] = "mol-" + str(id_num)
                        d["task_id_deprecated"] = id_num
                    logger.info("Inserting {} with taskid = {}"
                                .format(d["path"], d["task_id"]))
                elif self.update_duplicates:
                    d["task_id"] = result["task_id"]
                    logger.info("Updating {} with taskid = {}"
                                .format(d["path"], d["task_id"]))

                coll.update({"path": d["path"]}, {"$set": d}, upsert=True)
                return d["task_id"]
            else:
                logger.info("Skipping duplicate {}".format(d["dir_name"]))
            conn.close()
        else:
            d["task_id"] = 0
            logger.info("Simulated insert into database for {} with task_id {}"
                        .format(d["path"], d["task_id"]))
            return d

    def get_valid_paths(self, path):
        """
        There are some restrictions on the valid directory structures:

        1. There can be only one vasp run in each directory. Nested directories
           are fine.
        2. Directories designated "relax1", "relax2" are considered to be 2
           parts of an aflow style run.
        3. Directories containing vasp output with ".relax1" and ".relax2" are
           also considered as 2 parts of an aflow style run.
        """
        (parent, subdirs, files) = path

        return [os.path.join(parent, f) for f in files
                if f.endswith(".nwout")]

    def convert(self, d):
        return d

    def __str__(self):
        return self.__class__.__name__

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])

    @property
    def to_dict(self):
        init_args = {"host": self.host, "port": self.port,
                     "database": self.database, "user": self.user,
                     "password": self.password,
                     "collection": self.collection,
                     "simulate_mode": self.simulate}
        output = {"name": self.__class__.__name__,
                  "init_args": init_args, "version": __version__}
        return output


import unittest

from pymatgen.apps.borg.queen import BorgQueen


test_dir = os.path.join(os.path.dirname(__file__), "..", "..",
                        'test_files')


class DeltaSCFNwChemToDbTaskDroneTest(unittest.TestCase):

    def test_assimilate(self):
        import logging
        logging.basicConfig(level=logging.INFO)
        drone = DeltaSCFNwChemToDbTaskDrone(
            collection="mol_task_test")
        q = BorgQueen(drone)
        q.serial_assimilate(test_dir)



if __name__ == "__main__":
    unittest.main()