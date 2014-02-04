#!/usr/bin/env python

"""
TODO: Modify module doc.
"""

from __future__ import division
from pymatgen import Molecule
from pymatgen.analysis.molecule_structure_comparator import MoleculeStructureComparator
from pymatgen.io.qchemio import QcOutput
from rubicon.utils.snl.egsnl import EGStructureNL
from rubicon.utils.snl.egsnl_mongo import EGSNLMongoAdapter

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
from pymatgen.util.io_utils import clean_json
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.xyzio import XYZ

from rubicon.testset.parse_mol import get_nih_names

import json
import re

logger = logging.getLogger(__name__)


class DeltaSCFQChemToDbTaskDrone(AbstractDrone):
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

    def assimilate(self, path, fw_spec=None):
        """
        Parses nwchem runs. Then insert the result into the db. and return the
        mol_id or doc of the insertion.

        Returns:
            If in simulate_mode, the entire doc is returned for debugging
            purposes. Else, only the task_id of the inserted doc is returned.
        """
        try:
            d = self.get_task_doc(path, fw_spec)
            tid = self._insert_doc(d, fw_spec)
            return tid, d
        except Exception as ex:
            import traceback
            print traceback.format_exc(ex)
            logger.error(traceback.format_exc(ex))
            return False

    @classmethod
    def modify_svg(cls, svg):
        """
        Hack to the svg code to enhance the molecule.
        Because Xiaohui have no aeshetic cell, please change it a more
        beautiful color scheme
        """
        tokens = svg.split('\n')
        new_tokens = []
        for line in tokens:
            if "line" in line:
                line = re.sub('stroke-width="\d+.?\d*"',
                              'stroke-width="3.0"', line)
            if "rect" in line:
                line = re.sub(r'fill=".+?"', 'fill="Beige"', line)
            if ">H</text>" in line:
                line = re.sub(r'fill=".+?"', 'fill="DimGray"', line)
                line = re.sub(r'stroke=".+?"', 'stroke="DimGray"', line)
            new_tokens.append(line)
        new_svg = '\n'.join(new_tokens)
        return new_svg

    @staticmethod
    def _check_structure_change(mol1, mol2):
        """
        Check whether structure is changed:

        Return:
            True: structure changed, False: unchanged
        """
        with open("FW.json") as f:
            fw = json.load(f)
        if 'egsnl' not in fw['spec']:
            raise ValueError("Can't find initial SNL")
        if 'known_bonds' not in fw['spec']['egsnl']:
            raise ValueError("Can't find known bonds information")
        bonds = fw['spec']['egsnl']['known_bonds']
        msc = MoleculeStructureComparator(priority_bonds=bonds)
        return not msc.are_equal(mol1, mol2)

    @classmethod
    def get_task_doc(cls, path, fw_spec=None):
        """
        Get the entire task doc for a path, including any post-processing.
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        qcout = QcOutput(path)
        data = qcout.data
        mol = data[0]["molecules"][-1]
        bb = BabelMolAdaptor(mol)
        pbmol = bb.pybel_mol
        xyz = XYZ(mol)
        smiles = pbmol.write("smi").split()[0]
        can = pbmol.write("can").split()[0]
        inchi_final = pbmol.write("inchi").strip()
        svg = cls.modify_svg(pbmol.write("svg"))
        comp = mol.composition
        initial_mol = data[0]["molecules"][0]
        charge = mol.charge
        spin_mult = mol.spin_multiplicity
        data_dict = {}
        from pymatgen.symmetry.pointgroup import PointGroupAnalyzer

        pga = PointGroupAnalyzer(mol)
        sch_symbol = pga.sch_symbol
        stationary_type = None
        for d in data:
            # noinspection PyTypeChecker
            if d["jobtype"] == "opt":
                data_dict["geom_opt"] = d
            elif d["jobtype"] == "freq":
                data_dict["freq"] = d
                # noinspection PyTypeChecker
                if not d["has_error"]:
                    if d['frequencies'][0]["frequency"] < -0.00:
                        # it is stupied that -0.00 is less than 0.00
                        stationary_type = "non-minimum"
                    else:
                        stationary_type = "minimum"
                else:
                    stationary_type = "unknown"
            elif d["jobtype"] == "sp":
                # noinspection PyTypeChecker
                suffix = "" if d["solvent_method"] == "NA" \
                    else "_" + d["solvent_method"]
                data_dict["scf" + suffix] = d

        data = data_dict

        d = {"path": os.path.abspath(path),
             "folder": os.path.basename(os.path.dirname(os.path.abspath(
                 path))),
             "calculations": data,
             "molecule_initial": initial_mol.to_dict,
             "molecule_final": mol.to_dict,
             "pointgroup": sch_symbol,
             "pretty_formula": comp.reduced_formula,
             "reduced_cell_formula_abc": comp.alphabetical_formula,
             "formula": comp.formula,
             "charge": charge,
             "spin_multiplicity": spin_mult,
             "composition": comp.to_dict,
             "elements": list(comp.to_dict.keys()),
             "nelements": len(comp),
             "smiles": smiles, "can": can,
             "inchi_final": inchi_final,
             "svg": svg,
             "xyz": str(xyz),
             "names": get_nih_names(smiles)}

        if stationary_type:
            d['stationary_type'] = stationary_type
        if fw_spec:
            inchi_initial = fw_spec['inchi']
            if inchi_initial != d['inchi_final']:
                d['inchi_changed'] = True
            else:
                d['inchi_changed'] = False
        d['structure_changed'] = cls._check_structure_change(initial_mol, mol)
        if d['structure_changed']:
            d['state'] = 'rejected'
            d['reject_reason'] = 'structural change'
        if "state" not in d:
            for v in data_dict.values():
                if v['has_error']:
                    d['state'] = "error"
                    errors = d.get("errors", [])
                    errors += v["errors"]
                    d["errors"] = errors
        if "state" not in d:
            d["state"] = "successful"

        return clean_json(d)

    @staticmethod
    def update_tags(fw_spec, d, task_id):
        if fw_spec:
                d['task_type'] = fw_spec['task_type']
                if '_fizzled_parents' in fw_spec and not 'prev_task_type' in \
                        fw_spec:
                    d['task_type'] = fw_spec['_fizzled_parents'][0][
                        'task_type']
                else:
                    d['task_type'] = fw_spec['prev_task_type']
                d['run_tags'] = fw_spec['run_tags']
                d['implicit_solvent'] = fw_spec['implicit_solvent']
                d['user_tags'] = fw_spec["user_tags"]
                d['snl_initial'] = fw_spec['egsnl']
                d['snlgroup_id_initial'] = fw_spec['snlgroup_id']
                d['inchi_root'] = fw_spec['inchi_root']
                d['inchi_initial'] = fw_spec['inchi']
                if "geometry optimization" in d['task_type']:
                    new_s = Molecule.from_dict(d["molecule_final"])
                    old_snl = EGStructureNL.from_dict(d['snl_initial'])
                    history = old_snl.history
                    history.append(
                        {'name': 'Electrolyte Genome Project structure optimization',
                         'url': 'http://www.materialsproject.org',
                         'description': {'task_type': d['task_type'],
                                         'task_id': task_id,
                                         'when': datetime.datetime.utcnow()}})
                    new_snl = EGStructureNL(new_s, old_snl.authors, old_snl.projects,
                                            old_snl.references, old_snl.remarks,
                                            old_snl.data, history)

                    # enter new SNL into SNL db
                    # get the SNL mongo adapter
                    sma = EGSNLMongoAdapter.auto_load()

                    # add snl
                    egsnl, snlgroup_id = sma.add_snl(
                        new_snl, snlgroup_guess=d['snlgroup_id_initial'])
                    d['snl_final'] = egsnl.to_dict
                    d['snlgroup_id_final'] = snlgroup_id
                else:
                    d['snl_final'] = fw_spec['egsnl']
                    d['snlgroup_id_final'] = fw_spec['snlgroup_id']
                d['snlgroup_changed'] = (d['snlgroup_id_initial'] !=
                                             d['snlgroup_id_final'])

    def _insert_doc(self, d, fw_spec=None):
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
                self.update_tags(fw_spec, d, d["task_id"])
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
                if f.endswith(".qcout")]

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