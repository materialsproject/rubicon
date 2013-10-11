import copy
import json
import logging
import os
from fireworks.core.firework import FireTaskBase, FWAction, Workflow, FireWork
from fireworks.utilities.fw_serializers import FWSerializable
import sys
from pymatgen.core.structure import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.nwchemio import NwTask, NwInput
from pymatgen.symmetry.pointgroup import PointGroupAnalyzer
from pymongo import MongoClient
from rubicon.borg.hive import DeltaSCFNwChemToDbTaskDrone
from rubicon.firetasks.nwchem_task import NWChemTask


__author__ = 'xiaohuiqu'


class NWChemFireWorkCreator():
    def __init__(self, mol, molname, mission, additional_user_tags=None):
        theory_directive = {"iterations": 300, "vectors": "atomic"}
        symmetry_options = None
        pga = PointGroupAnalyzer(mol)
        if pga.sch_symbol == 'D*h' or pga.sch_symbol == 'C*v':
            # linear molecule, turn off symmetry
            symmetry_options = ['c1']
        initial_inchi = self.get_inchi(mol)
        user_tags = {'mission': mission,
                     "initial_inchi": initial_inchi,
                     "molname": molname}
        if additional_user_tags:
            user_tags.update(additional_user_tags)

        self.td = lambda: copy.deepcopy(theory_directive)
        self.sym = lambda: copy.deepcopy(symmetry_options)
        self.bs = '6-31+G*'
        self.ut = lambda: copy.deepcopy(user_tags)
        self.mol = mol


    def get_inchi(self, mol):
        bb = BabelMolAdaptor(mol)
        pbmol = bb.pybel_mol
        return pbmol.write("inchi")


    def geom_fw(self, charge_shift, fw_id):
        charge = self.mol.charge + charge_shift
        tasks_geom = [NwTask.dft_task(self.mol, charge=charge,
                                      operation="optimize",
                                      xc="b3lyp", basis_set=self.bs,
                                      theory_directives=self.td(),
                                      alternate_directives={"driver":
                                                            {"maxiter": 300}})]
        nwi = NwInput(self.mol, tasks_geom, symmetry_options=self.sym())
        spec = nwi.to_dict
        spec['user_tags'] = self.ut()
        charge_state_name = {0: "original", 1: "cation", -1: "anion"}
        spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
        spec['user_tags']['charge_shift'] = charge_shift
        task_name = charge_state_name[charge_shift] + ' geom opt'
        fw_geom = FireWork([NWChemTask(), NWChemGeomOptDBInsertionTask()],
                           spec=spec, name=task_name, fw_id=fw_id)
        return fw_geom


    def freq_fw(self, charge_shift, fw_id):
        charge = self.mol.charge + charge_shift
        tasks_geom = [NwTask.dft_task(self.mol, charge=charge,
                                      operation="freq",
                                      xc="b3lyp", basis_set=self.bs,
                                      theory_directives=self.td())]
        nwi = NwInput(self.mol, tasks_geom, symmetry_options=self.sym())
        spec = nwi.to_dict
        spec['user_tags'] = self.ut()
        charge_state_name = {0: "original", 1: "cation", -1: "anion"}
        spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
        spec['user_tags']['charge_shift'] = charge_shift
        task_name = charge_state_name[charge_shift] + ' freq'
        fw_freq = FireWork([NWChemTask(), NWChemFrequencyDBInsertionTask()],
                           spec=spec, name=task_name, fw_id=fw_id)
        return fw_freq


    def sp_fw(self, charge_shift, fw_id):
        charge = self.mol.charge + charge_shift
        tasks_geom = [NwTask.dft_task(self.mol, charge=charge,
                                      operation="energy",
                                      xc="b3lyp", basis_set=self.bs,
                                      theory_directives=self.td()),
                      NwTask.dft_task(self.mol, charge=charge,
                                      operation="energy",
                                      xc="b3lyp", basis_set=self.bs,
                                      theory_directives=self.td(),
                                      alternate_directives={'cosmo':
                                                            {"dielec": 78.0}})]
        nwi = NwInput(self.mol, tasks_geom, symmetry_options=self.sym())
        spec = nwi.to_dict
        spec['user_tags'] = self.ut()
        charge_state_name = {0: "original", 1: "cation", -1: "anion"}
        spec['user_tags']['charge_state'] = charge_state_name[charge_shift]
        spec['user_tags']['charge_shift'] = charge_shift
        task_name = charge_state_name[charge_shift] + ' single point energy'
        fw_sp = FireWork([NWChemTask(),
                          NWChemSinglePointEnergyDBInsertionTask()],
                         spec=spec, name=task_name, fw_id=fw_id)
        return fw_sp


class NWChemGeomOptDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem Geometry Optimization DB Insertion Task"

    def run_task(self, fw_spec):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'nwchem_calc_db.json')

        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('NWChemDrone')
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        logger.addHandler(sh)

        with open(db_path) as f:
            db_creds = json.load(f)
        drone = DeltaSCFNwChemToDbTaskDrone(
            host=db_creds['host'], port=db_creds['port'],
            database=db_creds['database'], user=db_creds['admin_user'],
            password=db_creds['admin_password'],
            collection=db_creds['collection'])
        assi_result = drone.assimilate(os.path.abspath(
                                os.path.join(os.getcwd(), "mol.nwout")))

        t_id = None
        d = None
        if assi_result:
            t_id, d = assi_result
        if t_id:
            if d["state"] == "successful":
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"]})
            else:
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"]},
                                defuse_children=True)
        else:
            return FWAction(defuse_children=True)



class NWChemFrequencyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem Frequency DB Insertion Task"

    def run_task(self, fw_spec):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'nwchem_calc_db.json')

        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('NWChemDrone')
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        logger.addHandler(sh)

        with open(db_path) as f:
            db_creds = json.load(f)
        drone = DeltaSCFNwChemToDbTaskDrone(
            host=db_creds['host'], port=db_creds['port'],
            database=db_creds['database'], user=db_creds['admin_user'],
            password=db_creds['admin_password'],
            collection=db_creds['collection'])
        assi_result = drone.assimilate(os.path.abspath(
                                os.path.join(os.getcwd(), "mol.nwout")))

        t_id = None
        d = None
        if assi_result:
            t_id, d = assi_result
        if t_id:
            if d["state"] == "successful":
                if d['stationary_type'] == 'minimum':
                    return FWAction(stored_data={'task_id': t_id},
                                    update_spec={"mol": d["final_molecule"]})
                else:
                    return self.img_freq_action(fw_spec, d, t_id)
            else:
                return FWAction(stored_data={'task_id': t_id},
                                defuse_children=True,
                                update_spec={"mol": d["final_molecule"]})
        else:
            return FWAction(defuse_children=True)

    def img_freq_action(self, fw_spec, d, t_id):
        if 'img_freq_max_trial' in fw_spec:
            max_fix_time = fw_spec['img_freq_max_trial']
        else:
            max_fix_time = 3
        if 'freq_fix_time' in d['user_tags']:
            fix_time = d['user_tags']['freq_fix_time'] + 1
        else:
            fix_time = 1
            db_dir = os.environ['DB_LOC']
            db_path = os.path.join(db_dir, 'nwchem_calc_db.json')
            with open(db_path) as f:
                db_creds = json.load(f)
            conn = MongoClient(db_creds['host'], db_creds['port'])
            db = conn[db_creds['database']]
            if db_creds['admin_user']:
                db.authenticate(db_creds['admin_user'],
                                db_creds['admin_password'])
            coll = db[db_creds['collection']]
            coll.update({"task_id": t_id},
                        {"$set": {'user_tags.freq_fix_time': 0}}, upsert=True)

        if fix_time > max_fix_time:
            return FWAction(stored_data={'task_id': t_id},
                            defuse_children=True,
                            update_spec={"mol": d["final_molecule"]})
        else:
            old_mol = Molecule.from_dict(d['final_molecule'])
            vib_mode = d['calculations']['freq']['frequencies']
            new_coords = [[c+v for c, v in zip(site.coords, mode[1])]
                          for site, mode in zip(old_mol.sites, vib_mode)]
            species = [site.specie.symbol
                       for site in old_mol.sites]
            new_mol = Molecule(species, new_coords)
            fw_creator = NWChemFireWorkCreator(new_mol,
                             d['user_tags']['molname'],
                             d['user_tags']['mission'],
                             additional_user_tags={'freq_fix_time': fix_time})
            geom_fwid, freq_fwid = 0, 1
            geom_fw = fw_creator.geom_fw(d['user_tags']['charge_shift'], geom_fwid)
            freq_fw = fw_creator.freq_fw(d['user_tags']['charge_shift'], freq_fwid)
            wf = Workflow([geom_fw, freq_fw],
                          links_dict={geom_fwid: freq_fwid})
            return FWAction(stored_data={'task_id': t_id},detours=wf)



class NWChemSinglePointEnergyDBInsertionTask(FireTaskBase, FWSerializable):
    _fw_name = "NWChem Single Point Energy DB Insertion Task"

    def run_task(self, fw_spec):
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'nwchem_calc_db.json')

        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('NWChemDrone')
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        logger.addHandler(sh)

        with open(db_path) as f:
            db_creds = json.load(f)
        drone = DeltaSCFNwChemToDbTaskDrone(
            host=db_creds['host'], port=db_creds['port'],
            database=db_creds['database'], user=db_creds['admin_user'],
            password=db_creds['admin_password'],
            collection=db_creds['collection'])
        assi_result = drone.assimilate(os.path.abspath(
                                os.path.join(os.getcwd(), "mol.nwout")))

        t_id = None
        d = None
        if assi_result:
            t_id, d = assi_result
        if t_id:
            if d["state"] == "successful":
                return FWAction(stored_data={'task_id': t_id},
                                update_spec={"mol": d["final_molecule"]})
            else:
                return FWAction(stored_data={'task_id': t_id},
                                defuse_children=True,
                                update_spec={"mol": d["final_molecule"]})
        else:
            return FWAction(defuse_children=True)

