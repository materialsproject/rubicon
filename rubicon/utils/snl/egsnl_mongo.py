# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import datetime
import os

from fireworks.utilities.fw_serializers import FWSerializable
from pymongo import MongoClient, DESCENDING

from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from rubicon.utils.snl.egsnl import EGStructureNL, SNLGroup


class EGSNLMongoAdapter(FWSerializable):
    def __init__(self, host='localhost', port=27017, db='snl', username=None,
                 password=None):
        self.host = host
        self.port = port
        self.db = db
        self.username = username
        self.password = password

        self.connection = MongoClient(host, port, j=False,
                                      connect=False)
        self.database = self.connection[db]
        if self.username:
            self.database.authenticate(username, password)

        self.snl = self.database.snl
        self.snlgroups = self.database.snlgroups
        self.id_assigner = self.database.id_assigner

        self._update_indices()

    def _reset(self):
        if "prod" in self.database.name:
            print("PROD database is not supposed to reset, please changed " \
                  "the code to reset")
            exit()
        self.restart_id_assigner_at(1, 1)
        self.snl.remove()
        self.snlgroups.remove()

    def _update_indices(self):
        self.snl.ensure_index('snl_id', unique=True)
        self.snl.ensure_index('autometa.natoms')
        self.snl.ensure_index('autometa.nelements')
        self.snl.ensure_index('autometa.formula')
        self.snl.ensure_index('autometa.reduced_cell_formula')
        self.snl.ensure_index('autometa.reduced_cell_formula_abc')
        self.snl.ensure_index('autometa.inchi')

        self.snlgroups.ensure_index('snlgroup_id', unique=True)
        self.snlgroups.ensure_index('all_snl_ids')
        self.snlgroups.ensure_index('canonical_snl.snl_id')
        self.snlgroups.ensure_index('autometa.atoms')
        self.snlgroups.ensure_index('autometa.nelements')
        self.snlgroups.ensure_index('autometa.formula')
        self.snlgroups.ensure_index('autometa.reduced_cell_formula')
        self.snlgroups.ensure_index('autometa.reduced_cell_formula_abc')

    def _get_next_snl_id(self):
        snl_id = self.id_assigner.find_and_modify(
            query={}, update={'$inc': {'next_snl_id': 1}})['next_snl_id']
        return snl_id

    def _get_next_snlgroup_id(self):
        snlgroup_id = self.id_assigner.find_and_modify(
            query={},
            update={'$inc': {'next_snlgroup_id': 1}})['next_snlgroup_id']
        return snlgroup_id

    def restart_id_assigner_at(self, next_snl_id, next_snlgroup_id):
        self.id_assigner.remove()
        self.id_assigner.insert(
            {"next_snl_id": next_snl_id, "next_snlgroup_id": next_snlgroup_id})

    def add_snl(self, snl, force_new=False, snlgroup_guess=None):
        snl_id = self._get_next_snl_id()
        pointgroup = PointGroupAnalyzer(snl.structure).sch_symbol
        egsnl = EGStructureNL.from_snl(snl, snl_id, pointgroup)
        snlgroup, add_new = self.add_egsnl(egsnl, force_new, snlgroup_guess)
        return egsnl, snlgroup.snlgroup_id

    def add_egsnl(self, egsnl, force_new=False, snlgroup_guess=None):
        snl_d = egsnl.as_dict()
        snl_d['snl_timestamp'] = datetime.datetime.utcnow().isoformat()
        self.snl.insert(snl_d)
        return self.build_groups(egsnl, force_new, snlgroup_guess)

    def _add_if_belongs(self, snlgroup, egsnl, testing_mode):
        if snlgroup.add_if_belongs(egsnl):
            print('MATCH FOUND, grouping (snl_id, snlgroup): {}'. \
                  format((egsnl.snl_id, snlgroup.snlgroup_id)))
            if not testing_mode:
                self.snlgroups.update({'snlgroup_id': snlgroup.snlgroup_id},
                                      snlgroup.as_dict())
            return True
        return False

    def build_groups(self, egsnl, force_new=False, snlgroup_guess=None,
                     testing_mode=False):
        # testing mode is used to see if something already exists in DB w/o
        # adding it to the db
        match_found = False
        if not force_new:
            if snlgroup_guess:
                sgp = self.snlgroups.find_one(
                    filter={'snlgroup_id': snlgroup_guess})
                snlgroup = SNLGroup.from_dict(sgp)
                match_found = self._add_if_belongs(snlgroup, egsnl,
                                                   testing_mode)

            if not match_found:
                # look at all potential matches
                for entry in self.snlgroups.find(
                        filter={'snlgroup_key': egsnl.snlgroup_key},
                        sort=[("num_snl", DESCENDING)]):
                    snlgroup = SNLGroup.from_dict(entry)
                    match_found = self._add_if_belongs(snlgroup, egsnl,
                                                       testing_mode)
                    if match_found:
                        break

        if not match_found:
            # add a new SNLGroup
            snlgroup_id = self._get_next_snlgroup_id()
            snlgroup = SNLGroup(snlgroup_id, egsnl)
            if not testing_mode:
                self.snlgroups.insert(snlgroup.as_dict())

        return snlgroup, not match_found

    def switch_canonical_snl(self, snlgroup_id, canonical_egsnl):
        sgp = self.snlgroups.find_one(filter={'snlgroup_id': snlgroup_id})
        snlgroup = SNLGroup.from_dict(sgp)

        all_snl_ids = [sid for sid in snlgroup.all_snl_ids]
        if canonical_egsnl.snl_id not in all_snl_ids:
            raise ValueError('Canonical SNL must already be in snlgroup to '
                             'switch!')

        new_group = SNLGroup(snlgroup_id, canonical_egsnl, all_snl_ids)
        self.snlgroups.update({'snlgroup_id': snlgroup_id},
                              new_group.as_dict())

    def to_dict(self):
        """
        Note: usernames/passwords are exported as unencrypted Strings!
        """
        return {'host': self.host, 'port': self.port, 'db': self.db,
                'username': self.username, 'password': self.password}

    @classmethod
    def from_dict(cls, d):
        return EGSNLMongoAdapter(d['host'], d['port'], d['db'], d['username'],
                                 d['password'])

    @classmethod
    def auto_load(cls):
        s_dir = os.environ['DB_LOC']
        s_file = os.path.join(s_dir, 'snl_db.yaml')
        return EGSNLMongoAdapter.from_file(s_file)
