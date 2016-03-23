# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import json
import shlex

import numpy
from pymongo import MongoClient

# TODO: Refactor PackMol module to python code and merge the ff_dev branch in pymatgen to master
from rubicon.analysis.lammps.MSD import MSD
from rubicon.analysis.lammps.calcCOM import calcCOM
from rubicon.analysis.lammps.gettimedata import gettimedata
from rubicon.analysis.lammps.getmoldata import getmoldata
from rubicon.analysis.lammps.COMradialnofort import COMradialdistribution
from rubicon.analysis.lammps.getatomcharges import getatomcharges
from rubicon.analysis.lammps.calcNEconductivity import calcNEconductivity

from fireworks import FireTaskBase, explicit_serialize

__author__ = 'navnidhirajput'


@explicit_serialize
class ParselammpsProperties(FireTaskBase):
    """
    Parse LAMMPS properties.

    Required params:


    Optional params:

    """

    _fw_name = "Lammps Properties Parser"

    def lampps_properties(self, trjfile, datafile, logfile):
        """
        calculate the MSD and diffusivity for all
        molecules in a system as well as the center of mass radial distribution
        function for all pairs of molecules in the system

        Requires the following comments in the lammps data file starting
        at the third line

        # "number" "molecule" molecules

        where "number" is the number of that molecule type and
        "molecule" is a name for that molecule

        Do not include blank lines in between the molecule types

        Outputs are stored in a dictionary called output to later be stored
        in JSON format
        :return: Output
        """
        c = calcCOM()
        m = MSD()
        gt = gettimedata()
        gm = getmoldata()
        crd = COMradialdistribution()
        gc = getatomcharges()
        ne = calcNEconductivity()

        output = {}
        output['RDF'] = {}
        output['RDF']['units'] = 'unitless and angstroms'
        output['Conductivity'] = {}
        output['Conductivity']['units'] = 'S/m'
        T = 298  # get from lammpsio

        tsjump = gt.getjump(trjfile)
        (nummoltype, moltypel, moltype) = gm.getmoltype(datafile)
        dt = gt.getdt(logfile)
        n = gc.findnumatoms(datafile)
        (molcharges, atomcharges, n) = gc.getmolcharges(datafile, n)
        molcharge = gc.molchargedict(molcharges, moltypel, moltype)
        (comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2) = c.calcCOM([trjfile],
                                                                  datafile)

        output = m.runMSD(comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, moltype,
                          moltypel, dt, tsjump, output)
        output = ne.calcNEconductivity(output, molcharge, Lx, Ly, Lz,
                                       nummoltype, moltypel, T)
        output = crd.runradial(datafile, comx, comy, comz, Lx, Ly, Lz, Lx2,
                               Ly2, Lz2, output, nummoltype, moltypel, moltype,
                               firststep=1)
        return output

    def _insert_doc(self, fw_spec=None, trjfile=None, datafile=None,
                    logfile=None):
        db_dir = shlex.os.environ['DB_LOC']
        db_path = shlex.os.path.join(db_dir, 'tasks_db.json')
        with open(db_path) as f:
            db_creds = json.load(f)
        conn = MongoClient(db_creds['host'], db_creds['port'], )
        db = conn[db_creds['database']]
        if db_creds['admin_user']:
            db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
            coll = db['lammps_properties']
        parse_lammps_prop = self.lampps_properties(trjfile, datafile, logfile)
        docs = parse_lammps_prop
        docs = {k: list(v) if isinstance(v, numpy.ndarray) else v for k, v in
                docs.items()}
        coll.insert(docs)

    def run_task(self, fw_spec):
        mol_traj_file = fw_spec["prev_lammps_trj"]
        mol_data_file = fw_spec["prev_lammps_data"]
        mol_log_file = fw_spec["prev_lammps_log"]
        self._insert_doc(trjfile=mol_traj_file, datafile=mol_data_file,
                         logfile=mol_log_file)
