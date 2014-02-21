import copy
import logging
import re
import shlex
import os
import socket
import datetime
import sys
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QchemJob
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.core.fw_config import FWData
from fireworks.utilities.fw_serializers import FWSerializable

from custodian.custodian import Custodian
from pymatgen.core.structure import Molecule
from pymatgen.io.qchemio import QcInput

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Jun 07, 2013'


DATETIME_HANDLER = lambda obj: obj.isoformat() \
    if isinstance(obj, datetime.datetime) else None


class QChemTask(FireTaskBase, FWSerializable):
    """
    Write QChem task and run QChem
    """

    _fw_name = "QChem Task"

    def run_task(self, fw_spec):
        qcinp = QcInput.from_dict(fw_spec["qcinp"])
        if 'mol' in fw_spec:
            mol = Molecule.from_dict(fw_spec["mol"])
            for qj in qcinp.jobs:
                if isinstance(qj.mol, Molecule):
                    qj.mol = copy.deepcopy(mol)
        mol = qcinp.jobs[0].mol

        carver_name_pattern = re.compile("c[0-9]{4}-ib")
        fw_data = FWData()
        half_cpus_cmd = shlex.split("qchem -np 12")
        if "PBS_JOBID" in os.environ and \
                ("hopque" in os.environ["PBS_JOBID"] or
                         "edique" in os.environ["PBS_JOBID"]):
        # hopper or edison compute nodes
            if (not fw_data.MULTIPROCESSING) or (fw_data.SUB_NPROCS is None):
                qc_exe = shlex.split("qchem -np {}".format(min(24, len(mol))))
                half_cpus_cmd = shlex.split("qchem -np {}".format(
                    min(12, len(mol))))
            else:
                qc_exe = shlex.split("qchem -np {}".format(
                    min(fw_data.SUB_NPROCS, len(mol))))
                half_cpus_cmd = shlex.split("qchem -np {}".format(
                    min(fw_data.SUB_NPROCS/2, len(mol))))
        elif carver_name_pattern.match(socket.gethostname()):
        # mendel compute nodes
            qc_exe = shlex.split("qchem -np {}".format(min(8, len(mol))))
        else:
            qc_exe = ["qchem"]

        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('QChemDrone')
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        logger.addHandler(sh)

        scf_max_cycles = 200
        geom_max_cycles = 200
        alt_cmd = {"half_cpus": half_cpus_cmd,
                   "openmp": shlex.split("qchem -seq -nt 24")}
        if fw_spec['num_atoms'] > 50:
            scf_max_cycles = 300
            geom_max_cycles = 500
            if "PBS_JOBID" in os.environ and \
                    "hopque" in os.environ["PBS_JOBID"]:
                # on Hopper use half core for large molecule
                qc_exe = half_cpus_cmd
            alt_cmd = None

        qcinp.write_file("mol.qcinp")
        job = QchemJob(qc_exe, input_file="mol.qcinp", output_file="mol.qcout",
                       qclog_file="mol.qclog", alt_cmd=alt_cmd, gzipped=True)
        handler = QChemErrorHandler(qchem_job=job,
                                    scf_max_cycles=scf_max_cycles,
                                    geom_max_cycles=geom_max_cycles)
        c = Custodian(handlers=[handler], jobs=[job], max_errors=50)
        custodian_out = c.run()

        all_errors = set()
        for run in custodian_out:
            for correction in run['corrections']:
                all_errors.update(correction['errors'])

        stored_data = {'error_list': list(all_errors)}
        update_spec = {'prev_qchem_dir': os.getcwd(),
                       'prev_task_type': fw_spec['task_type'],
                       'egsnl': fw_spec['egsnl'],
                       'snlgroup_id': fw_spec['snlgroup_id'],
                       'run_tags': fw_spec['run_tags'],
                       'inchi_root': fw_spec['inchi_root']}
        if 'mol' in fw_spec:
            update_spec['mol'] = fw_spec['mol']

        return FWAction(stored_data=stored_data, update_spec=update_spec)
