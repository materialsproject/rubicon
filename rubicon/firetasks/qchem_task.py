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
        qcinp.write_file("mol.qcinp")
        hopper_name_pattern = re.compile("nid\d+")
        carver_name_pattern = re.compile("c[0-9]{4}-ib")
        fw_data = FWData()
        if hopper_name_pattern.match(socket.gethostname()):
        # hopper compute nodes
            if (not fw_data.MULTIPROCESSING) or (fw_data.SUB_NPROCS is None):
                qc_exe = shlex.split("qchem -np {}".format(min(24, len(mol))))
            else:
                qc_exe = shlex.split("qchem -np {}".format(
                    min(fw_data.SUB_NPROCS, len(mol))))
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

        job = QchemJob(qc_exe, input_file="mol.qcinp", output_file="mol.qcout",
                       qclog_file="mol.qclog",
                       alt_cmd={"half_cpus": shlex.split("qchem -np 12"),
                                "openmp": shlex.split("qchem -nt 24")})
        handler = QChemErrorHandler(qchem_job=job)
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
