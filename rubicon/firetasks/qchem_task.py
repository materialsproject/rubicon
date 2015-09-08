import copy
import logging
import re
import shlex
import os
import shutil
import socket
import datetime
import sys
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QchemJob
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.fw_config import FWData
from fireworks.utilities.fw_serializers import FWSerializable

from custodian.custodian import Custodian
from pymatgen.core.structure import Molecule
from pymatgen.io.qchem import QcInput
from rubicon.utils.eg_wf_utils import move_to_eg_garden
from rubicon.workflows.wf_settings import MOVE_TO_EG_GARDEN

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

    @staticmethod
    def _calibrate_alcf_cmd(input_file="mol.qcinp", max_minutes=60, num_nodes=8, ranks_per_node=1,
                            num_threads=64, scr_size_GB=4, use_runjob=False):
        qc_package_path = "/projects/JCESR/pkcoff/qchem422.mod2"
        qcaux_path = "/projects/JCESR/qcaux"
        qc_exe_path = "/projects/JCESR/pkcoff/public/qcprog-optv1.exe"
        qc_scr_dir = "/dev/local/qchem"
        #scr_size_bytes = scr_size_GB * (2 ** 30)
        scr_size_bytes = scr_size_GB * (10**9)
        qc_envs = {
            "HOME": os.environ["HOME"],
            "QC": qc_package_path,
            "QCAUX": qcaux_path,
            "QCSCRATCH": qc_scr_dir,
            "QCFILEPREF": os.path.join(qc_scr_dir, "qc_job"),
            "QCTMPDIR": qc_scr_dir,
            "QCTHREADS": num_threads,
            "OMP_NUM_THREADS": num_threads,
            "INT_OMP_MIN_LENSTEP": 10,
            "INT_OMP_MAX_LENSTEP": 50,
            "INT_OMP_MIN_LENSTEP_PATH1": 25,
            "INT_OMP_MAX_LENSTEP_PATH1": 100,
            "QCLOCALFSSIZE": scr_size_bytes
        }
        qsub_env_text = ":".join(["{}={}".format(k, v) for k, v in qc_envs.items()])
        runjob_env_text = " ".join(["--envs {}={}".format(k, v) for k, v in qc_envs.items()])
        block_name = os.environ["COBALT_PARTNAME"]
        cur_dir = os.getcwd()
        qc_runjob_cmd = "runjob --block {block_name} -n {num_processes} --ranks-per-node " \
                        "{ranks_per_node} --verbose 2 {envs} --cwd={cur_dir} : {qc_exe_path} {qc_input_file} " \
                        "{scr_dir}".format(
                            block_name=block_name, num_processes=num_nodes * ranks_per_node,
                            ranks_per_node=ranks_per_node, envs=runjob_env_text, cur_dir=cur_dir,
                            qc_exe_path=qc_exe_path, qc_input_file=input_file, scr_dir=qc_scr_dir)
        qc_qsub_cmd = "qsub -A JCESR -t {max_minutes} -n {num_nodes} --mode c{ranks_per_node} --env {envs} " \
                      "--cwd={cur_dir} {qc_exe_path} {qc_input_file} {scr_dir}".format(
                          max_minutes=max_minutes, num_nodes=num_nodes, ranks_per_node=ranks_per_node,
                          envs=qsub_env_text, cur_dir=cur_dir, qc_exe_path=qc_exe_path, qc_input_file=input_file,
                          scr_dir=qc_scr_dir)
        if use_runjob:
            qc_cmd = qc_runjob_cmd
        else:
            qc_cmd = qc_qsub_cmd
        return qc_cmd

    @staticmethod
    def _customize_alcf_qcinp(qcinp, num_nodes=8):
        for qj in qcinp.jobs:
            qj.params["rem"]["parallel_tasks"] = num_nodes
            qj.params["rem"]["BLAS3_DFT"] = 1
            qj.params["rem"]["PDIAG_ON"] = 1
        # use Paul Coffman's version of pathtable
        pathtable_src = "/projects/JCESR/pkcoff/public/pathtable"
        shutil.copy(pathtable_src, "pathtable")

    def run_task(self, fw_spec):
        qcinp = QcInput.from_dict(fw_spec["qcinp"])
        if 'mol' in fw_spec:
            mol = Molecule.from_dict(fw_spec["mol"])
            for qj in qcinp.jobs:
                if isinstance(qj.mol, Molecule):
                    qj.mol = copy.deepcopy(mol)
        mol = qcinp.jobs[0].mol
        if "mixed_basis" in fw_spec:
            for qj in qcinp.jobs:
                if qj.params["rem"]["jobtype"] != "sp":
                    qj.set_basis_set(fw_spec["mixed_basis"])
        if "mixed_aux_basis" in fw_spec:
            for qj in qcinp.jobs:
                if qj.params["rem"]["jobtype"] != "sp":
                    qj.set_aux_basis_set(fw_spec["mixed_aux_basis"])

        carver_name_pattern = re.compile("c[0-9]{4}-ib")
        fw_data = FWData()
        half_cpus_cmd = shlex.split("qchem -np 12")
        if "PBS_JOBID" in os.environ and "edique" in os.environ["PBS_JOBID"]:
            # edison compute nodes
            if (not fw_data.MULTIPROCESSING) or (fw_data.SUB_NPROCS is None):
                qc_exe = shlex.split("qchem -np {}".format(min(24, len(mol))))
                half_cpus_cmd = shlex.split("qchem -np {}".format(
                    min(12, len(mol))))
            else:
                nodelist = ",".join(fw_data.NODE_LIST)
                os.environ["QCNODE"] = nodelist
                qc_exe = shlex.split("qchem -np {}".format(
                    min(fw_data.SUB_NPROCS, len(mol))))
                half_cpus_cmd = shlex.split("qchem -np {}".format(
                    min(fw_data.SUB_NPROCS / 2, len(mol))))
        elif "PBS_JOBID" in os.environ and "hopque" in os.environ["PBS_JOBID"]:
            if (not fw_data.MULTIPROCESSING) or (fw_data.SUB_NPROCS is None):
                qc_exe = shlex.split("qchem -np {}".format(min(24, len(mol))))
                half_cpus_cmd = shlex.split("qchem -np {}".format(
                    min(12, len(mol))))
            else:
                nodelist = ",".join(fw_data.NODE_LIST)
                os.environ["QCNODE"] = nodelist
                qc_exe = shlex.split("qchem -np {}".format(
                    min(fw_data.SUB_NPROCS, len(mol))))
                half_cpus_cmd = shlex.split("qchem -np {}".format(
                    min(fw_data.SUB_NPROCS / 2, len(mol))))
        elif carver_name_pattern.match(socket.gethostname()):
            # mendel compute nodes
            qc_exe = shlex.split("qchem -np {}".format(min(8, len(mol))))
        elif 'vesta' in socket.gethostname():
            # ALCF, Blue Gene
            num_nodes = 8
            qc_exe = shlex.split(self._calibrate_alcf_cmd())
            half_cpus_cmd = shlex.split(self._calibrate_alcf_cmd(num_nodes=num_nodes, scr_size_GB=8))
            self._customize_alcf_qcinp(qcinp, num_nodes=num_nodes)
        elif "macqu" in socket.gethostname().lower():
            qc_exe = shlex.split("qchem -nt 2")
        else:
            qc_exe = ["qchem"]

        logging.basicConfig(level=logging.INFO)
        qchem_logger = logging.getLogger('QChemDrone')
        qchem_logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        qchem_logger.addHandler(sh)

        scf_max_cycles = 200
        geom_max_cycles = 200
        alt_cmd = {"half_cpus": half_cpus_cmd,
                   "openmp": shlex.split("qchem -seq -nt 24")}
        if 'vesta' in socket.gethostname():
            alt_cmd.pop("openmp")
        if fw_spec['num_atoms'] > 50:
            scf_max_cycles = 300
            geom_max_cycles = 500

        qcinp.write_file("mol.qcinp")
        if 'implicit_solvent' in fw_spec and \
                        'solvent_data' in fw_spec['implicit_solvent']:
            solvent_data = fw_spec['implicit_solvent']['solvent_data']
            values = ['{:.4f}'.format(solvent_data[t]) for t in
                      ['Dielec', 'SolN', 'SolA', 'SolB', 'SolG', 'SolC',
                       'SolH']]
            solvent_text = ' '.join(values)
            with open('solvent_data', 'w') as f:
                f.write(solvent_text)
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

        prev_qchem_dir = os.getcwd()
        if MOVE_TO_EG_GARDEN:
            prev_qchem_dir = move_to_eg_garden(prev_qchem_dir)

        stored_data = {'error_list': list(all_errors)}
        update_spec = {'prev_qchem_dir': prev_qchem_dir,
                       'prev_task_type': fw_spec['task_type']}
        propagate_keys = ['egsnl', 'snlgroup_id', 'inchi_root',
                          'mixed_basis', 'mixed_aux_basis', 'mol']
        for k in propagate_keys:
            if k in fw_spec:
                update_spec[k] = fw_spec[k]

        return FWAction(stored_data=stored_data, update_spec=update_spec)
