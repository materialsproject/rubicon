import json

__author__ = 'xiaohuiqu'

import copy
import logging
import re
import shlex
import os
import socket
from custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QchemJob
import math
from pymatgen import Molecule
from pymatgen.io.qchemio import QcInput, QcOutput
import sys

__author__ = 'xiaohuiqu'


def run_qchem(filename):
    base_filename = os.path.splitext(filename)[0]
    output_filename = base_filename + ".qcout"
    log_filename = base_filename + ".qclog"
    qcinp = QcInput.from_file(filename)
    mol = qcinp.jobs[0].mol
    carver_name_pattern = re.compile("c[0-9]{4}-ib")
    half_cpus_cmd = None
    openmp_cmd = None
    if "PBS_JOBID" in os.environ and \
            ("hopque" in os.environ["PBS_JOBID"] or
                     "edique" in os.environ["PBS_JOBID"]):
        # hopper or edison compute nodes
        qc_exe = shlex.split("qchem -np {}".format(min(24, len(mol))))
        half_cpus_cmd = shlex.split("qchem -np {}".format(min(12, len(mol))))
        openmp_cmd = shlex.split("qchem -seq -nt 24")
    elif carver_name_pattern.match(socket.gethostname()):
        # mendel compute nodes
        qc_exe = shlex.split("qchem -np {}".format(min(8, len(mol))))
        half_cpus_cmd = shlex.split("qchem -np {}".format(min(4, len(mol))))
        openmp_cmd = shlex.split("qchem -seq -nt 8")
    else:
        qc_exe = ["qchem"]

    logging.basicConfig(level=logging.INFO)
    qchem_logger = logging.getLogger('QChemIndependentRun')
    qchem_logger.setLevel(logging.INFO)
    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setLevel(getattr(logging, 'INFO'))
    qchem_logger.addHandler(sh)

    scf_max_cycles = 200
    geom_max_cycles = 200
    alt_cmd = {"half_cpus": half_cpus_cmd,
               "openmp": openmp_cmd}
    if len(mol) > 50:
        scf_max_cycles = 300
        geom_max_cycles = 500

    job = QchemJob(qc_exe, input_file=filename, output_file=output_filename,
                   qclog_file=log_filename, alt_cmd=alt_cmd, gzipped=False)
    handler = QChemErrorHandler(qchem_job=job,
                                scf_max_cycles=scf_max_cycles,
                                geom_max_cycles=geom_max_cycles,
                                input_file=filename,
                                output_file=output_filename)
    c = Custodian(handlers=[handler], jobs=[job], max_errors=50)
    c.run()

def write_smx_solvent_data(solvent):
    smx_data_file = os.path.join(os.path.dirname(__file__),
                                 "../utils/data", "smx_data.json")
    with open(smx_data_file) as f:
        smx_data = json.load(f)
    if solvent not in smx_data["builtin_solvent"]:
        if solvent not in smx_data["custom_solvent"]:
            raise Exception("Don't know the SMx parameters for "
                            "solvent '{}'".format(solvent))
        solvent_data = smx_data["custom_solvent"][solvent]
    values = ['{:.4f}'.format(solvent_data[t]) for t in
              ['Dielec', 'SolN', 'SolA', 'SolB', 'SolG', 'SolC',
               'SolH']]
    solvent_text = ' '.join(values)
    with open('solvent_data', 'w') as f:
        f.write(solvent_text)

def perturb_molecule(old_mol, vib_mode, reversed_direction=False, perturb_scale=0.3):
    max_dis = max([math.sqrt(sum([x ** 2 for x in mode]))
                   for mode in vib_mode])
    scale = perturb_scale / max_dis
    normalized_mode = [[x * scale for x in mode]
                       for mode in vib_mode]
    direction = 1.0
    if reversed_direction:
        direction = -1.0
    new_coords = [[c+v*direction for c, v in zip(site.coords, mode)]
                  for site, mode in zip(old_mol.sites, normalized_mode)]
    species = [site.specie.symbol
               for site in old_mol.sites]
    charge = old_mol.charge
    spin_multiplicity = old_mol.spin_multiplicity
    new_mol = Molecule(species, new_coords, charge=charge,
                       spin_multiplicity=spin_multiplicity)
    return new_mol

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Run A QChem Job for a QChem Input File")
    parser.add_argument("-f", "--filename", dest="filename", type=str,
                        required=True,
                        help="the QChem input filename")
    parser.add_argument("-e", "--eliminate", dest="eli_img", action="store_true",
                        help="whether to eliminate imaginary frequency")
    parser.add_argument("-s", "--solvent", dest="solvent", type=str,
                        required=False,
                        help="the implicit solvent")
    options = parser.parse_args()
    if options.solvent:
        write_smx_solvent_data(options.solvent)
    run_qchem(options.filename)
    if options.eli_img:
        base_filename = os.path.splitext(options.filename)[0]
        output_filename = base_filename + ".qcout"
        qcout = QcOutput(output_filename)
        if qcout.data[-1]['frequencies'][0]["frequency"] < -0.00:
            os.system("tar czvf img_freq_1.tar.gz *")
            old_mol = qcout.data[-1]["molecules"][-1]
            vib_mode = ['frequencies'][0]["vib_mode"]
            new_mol = perturb_molecule(old_mol, vib_mode)
            qctask_freq = qcout.data[-1]["input"]
            qctask_freq.mol = "read"
            qctask_opt = copy.deepcopy(qctask_freq)
            qctask_opt.params["rem"]["jobtype"] = "opt"
            qctask_opt.mol = new_mol
            qcinp = QcInput([qctask_opt, qctask_freq])
            eli_file_1 = base_filename + "_eli_img_1.qcinp"
            qcinp.write_file(eli_file_1)
            run_qchem(options.filename)

            output_filename = base_filename + "_eli_img_1.qcout"
            qcout = QcOutput(output_filename)
            if qcout.data[-1]['frequencies'][0]["frequency"] < -0.00:
                os.system("tar czvf img_freq_2.tar.gz *")
                old_mol = qcout.data[-1]["molecules"][-1]
                vib_mode = ['frequencies'][0]["vib_mode"]
                new_mol = perturb_molecule(old_mol, vib_mode)
                qctask_freq = qcout.data[-1]["input"]
                qctask_freq.mol = "read"
                qctask_opt = copy.deepcopy(qctask_freq)
                qctask_opt.params["rem"]["jobtype"] = "opt"
                qctask_opt.mol = new_mol
                qcinp = QcInput([qctask_opt, qctask_freq])
                for j in qcinp.jobs:
                    j.set_dft_grid(128, 302)
                    if j.params["rem"]["jobtype"] == "opt":
                        j.scale_geom_opt_threshold(0.1, 0.1, 0.1)
                        j.set_geom_max_iterations(100)
                eli_file_2 = base_filename + "_eli_img_2.qcinp"
                qcinp.write_file(eli_file_2)
                run_qchem(options.filename)

            output_filename = base_filename + "_eli_img_2.qcout"
            qcout = QcOutput(output_filename)
            if qcout.data[-1]['frequencies'][0]["frequency"] < -0.00:
                os.system("tar czvf img_freq_3.tar.gz *")
                old_mol = qcout.data[-1]["molecules"][-1]
                vib_mode = ['frequencies'][0]["vib_mode"]
                new_mol = perturb_molecule(old_mol, vib_mode)
                qctask_freq = qcout.data[-1]["input"]
                qctask_freq.mol = "read"
                qctask_opt = copy.deepcopy(qctask_freq)
                qctask_opt.params["rem"]["jobtype"] = "opt"
                qctask_opt.mol = new_mol
                qcinp = QcInput([qctask_opt, qctask_freq])
                for j in qcinp.jobs:
                    j.set_dft_grid(90, 590)
                    if j.params["rem"]["jobtype"] == "opt":
                        j.scale_geom_opt_threshold(0.1, 0.1, 0.1)
                        j.set_geom_max_iterations(100)
                eli_file_3 = base_filename + "_eli_img_3.qcinp"
                qcinp.write_file(eli_file_3)
                run_qchem(options.filename)

if __name__ == '__main__':
    main()

