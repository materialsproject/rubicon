# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import copy
import json
import logging
import math
import os
import re
import shlex
import socket
import sys

from custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QchemJob
from six.moves import zip

from pymatgen import Molecule
from pymatgen.io.qchem import QcInput, QcOutput

from rubicon.firetasks.qchem.qchem_task import QChemTask
from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator

__author__ = 'xiaohuiqu'


def call_qchem_task(filename, solvent=None, mixed_basis=None, mixed_aux_basis=None):
    base_filename = os.path.splitext(filename)[0]
    output_filename = base_filename + ".qcout"
    log_filename = base_filename + ".qclog"
    qcinp = QcInput.from_file(filename)
    solvent_token = set_solvent_data(qcinp, solvent)
    QChemTask.run_qchem(qcinp, solvent_token, mixed_aux_basis, mixed_basis)

def set_solvent_data(qcinp, solvent, vdw_surface=True):
    if not solvent:
        return None
    mol = qcinp.jobs[0].mol
    qctask_creator = QChemFireWorkCreator(mol, molname=mol.formula,
                                          mission="Standalone calculation")
    solvent_token = None
    for qctask in qcinp.jobs:
        if "solvent_method" in qctask.params["rem"]:
            # Is a solution phase calculation
            solvent_method = qctask.params["rem"].get("solvent_method")
            solvent_token = qctask_creator.set_solvent_method(
                qctask, solvent, solvent_method, use_vdw_surface=vdw_surface)
    return solvent_token



def perturb_molecule(old_mol, vib_mode, reversed_direction=False,
                     perturb_scale=0.3):
    max_dis = max([math.sqrt(sum([x ** 2 for x in mode]))
                   for mode in vib_mode])
    scale = perturb_scale / max_dis
    normalized_mode = [[x * scale for x in mode]
                       for mode in vib_mode]
    direction = 1.0
    if reversed_direction:
        direction = -1.0
    new_coords = [[c + v * direction for c, v in zip(site.coords, mode)]
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
    parser.add_argument("-e", "--eliminate", dest="eli_img",
                        action="store_true",
                        help="whether to eliminate imaginary frequency")
    parser.add_argument("-b", "--mixed_basis", dest="mixed_basis", type=json.loads,
                        required=False,
                        help="The mixed basis as a dict")
    parser.add_argument("-a", "--mixed_aux_basis", dest="mixed_aux_basis", type=json.loads,
                        required=False,
                        help="The mixed auxiliary basis as a dict")
    parser.add_argument("-s", "--solvent", dest="solvent", type=str,
                        required=False,
                        help="the implicit solvent")
    options = parser.parse_args()
    call_qchem_task(options.filename, options.solvent, options.mixed_basis, options.mixed_aux_basis)
    if options.eli_img:
        base_filename = os.path.splitext(options.filename)[0]
        output_filename = base_filename + ".qcout"
        qcout = QcOutput(output_filename)
        charge = None
        spin_multiplicity = None
        for d in qcout.data:
            j = d["input"]
            if j.charge is not None:
                charge = j.charge
            if j.spin_multiplicity is not None:
                spin_multiplicity = j.spin_multiplicity
        if qcout.data[-1]['frequencies'][0]["frequency"] < -0.00:
            os.system("tar czvf img_freq_1.tar.gz *")
            old_mol = qcout.data[-1]["molecules"][-1]
            vib_mode = qcout.data[-1]['frequencies'][0]["vib_mode"]
            new_mol = perturb_molecule(old_mol, vib_mode)
            qctask_freq = qcout.data[-1]["input"]
            qctask_freq.mol = "read"
            qctask_freq.charge = charge
            qctask_freq.spin_multiplicity = spin_multiplicity
            qctask_opt = copy.deepcopy(qctask_freq)
            qctask_opt.params["rem"]["jobtype"] = "opt"
            qctask_opt.params["rem"].pop("scf_guess", None)
            qctask_opt.mol = new_mol
            qcinp = QcInput([qctask_opt, qctask_freq])
            eli_file_1 = base_filename + "_eli_img_1.qcinp"
            qcinp.write_file(eli_file_1)
            call_qchem_task(eli_file_1)

            output_filename = base_filename + "_eli_img_1.qcout"
            qcout = QcOutput(output_filename)
            if qcout.data[-1]['frequencies'][0]["frequency"] < -0.00:
                os.system("tar czvf img_freq_2.tar.gz *")
                old_mol = qcout.data[-1]["molecules"][-1]
                vib_mode = qcout.data[-1]['frequencies'][0]["vib_mode"]
                new_mol = perturb_molecule(old_mol, vib_mode)
                qctask_freq = qcout.data[-1]["input"]
                qctask_freq.mol = "read"
                qctask_freq.charge = charge
                qctask_freq.spin_multiplicity = spin_multiplicity
                qctask_opt = copy.deepcopy(qctask_freq)
                qctask_opt.params["rem"]["jobtype"] = "opt"
                qctask_opt.params["rem"].pop("scf_guess", None)
                qctask_opt.mol = new_mol
                qcinp = QcInput([qctask_opt, qctask_freq])
                for j in qcinp.jobs:
                    j.set_dft_grid(128, 302)
                    j.set_integral_threshold(12)
                    if j.params["rem"]["jobtype"] == "opt":
                        j.scale_geom_opt_threshold(0.1, 0.1, 0.1)
                        j.set_geom_max_iterations(100)
                eli_file_2 = base_filename + "_eli_img_2.qcinp"
                qcinp.write_file(eli_file_2)
                call_qchem_task(eli_file_2)

                output_filename = base_filename + "_eli_img_2.qcout"
                qcout = QcOutput(output_filename)
                if qcout.data[-1]['frequencies'][0]["frequency"] < -0.00:
                    os.system("tar czvf img_freq_3.tar.gz *")
                    old_mol = qcout.data[-1]["molecules"][-1]
                    vib_mode = qcout.data[-1]['frequencies'][0]["vib_mode"]
                    new_mol = perturb_molecule(old_mol, vib_mode)
                    qctask_freq = qcout.data[-1]["input"]
                    qctask_freq.mol = "read"
                    qctask_freq.charge = charge
                    qctask_freq.spin_multiplicity = spin_multiplicity
                    qctask_opt = copy.deepcopy(qctask_freq)
                    qctask_opt.params["rem"]["jobtype"] = "opt"
                    qctask_opt.params["rem"].pop("scf_guess", None)
                    qctask_opt.mol = new_mol
                    qcinp = QcInput([qctask_opt, qctask_freq])
                    for j in qcinp.jobs:
                        j.set_dft_grid(90, 590)
                        j.set_integral_threshold(12)
                        if j.params["rem"]["jobtype"] == "opt":
                            j.scale_geom_opt_threshold(0.1, 0.1, 0.1)
                            j.set_geom_max_iterations(100)
                    eli_file_3 = base_filename + "_eli_img_3.qcinp"
                    qcinp.write_file(eli_file_3)
                    call_qchem_task(eli_file_3, options.solvent, options.mixed_basis, options.mixed_aux_basis)


if __name__ == '__main__':
    main()
