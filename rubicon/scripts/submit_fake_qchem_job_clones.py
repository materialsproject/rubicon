# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import glob
import json
import os

from monty.io import zopen
from monty.os.path import zpath

from fireworks import Firework, Workflow, LaunchPad
from fireworks.utilities.fw_utilities import get_slug
from pymatgen import Molecule
from pymatgen.matproj.snl import StructureNL
from rubicon.firetasks.qchem.egsnl_tasks import AddEGSNLTask
from rubicon.firetasks.qchem.fake_run_qchem_task import FakeRunQChemTask

__author__ = 'xiaohuiqu'


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Run A QChem Job for a QChem Input File")
    parser.add_argument("-d", "--directory", dest="directory", type=str,
                        required=True,
                        help="the directory contains all the QChem jobs to be pretended to run again")
    parser.add_argument("-p", "--priority", dest="priority", type=int,
                        default=100,
                        help="the FireWorks priority")
    parser.add_argument("-b", "--batch_size", dest="batch_size", type=int,
                        default=100,
                        help="the number of FireWorks in a Workflow")
    options = parser.parse_args()

    fw_priority = options.priority
    batch_size = options.batch_size

    lp = LaunchPad.auto_load()

    src_dir = os.path.abspath(options.directory)
    src_dir_sub_dirs = glob.glob(os.path.join(src_dir, "*"))
    num_dirs = len(src_dir_sub_dirs)
    current_fwid = 1
    links_dict = dict()
    fws_all = []
    num_fw_in_current_batch = 0
    batch_num = 1
    for i, sd in enumerate(src_dir_sub_dirs):
        if not os.path.isdir(sd):
            continue
        fw_json_filename = os.path.join(sd, "FW.json")
        if not (os.path.exists(fw_json_filename) or os.path.exists(
                    fw_json_filename + ".gz")):
            continue
        with zopen(zpath(fw_json_filename)) as f:
            fw_dict = json.load(f)
        print("{percent:4.2%} completed, processing directory {d:s}, "
              "molecule name {molname:s}," \
              " mission {mission:s}".format(percent=i / float(num_dirs), d=sd,
                                            molname=
                                            fw_dict['spec']['user_tags'][
                                                'molname'],
                                            mission=
                                            fw_dict['spec']['user_tags'][
                                                'mission']))

        molname = fw_dict['spec']['user_tags']['molname']
        egsnl_tasks = [AddEGSNLTask()]
        if 'mol' in fw_dict:
            mol = Molecule.from_dict(fw_dict['spec']['mol'])
        else:
            mol = Molecule.from_dict(
                fw_dict['spec']['qcinp']['jobs'][0]['molecule'])
        snl = StructureNL(mol, "Xiaohui Qu <xqu@lbl.gov>",
                          "Electrolyte Genome")
        egsnl_task_spec = {'task_type': 'Add to SNL database',
                           'snl': snl.as_dict(),
                           '_category': 'Parse Previous QChem Job',
                           '_priority': fw_priority}
        snl_fw_id = current_fwid
        current_fwid += 1
        fws_all.append(Firework(egsnl_tasks, egsnl_task_spec,
                                name=get_slug(
                                    molname + ' -- Add to SNL database For fake QChem Task'),
                                fw_id=snl_fw_id))

        fake_qchem_tasks = [FakeRunQChemTask()]
        src_qchem_dir = sd
        fake_qchem_spec = {'_priority': fw_priority * 2,
                           'src_qchem_dir': src_qchem_dir,
                           '_category': 'Parse Previous QChem Job',
                           'run_tags': fw_dict['spec']['run_tags'],
                           'implicit_solvent': fw_dict['spec'][
                               'implicit_solvent'],
                           'task_type': fw_dict['spec']['task_type'],
                           'charge': fw_dict['spec']['charge'],
                           'spin_multiplicity': fw_dict['spec'][
                               'spin_multiplicity'],
                           'num_atoms': fw_dict['spec']['num_atoms'],
                           'user_tags': fw_dict['spec']['user_tags'],
                           'mol': mol.as_dict(),
                           'inchi': fw_dict['spec']['inchi'],
                           '_dupefinder': fw_dict['spec']['_dupefinder'],
                           'qcinp': fw_dict['spec']['qcinp'],
                           'qm_method': fw_dict['spec']['qm_method'],
                           'inchi_root': fw_dict['spec']['inchi_root']}
        for k in ['mixed_basis', 'mixed_aux_basis']:
            if k in fw_dict['spec']:
                fake_qchem_spec[k] = fw_dict['spec'][k]
        fake_qchem_fw_id = current_fwid
        current_fwid += 1
        fws_all.append(Firework(fake_qchem_tasks, fake_qchem_spec,
                                name='Fake' + fw_dict['name'],
                                fw_id=fake_qchem_fw_id))
        links_dict[snl_fw_id] = fake_qchem_fw_id

        num_fw_in_current_batch += 1
        if num_fw_in_current_batch >= 100:
            wf = Workflow(fws_all, links_dict,
                          "Read Previous QChem Jobs Id-{}".format(batch_num))
            lp.add_wf(wf)
            batch_num += 1
            links_dict = dict()
            fws_all = []
            num_fw_in_current_batch = 0

    if num_fw_in_current_batch > 0:
        wf = Workflow(fws_all, links_dict, "Read Previous QChem Jobs")
        lp.add_wf(wf)


if __name__ == '__main__':
    main()
