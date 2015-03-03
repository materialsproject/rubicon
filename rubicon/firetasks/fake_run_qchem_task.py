import glob
import json
import logging
import os
from fireworks import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import FWSerializable
import sys
from monty.io import zopen
from monty.os.path import zpath
import shutil
from rubicon.utils.eg_wf_utils import move_to_eg_garden
from rubicon.workflows.wf_settings import MOVE_TO_EG_GARDEN

__author__ = 'xiaohuiqu'


class FakeRunQChemTask(FireTaskBase, FWSerializable):
    """
    Read Previous QChem job files, and pretend to be a QChem Task
    """

    _fw_name = "Fake Run QChem Task"

    def run_task(self, fw_spec):
        logging.basicConfig(level=logging.INFO)
        qchem_logger = logging.getLogger('QChemDrone')
        qchem_logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        #sh.setLevel(getattr(logging, 'INFO'))
        qchem_logger.addHandler(sh)

        cur_dir = os.getcwd()
        src_qchem_dir = fw_spec['src_qchem_dir']
        for filename in glob.glob(os.path.join(src_qchem_dir, '*')):
            if os.path.isfile(filename):
                shutil.copy(filename, cur_dir)

        if os.path.exists("custodian.json") or os.path.exists("custodian.json"+".gz"):
            with zopen(zpath("custodian.json")) as f:
                custodian_out = json.load(f)
        else:
            custodian_out = []

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