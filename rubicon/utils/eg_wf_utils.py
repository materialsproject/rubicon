import json
import os
import shutil
from rubicon.workflows.wf_settings import EG_RUN_LOCS

__author__ = 'xiaohuiqu'


def get_eg_block_part(m_dir):
    grandpa_dir, parent_dir = os.path.split(m_dir)[-2:]
    return grandpa_dir, parent_dir


def get_eg_dir_loc(m_dir):
    if os.path.exists(m_dir):
        return m_dir

    all_locs = []
    if "GARDEN_LOC" in os.environ:
        all_locs.append(os.environ["GARDEN_LOC"])
    for scr in ["GSCRATCH", "SCRATCH", "SCRATCH2"]:
        for loc in EG_RUN_LOCS:
            all_locs.append(os.path.join(scr, loc))

    grandpa_dir, parent_dir = get_eg_block_part(m_dir)

    for preamble in all_locs:
        new_loc = os.path.join(preamble, grandpa_dir, parent_dir)
        if os.path.exists(new_loc):
            return new_loc

    raise ValueError('get_loc() -- dir does not exist!! Make sure your base '
                     'directory is listed in RUN_LOCS of wf_settings.py\n'
                     'and your garden location is set in the environment'
                     ' variable GARDEN_LOC')


def get_eg_file_loc(m_file):
    m_dir = os.path.dirname(m_file)
    filename = os.path.basename(m_file)
    dir_loc = get_eg_dir_loc(m_dir)
    return os.path.join(dir_loc, filename)


def move_to_eg_garden(m_dir):
    if "GARDEN_LOC" not in os.environ:
        return m_dir
    garden_part = os.path.abspath(os.environ["GARDEN_LOC"])
    if os.path.exists(m_dir) or os.path.exists(m_dir + ".gz"):
        if not os.path.isdir(m_dir):
            m_dir = m_dir
    else:
        raise ValueError("The folder \"{}\" doesn't exist".format(m_dir))
    grandpa_dir, parent_dir = get_eg_block_part(m_dir)
    dest_dir = os.path.join(garden_part, grandpa_dir, parent_dir)
    dest_parent_dir = os.path.join(garden_part, grandpa_dir)
    if not os.path.exists(dest_parent_dir):
        # noinspection PyBroadException
        try:
            os.makedirs(dest_parent_dir)
        except:
            if not os.path.exists(dest_parent_dir):
                raise ValueError("Couldn't create parent folder \"{}\" in "
                                 "garden".format(dest_parent_dir))
    if os.path.exists(dest_dir):
        raise ValueError("A same name folder \"{}\" already exists in garden"
                         .format(dest_dir))
    shutil.move(m_dir, dest_dir)
    return dest_dir


def get_defuse_causing_qchem_fwid(qcout_path):
    dirname = os.path.dirname(qcout_path)
    fw_spec_path = os.path.join(dirname, "FW.json")
    with zopen(zpath(fw_spec_path)) as f:
        fw_dict = json.load(f)
    fw_id = fw_dict["fw_id"]
    return fw_id