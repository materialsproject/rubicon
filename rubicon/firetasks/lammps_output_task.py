

__author__ = 'navnidhirajput'

import shlex
import subprocess
from monty import logging
from pymatgen import Molecule
try:
    # just a walkaround before the packmol is merged to master branch
    # after packmol is merged to master branch, the try...catch block
    # should be removed
    from pymatgen.packmol.packmol import PackmolRunner
    from pymatgen.packmol.lammpsio import LammpsLog
except:
    pass
from rubicon.gff.boxmol import BoxMol
from rubicon.gff.lammpsin import DictLammpsInputSet
from rubicon.gff.lamppsio import LmpInput
from rubicon.gff.antechamberio import AntechamberRunner


__author__ = 'navnidhirajput'


from fireworks import FireTaskBase, explicit_serialize, Firework, Workflow


@explicit_serialize
class WritelammpsOutputTask(FireTaskBase):
    """
    Writes LAMMPS Output files.

    Required params:


    Optional params:

    """

    _fw_name = "Lammps Output Writer"


    def run_task(self, fw_spec):

        # parsing code
        #qcout = QcOutput(zpath(path))

                    #dn update/insertion code
                    # coll.update({"super_mol_snlgroup_id": fw_spec["snlgroup_id"],
                    #      "fragments_def": fw_spec["fragments"]},
                    #     {"$set": d},
                    #     upsert=True)
        llog=LammpsLog()
        print llog.avgs

# task2 = WritelammpsOutputTask()
# coords_n1c=[[4.522,   8.999,   5.512],
#            [6.666,   9.383,   5.971],
#            [7.522,   9.623,   6.578],
#            [6.567,   9.114,   4.645],
#            [7.322,   9.072,   3.879],
#            [5.388,   9.312,   6.496],
#            [5.231,   8.881,   4.373],
#            [5.041,   9.544,   7.906],
#            [5.944,   9.842,   8.432],
#            [4.305,  10.344,   7.985],
#            [4.653,   8.630,   8.355],
#            [4.699,   8.519,   3.03],
#            [3.676,   8.890,   2.982],
#            [5.285,   9.084,   2.312],
#            [4.774,   7.019,   2.765],
#            [4.386,   6.818,   1.765],
#            [5.802,   6.657,   2.807],
#            [4.176,   6.448,   3.480],
#            [3.050,   8.855,   5.663],
#            [2.657,   8.107,   4.974],
#            [2.796,   8.543,   6.676],
#            [2.542,   9.803,   5.463]]
#
#     coords_pc=[[2.516,   1.762,   7.803],
#                [3.382,   2.859,   7.178],
#                [4.673,   2.240,   7.066],
#                [3.476,   3.751,   7.797],
#                [3.037,   3.136,   6.177],
#                [4.526,   0.887,   7.121],
#                [3.246,   0.565,   7.439],
#                [5.398,   0.103,   6.918],
#                [1.090,   1.684,   7.302],
#                [2.539,   1.829,  8.896],
#                [1.069,   1.564,   6.216],
#                [0.552,   2.599,   7.564],
#                [0.568,   0.839,   7.755]]
#
#
#
#     coords_tfn=[[3.734,   8.436,  10.848],
#                [5.713,   9.869,   7.928],
#                [4.712,   9.816,   8.816],
#                [4.981,   8.437,  10.088],
#                [6.222,   8.807,  10.764],
#                [4.651,  11.005,  9.445],
#                [5.009,   7.118,   9.161],
#                [5.914,   5.626,   7.385],
#                [3.564,   9.647,   8.144],
#                [6.318,   6.404,   8.552],
#                [7.545,   7.196,   8.484],
#                [7.774,   4.365,   9.475],
#                [5.668,   4.219,   9.986],
#                [6.692,   5.074,   9.850],
#                [6.947,   5.614,  11.049]]
#     tfn = Molecule(["O", "F", "C", "S", "O", "F", "N", "O", "F", "S", "O", "F", "F", "C", "F"], coords_tfn,site_properties={"mol_name":["TFN"]*len(coords_tfn)})
#     n1c = Molecule(["C", "C", "H", "C", "H", "N", "N", "C", "H", "H", "H", "C", "H", "H", "C", "H", "H", "H", "C", "H", "H", "H"], coords_n1c,site_properties={"mol_name":["N1C"]*len(coords_n1c)})
#     pc = Molecule(["C", "C", "O", "H", "H", "C", "O", "O", "C", "H", "H", "H", "H"], coords_pc,site_properties={"mol_name":["PC"]*len(coords_pc)})
    #
    # fw = Firework([task2], spec={"molecules": [tfn.as_dict(),n1c.as_dict(),pc.as_dict()]})
    # wf = Workflow([fw])
    # print task2.run_task(spec={"molecules": [tfn.as_dict(),n1c.as_dict(),pc.as_dict()]})


test_yong = LammpsLog.from_file('mol.log')
docs = test_yong.llog

db_dir = shlex.os.environ['DB_LOC']
db_path = shlex.os.path.join(db_dir, 'tasks_db.json')

# if db.counter.find({"_id": "mol_taskid"}).count() == 0:
#         db.counter.insert({"_id": "mol_taskid", "c": 1})
#             conn.close()

docs["task_id"] = 348
