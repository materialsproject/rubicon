from fireworks import FireWork, Workflow, LaunchPad
from pymatgen import Molecule
from rubicon.firetasks.lammps_task import WritelammpsInputTask

__author__ = 'navnidhirajput'


if __name__ == '__main__':
    task1 = WritelammpsInputTask()
    #task2 = LamppsCustodianTask()
    # task3 = VaspAnalyzeTask()
    coords_pc=[[2.516,   1.762,   7.803],
               [3.382,   2.859,   7.178],
               [4.673,   2.240,   7.066],
               [3.476,   3.751,   7.797],
               [3.037,   3.136,   6.177],
               [4.526,   0.887,   7.121],
               [3.246,   0.565,   7.439],
               [5.398,   0.103,   6.918],
               [1.090,   1.684,   7.302],
               [2.539,   1.829,  8.896],
               [1.069,   1.564,   6.216],
               [0.552,   2.599,   7.564],
               [0.568,   0.839,   7.755]]

    pc = Molecule(["C", "C", "O", "H", "H", "C", "O", "O", "C", "H", "H", "H", "H"], coords_pc)

    fw = FireWork([task1], spec={"molecules": [pc.to_dict]})
    wf = Workflow([fw])

    lp = LaunchPad.auto_load()
    lp.add_wf(wf)

