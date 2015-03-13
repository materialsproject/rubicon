from fireworks import  Workflow, LaunchPad
from fireworks.core.firework import Firework
from pymatgen import Molecule
from rubicon.firetasks.lammps_input_task import WritelammpsInputTask
from rubicon.firetasks.lammps_output_task import WritelammpsOutputTask

import glob
from pymatgen import Molecule

__author__ = 'navnidhirajput'


if __name__ == '__main__':
    task1 = WritelammpsInputTask()
    #task2 = WritelammpsOutputTask()


    coords_n1c=[[4.522,   8.999,   5.512],
           [6.666,   9.383,   5.971],
           [7.522,   9.623,   6.578],
           [6.567,   9.114,   4.645],
           [7.322,   9.072,   3.879],
           [5.388,   9.312,   6.496],
           [5.231,   8.881,   4.373],
           [5.041,   9.544,   7.906],
           [5.944,   9.842,   8.432],
           [4.305,  10.344,   7.985],
           [4.653,   8.630,   8.355],
           [4.699,   8.519,   3.03],
           [3.676,   8.890,   2.982],
           [5.285,   9.084,   2.312],
           [4.774,   7.019,   2.765],
           [4.386,   6.818,   1.765],
           [5.802,   6.657,   2.807],
           [4.176,   6.448,   3.480],
           [3.050,   8.855,   5.663],
           [2.657,   8.107,   4.974],
           [2.796,   8.543,   6.676],
           [2.542,   9.803,   5.463]]

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



    coords_tfn=[[3.734,   8.436,  10.848],
               [5.713,   9.869,   7.928],
               [4.712,   9.816,   8.816],
               [4.981,   8.437,  10.088],
               [6.222,   8.807,  10.764],
               [4.651,  11.005,  9.445],
               [5.009,   7.118,   9.161],
               [5.914,   5.626,   7.385],
               [3.564,   9.647,   8.144],
               [6.318,   6.404,   8.552],
               [7.545,   7.196,   8.484],
               [7.774,   4.365,   9.475],
               [5.668,   4.219,   9.986],
               [6.692,   5.074,   9.850],
               [6.947,   5.614,  11.049]]

    coords_acn=[[-4.452,  -0.223,  -0.000],
               [-3.389,   0.780,   0.000],
               [-4.367,  -0.861,   0.889],
               [-4.375,  -0.852,  -0.896],
               [-5.433,   0.268,   0.007],
               [-2.545,   1.577,   0.000]]

    coords_mg=[[0.000,  0.000,  0.000]]
    coords_al=[[0.000,  0.000,  0.000]]
    coords_cl=[[0.000,  0.000,  0.000]]
    coords_thf=[[-1.088,   2.692  -1.000],
                [-0.728,  -0.308,  -1.000],
                [0.692,  -0.308,   4.000],
                [-1.378,   1.692,  -3.000],
                [-0.748,  -0.308,   3.000],
                [1.122,   0.692,  -2.000],
                [0.712,   0.692,  -3.000],
                [1.302,  -3.308,  -3.000],
                [0.032,  -3.308,   0.000],
                [2.042,  -2.308,   3.000],
                [1.312,   2.692,   2.000],
                [-1.998,  -0.308,  -2.000],
                [-1.278,   1.692,   3.000]]




    coords_tfsi=[[8.570,  13.900,  10.430],
            [ 8.020,  15.350,  10.640],
            [ 6.840,  15.410,  11.470],
            [ 9.020,  16.380,  10.790],
            [ 7.330,  15.680,   8.940],
            [ 6.600,  14.670,   8.470],
            [ 8.280,  15.890,   8.030],
            [ 6.540,  16.750,   8.900],
            [ 9.440,  13.090,  11.410],
            [ 9.330,  11.670,  11.100],
            [ 9.350,  13.420,  12.810],
            [11.180,  13.520,  10.920],
            [11.490,  13.050,   9.710],
            [11.450,  14.810,  10.890],
            [12.080,  12.980,  11.740]]

    coords_diglyme=[[  1.023,   0.136,  -0.063],
                [  0.506,  1.443,   0.126],
                [  2.441,   0.089,  -0.070],
                [ -1.013,   1.345,   0.112],
                [  2.854,  -1.361,  -0.280],
                [ -1.529,   2.651,   0.301],
                [  4.270,  -1.407,  -0.287],
                [ -2.943,   2.693,   0.307],
                [  4.781,  -2.713,  -0.476],
                [  0.846,   2.116,  -0.677],
                [  0.846,   1.860,   1.087],
                [  2.847,   0.460,   0.884],
                [  2.847,   0.715,  -0.880],
                [ -1.353,  0.671,   0.915],
                [ -1.353,   0.927,  -0.849],
                [ 2.447,  -1.733,  -1.235],
                [ 2.447,  -1.988,   0.530],
                [-3.232,   3.736,   0.459],
                [-3.359,  2.080,   1.123 ],
                [-3.359,   2.337,  -0.649 ],
                [5.871,  -2.635,  -0.465 ],
                [4.458,  -3.133,  -1.442 ],
                [4.459,  -3.390,   0.331 ]]




    tfn = Molecule(["O", "F", "C", "S", "O", "F", "N", "O", "F", "S", "O", "F", "F", "C", "F"], coords_tfn,site_properties={"mol_name":["TFN"]*len(coords_tfn)})
    n1c = Molecule(["C", "C", "H", "C", "H", "N", "N", "C", "H", "H", "H", "C", "H", "H", "C", "H", "H", "H", "C", "H", "H", "H"], coords_n1c,site_properties={"mol_name":["N1C"]*len(coords_n1c)})
    pc = Molecule(["C", "C", "O", "H", "H", "C", "O", "O", "C", "H", "H", "H", "H"], coords_pc,site_properties={"mol_name":["PC"]*len(coords_pc)})
    acn = Molecule(["C", "C", "H", "H", "H", "N"], coords_acn,site_properties={"mol_name":["ACN"]*len(coords_acn)})
    mg = Molecule(["Mg"], coords_mg,site_properties={"mol_name":["MAG"]*len(coords_mg)})
    al = Molecule(["Al"], coords_al,site_properties={"mol_name":["ALU"]*len(coords_al)})
    cl = Molecule(["Cl"], coords_cl,site_properties={"mol_name":["CHL"]*len(coords_cl)})
    thf = Molecule(["C","C","C","H","H","C","H","H","O","H","H","H","H"], coords_thf,site_properties={"mol_name":["THF"]*len(coords_thf)})
    tfsi = Molecule(["N","S","O","O","C","F","F","F","S","O","O","C","F","F","F"], coords_tfsi,site_properties={"mol_name":["tfsi"]*len(coords_tfsi)})
    diglyme = Molecule(["O","C","C","C","C","O","O","C","C","H","H","H","H","H","H","H","H","H","H","H","H","H","H"], coords_diglyme,site_properties={"mol_name":["diglyme"]*len(coords_diglyme)})



    #fw = Firework([task1], spec={"molecules": [tfn.as_dict(),n1c.as_dict(),pc.as_dict()]})
    #fw = Firework([task1], spec={"molecules": [mg.as_dict(),tfsi.as_dict(),diglyme.as_dict()]})
    #fw = Firework([task1], spec={"molecules": [diglyme.as_dict()]})
    filelist = glob.glob("/Users/navnidhirajput/Dropbox/solvent_molecules/*")
    for mol in filelist:
        mol = Molecule.from_file(mol)

        #fw = Firework([task1], spec={"molecules": [mg.as_dict()]})
        fw = Firework([task1], spec={"molecules": [mol.as_dict()]})

        wf = Workflow([fw])

    lp = LaunchPad.auto_load()
    lp.add_wf(wf)

