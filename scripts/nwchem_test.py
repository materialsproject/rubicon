from pymatgen import Molecule
from pymatgen.io.nwchem import NwInput, NwTask

coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]
mol = Molecule(["C", "H", "H", "H", "H"], coords)

tasks = [
    NwTask.dft_task(mol, operation="optimize", xc="b3lyp",
                    basis_set="6-31++G*"),
    NwTask.dft_task(mol, operation="freq", xc="b3lyp",
                    basis_set="6-31++G*")

]
nwi = NwInput(mol, tasks)

print str(nwi)
