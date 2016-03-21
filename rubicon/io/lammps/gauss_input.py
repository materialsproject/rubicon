__author__ = 'navnidhirajput'

from pymatgen.core import Molecule

coords_mg = [[-1.695, 0.443, 0.000]]
coords_thf = [[-3.682, 1.377, -0.211],
              [-3.304, -0.066, 0.141],
              [-1.806, -0.095, -0.204],
              [-3.893, -0.804, -0.419],
              [-3.452, -0.252, 1.215],
              [-1.362, 1.314, 0.202],
              [-1.666, -0.247, -1.284],
              [-1.254, -0.881, 0.327],
              [-2.501, 2.169, 0.007],
              [-0.532, 1.700, -0.406],
              [-1.063, 1.350, 1.264],
              [-4.490, 1.780, 0.416],
              [-3.984, 1.466, -1.268]]

mg = Molecule(["Mg"], coords_mg,
              site_properties={"mol_name": ["MAG"] * len(coords_mg)})
thf = Molecule(
    ["C", "C", "C", "H", "H", "C", "H", "H", "O", "H", "H", "H", "H"],
    coords_thf, site_properties={"mol_name": ["THF"] * len(coords_thf)})




class GaussianInput():
    def gaus_input(self, mol, multiplicity=None, charge=None):
        mol_name = mol.site_properties['mol_name'][0]
        gaus_lines = []
        gaus_lines.append(
            '{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}'.format('--Link1--', '\n',
                                                      '%mem=256MW', '\n',
                                                      '%NProcShared=4', '\n',
                                                      '%LindaWorker=localhost',
                                                      '\n', '%chk=', mol_name,
                                                      '.chk'
                                                      , '\n',
                                                      '# b3lyp/aug-cc-pvdz opt=(calcfc,tight) int=ultrafine',
                                                      '\n',
                                                      '# SCF=tight nosymm test',
                                                      '\n', ))

        gaus_lines.append(
            '{}{}{}{}{}'.format('\n', 'created by gaus_input.py from ',
                                mol_name, '.pdb', '\n'))
        gaus_lines.append('{}{}   {}'.format('\n', charge, multiplicity))

        for k, v in enumerate(mol.cart_coords):
            gaus_lines.append(
                '{}     {}   {}    {}    {}'.format('\n', mol.species[k],
                                                    mol.cart_coords[k][0],
                                                    mol.cart_coords[k][1],
                                                    mol.cart_coords[k][2]))

        gaus_lines.append(
            '{} {}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}'.format('\n', '\n',
                                                                   '\n',
                                                                   '--Link1',
                                                                   '\n',
                                                                   '%mem=256MW',
                                                                   '\n',
                                                                   '%NProcShared=4',
                                                                   '\n',
                                                                   '%LindaWorker=localhost',
                                                                   '\n',
                                                                   '%chk=',
                                                                   mol_name,
                                                                   '.chk',
                                                                   '\n',
                                                                   '# b3lyp/aug-cc-pvdz freq',
                                                                   '\n',
                                                                   '# geom=allcheck guess=read',
                                                                   '\n',
                                                                   '# SCF=tight nosymm test',
                                                                   '\n',
                                                                   '# pop=MK iop(6/33=2,6/41=10,6/42=10,7/33=1)',
                                                                   '\n'))

        with open('mol.gau', 'w') as f:
            f.writelines(gaus_lines)


mol = thf

gaus = GaussianInput()
gaus.gaus_input(mol, '1', '-2')
