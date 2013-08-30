__author__ = 'xiaohuiqu'

import json
import os

def singleton(class_):
    instances = {}

    def getinstance(*args, **kwargs):
        if class_ not in instances:
            instances[class_] = class_(*args, **kwargs)
        return instances[class_]

    return getinstance


@singleton
class NWChemBasisSetFinder():

    def __init__(self):
        data_file = os.path.join(os.path.dirname(__file__), 'data', 'available_nwchem_basis_set.json')
        with open(data_file) as f:
            self.basis_set = json.load(f)

    def basis_set_exist(self, element, b):
        if b.lower() not in self.basis_set.keys():
            return False
        return element in self.basis_set[b.lower()]

    def get_first_available_basis_set(self, element, preferred_basis_list):
        for b in preferred_basis_list:
            if self.basis_set_exist(element, b):
                return b
        else:
            return None

    def _write_nwchem_basis_set_json(self, library_path=None):
        import os
        import glob
        known_elements = set(['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                          'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
                          'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                          'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
                          'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
                          'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
                          'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
                          'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
                          'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
                          'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
                          'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
                          'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo', 'Uuu',
                          'Uub', 'Uuq', 'Uup', 'Uuh'])
        if library_path is None:
            nwchem_config_file = os.path.join(os.path.expanduser('~'), '.nwchemrc')
            if os.path.isfile(nwchem_config_file):
                with open(nwchem_config_file) as f:
                    path_dict = {line.split()[0]: line.strip().split()[1] for line in f.readlines()}
                    if 'nwchem_basis_library' in path_dict:
                        library_path =path_dict['nwchem_basis_library']
        basis_dict = {}
        for filename in glob.glob(os.path.join(library_path, '*')):
            with open(filename) as f:
                text = f.readlines()
            first_token = set([line.strip()[:3].strip() for line in text])
            basis_set = None
            for line in text:
                if line.startswith('basis'):
                    t1 = line[line.find('"')+1:line.rfind('"')].lower()
                    basis_set = t1[t1.find('_')+1:]
                    break
            if basis_set is None:
                continue
            elements = first_token.intersection(known_elements)
            basis_dict.update({basis_set: list(elements)})
        data_file = os.path.join(os.path.dirname(__file__), 'data', 'available_nwchem_basis_set.json')
        with open(data_file, 'w') as f:
            json.dump(basis_dict, f, indent=4)


if __name__ == '__main__':
    finder = NWChemBasisSetFinder()
    finder._write_nwchem_basis_set_json()
