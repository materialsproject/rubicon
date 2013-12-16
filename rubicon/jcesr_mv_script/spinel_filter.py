from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen import MPRester, Element, Structure

__author__ = 'Miao Liu'

import os
import numpy as np
from openpyxl import Workbook
import pymongo
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

__author__ = 'miaoliu'


MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
API_KEY = "Neo5xcHNOTE41tZi"

connection= pymongo.Connection('mongodb://admin_ml_prod:ydxr!oX^@mongodb03.nersc.gov:27017/vasp_ml_prod')

db = connection['vasp_ml_prod']

materials = db.materials

Working_ion_pool = ["Y"]

Redox_pool = ["Ti","V","Cr","Mn","Fe","Co","Ni"]

#Redox_pool = ["Ti","Ni"]

mpr = MPRester(api_key=API_KEY, host="www.materialsproject.org")


sm = StructureMatcher(ltol=0.2, stol=0.8, angle_tol=10, primitive_cell=True,
                 scale=True, attempt_supercell=True, allow_subset=True,
                 comparator=ElementComparator(), supercell_size='num_sites')

# mp-25545 MnO2, mp-619628 MgMn2O4
charged_prototype = "mp-25545"
discharged_prototype = "mp-619628"

formula_list = []

spinel_list = []

for working_ion in Working_ion_pool:
    for redox_ion in Redox_pool:
        formula_list.append([redox_ion,working_ion])


for formula in formula_list:
    print formula
    charged_entry = mpr.get_entries(charged_prototype,inc_structure="final")[0]
    charged_structure = charged_entry.structure
    discharged_entry = mpr.get_entries(discharged_prototype,inc_structure="final")[0]
    discharged_structure = discharged_entry.structure


    print charged_structure.formula,":",discharged_structure.formula
    charged_structure.replace_species({Element("Mn"):Element(formula[0])})
    discharged_structure.replace_species({Element("Mn"):Element(formula[0]),Element("Mg"):Element(formula[1])})
    print charged_structure.formula,":",discharged_structure.formula


    cursor = materials.find({"pretty_formula":formula[0]+"O2"},{"_id":0})
    counter = 1
    cl = []
    for c in cursor:
        #print counter,c["pretty_formula"],c["spacegroup"]["number"],c["task_id"]
        c_structure = Structure.from_dict(c["structure"])
        if sm.fit(charged_structure,c_structure):
            print counter,c["pretty_formula"],c["spacegroup"]["number"],c["task_id"],c["final_energy_per_atom"]
            cl.append(c["task_id"])
        counter += 1


    cursor = materials.find({"pretty_formula":formula[0]+"2"+formula[1]+"O4"},{"_id":0})
    if not cursor.count():
        cursor = materials.find({"pretty_formula":formula[1]+formula[0]+"2O4"},{"_id":0})
        if not cursor.count():
            cursor = materials.find({"pretty_formula":formula[1]+"("+formula[0]+"O2)2"},{"_id":0})
    counter = 1
    dcl = []
    for c in cursor:
        #print counter,c["pretty_formula"],c["spacegroup"]["number"],c["task_id"]
        c_structure = Structure.from_dict(c["structure"])
        print counter,c["pretty_formula"],c["spacegroup"]["number"],c["task_id"],c["final_energy_per_atom"]
        if sm.fit(discharged_structure,c_structure):
            print counter,c["pretty_formula"],c["spacegroup"]["number"],c["task_id"],c["final_energy_per_atom"]
            dcl.append(c["task_id"])
        counter += 1

    spinel_list.append([formula,cl,dcl])

for spinel in spinel_list:
    print spinel