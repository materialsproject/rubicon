from pymatgen.apps.battery.insertion_battery import InsertionElectrode
from pymatgen import Structure, MPRester
from pymatgen.entries.computed_entries import ComputedStructureEntry
import pymongo
from pymatpro.safety_paper.views_plot import HeatMapPlot

API_KEY = "x"


working_ion_pool = ["Al","Y","Mg","Ca","Zn"]
redox_ion_pool = ["Ti","V","Cr","Mn","Fe","Co","Ni"]

connection= pymongo.Connection('mongodb://USERNAME:PASSWORD@mongodb03.nersc.gov:27017/DBNAME')
db = connection['vasp_ml_prod']
materials = db.materials


results = [
[['Ti', 'Zn'], [u'mp-25262'], [u'sandbox-mvc-1514']],
[['V', 'Zn'], [u'mp-634430'], [u'mp-19032']],
[['Cr', 'Zn'], [u'sandbox-mvc-832'], [u'mp-647558']],
[['Mn', 'Zn'], [u'mp-25545'], [u'mp-649980']],
[['Fe', 'Zn'], [u'sandbox-mvc-819'], [u'mp-606972']],
[['Co', 'Zn'], [u'sandbox-mvc-1596', u'mp-632560'], []],
[['Ni', 'Zn'], [u'sandbox-mvc-1576', u'mp-25595'], []],

[['Ti', 'Mg'], [u'mp-25262'], [u'mp-504721']],
[['V', 'Mg'], [u'mp-634430'], [u'mp-18900']],
[['Cr', 'Mg'], [u'sandbox-mvc-832'], [u'mp-636411']],
[['Mn', 'Mg'], [u'mp-25545'], [u'mp-619628']],
[['Fe', 'Mg'], [u'sandbox-mvc-819'], [u'sandbox-mvc-1573']],
[['Co', 'Mg'], [u'sandbox-mvc-1596'], [u'sandbox-mvc-1774']],
[['Ni', 'Mg'], [u'mp-25595'], [u'sandbox-mvc-1597']],

[['Ti', 'Ca'], [u'mp-25262'], [u'sandbox-mvc-1437']],
[['V', 'Ca'], [u'mp-634430'], [u'sandbox-mvc-1570']],
[['Cr', 'Ca'], [u'sandbox-mvc-832'], [u'sandbox-mvc-1577']],
[['Mn', 'Ca'], [u'mp-25545'], [u'sandbox-mvc-1788']],
[['Fe', 'Ca'], [u'sandbox-mvc-819'], [u'sandbox-mvc-1608']],
[['Co', 'Ca'], [u'sandbox-mvc-1596'], [u'sandbox-mvc-1786']],
[['Ni', 'Ca'], [u'mp-25595'], [u'sandbox-mvc-1606']],

[['Ti', 'Al'], [u'mp-25262'], [u'sandbox-mvc-1531']],
[['V', 'Al'], [u'mp-634430'], [u'mp-565360']],
[['Cr', 'Al'], [u'sandbox-mvc-832'], [u'sandbox-mvc-1792']],
[['Mn', 'Al'], [u'mp-25545'], [u'sandbox-mvc-1800']],
[['Fe', 'Al'], [u'sandbox-mvc-819'], [u'sandbox-mvc-1803']],
[['Co', 'Al'], [u'sandbox-mvc-1596'], [u'sandbox-mvc-1768']],
[['Ni', 'Al'], [u'mp-25595'], [u'sandbox-mvc-1795']],

[['Ti', 'Y'], [u'mp-25262'], [u'sandbox-mvc-1565']],
[['V', 'Y'], [u'mp-634430'], []],
[['Cr', 'Y'], [u'sandbox-mvc-832'], []],
[['Mn', 'Y'], [u'mp-25545'], [u'sandbox-mvc-1781']],
[['Fe', 'Y'], [u'sandbox-mvc-819'], [u'sandbox-mvc-1830']],
[['Co', 'Y'], [u'sandbox-mvc-1596'], []],
[['Ni', 'Y'], [u'mp-25595'], [u'sandbox-mvc-1778']]
]


batt_list = []

for result in results:
    if [] not in result:
        id_charge = result[1][0]
        id_discharge = result[2][0]

        dict_charge = materials.find({"task_id":id_charge})[0]
        dict_discharge = materials.find({"task_id":id_discharge})[0]

        charge_struc = Structure.from_dict(dict_charge["structure"])
        charge_energy = dict_charge["final_energy"]
        entry_charge = ComputedStructureEntry(charge_struc,charge_energy)
        charge_ehull = dict_charge["e_above_hull"]

        discharge_struc = Structure.from_dict(dict_discharge["structure"])
        discharge_energy = dict_discharge["final_energy"]
        entry_discharge = ComputedStructureEntry(discharge_struc,discharge_energy)
        discharge_ehull = dict_discharge["e_above_hull"]


        mpr = MPRester(api_key=API_KEY, host="www.materialsproject.org")
        entries = mpr.get_entries(result[0][1],inc_structure="final")

        counter = 0
        energies = []
        for entry in entries:
            energies.append(entry.energy_per_atom)
        working_ion_entry = entries[(energies.index(min(energies)))]


        cathode = InsertionElectrode([entry_charge,entry_discharge],working_ion_entry)

        batt_list.append([result[0],result[1],result[2],charge_struc.formula,charge_ehull,discharge_struc.formula,discharge_ehull,cathode.max_voltage_step,cathode.get_average_voltage(),cathode.get_capacity_vol(),cathode.get_capacity_grav(),cathode.max_delta_volume])

voltage = []
capacity_grav = []
discharge_ehull = []
charge_ehull = []

for working_ion in working_ion_pool:
    ion_row = []
    v_row = []
    cap_g_row = []
    dc_ehull_row = []
    c_ehull_row = []

    for redox_ion in redox_ion_pool:
        ion_comb = [redox_ion,working_ion]
        battinfo = []
        for batt in batt_list:
            if ion_comb in batt:
                battinfo = batt
        if battinfo:
            ion_row.append(battinfo[0])
            v_row.append(battinfo[8])
            cap_g_row.append(battinfo[10])
            dc_ehull_row.append(battinfo[6])
            c_ehull_row.append(battinfo[4])

        else:
            v_row.append(0)
            cap_g_row.append(0)
            dc_ehull_row.append(0)
            c_ehull_row.append(0)
    voltage.append(v_row)
    capacity_grav.append(cap_g_row)
    discharge_ehull.append(dc_ehull_row)
    charge_ehull.append(c_ehull_row)

for vrow in voltage:
    print vrow

for crow in capacity_grav:
    print crow

for drow in discharge_ehull:
    print drow


plot_params = {}
hmp4 = HeatMapPlot([charge_ehull[0]], redox_ion_pool, [" "], m_props=plot_params)
hmp4.plot(export_filename="ehull_charge")

