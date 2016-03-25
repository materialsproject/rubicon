# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from rubicon.testset.get_g3_benchmark import get_g3_bench_collection

__author__ = 'Xiaohui Qu'

import urllib2
from BeautifulSoup import BeautifulSoup
from collections import defaultdict
import json

ip_url = "http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g3energies/G3IP.htm"
ea_url = "http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g3energies/G3ea.htm"

bench_dict = defaultdict(dict)

bench_dict['unit'] = 'kcal/mol'


def get_table(url, title, keyname, result):
    global soup, row, tds
    soup = BeautifulSoup(urllib2.urlopen(url).read())
    if soup.title.string != title:
        print("Warning! webpage", url, "changed")
    for row in soup.table.findAll('tr')[1:]:
        tds = row.findAll('td')
        formula = str(tds[0].p.text).strip()
        d = result[formula]
        d[keyname] = dict(G3=float(tds[1].p.string),
                          Expt=float(tds[2].p.string))


def refname2inchi(mission_tag):
    collection = get_g3_bench_collection()
    result_cursor = collection.find(filter={"user_tags.mission": mission_tag},
                                    projection=['inchi', 'user_tags.fw_name'])
    calc_result = list(result_cursor)
    fw_name_2_inchi = {m['user_tags']['fw_name'].strip(): m['inchi'] for m in
                       calc_result}

    with open("gauname2refname.json") as f:
        gau2ref = json.load(f)
    gau2ref.pop('NA')
    ref2gau = {v: k for k, v in gau2ref.items()}

    ref2inchi = {ref_name: fw_name_2_inchi[fw_name] for ref_name, fw_name in
                 ref2gau.items()
                 if fw_name in fw_name_2_inchi}

    return ref2inchi


if __name__ == '__main__':
    get_table(ip_url, "G3 Ionization Potentials", "IP", bench_dict)
    get_table(ea_url, "G3 Electron Affinities", "EA", bench_dict)

    bench_dict["CH3O CS (2A')"]["EA"] = bench_dict["CH3O"]["EA"]
    bench_dict.pop('CH3O')

    bench_dict["HS"]['IP'] = bench_dict['SH']['IP']
    bench_dict.pop('SH')

    r2i = refname2inchi("G2-97 Test Set Benchmark (Larry Scheme)")
    r2i.update(refname2inchi("G2-97 Test Set Benchmark (Shyue Scheme)"))

    for mol in bench_dict.keys():
        if mol in r2i:
            bench_dict[mol]["inchi"] = r2i[mol]

    with open("G3_ref_with_inchi.json", 'w') as f:
        json.dump(bench_dict, f, indent=4, sort_keys=True)
