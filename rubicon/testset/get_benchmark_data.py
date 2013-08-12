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
        print "Warning! webpage", url, "changed"
    for row in soup.table.findAll('tr')[1:]:
        tds = row.findAll('td')
        formula = str(tds[0].p.text).strip()
        d = result[formula]
        d[keyname] = dict(G3=float(tds[1].p.string), Expt=float(tds[2].p.string))


get_table(ip_url, "G3 Ionization Potentials", "IP", bench_dict)
get_table(ea_url, "G3 Electron Affinities", "EA", bench_dict)

with open("G3_Bench.json", 'w') as f:
    json.dump(bench_dict, f, indent=4, sort_keys=True)