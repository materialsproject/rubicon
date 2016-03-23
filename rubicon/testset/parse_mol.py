# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import glob
import re
import sys
import traceback

import requests
from six.moves import range

try:
    import pybel as pb
except ImportError:
    print("WARNING: Error importing pybel, setting pb to None!")
    pb = None

from pymatgen import Element, Molecule
from pymatgen.io.gaussian import GaussianInput
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.io.xyz import XYZ


def parse_file(filename):
    with open(filename) as f:
        txt = f.read()
    toks = txt.split("--link1--")
    for t in toks:
        try:
            lines = t.strip().split("\n")
            lines = [l.strip() for l in lines]
            gau = GaussianInput.from_string("\n".join(lines))
            yield (gau.molecule, gau.charge, gau.spin_multiplicity)
        except:
            print("error in {}".format(t))


session = requests.Session()


def get_nih_names(smiles):
    # noinspection PyBroadException
    try:
        smi = re.sub("#", "%23", smiles)
        response = session.get(
            "http://cactus.nci.nih.gov/chemical/structure/{}/names".format(
                smi))
        if response.status_code == 200:
            names = response.text.split("\n")
            return [n.strip() for n in names if n.strip() != ""]
        else:
            print("{} not found.\n".format(smiles))
            return []
    except:
        return []


def insert_g3testset(coll):
    for f in glob.glob("g*.txt"):
        print("Parsing " + f)
        for (m, charge, spin) in parse_file(f):
            try:
                clean_sites = []
                for site in m:
                    if Element.is_valid_symbol(site.specie.symbol):
                        clean_sites.append(site)
                clean_mol = Molecule.from_sites(clean_sites, charge=charge,
                                                spin_multiplicity=spin)
                xyz = XYZ(clean_mol)
                bb = BabelMolAdaptor.from_string(str(xyz), "xyz")
                pbmol = pb.Molecule(bb.openbabel_mol)
                smiles = pbmol.write("smi").split()[0]
                can = pbmol.write("can").split()[0]
                inchi = pbmol.write("inchi")
                svg = pbmol.write("svg")
                d = {"molecule": clean_mol.as_dict()}
                comp = clean_mol.composition
                d["pretty_formula"] = comp.reduced_formula
                d["formula"] = comp.formula
                d["composition"] = comp.as_dict()
                d["elements"] = list(comp.as_dict().keys())
                d["nelements"] = len(comp)
                d["charge"] = charge
                d["spin_multiplicity"] = spin
                d["smiles"] = smiles
                d["can"] = can
                d["inchi"] = inchi
                # d["names"] = get_nih_names(smiles)
                d["svg"] = svg
                d["xyz"] = str(xyz)
                d["tags"] = ["G305 test set"]
                coll.update({"inchi": inchi, "charge": charge,
                             "spin_multiplicity": spin}, {"$set": d},
                            upsert=True)
            except Exception as ex:
                print("Error in {}".format(f))
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_traceback,
                                          limit=2, file=sys.stdout)
        print("{} parsed!".format(f))


def insert_solvents(coll):
    names = """THF
    monoglyme
    Dimethylacetamide
    propylene carbonate
    ethylene carbonate
    dimethylcarbonate
    DMSO
    ACN
    Diethylcarbonate
    propyl glyme
    ethyl glyme
    ethyl diglyme
    diglyme
    Butyldiglyme
    tetraglyme
    Pentaethylene glycol diethyl ether
    Tetraethylene glycol diethyl ether
    Tetraethylene glycol dibutyl ether
    Butyldiglyme"""

    for n in names.split("\n"):
        response = requests.get(
            "http://cactus.nci.nih.gov/chemical/structure/{}/file?format=xyz".format(
                n))
        if response.status_code == 200:
            xyz = XYZ.from_string(response.text)
            clean_mol = xyz.molecule
            bb = BabelMolAdaptor(clean_mol)
            pbmol = pb.Molecule(bb.openbabel_mol)
            smiles = pbmol.write("smi").split()[0]
            can = pbmol.write("can").split()[0]
            inchi = pbmol.write("inchi")
            svg = pbmol.write("svg")

            d = {"molecule": clean_mol.as_dict()}
            comp = clean_mol.composition
            d["pretty_formula"] = comp.reduced_formula
            d["formula"] = comp.formula
            d["composition"] = comp.as_dict()
            d["elements"] = list(comp.as_dict().keys())
            d["nelements"] = len(comp)
            d["charge"] = clean_mol.charge
            d["spin_multiplicity"] = clean_mol.spin_multiplicity
            d["smiles"] = smiles
            d["can"] = can
            d["inchi"] = inchi
            # d["names"] = get_nih_names(smiles)
            d["svg"] = svg
            d["xyz"] = str(xyz)
            d["tags"] = ["Solvents"]
            coll.update({"inchi": inchi, "charge": clean_mol.charge,
                         "spin_multiplicity": clean_mol.spin_multiplicity},
                        {"$set": d}, upsert=True)
        else:
            print("{} not found.\n".format(n))


def insert_elements(coll):
    print("adding missing elements.")
    for z in range(1, 19):
        el = Element.from_Z(z)
        r = coll.find(filter={"formula": "{}1".format(el.symbol)})
        if r.count() == 0:
            try:
                clean_mol = Molecule([el], [[0, 0, 0]])
                xyz = XYZ(clean_mol)
                bb = BabelMolAdaptor.from_string(str(xyz), "xyz")
                pbmol = pb.Molecule(bb.openbabel_mol)
                smiles = pbmol.write("smi").split()[0]
                can = pbmol.write("can").split()[0]
                inchi = pbmol.write("inchi")
                svg = pbmol.write("svg")
                d = {"molecule": clean_mol.as_dict()}
                comp = clean_mol.composition
                d["pretty_formula"] = comp.reduced_formula
                d["formula"] = comp.formula
                d["composition"] = comp.as_dict()
                d["elements"] = list(comp.as_dict().keys())
                d["nelements"] = len(comp)
                d["charge"] = 0
                d["spin_multiplicity"] = clean_mol.spin_multiplicity
                d["smiles"] = smiles
                d["can"] = can
                d["inchi"] = inchi
                # d["names"] = get_nih_names(smiles)
                d["svg"] = svg
                d["xyz"] = str(xyz)
                d["tags"] = ["G305 test set"]
                coll.insert(d)
            except Exception as ex:
                print("Error in {}".format(el))
        elif r.count() > 1:
            print("More than 1 {} found. Removing...".format(el))
            results = list(r)
            for r in results[1:]:
                print(r["_id"])
                coll.remove({"_id": r["_id"]})


if __name__ == "__main__":
    try:
        from pymatpro.db.mongo.query_engine_mongo import MongoQueryEngine
    except ImportError:
        print("Install pymatpro")
    qe = MongoQueryEngine()
    db = qe.db
    coll = db["molecules"]
    # coll.remove({})
    insert_g3testset(coll)
    # insert_solvents(coll)
    # insert_elements(coll)
