__author__ = 'Xiaohui Qu'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Xiaohui Qu'
__email__ = 'xqu@lbl.gov'
__date__ = 'Aug 08, 2013'

import glob
import json
import os
import re

import openbabel as ob

from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.matproj.snl import StructureNL


def read_mol_text(text):
    obconv = ob.OBConversion()
    obconv.SetInFormat('sdf')
    obmol = ob.OBMol()
    obconv.ReadString(obmol, text)
    return obmol


def build3d(obmol):
    gen3d = ob.OBOp.FindType("Gen3D")
    gen3d.Do(obmol)


def get_cas_number(text):
    m = re.search(r"> <cas.rn> \((.*)\)", text)
    if m is None:
        return "CAS Number Unknown"
    return m.group(1)


def str_to_obmols(text):
    mol_texts = [s + '$$$$\n' for s in text.split('$$$$\n')]
    return [(read_mol_text(s), get_cas_number(s)) for s in mol_texts if
            len(s.split('\n')) > 5]


if __name__ == '__main__':
    filenames = sorted(glob.glob("cas_files/*.sdf"), reverse=True)
    for filename in filenames:
        dirname = filename[:-4]
        if os.path.exists(dirname):
            print "directory " + dirname + " already exists"
            print "please delete it before use this script"
            exit(0)
    for filename in filenames:
        dirname = filename[:-4]
        os.mkdir(dirname)
        print "reading", filename
        text = None
        with open(filename) as f:
            text = f.read()
        mol_tokens = str_to_obmols(text)
        for (i, (mol, cas)) in enumerate(mol_tokens):
            print "processing molecule", i + 1, cas, "of", len(
                mol_tokens), "molecules"
            try:
                build3d(mol)
            except:
                os.system("echo " + cas + " >> failed_mols.txt")
        pmg_mols = [(BabelMolAdaptor(obmol).pymatgen_mol, cas) for (obmol, cas)
                    in mol_tokens]
        snl_texts = [StructureNL(mol, "Xiaohui Qu <xqu@lbl.gov>", remarks=cas)
                     for (mol, cas) in pmg_mols]
        for snl in snl_texts:
            with open(dirname + "/" + snl.remarks[0] + ".snl", 'w') as f:
                json.dump(snl.as_dict(), f, indent=4)
    print "Done"
