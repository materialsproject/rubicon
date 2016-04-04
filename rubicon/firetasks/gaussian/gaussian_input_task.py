# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import shlex
import subprocess

__author__ = 'navnidhirajput'

from fireworks import FireTaskBase, explicit_serialize, FWAction
from pymatgen.io import gaussian


@explicit_serialize
class WritegaussianGeoTask(FireTaskBase):
    """
    Writes Gaussian Input files for geometry optimization.

    """

    _fw_name = "Gaussian Geometry Writer"

    def run_task(self, fw_spec):
        mol = fw_spec["molecule"]
        mol_name = fw_spec["mol_name"]
        charge = fw_spec["charge"]
        spin_multiplicity = fw_spec["spin_multiplicity"]

        gaus_lines = gaussian.GaussianInput(mol, charge=charge,
                                            spin_multiplicity=spin_multiplicity,
                                            title='created by gaussian_geo_task from' + ' ' + mol_name,
                                            functional="b3lyp",
                                            basis_set="aug-cc-pvdz",
                                            route_parameters={
                                                'opt': "(calcfc,tight)",
                                                'int': "ultrafine",
                                                "\n# SCF": "tight"},
                                            input_parameters=None,
                                            link0_parameters={
                                                "%mem": "256MW",
                                                "%NProcShared": 4,
                                                "%LindaWorker": "localhost",
                                                "%chk": mol_name + ".chk"},
                                            dieze_tag="#",
                                            gen_basis=None)

        gaus_lines.write_file('mol_geo.gau', cart_coords=True)

        with open('mol_geo.gau') as f, open("mol_geo.out", 'w') as fo:
            subprocess.call(shlex.split("g09launch"), stdin=f, stdout=fo)

        prev_gaussian_geo = shlex.os.path.join(shlex.os.getcwd(),
                                               'mol_geo.out')
        update_spec = {'prev_gaussian_geo': prev_gaussian_geo}

        return FWAction(update_spec=update_spec)


@explicit_serialize
class WritegaussianFreqESPTask(FireTaskBase):
    """
    Writes Gaussian Input files for frequency and charge calculation.

    """

    _fw_name = "Gaussian frequency and charge Writer"

    def run_task(self, fw_spec):
        filename = fw_spec['prev_gaussian_geo']
        gaus_geo = gaussian.GaussianOutput(filename)
        mol_opt = gaus_geo.final_structure
        mol_name = fw_spec["mol_name"]
        charge = fw_spec["charge"]
        spin_multiplicity = fw_spec["spin_multiplicity"]

        gaus_freq_charge = gaussian.GaussianInput(mol_opt, charge=charge,
                                                  spin_multiplicity=spin_multiplicity,
                                                  title='created by gaussian_frq_task from' + ' ' + mol_name,
                                                  functional="b3lyp",
                                                  basis_set="aug-cc-pvdz  freq",
                                                  route_parameters={
                                                      "SCF": "tight",
                                                      "pop": "MK iop(6/33=2,6/41=10,6/42=10,7/33=1)"},
                                                  input_parameters=None,
                                                  link0_parameters={
                                                      "%mem": "30GB",
                                                      "%NProcShared": 4,
                                                      "%LindaWorker": "localhost",
                                                      "%chk": mol_name + ".chk"},
                                                  dieze_tag="#",
                                                  gen_basis=None)
        gaus_freq_charge.write_file('mol_freq.gau', cart_coords=True)

        with open('mol_freq.gau') as f, open("mol_freq.out", 'w') as fo:
            subprocess.call(shlex.split("g09launch"), stdin=f, stdout=fo)

        prev_gaussian_freq = shlex.os.path.join(shlex.os.getcwd(),
                                                'mol_freq.out')
        update_spec = {'prev_gaussian_freq': prev_gaussian_freq}

        return FWAction(update_spec=update_spec)
