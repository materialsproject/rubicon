# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
    A wrapper for AntechamberRunner which generates force field files
    for a specified molecule using gaussian output file as input
"""

import shlex
import subprocess
import tempfile
from collections import namedtuple

from monty.dev import requires
from monty.os.path import which
from monty.tempfile import ScratchDir

from rubicon.io.amber.topology import Topology, TopCorruptionException
from rubicon.io.amber.topology import correct_corrupted_top_files
from rubicon.io.amber.generalized_force_field import GeneralizedForceField, FFCorruptionException
from rubicon.io.amber.generalized_force_field import correct_corrupted_frcmod_files


__author__ = 'Navnidhi Rajput, Kiran Mathew'


class AntechamberRunner(object):
    """
    A wrapper for AntechamberRunner software

    """

    @requires(which('parmchk'), "Requires the binary parmchk."
                                "Install AmberTools from http://ambermd.org/#AmberTools")
    @requires(which('antechamber'), "Requires the binary antechamber."
                                    "Install AmberTools from http://ambermd.org/#AmberTools")
    def __init__(self, mols):
        """
        Args:
            mols: List of molecules
        """

        self.mols = mols

    def _run_parmchk(self, filename='ANTECHAMBER_AC.AC'):
        """
        run parmchk using ANTECHAMBER_AC.AC file

        Args:
            filename: pdb file of the molecule
        """
        command_parmchk = (
            'parmchk -i ' + filename + ' -f ac -o mol.frcmod -a Y')
        exit_code = subprocess.call(shlex.split(command_parmchk))
        return exit_code

    def _run_antechamber(self, filename):
        """
        run antechmaber using the provided gaussian output file
        """
        command = 'antechamber -i ' + filename + " -fi gout -o mol -fo charmm -c resp -s 2 runinfo"
        exit_code = subprocess.call(shlex.split(command))
        return exit_code

    def _get_gaussian_ff_top_single(self, filename=None):
        """
        run antechamber using gaussian output file, then run parmchk
        to generate missing force field parameters. Store and return
        the force field and topology information in ff_mol.

        Args:
            filename: gaussian output file of the molecule

        Returns:
            Amberff namedtuple object that contains information on force field and
            topology
        """
        scratch = tempfile.gettempdir()
        Amberff = namedtuple("Amberff", ["force_field", "topology"])
        with ScratchDir(scratch, copy_from_current_on_enter=True,
                        copy_to_current_on_exit=True) as d:
            # self._convert_to_pdb(mol, 'mol.pdb')
            # self.molname = filename.split('.')[0]
            self._run_antechamber(filename)
            self._run_parmchk()
            # if antechamber can't find parameters go to gaff_example.dat
            try:
                top = Topology.from_file('mol.rtf')
            except TopCorruptionException:
                correct_corrupted_top_files('mol.rtf', 'gaff_example.txt')
                top = Topology.from_file('mol.rtf')
            try:
                gff = GeneralizedForceField.from_file('mol.frcmod')
            except FFCorruptionException:
                correct_corrupted_frcmod_files('ANTECHAMBER.FRCMOD', 'gaff_example.txt')
                gff = GeneralizedForceField.from_file('ANTECHAMBER.FRCMOD')
            # gff.set_atom_mappings('ANTECHAMBER_AC.AC')
            # gff.read_charges()
            # decorate the molecule with the sire property "atomname"
            #mol.add_site_property("atomname", (list(gff.atom_index.values())))
        return Amberff(gff, top)

    def get_gaussian_ff_top(self, filenames):
        """
        return a list of amber force field and topology for the list of
        gaussian output filenames corresponding to each molecule in mols list.

        Args:
            filenames (list): list of gaussian output files for each type of molecule

        Returns:
            list of Amberff namedtuples
        """
        amber_ffs = []
        for fname in filenames:
            amber_ffs.append(self._get_gaussian_ff_top_single(filename=fname))
        return amber_ffs
