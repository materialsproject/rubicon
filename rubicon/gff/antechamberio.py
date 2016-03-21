""" A wrapper for AntechamberRunner which generates force field files
    for a specified molecule using gaussian output file as input
"""

__author__ = 'navnidhirajput'

import subprocess
import shlex
import tempfile

from monty.os.path import which
from monty.dev import requires
from monty.tempfile import ScratchDir

from rubicon.gff.ffmol import FFmol
from rubicon.gff import Gff, FFCorruptionException, \
    correct_corrupted_frcmod_files
from rubicon.gff.topology import TopMol, TopCorruptionException, \
    correct_corrupted_top_files


class AntechamberRunner():
    """
    A wrapper for AntechamberRunner software

    """
    @requires(which('parmchk'), "Requires the binary parmchk")
    @requires(which('antechamber'), "Requires the binary antechamber")
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
        return_cmd = subprocess.call(shlex.split(command_parmchk))
        return return_cmd

    def get_ff_top_mol(self, mol, filename=None):
        """
        run antechamber using gaussian output file then run parmchk
        to generate missing force field parameters. Store and return
        the force field and topology information in ff_mol.

        Args:
            filename: gaussian output file of the molecule
            mol: molecule
        Returns:
            ff_mol : Object of FFmol which contains information of
                     force field and topology
        """
        scratch = tempfile.gettempdir()

        with ScratchDir(scratch, copy_from_current_on_enter=True,
                        copy_to_current_on_exit=True) as d:

            # self._convert_to_pdb(mol, 'mol.pdb')
            command = 'antechamber -i ' + filename + " -fi gout -o mol -fo charmm -c resp -s 2 runinfo"
            return_cmd = subprocess.call(shlex.split(command))
            self.molname = filename.split('.')[0]
            self._run_parmchk()
            # if antechamber can't find parameters go to gaff_nidhi.dat
            try:
                top = TopMol.from_file('mol.rtf')


            except TopCorruptionException:
                correct_corrupted_top_files('mol.rtf', 'gaff_nidhi.txt')
                top = TopMol.from_file('mol.rtf')

            try:
                gff = Gff.from_forcefield_para('mol.frcmod')

            except FFCorruptionException:
                correct_corrupted_frcmod_files('ANTECHAMBER.FRCMOD',
                                               'gaff_nidhi.txt')
                gff = Gff.from_forcefield_para('ANTECHAMBER.FRCMOD')
            gff.read_atom_index(mol, 'ANTECHAMBER_AC.AC')
            # gff.read_charges()

            mol.add_site_property("atomname", (gff.atom_index.values()))
        ffmol = FFmol(gff, top)
        return ffmol
