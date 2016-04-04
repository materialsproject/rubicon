# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module implements classes for processing Lammps output files.
"""

import re
from io import open
from collections import defaultdict

import numpy as np

from monty.json import MSONable


__author__ = "Kiran Mathew, Navnidhi Rajput, Michael Humbert"
__email__ = "kmathew@lbl.gov"


class LammpsRun(object):
    """
    Parse the lammps data file, trajectory file and the log file to extract
    useful info about the system.

    Args:
        data_file (str): path to the data file
        trajectory_file (str): path to the trajectory file
        log_file (str): path to the log file
    """
    def __init__(self, data_file, trajectory_file, log_file):
        self.data_file = data_file
        self.trajectory_file = trajectory_file
        self.log_file = log_file
        self.lammpslog = LammpsLog.from_file(log_file)
        self._parse_data()
        self._set_mol_masses_and_charges()
        self._parse_trajectory()

    def _parse_data(self):
        """
        parse lammps data file
        """
        self.atomic_masses = {} # atom_type: mass
        self.molecules = defaultdict(list) # mol_num:[[atom,atom_type]]
        self.atom_to_mol = {} # atom_num: mol_num
        self.mol_types = {} # 'H2O': 20
        mols_pattern = re.compile("^#\s+(\d+)\s+([0-9a-zA-Z]+)\s+molecules$")
        atoms_pattern = re.compile("^(\d+)\s+(\d+)\s+(\d+)\s+([0-9eE\.+-]+)\s+("
                                 "[0-9eE\.+-]+)\s+([0-9eE\.+-]+)\w*")
        masses_pattern = re.compile("^(\d+)\s+([0-9\.]+)$")
        with open(self.data_file) as df:
            for line in df:
                if mols_pattern.search(line):
                    m = mols_pattern.search(line)
                    self.mol_types[m.group(2)] = int(m.group(1))
                if masses_pattern.search(line):
                    m = masses_pattern.search(line)
                    self.atomic_masses[int(m.group(1))] = float(m.group(
                        2))
                if "atoms" in line:
                    self.natoms = int(line.split()[0])
                    self.atomic_charges = np.zeros(self.natoms)
                m = atoms_pattern.search(line)
                if m:
                    self.atomic_charges[int(m.group(1))-1] = float(m.group(4))
                    self.atom_to_mol[int(m.group(1))] = int(m.group(2))
                    self.molecules[int(m.group(2))-1].append([
                        int(m.group(1))-1,
                        int(m.group(3))])
        self.nmols = len(self.molecules.keys())
        for k, v in self.molecules.items():
            self.molecules[k] = np.array(v)

    def _parse_trajectory(self):
        """
        parse the trajectory file
        """
        self.box_size = []
        traj_timesteps = []
        trajectory = []
        timestep_label = "ITEM: TIMESTEP"
        # box_label "ITEM: BOX BOUNDS"
        # traj_label, "ITEM: ATOMS id type x y z vx vy vz mol"
        box_pattern = re.compile("^([0-9\.-]+)\s+([0-9\.-]+)$")
        traj_pattern = re.compile("\s*(\d+)\s+(\d+)\s+([0-9eE\.+-]+)\s+([0-9eE\.+-]+)\s+"
                                  "([0-9eE\.+-]+)\s+"
                                  "([0-9eE\.+-]+)\s+"
                                "([0-9eE\.+-]+)\s+([0-9eE\.+-]+)\s+(\d+)\s*")
        parse_timestep = False
        with open(self.trajectory_file) as tf:
            for line in tf:
                if timestep_label in line:
                    parse_timestep = True
                    continue
                if parse_timestep:
                    traj_timesteps.append(float(line))
                    parse_timestep = False
                if box_pattern.search(line) and len(traj_timesteps) == 1:
                        m = box_pattern.search(line)
                        #self.box_size.append(m.group(1))
                        self.box_size.append(float(m.group(2)))
                if traj_pattern.search(line):
                    m = traj_pattern.search(line)
                    line_data = []
                    line_data.append(int(m.group(1))-1)
                    line_data.append(self.atomic_masses[int(m.group(2))])
                    line_data.extend([float(x) for i, x in enumerate(
                        m.groups()) if i+1 > 2 and i+1 < 9])
                    line_data.append(int(m.group(9)))
                    trajectory.append(tuple(line_data))
        traj_dtype = np.dtype([('Atoms_id'.encode("ascii"), np.int64),
                               ('mass'.encode("ascii"), np.float64),
                               ('x'.encode("ascii"), np.float64),
                               ('y'.encode("ascii"), np.float64),
                               ('z'.encode("ascii"), np.float64),
                               ('vx'.encode("ascii"), np.float64),
                               ('vy'.encode("ascii"), np.float64),
                               ('vz'.encode("ascii"), np.float64),
                               ('mol'.encode("ascii"), np.int64)])
        self.trajectory = np.array(trajectory, dtype=traj_dtype)
        self.traj_timesteps = np.array(traj_timesteps, dtype=np.float64)
        for step in range(self.traj_timesteps.size):
            begin = step * self.natoms
            end = (step + 1) * self.natoms
            self.trajectory[begin:end] = np.sort(self.trajectory[begin:end],
                                                 order="Atoms_id".encode("ascii"))

    def _set_mol_masses_and_charges(self):
        """
        set the charge and mass for each molecule
        """
        self.mol_charges = np.zeros(self.nmols)
        self.mol_masses = np.zeros(self.nmols)
        for mol_id, val in self.molecules.items():
            self.mol_masses[mol_id-1] = sum([self.atomic_masses[atom[1]] for
                                        atom in val])
            self.mol_charges[mol_id-1] = sum([self.atomic_charges[atom[0]-1]
                                              for atom in val])

    @property
    def timestep(self):
        return self.lammpslog.timestep

    def _weighted_average(self, step, mol_id, mol_vector):
        """
        Calculate the weighted average of the array comprising of
        atomic vectors corresponding to the molecule with id mol_id for
        each time step.

        Args:
            step (int): time step
            mol_id (int): molecule id
            mol_vector (numpy array): array of shape,
                natoms_in_molecule with id mol_id x 3

        Returns:
            1D numpy array(3 x 1) of weighted averages in x, y, z directions
        """
        mol_masses = self._get_mol_masses(step, mol_id)
        return np.array([np.dot(mol_vector[:, dim],
                                mol_masses) / np.sum(mol_masses) for dim in range(3)])

    def _get_mol_vector(self, step, mol_id, param=["x", "y", "z"]):
        """
        Returns numpy array corresponding to atomic vectors of paramter
        "param" for the given time step and molecule id

        Args:
            step (int): time step
            mol_id (int): molecule id
            param (list): the atomic parameter for which the weighted
                average is to be computed

        Returns:
            2D numpy array(natoms_in_molecule x 3) of atomic vectors
        """
        begin = step * self.natoms
        end = (step + 1) * self.natoms
        mol_vector_structured = \
            self.trajectory[begin:end][self.molecules[mol_id][:,0]][param]
        new_shape = mol_vector_structured.shape + (-1,)
        mol_vector = mol_vector_structured.view(np.float64).reshape(new_shape)
        return mol_vector

    def _get_mol_masses(self, step, mol_id):
        """
        Return the array comprising of atomic masses for the given
        molecule and time step.

        Args:
            step (int): time step
            mol_id (int): molecule id

        Returns:
            1D numpy array(natoms_in_molecule x 1).
        """
        begin = step * self.natoms
        end = (step + 1) * self.natoms
        mol_masses_structured = \
        self.trajectory[begin:end][self.molecules[mol_id][:, 0]]["mass"]
        new_shape = mol_masses_structured.shape
        mol_masses = mol_masses_structured.view(np.float64).reshape(new_shape)
        return mol_masses

    @property
    def mol_trajectory(self):
        """
        Compute the weighted average trajectory of each molecule at each
        timestep

        Returns:
            2D numpy array ((n_timesteps*mols_number) x 3)
        """
        traj = []
        for step in range(self.traj_timesteps.size):
            tmp_mol = []
            for mol_id in range(self.nmols):
                mol_coords = self._get_mol_vector(step, mol_id,
                                                  param=["x", "y", "z"])
                # take care of periodic boundary conditions
                pbc_wrap(mol_coords, self.box_size)
                tmp_mol.append(self._weighted_average(step, mol_id, mol_coords))
            traj.append(tmp_mol)
        return np.array(traj)

    @property
    def mol_velocity(self):
        """
         Compute the weighted average velcoity of each molecule at each
         timestep.

         Returns:
            2D numpy array ((n_timesteps*mols_number) x 3)
        """
        velocity = []
        for step in range(self.traj_timesteps.size):
            tmp_mol = []
            for mol_id in range(self.nmols):
                mol_velocities = self._get_mol_vector(step, mol_id,
                                                      param=["vx", "vy", "vz"])
                tmp_mol.append(self._weighted_average(step, mol_id, mol_velocities))
            velocity.append(tmp_mol)
        return np.array(velocity)


class LammpsLog(MSONable):
    """
    Parser for LAMMPS log file (parse function).
    Saves the output properties (log file) in the form of a dictionary (LOG)
    with the key being the LAMMPS output property (see 'thermo_style custom'
    command in the LAMMPS documentation).

    For example, LOG['temp'] will return the temperature data array in the log
    file.
    """

    def __init__(self, llog, avgs=None):
        """
        Args:
            llog (dict):
                Dictionary LOG has all the output property data as numpy 1D a
                rrays with the property name as the key
            avgs:
                Dictionary of averages, will be generated automatically if
                unspecified
        """
        self.llog = llog
        if avgs:
            self.avgs = avgs
        else:
            self.avgs = {}
            # calculate the average
            for key in self.llog.keys():
                self.avgs[str(key)] = np.mean(self.llog[key])

    @classmethod
    def from_file(cls, filename):
        """
        Parse the log file.
        """
        md = 0  # To avoid reading the minimization data steps
        header = 0
        footer_blank_line = 0
        llog = {}
        with open(filename, 'r') as logfile:
            total_lines = len(logfile.readlines())
            logfile.seek(0)
            for line in logfile:
                # timestep
                time = re.search('timestep\s+([0-9]+)', line)
                if time:
                    timestep = float(time.group(1))
                    llog['timestep'] = timestep
                # total steps of MD
                steps = re.search('run\s+([0-9]+)', line)
                if steps:
                    md_step = float(steps.group(1))
                    md = 1
                # save freq to log
                thermo = re.search('thermo\s+([0-9]+)', line)
                if thermo:
                    log_save_freq = float(thermo.group(1))
                # log format
                format = re.search('thermo_style.+', line)
                if format:
                    data_format = format.group().split()[2:]
                if all(isinstance(x, float) for x in
                       list(_list2float(line.split()))) and md == 1: break
                header += 1
            # note: we are starting from the "break" above
            for line in logfile:
                if line == '\n':
                    footer_blank_line += 1
            if total_lines >= header + md_step / log_save_freq:
                rawdata = np.genfromtxt(fname=filename, dtype=float,
                                        skip_header=header, skip_footer=int(
                        total_lines - header - md_step / log_save_freq - 1) - footer_blank_line)
            else:
                rawdata = np.genfromtxt(fname=filename, dtype=float,
                                        skip_header=header, skip_footer=1)
            for column, property in enumerate(data_format):
                llog[property] = rawdata[:, column]
            return LammpsLog(llog)

    @property
    def timestep(self):
        return self.llog['timestep']

    @property
    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "lammpslog": self.llog,
                "avgs": self.avgs}

    @classmethod
    def from_dict(cls, d):
        return LammpsLog(d['lammpslog'], d['avgs'])


def pbc_wrap(array, box_size):
    """
    wrap the array for molecule coordinates around the periodic boundary.

    Args:
        array (numpy.ndarray): molecule coordinates, [[x1,y1,z1],[x2,y2,z2],..]
        box_size (list): [x_length, y_length, z_length]
    """
    ref = array[0,0]
    for i in range(3):
        array[:, i] = np.where((array[:, i] - ref) >= box_size[i] / 2,
                               array[:, i] - box_size[i], array[:, i])
        array[:, i] = np.where((array[:, i] - ref) < -box_size[i] / 2,
                               array[:, i] + box_size[i], array[:, i])


def _list2float(seq):
    for x in seq:
        try:
            yield float(x)
        except ValueError:
            yield x


if __name__ == '__main__':
    filename = 'visc.log'
    log = LammpsLog.from_file(filename)
    # log.viscosity(100001)
