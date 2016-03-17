#!/usr/bin/env python

import re

import numpy as np
# from json_coders import MSONable
import scipy.integrate
from multiprocessing import Pool
import time


def _list2float(seq):
    for x in seq:
        try:
            yield float(x)
        except ValueError:
            yield x


def autocorrelate(a):
    b = np.concatenate((a, np.zeros(len(a))), axis=1)
    c = np.fft.ifft(np.fft.fft(b) * np.conjugate(np.fft.fft(b))).real
    d = c[:len(c) / 2]
    d = d / (np.array(range(len(a))) + 1)[::-1]
    return d


class LammpsLog():
    """
    Parser for LAMMPS log file (parse function).
    Saves the output properties (log file) in the form of a dictionary (LOG) with the key being
    the LAMMPS output property (see 'thermo_style custom' command in the LAMMPS documentation).
    For example, LOG['temp'] will return the temperature data array in the log file.
    """

    def __init__(self, llog, avgs=None):
        """
        Args:
            llog:
                Dictionary of lamps log
            avgs:
                Dictionary of averages, will be generated automatically if unspecified
        """

        self.llog = llog  # Dictionary LOG has all the output property data as numpy 1D arrays with the property name as the key

        if avgs:
            self.avgs = avgs  # Dictionary of averages for storage / query
        else:
            self.avgs = {}
            # calculate the average
            for key in self.llog.keys():
                self.avgs[str(key)] = np.mean(self.llog[key])

    @classmethod
    def from_file(cls, filename):
        """
        Parses the log file. 
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
            print int(md_step / log_save_freq)

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

    def list_properties(self):
        """
        print the list of properties
        """
        print log.llog.keys()

    # viscosity
    def viscosity(self, cutoff):

        """
            cutoff: initial lines ignored during the calculation
            output: a file named viscosity_parallel.txt, 
            which saves the correlation function and its integration which is teh viscosity in cP
        """

        NCORES = 4
        p = Pool(NCORES)

        a1 = self.llog['pxy'][cutoff:]
        a2 = self.llog['pxz'][cutoff:]
        a3 = self.llog['pyz'][cutoff:]
        a4 = self.llog['pxx'][cutoff:] - self.llog['pyy'][cutoff:]
        a5 = self.llog['pyy'][cutoff:] - self.llog['pzz'][cutoff:]
        a6 = self.llog['pxx'][cutoff:] - self.llog['pzz'][cutoff:]
        array_array = [a1, a2, a3, a4, a5, a6]
        pv = p.map(autocorrelate, array_array)
        pcorr = (pv[0] + pv[1] + pv[2]) / 6 + (pv[3] + pv[4] + pv[5]) / 24

        temp = np.mean(self.llog['temp'][cutoff:])

        visco = (scipy.integrate.cumtrapz(pcorr,
                                          self.llog['step'][:len(pcorr)])) * \
                self.llog['timestep'] * 10 ** -15 * 1000 * 101325. ** 2 * \
                self.llog['vol'][-1] * 10 ** -30 / (1.38 * 10 ** -23 * temp)
        output = open('viscosity_parallel.txt', 'w')
        output.write(
            '#Time (fs), Average Pressure Correlation (atm^2), Viscosity (cp)\n')
        for line in zip(np.array(self.llog['step'][:len(pcorr) - 1]) *
                                self.llog['timestep'], pcorr, visco):
            output.write(' '.join(str(x) for x in line) + '\n')
        output.close()
        p.close()

    @property
    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "llog": self.llog,
                "avgs": self.avgs}

    @classmethod
    def from_dict(cls, d):
        return LammpsLog(d['llog'], d['avgs'])


if __name__ == '__main__':
    filename = 'visc.log'
    log = LammpsLog.from_file(filename)
    # print log.LOG.keys()
    # print log.llog
    # log.list_properties()
    log.viscosity(100001)
    # print np.mean(log.llog['density'])
    # print log.ave['step']
