# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import re

import numpy as np

from six.moves import range

from monty.json import MSONable


def _list2float(seq):
    for x in seq:
        try:
            yield float(x)
        except ValueError:
            yield x


class LammpsLog(MSONable):
    """
    Parser for LAMMPS log file (parse function).
    Saves the output properties (log file) in the form of a dictionary (LOG)
    with the key being the LAMMPS output property (see 'thermo_style custom'
    command in the LAMMPS documentation).
    For example, LOG['temp'] will return the temperature data array in the log file.
    """

    def __init__(self, llog, avgs=None):
        """
        Args:
            llog (dict):
                Dictionary LOG has all the output property data as numpy 1D arrays with the property name as the key
            avgs:
                Dictionary of averages, will be generated automatically if unspecified
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
            print(int(md_step / log_save_freq))

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
    #log.viscosity(100001)

