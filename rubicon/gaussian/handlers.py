#!/usr/bin/env python

"""
TODO: Change the module doc.
"""

from __future__ import division

from monty.json import MSONable

__author__ = "Shyue Ping Ong"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Beta"
__date__ = "5/20/13"

from custodian.custodian import ErrorHandler


class GaussianErrorHandler(ErrorHandler, MSONable):
    """
    Error handler for Gaussian Jobs.
    """

    def __init__(self, output_filename="gau.out"):
        self.output_filename = output_filename

    def check(self):
        # TODO: Implement actual checks.
        return False

    def correct(self):
        # TODO: Implement corrections. Right now, this should never be reached
        # since check always return False.
        raise NotImplementedError("Not yet implemented corrections.")
        return {"errors": [], "actions": []}

    @property
    def is_monitor(self):
        return False

    def __str__(self):
        return "GaussianErrorHandler"

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "output_filename": self.output_filename}

    @classmethod
    def from_dict(cls, d):
        return cls(d["output_filename"])
