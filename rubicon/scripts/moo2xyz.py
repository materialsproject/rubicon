# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import argparse

from rubicon.io.mopac.mopacio import MopOutput

__author__ = 'xiaohuiqu'


def main():
    parser = argparse.ArgumentParser(
        description="Convert a MOPAC output file to an XYZ file")
    parser.add_argument("-i", "--moo", dest="moo_file", type=str,
                        required=True,
                        help="the MOPAC output file")
    parser.add_argument("-o", "--xyz", dest="xyz_file", type=str,
                        required=True,
                        help="The XYZ file to write to")
    options = parser.parse_args()
    moo = MopOutput(options.moo_file)
    moo.data["molecules"][-1].to("xyz", options.xyz_file)
