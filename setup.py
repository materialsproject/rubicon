#!/usr/bin/env python

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Apr 10, 2013"

import os
import sys
import glob

from setuptools import find_packages

try:
    from numpy.distutils.core import setup, Extension
except ImportError:
    print "numpy.distutils.core cannot be imported. Install numpy"
    sys.exit(-1)

from rubicon import __version__

module_dir = os.path.dirname(os.path.abspath(__file__))
f90_sources = glob.glob(os.path.join(module_dir, "rubicon", "analysis",
                                     "lammps", "*.f90"))
md_analyzer = Extension(name = 'rubicon.analysis.lammps._md_analyzer',
                           sources = f90_sources)

if __name__ == "__main__":
    setup(name='rubicon',
          version=__version__,
          description='codes for JCESR / Electrolyte Genome',
          long_description=open(os.path.join(module_dir, 'README.rst')).read(),
          url='https://github.com/materialsproject/rubicon',
          author='Anubhav Jain',
          author_email='anubhavster@gmail.com',
          license='modified BSD',
          packages=find_packages(),
          zip_safe=False,
          install_requires=['pymatgen>=3.0', 'fireworks>=0.9',
                            'custodian>=0.7', 'monty>=0.7.0'],
          extras_require={"plotting": ["matplotlib>=1.1"]},
          classifiers=["Programming Language :: Python :: 2.7",
                       "Development Status :: 2 - Pre-Alpha",
                       "Intended Audience :: Science/Research",
                       "Intended Audience :: System Administrators",
                       "Intended Audience :: Information Technology",
                       "Operating System :: OS Independent",
                       "Topic :: Other/Nonlisted Topic",
                       "Topic :: Scientific/Engineering"],
          test_suite='nose.collector',
          tests_require=['nose'],
          ext_modules=[md_analyzer],
          scripts=[os.path.join(os.path.join(module_dir, "scripts", f))
                   for f in os.listdir(os.path.join(module_dir, "scripts"))])
