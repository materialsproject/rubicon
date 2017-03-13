# coding: utf-8

import glob
import os
import sys

from setuptools import find_packages

try:
    from numpy.distutils.core import setup, Extension
except ImportError:
    print("numpy.distutils.core cannot be imported. Install numpy")
    sys.exit(-9)

from rubicon import __version__

module_dir = os.path.dirname(os.path.abspath(__file__))
f90_sources = glob.glob(os.path.join(module_dir, "rubicon", "analysis",
                                     "lammps", "fortran_files", "*.f90"))
md_analyzer = Extension(name='rubicon.analysis.lammps._md_analyzer',
                        sources=f90_sources)

if __name__ == "__main__":
    setup(name='rubicon',
          version=__version__,
          description='Codes for JCESR / Electrolyte Genome',
          long_description=open(os.path.join(module_dir, 'README.rst')).read(),
          url='https://github.com/materialsproject/rubicon',
          author='Xiaohui Qu, Navnidhi Rajput, Kiran Mathew, Michael Humbert, '
                 'Zhang Yong, Shyue Ping Ong, Anubhav Jain',
          author_email='xqu@lbl.gov',
          license='modified BSD',
          packages=find_packages(),
          install_requires=['pymatgen>=3.6.0', 'fireworks>=1.3.0',
                            'custodian>=1.0', 'monty>=0.8.0'],
          extras_require={"plotting": ["matplotlib>=1.1"],
                          "molecules": ["openbabel"],
                          "molecular dynamics": ["Packmol", "AmberTools"]},
          keywords=["materials", "project", "electrolyte",
                    "molecular dynamics",
                    "lammps", "qchem", "analysis"],
          classifiers=["Programming Language :: Python :: 2.7",
                       "Development Status :: 2 - Pre-Alpha",
                       "Intended Audience :: Science/Research",
                       "Intended Audience :: System Administrators",
                       "Intended Audience :: Information Technology",
                       "Operating System :: OS Independent",
                       "Topic :: Other/Nonlisted Topic",
                       "Topic :: Scientific/Engineering"],
          ext_modules=[md_analyzer],
          scripts=[os.path.join(os.path.join(module_dir, "scripts", f))
                   for f in os.listdir(os.path.join(module_dir, "scripts"))])
