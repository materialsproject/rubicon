from __future__ import print_function, unicode_literals

import os
import warnings


MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


def _warning(message, category=UserWarning, filename='', lineno=-1):
    print("{} from {}: {}".format(category.__name__, MODULE_DIR, message))

warnings.showwarning = _warning


def lammps_warning():
    warnings.warn("io.lammps has been merged with pymatgen. "
                  "The module will be removed from rubicon after the next "
                  "release of pymatgen.", FutureWarning)
