# coding=utf-8

"""
This module implements input and output processing from MOPAC
"""
import copy
from textwrap import TextWrapper
from monty.io import zopen
from pymatgen.core.structure import Molecule
from pymatgen.serializers.json_coders import MSONable

__author__ = "Xiaohui Qu"
__copyright__ = "Copyright 2014, The Electrolyte Genome Project"
__version__ = "0.1"
__maintainer__ = "Xiaohui Qu"
__email__ = "xhqu1981@gmail.com"
__date__ = "7/28/14"


class MopTask(MSONable):
    """
    An object representinbg a MOPAC input file.
    Currently, defaults to XYZ coordinates.

    Args:
        molecule (pymatgen Molecule object): The input molecule.
        charge (int): the charge of the molecule.
        jobtype (str): The type of MOPAC job. "SP" for single point energy,
            "opt" for geometry optimization, "freq" for vibration frequency.
        title (str): comments for the job. Limited to two lines.
        sqm_method (str): theoretical level. e. g. PM6, PM7, AM1
        optional_params (dict): addictional parameters. If the keyword isn't
            associated with any values, set the value to None.
            e. g. {"CYCLES": 500, "XYZ": None}
    """

    available_sqm_methods = {"PM7", "PM7-TS", "PM6", "PM3", "AM1", "PM6-D3",
                             "PM6-DH+", "PM6-DH2", "PM6-DH2X"}
    jobtype2text = {"SP": "1SCF", "OPT": "EF", "OPT_BFGS": "BFGS", "FREQ": "THERMO"}
    jobtext2type = {v: k for k, v in jobtype2text.iteritems()}
    available_sqm_tasktext = set(jobtype2text.values())

    def __init__(self, molecule, charge, jobtype, title, sqm_method="PM7",
                 optional_params=None):
        if not isinstance(molecule, Molecule):
            raise Exception("only pymatgen Molecule object is accepted")
        self.mol = copy.deepcopy(molecule)
        if not isinstance(charge, int):
            raise Exception("charge must an integer")
        self.keywords = dict()

        explicit_keywords = {"CHARGE"} | self.available_sqm_tasktext | self.available_sqm_methods
        optional_keywords = set([k.upper() for k in optional_params.keys()]) if optional_params else set()
        overlap_keywords = explicit_keywords & optional_keywords
        if len(overlap_keywords) > 0:
            raise Exception(" ".join(overlap_keywords) + " are duplicated"
                            "in optional_params")
        self.keywords["CHARGE"] = charge
        if jobtype.upper() not in ["SP", "OPT", "FREQ"]:
            raise Exception('Job type "{}" is not supported currently'.format(jobtype))
        self.keywords[self.jobtype2text[jobtype.upper()]] = None
        self.title = title
        if sqm_method.upper() not in self.available_sqm_methods:
            raise Exception('Semi-empirical methods "{}" is not supported currently'.format(sqm_method))
        self.keywords[sqm_method.upper()] = None
        if optional_params:
            for k, v in optional_params.iteritems():
                if isinstance(v, str) or isinstance(v, unicode):
                    self.keywords[k.upper()] = v.upper()
                else:
                    self.keywords[k.upper()] = v
        if "INT" in self.keywords.keys():
            raise Exception("Internal coordinates is not supported yet")
        self.keywords["XYZ"] = None

    @property
    def molecule(self):
        return self.mol

    def set_geom_max_iterations(self, iterations):
        """
        Set the max iterations of geometry optimization.

        Args:
            iterations: the maximum iterations of geometry optimization.
            (Integer)
        """
        if not isinstance(iterations, int):
            raise Exception("max iterations must be an integer")
        self.keywords["CYCLES"] = iterations

    def __str__(self):
        lines = []
        lines.extend(self._format_keywords())
        lines.extend(self._format_title())
        lines.extend(self._format_molecule())
        lines.append('\n')
        return '\n'.join(lines)

    def _format_keywords(self):
        all_keys = set(self.keywords.keys())
        sqm_method = (all_keys & self.available_sqm_methods).pop()
        task = (all_keys & self.available_sqm_tasktext).pop()
        priority_key_list = [sqm_method, task, "CHARGE"]
        full_key_list = priority_key_list + list(all_keys - set(priority_key_list))
        tokens = []
        for k in full_key_list:
            v = self.keywords[k]
            if v is None:
                tokens.append(k)
            elif isinstance(v, int):
                tokens.append("{k:s}={v:d}".format(k=k, v=v))
            elif isinstance(v, str) or isinstance(v, unicode):
                tokens.append("{k:s}={v:s}".format(k=k, v=v))
        lines = []
        lines.append(' '.join(tokens))
        return lines

    def _format_title(self):
        wrapper = TextWrapper(width=80, expand_tabs=True, replace_whitespace=True, break_long_words=True)
        raw_lines = wrapper.wrap(self.title)[:2]
        lines = raw_lines + [""] * max([0, 2 - len(raw_lines)])
        return lines

    def _format_molecule(self):
        lines = []
        for site in self.mol.sites:
            lines.append(" {element:<4} {x:>17.8f} {y:>17.8f} "
                         "{z:>17.8f}".format(element=site.species_string,
                                             x=site.x, y=site.y, z=site.z))
        return lines

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "molecule": self.mol.to_dict,
                "keywords": self.keywords,
                "title": self.title}

    @classmethod
    def from_dict(cls, d):
        mol = Molecule.from_dict(d["molecule"])
        charge = d["keywords"]["CHARGE"]
        all_keys = set(d["keywords"].keys())
        sqm_method = (all_keys & cls.available_sqm_methods).pop()
        jobtext = (all_keys & cls.available_sqm_tasktext).pop()
        jobtype = cls.jobtext2type[jobtext]
        title = d["title"]
        used_key = ["CHARGE", sqm_method, jobtext]
        optional_key = list(all_keys - set(used_key))
        optional_params = {k: d["keywords"][k] for k in optional_key}
        mop = MopTask(mol, charge, jobtype, title, sqm_method, optional_params)
        return mop

    def write_file(self, filename):
        with zopen(filename, "w") as f:
            f.write(self.__str__())

    @classmethod
    def from_file(cls, filename):
        with zopen(filename) as f:
            return cls.from_string(f.read())


    @classmethod
    def from_string(cls, contents):
        pass
