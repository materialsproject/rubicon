# coding=utf-8

"""
This module implements input and output processing from MOPAC
"""
import copy
from textwrap import TextWrapper
from monty.io import zopen
from pymatgen import SymmOp
from pymatgen.core.structure import Molecule
from pymatgen.serializers.json_coders import MSONable
import re
import numpy as np
from pymatgen.util.coord_utils import get_angle

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
        title (str): comments for the job. Limited to two lines. The program
            will wrap the title into two lines automatically.
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
    zmat_patt = re.compile("^(\w+)*([\s,]+(\w+)[\s,]+(\w+))*[\-\.\s,\w]*$")
    xyz_patt = re.compile("^(\w+)[\s,]+([\d\.eE\-]+)[\s,]+([\d\.eE\-]+)[\s,]+"
                          "([\d\.eE\-]+)[\-\.\s,\w.]*$")

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

        title_wrapper = TextWrapper(width=80, expand_tabs=True, replace_whitespace=True, break_long_words=True,
                                    drop_whitespace=True)
        raw_title_lines = title_wrapper.wrap(title)[:2]
        title_lines = raw_title_lines + [""] * max([0, 2 - len(raw_title_lines)])
        self.title = title_lines
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

    def use_bfgs(self):
        self.keywords.pop("EF", None)
        self.keywords["BFGS"] = None

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
        lines = [' '.join(tokens)]
        return lines

    def _format_title(self):
        return self.title

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
        title = ' '.join(d["title"]) if len(d["title"][1]) > 0 else d["title"][0]
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
            return MopTask.from_string(f.read())

    @classmethod
    def from_string(cls, contents):
        """
        Creates a MopTask object from a string

        Args:
            contents: String representing a MOPAC input file

        Returns:
            MopTask object
        """
        lines = contents.split('\n')
        keywords = cls._parse_keywords(lines[0:1])
        title = lines[1: 3]
        mol = cls._parse_molecule(lines[3:])
        d = {"keywords": keywords, "title": title, "molecule": mol.to_dict,
             "@module": cls.__module__, "@class": cls.__name__}
        return MopTask.from_dict(d)

    @classmethod
    def _parse_keywords(cls, contents):
        line = contents[0]
        d = dict()
        int_pattern = re.compile('^[-+]?\d+$')
        float_pattern = re.compile('^[-+]?\d+\.\d+([eE][-+]?\d+)?$')
        for t in line.split():
            if "=" not in t:
                d[t.upper()] = None
            else:
                k, v = t.split("=")
                if int_pattern.match(v):
                    d[k.upper()] = int(v)
                elif float_pattern.match(v):
                    d[k.upper()] = float(v)
                else:
                    d[k.k.upper()] = v.upper()
        return d

    @classmethod
    def _parse_molecule(cls, contents):
        """
        Helper method to parse coordinates of Molecule. Copied from GaussianInput class.
        """
        paras = {}
        var_pattern = re.compile("^([A-Za-z]+\S*)[\s=,]+([\d\-\.]+)$")
        for l in contents:
            m = var_pattern.match(l.strip())
            if m:
                paras[m.group(1)] = float(m.group(2))

        species = []
        coords = []
        # Stores whether a Zmatrix format is detected. Once a zmatrix format
        # is detected, it is assumed for the remaining of the parsing.
        zmode = False
        for l in contents:
            l = l.strip()
            if not l:
                break
            if (not zmode) and cls.xyz_patt.match(l):
                m = cls.xyz_patt.match(l)
                species.append(m.group(1))
                toks = re.split("[,\s]+", l.strip())
                if len(toks) > 4:
                    coords.append(map(float, toks[2:5]))
                else:
                    coords.append(map(float, toks[1:4]))
            elif cls.zmat_patt.match(l):
                zmode = True
                toks = re.split("[,\s]+", l.strip())
                species.append(toks[0])
                toks.pop(0)
                if len(toks) == 0:
                    coords.append(np.array([0.0, 0.0, 0.0]))
                else:
                    nn = []
                    parameters = []
                    while len(toks) > 1:
                        ind = toks.pop(0)
                        data = toks.pop(0)
                        try:
                            nn.append(int(ind))
                        except ValueError:
                            nn.append(species.index(ind) + 1)
                        try:
                            val = float(data)
                            parameters.append(val)
                        except ValueError:
                            if data.startswith("-"):
                                parameters.append(-paras[data[1:]])
                            else:
                                parameters.append(paras[data])
                    if len(nn) == 1:
                        coords.append(np.array(
                            [0.0, 0.0, float(parameters[0])]))
                    elif len(nn) == 2:
                        coords1 = coords[nn[0] - 1]
                        coords2 = coords[nn[1] - 1]
                        bl = parameters[0]
                        angle = parameters[1]
                        axis = [0, 1, 0]
                        op = SymmOp.from_origin_axis_angle(coords1, axis,
                                                           angle, False)
                        coord = op.operate(coords2)
                        vec = coord - coords1
                        coord = vec * bl / np.linalg.norm(vec) + coords1
                        coords.append(coord)
                    elif len(nn) == 3:
                        coords1 = coords[nn[0] - 1]
                        coords2 = coords[nn[1] - 1]
                        coords3 = coords[nn[2] - 1]
                        bl = parameters[0]
                        angle = parameters[1]
                        dih = parameters[2]
                        v1 = coords3 - coords2
                        v2 = coords1 - coords2
                        axis = np.cross(v1, v2)
                        op = SymmOp.from_origin_axis_angle(
                            coords1, axis, angle, False)
                        coord = op.operate(coords2)
                        v1 = coord - coords1
                        v2 = coords1 - coords2
                        v3 = np.cross(v1, v2)
                        adj = get_angle(v3, axis)
                        axis = coords1 - coords2
                        op = SymmOp.from_origin_axis_angle(
                            coords1, axis, dih - adj, False)
                        coord = op.operate(coord)
                        vec = coord - coords1
                        coord = vec * bl / np.linalg.norm(vec) + coords1
                        coords.append(coord)

        def parse_species(sp_str):
            """
            The species specification can take many forms. E.g.,
            simple integers representing atomic numbers ("8"),
            actual species string ("C") or a labelled species ("C1").
            Sometimes, the species string is also not properly capitalized,
            e.g, ("c1"). This method should take care of these known formats.
            """
            try:
                return int(sp_str)
            except ValueError:
                sp = re.sub("\d", "", sp_str)
                return sp.capitalize()

        species = map(parse_species, species)

        return Molecule(species, coords)

    def use_precise(self, use=True):
        """
        Tighten criteria of the all the calculations.
        """
        if use:
            self.keywords["PRECISE"] = None
        else:
            self.keywords.pop("PRECISE", None)


class MopOutput(object):

    kcal_per_mol_2_eV = 4.3363E-2

    def __init__(self, filename):
        self.filename = filename
        with zopen(filename) as f:
            chunk = f.read()
        self.data = self._parse_job(chunk)

    @classmethod
    def _expected_successful_pattern(cls, input_keywords):
        text = ["SCF FIELD WAS ACHIEVED"]
        all_keys = input_keywords.keys()
        if "EF" in all_keys or "BFGS" in all_keys:
            text.append("GEOMETRY OPTIMISED USING .*\(.*\)\.")
        return text

    @classmethod
    def _parse_job(cls, output):
        heat_pattern = re.compile("FINAL HEAT OF FORMATION =\s+(?P<energy>-?\d+\.\d+)\s+KCAL/MOL")
        total_energy_pattern = re.compile("TOTAL ENERGY\s+=\s+(?P<energy>-\d+\.\d+)\s+EV")
        coord_pattern = re.compile("\s*\d+\s+(?P<element>[A-Z][a-z]*)\s+"
                                   "(?P<x>\-?\d+\.\d+)\s+"
                                   "(?P<y>\-?\d+\.\d+)\s+"
                                   "(?P<z>\-?\d+\.\d+)")
        energies = []
        parse_keywords = None
        result_section = False
        star_line_count = 0
        parse_coords = False
        input_keywords = None
        jobtype = None
        gracefully_terminated = False
        errors = []
        coords = []
        species = []
        molecules = []
        for line in output.split('\n'):
            if parse_keywords and input_keywords is None:
                input_keywords = MopTask._parse_keywords([line])
                jobtext = (set(input_keywords.keys()) & MopTask.available_sqm_tasktext).pop()
                jobtype = MopTask.jobtext2type[jobtext]
                parse_keywords = False
            if result_section and "*" * 50 in line:
                star_line_count += 1
                if star_line_count == 2:
                    parse_keywords = True
            if "PM7 CALCULATION RESULTS" in line:
                result_section = True
            if parse_coords:
                if "ATOM" in line:
                    continue
                if len(line.strip()) == 0:
                    if len(coords) == 0:
                        continue
                    else:
                        parse_coords = False
                        molecules.append(Molecule(species, coords))
                        species = None
                        coords = None
                        continue
                m = coord_pattern.match(line)
                coords.append([float(m.group("x")), float(m.group("y")),
                               float(m.group("z"))])
                species.append(m.group("element"))
            if "CARTESIAN COORDINATES" in line:
                parse_coords = True
                coords = []
                species = []
            m = heat_pattern.search(line)
            if m:
                heat_of_formation = float(m.group("energy")) * cls.kcal_per_mol_2_eV
                energies.append(tuple(["Heat of Formation", heat_of_formation]))
            m = total_energy_pattern.search(line)
            if m:
                total_energy = float(m.group("energy"))
                energies.append(tuple(["Total Energy", total_energy]))
            if "== MOPAC DONE ==" in line:
                gracefully_terminated = True

        if len(errors) == 0:
            for text in cls._expected_successful_pattern(input_keywords):
                sucess_pattern = re.compile(text)
                if not sucess_pattern.search(output):
                    errors.append("Can't find text to indicate success")

        data = {
            "jobtype": jobtype,
            "energies": energies,
            "molecules": molecules,
            "errors": errors,
            "has_error": len(errors) > 0,
            "gracefully_terminated": gracefully_terminated
        }
        return data

