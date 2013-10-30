import datetime
from pymatgen import Composition, PMGJSONDecoder, Structure, Molecule
from pymatgen.analysis.molecule_matcher import MoleculeMatcher, InchiMolAtomMapper
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.matproj.snl import StructureNL


__author__ = 'xiaohuiqu'


def get_meta_from_structure(mol):
    '''
    set basis information for the molecule

    Args:
        mol: pymatgen Molecule object

    Returns:
        the meta information dict
    '''
    comp = mol.composition
    elsyms = sorted(set([e.symbol for e in comp.elements]))
    bb = BabelMolAdaptor(mol)
    pbmol = bb.pybel_mol
    inchi = pbmol.write("inchi").strip()
    meta = {'nsites': len(mol),
            'elements': elsyms,
            'nelements': len(elsyms),
            'formula': comp.formula,
            'reduced_cell_formula': comp.reduced_formula,
            'reduced_cell_formula_abc': Composition(comp.reduced_formula)
            .alphabetical_formula,
            'anonymized_formula': comp.anonymized_formula,
            'chemsystem': '-'.join(elsyms),
            'is_valid': mol.is_valid(),
            'inchi': inchi}
    return meta


class EGStructureNL(StructureNL):
    # adds snl_id, pointgroup, and autometa properties to StructureNL.

    def __init__(self, *args, **kwargs):
        super(EGStructureNL, self).__init__(*args, **kwargs)
        if not self.pointgroup:
            raise ValueError('An EGStructureNL must have a pointgroup '
                             'assigned!')
        self.snl_autometa = get_meta_from_structure(self.structure)

    @property
    def snl_id(self):
        return self.data['_electrolytegenome']['snl_id']

    @property
    def pointgroup(self):
        return self.data['_electrolytegenome']['pointgroup']

    @property
    def snlgroup_key(self):
        return self.snl_autometa['inchi']

    @property
    def to_dict(self):
        m_dict = super(EGStructureNL, self).to_dict
        m_dict.update(self.snl_autometa)
        m_dict['snl_id'] = self.snl_id
        m_dict['snlgroup_key'] = self.snlgroup_key
        return m_dict

    @staticmethod
    def from_dict(d):
        a = d["about"]
        dec = PMGJSONDecoder()

        created_at = dec.process_decoded(a["created_at"]) if "created_at" in a \
            else None
        data = {k: v for k, v in d["about"].items()
                if k.startswith("_")}
        data = dec.process_decoded(data)

        structure = Structure.from_dict(d) if "lattice" in d \
            else Molecule.from_dict(d)
        return EGStructureNL(structure, a["authors"],
                           projects=a.get("projects", None),
                           references=a.get("references", ""),
                           remarks=a.get("remarks", None), data=data,
                           history=a.get("history", None),
                           created_at=created_at)

    @staticmethod
    def from_snl(snl, snl_id, pointgroup):
        # make a copy of SNL
        snl2 = StructureNL.from_dict(snl.to_dict)
        if '_electrolytegenome' not in snl2.data:
            snl2.data['_electrolytegenome'] = {}

        snl2.data['_electrolytegenome']['snl_id'] = snl_id
        snl2.data['_electrolytegenome']['pointgroup'] = pointgroup

        return EGStructureNL.from_dict(snl2.to_dict)


class SNLGroup():
    def __init__(self, snlgroup_id, canonical_snl, all_snl_ids=None):
        # Auto fields
        self.created_at = datetime.datetime.utcnow()
        self.updated_at = datetime.datetime.utcnow()

        # User fields
        self.snlgroup_id = snlgroup_id
        self.canonical_snl = canonical_snl

        self.all_snl_ids = all_snl_ids if all_snl_ids else []
        if self.canonical_snl.snl_id not in self.all_snl_ids:
            self.all_snl_ids.append(canonical_snl.snl_id)

        # Convenience fields
        self.canonical_structure = canonical_snl.structure
        self.snl_autometa = get_meta_from_structure(self.canonical_structure)

    @property
    def to_dict(self):
        d = self.snl_autometa
        d['created_at'] = self.created_at
        d['updated_at'] = self.updated_at
        d['snlgroup_id'] = self.snlgroup_id
        d['canonical_snl'] = self.canonical_snl.to_dict
        d['all_snl_ids'] = self.all_snl_ids
        d['num_snl'] = len(self.all_snl_ids)

        d['snlgroup_key'] = self.canonical_snl.snlgroup_key
        return d

    @staticmethod
    def from_dict(d):
        return SNLGroup(d['snlgroup_id'],
                        EGStructureNL.from_dict(d['canonical_snl']),
                        d['all_snl_ids'])

    def add_if_belongs(self, cand_snl, exact_match=True):

        # no need to compare if structue is different
        if cand_snl.snlgroup_key != self.canonical_snl.snlgroup_key:
            return False

        # make sure the structure is not already in all_structures
        if cand_snl.snl_id in self.all_snl_ids:
            print 'WARNING: add_if_belongs() has detected that you are trying to add the same SNL id twice!'
            return False

        if exact_match:
            mm = MoleculeMatcher(tolerance=0.2,
                                 mapper=InchiMolAtomMapper(angle_tolerance=5.0))

            if not mm.fit(cand_snl.structure, self.canonical_structure):
                return False

        # everything checks out, add to the group
        self.all_snl_ids.append(cand_snl.snl_id)
        self.updated_at = datetime.datetime.utcnow()

        return True