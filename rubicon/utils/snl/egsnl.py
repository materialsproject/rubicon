from pymatgen import Composition, PMGJSONDecoder, Structure, Molecule
from pymatgen.matproj.snl import StructureNL

__author__ = 'xiaohuiqu'


def get_meta_from_structure(structure):
    comp = structure.composition
    elsyms = sorted(set([e.symbol for e in comp.elements]))
    meta = {'natoms': len(structure),
            'elements': elsyms,
            'nelements': len(elsyms),
            'formula': comp.formula,
            'reduced_cell_formula': comp.reduced_formula,
            'reduced_cell_formula_abc': Composition(comp.reduced_formula)
            .alphabetical_formula,
            'anonymized_formula': comp.anonymized_formula,
            'chemsystem': '-'.join(elsyms),
            'is_ordered': structure.is_ordered,
            'is_valid': structure.is_valid()}
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
        return self.snl_autometa['reduced_cell_formula_abc'] + "--" + str(self.sg_num)

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