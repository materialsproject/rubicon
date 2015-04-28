from fireworks.core.firework import Firework, Workflow
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen.core.structure import Molecule
from rubicon.firetasks.multistep_qchem_task import BasisSetSuperpositionErrorCalculationTask, \
    CounterpoiseCorrectionGenerationTask
from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator
from rubicon.utils.snl.egsnl import EGStructureNL

__author__ = 'xiaohuiqu'


class BSSEFragment(FWSerializable):
    OVERLAPPED = "overlapped"
    ISOLATED = "isolated"

    def __init__(self, charge, spin_multiplicity, ghost_atoms):
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.ghost_atoms = sorted(set(ghost_atoms))
        if not isinstance(self.ghost_atoms, list):
            raise ValueError("ghost atoms must be list of integers")
        for atom in self.ghost_atoms:
            if not isinstance(atom, int):
                raise ValueError("Each element of ghost atoms must be an integer")

    @staticmethod
    def get_host_atoms(frag_atoms, mol):
        all_atoms = set(range(len(mol)))
        return sorted(all_atoms - set(frag_atoms))

    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "charge": self.charge,
                "spin_multiplicity": self.spin_multiplicity,
                "ghost_atoms": self.ghost_atoms}

    @classmethod
    def from_dict(cls, m_dict):
        return BSSEFragment(m_dict["charge"], m_dict["spin_multiplicity"], m_dict["ghost_atoms"])


def get_sub_mol(mol, frag):
    species = [site.species_string for i, site in enumerate(mol.sites) if i not in frag.ghost_atoms]
    coords = [site.coords for i, site in enumerate(mol.sites) if i not in frag.ghost_atoms]
    sub_mol = Molecule(species, coords, frag.charge, frag.spin_multiplicity)
    return sub_mol

def counterpoise_correction_generation_fw(molname, charge, spin_multiplicity, qm_method, fragments,
                       mission, priority=1, parent_fwid=None, additional_user_tags=None, large=False):
    fw_spec = dict()
    fw_spec["user_tags"] = dict()
    fw_spec["user_tags"]["molname"] = molname
    fw_spec["user_tags"]["mission"] = mission
    fw_spec["qm_method"] = qm_method
    fw_spec["fragments"] = fragments
    fw_spec["charge"] = charge
    fw_spec["spin_multiplicity"] = spin_multiplicity
    fw_spec["large"] = large
    if priority:
        fw_spec['_priority'] = priority
    fw_spec["user_tags"].update(additional_user_tags)
    fwid_base = 1
    if parent_fwid:
        if not (isinstance(parent_fwid, int) or isinstance(parent_fwid, list)):
            raise ValueError("Parent FireWork ID must be integer or list")
        parent_fwid = parent_fwid if isinstance(parent_fwid, list) \
            else [parent_fwid]
        fwid_base = max(parent_fwid) + 1
    current_fwid = fwid_base
    links_dict = dict()
    fw_cp = Firework([CounterpoiseCorrectionGenerationTask()],
                     spec=fw_spec, name=molname + " Counterpoise Correction Generation", fw_id=current_fwid)
    for p_fwid in parent_fwid:
        links_dict[p_fwid] = current_fwid
    return [fw_cp], links_dict


def bsse_fws(super_mol_egsnl, name, super_mol_snlgroup_id, super_mol_charge, super_mol_spin_multiplicity,
             super_mol_inchi_root, qm_method, fragments, mission, dupefinder=None, priority=1,
             parent_fwid=None, additional_user_tags=None, is_spawnned=False, large=False):
    super_mol = EGStructureNL.from_dict(super_mol_egsnl).structure
    fwid_base = 1
    if parent_fwid:
        if not (isinstance(parent_fwid, int) or isinstance(parent_fwid, list)):
            raise ValueError("Parent FireWork ID must be integer or list")
        parent_fwid = parent_fwid if isinstance(parent_fwid, list) \
            else [parent_fwid]
        fwid_base = max(parent_fwid) + 1
    if is_spawnned:
        fwid_base = -1
    fwid_incr_factor = 1 if not is_spawnned else -1
    fws = []
    db_fwids = []
    links_dict = dict()
    current_fwid = fwid_base
    for frag in fragments:
        if len(frag.ghost_atoms) == 0:
            continue
        frag_name = name + "_" + BasisSetSuperpositionErrorCalculationTask.get_fragment_name(frag.ghost_atoms)
        fw_ov_creator = QChemFireWorkCreator(mol=super_mol, molname=frag_name, mission=mission, dupefinder=dupefinder,
                                             priority=priority, additional_user_tags=additional_user_tags, large=large)
        fw_ov_cal_id = current_fwid
        current_fwid += fwid_incr_factor
        fw_ov_db_id = current_fwid
        current_fwid += fwid_incr_factor
        fws_ov = fw_ov_creator.vacuum_only_sp_fw(frag.charge, frag.spin_multiplicity, fw_ov_cal_id, fw_ov_db_id,
                                                 priority=priority, qm_method=qm_method,
                                                 super_mol_snlgroup_id=super_mol_snlgroup_id,
                                                 super_mol_egsnl=super_mol_egsnl,
                                                 super_mol_inchi_root=super_mol_inchi_root,
                                                 ghost_atoms=frag.ghost_atoms, bs_overlap=True)
        fws.extend(fws_ov)
        db_fwids.append(fw_ov_db_id)
        links_dict[fw_ov_cal_id] = fw_ov_db_id
        if parent_fwid:
            for p_fwid in parent_fwid:
                links_dict[p_fwid] = fw_ov_cal_id

        sub_mol = get_sub_mol(super_mol, frag)
        fw_iso_creator = QChemFireWorkCreator(mol=sub_mol, molname=frag_name, mission=mission, dupefinder=dupefinder,
                                              priority=priority, additional_user_tags=additional_user_tags, large=large)
        fw_iso_cal_id = current_fwid
        current_fwid += fwid_incr_factor
        fw_iso_db_id = current_fwid
        current_fwid += fwid_incr_factor
        fws_iso = fw_iso_creator.vacuum_only_sp_fw(frag.charge, frag.spin_multiplicity, fw_iso_cal_id, fw_iso_db_id,
                                                   priority=priority, qm_method=qm_method,
                                                   super_mol_snlgroup_id=super_mol_snlgroup_id,
                                                   super_mol_egsnl=super_mol_egsnl,
                                                   super_mol_inchi_root=super_mol_inchi_root,
                                                   ghost_atoms=None, bs_overlap=False)
        fws.extend(fws_iso)
        db_fwids.append(fw_iso_db_id)
        links_dict[fw_iso_cal_id] = fw_iso_db_id
        if parent_fwid:
            for p_fwid in parent_fwid:
                links_dict[p_fwid] = fw_ov_cal_id
            for p_fwid in parent_fwid:
                links_dict[p_fwid] = fw_iso_cal_id
    user_tags = {"molname": name, "mission": mission}
    if additional_user_tags:
        user_tags.update(additional_user_tags)
    bsse_spec = {"mol": super_mol,
                 "fragments": fragments,
                 "snlgroup_id": super_mol_snlgroup_id,
                 "egsnl": super_mol_egsnl,
                 "user_tags": user_tags,
                 "qm_method": qm_method,
                 "charge": super_mol_charge,
                 "spin_multiplicity": super_mol_spin_multiplicity,
                 "inchi_root": super_mol_inchi_root}
    fw_bsse = Firework([BasisSetSuperpositionErrorCalculationTask()],
                       spec=bsse_spec, name=name+" BSSE Calculation", fw_id=current_fwid)
    for i in db_fwids:
        links_dict[i] = current_fwid
    fws.append(fw_bsse)
    return fws, links_dict

def bsse_wf(super_mol, name, **kwargs):
    fws, links_dict = bsse_fws(super_mol, name, **kwargs)
    return Workflow(fws, links_dict, name)