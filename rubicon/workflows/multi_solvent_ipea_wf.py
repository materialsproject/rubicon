from rubicon.workflows.multistep_ipea_wf import QChemFireWorkCreator


def multi_solvent_ipea_fws(mol, name, mission, dupefinder=None, priority=1,
                       parent_fwid=None):
    large = False
    if len(mol) > 50:
        large = True
    fw_creator = QChemFireWorkCreator(mol=mol, molname=name, mission=mission,
                                      dupefinder=dupefinder,
                                      priority=priority, large=large)
    fwid_base = 1
    if parent_fwid:
        if not (isinstance(parent_fwid, int) or isinstance(parent_fwid, list)):
            raise ValueError("Parent FireWork ID must be integer or list")
        parent_fwid = parent_fwid if isinstance(parent_fwid, list) \
            else [parent_fwid]
        fwid_base = max(parent_fwid) + 1
    fireworks = []
    links_dict = dict()
    # the task in the order of anion, neutral, cation
    cgi_cal, ngi_cal, agi_cal = (None, None, None)
    cfi_db, nfi_db, afi_db = (None, None, None)
    cgi_db, ngi_db, agi_db = (None, None, None)
    charge = (-1, 0, 1)
    spin_multiplicity = (2, 1, 2)
    if len(mol) > 1:
        fw_ids = zip(* [iter(range(fwid_base + 0, fwid_base + 6))] * 2)
        fws = (fw_creator.geom_fw(ch, spin, fwid_cal, fwid_db)
               for ch, spin, (fwid_cal, fwid_db)
               in zip(charge, spin_multiplicity, fw_ids))
        (cgi_cal, cgi_db), (ngi_cal, ngi_db), (agi_cal, agi_db) = fw_ids
        fireworks.extend(itertools.chain.from_iterable(fws))
        links_dict.update(dict(fw_ids))

        if not large:
            fw_ids = zip(* [iter(range(fwid_base + 6, fwid_base + 6 + 6))] * 2)
            fws = (fw_creator.freq_fw(ch, spin, fwid_cal, fwid_db)
                   for ch, spin, (fwid_cal, fwid_db)
                   in zip(charge, spin_multiplicity, fw_ids))
            (cfi_cal, cfi_db), (nfi_cal, nfi_db), (afi_cal, afi_db) = fw_ids
            fireworks.extend(itertools.chain.from_iterable(fws))
            links_dict.update(dict(fw_ids))
            links_dict.update({cgi_db: cfi_cal,
                               ngi_db: nfi_cal,
                               agi_db: afi_cal})

    fw_ids = zip(* [iter(range(fwid_base + 12, fwid_base + 12 + 6))] * 2)
    fws = (fw_creator.sp_fw(ch, spin, fwid_cal, fwid_db)
           for ch, spin, (fwid_cal, fwid_db)
           in zip(charge, spin_multiplicity, fw_ids))
    (cspi_cal, cspi_db), (nspi_cal, nspi_db), (aspi_cal, aspi_db) = fw_ids
    links_dict.update((dict(fw_ids)))
    fireworks.extend(itertools.chain.from_iterable(fws))
    if len(mol) > 1:
        if large:
            links_dict.update({cgi_db: cspi_cal, ngi_db: nspi_cal,
                               agi_db: aspi_cal})
        else:
            links_dict.update({cfi_db: cspi_cal, nfi_db: nspi_cal,
                               afi_db: aspi_cal})
        links_dict.update({nspi_db: [cgi_cal, agi_cal]})
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = ngi_cal
    else:
        links_dict.update({nspi_db: [cspi_cal, aspi_cal]})
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = nspi_cal
    return fireworks, links_dict


def mol_to_solvent_ipea_wf(mol, name, mission, dupefinder=None, priority=1,
                   parent_fwid=None):
    fireworks, links_dict = multi_solvent_ipea_fws(mol, name, mission,
                                                 dupefinder,
                                               priority, parent_fwid)
