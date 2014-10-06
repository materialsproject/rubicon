from fireworks.core.firework import Firework, Workflow
from fireworks.utilities.fw_utilities import get_slug
from pymatgen import Composition
from rubicon.dupefinders.dupefinder_eg import DupeFinderEG
from rubicon.firetasks.egsnl_tasks import AddEGSNLTask
from rubicon.utils.snl.egsnl import EGStructureNL, get_meta_from_structure
from rubicon.workflows.bsse_wf import counterpoise_correction_generation_fw
from rubicon.workflows.md_relax_workflow import md_relax_fws
from rubicon.workflows.multi_solvent_ipea_wf import multi_solvent_ipea_fws
from rubicon.workflows.multistep_ipea_wf import multistep_ipea_fws
from rubicon.workflows.single_point_energy_wf import single_point_energy_fws
from rubicon.workflows.solvation_energy_wf import solvation_energy_fws


def snl_to_eg_wf(snl, parameters=None):
    fws_all = []
    parameters = parameters if parameters else {}

    snl_priority = parameters.get('priority', 1)
    mission = parameters.get('mission', 'Electron Genome Production')
    priority = snl_priority * 2  # once we start a job, keep going!

    f = Composition(snl.structure.composition.reduced_formula).\
        alphabetical_formula
    molname = parameters.get("nick_name", f)

    # add the SNL to the SNL DB and figure out duplicate group
    tasks = [AddEGSNLTask()]
    spec = {'task_type': 'Add to SNL database',
            'snl': snl.as_dict(),
            '_priority': snl_priority}
    if 'snlgroup_id' in parameters and isinstance(snl, EGStructureNL):
        spec['force_egsnl'] = snl.to_dict
        spec['force_snlgroup_id'] = parameters['snlgroup_id']
        del spec['snl']
    fws_all.append(Firework(tasks, spec,
                   name=get_slug(molname + ' -- Add to SNL database'),
                   fw_id=1))

    default_solvents = ['diglym', 'acetonitrile', 'dmso', 'thf',
                        'dimethylamine', 'dimethoxyethane',
                        'dimethylaniline', 'tetraglyme']

    workflow_type = parameters.get('workflow', 'ipea')
    ref_charge = parameters.get('ref_charge', 0)
    spin_multiplicities = parameters.get('spin_multiplicities', (2, 1, 2))
    user_tags = {"initial_charge": ref_charge}
    solvent_method = parameters.get("solvent_method", "ief-pcm")
    use_vdW_surface = parameters.get("use_vdW_surface", False)
    qm_method = parameters.get("qm_method", None)
    population_method = parameters.get("population_method", None)
    check_large = parameters.get("check_large", True)
    if workflow_type == 'ipea':
        solvent = parameters.get('solvent', "water")
        fws_tasks, connections = multistep_ipea_fws(
            mol=snl.structure, name=molname, mission=mission, solvent=solvent, solvent_method=solvent_method,
            use_vdW_surface=use_vdW_surface,
            ref_charge=ref_charge, spin_multiplicities=spin_multiplicities, dupefinder=DupeFinderEG(),
            priority=priority, parent_fwid=1, additional_user_tags=user_tags, qm_method=qm_method,
            check_large=check_large)
    elif workflow_type == 'multiple solvent ipea':
        solvents = parameters.get('solvents', default_solvents)
        fws_tasks, connections = multi_solvent_ipea_fws(
            mol=snl.structure, name=molname, mission=mission, solvents=solvents, solvent_method=solvent_method,
            use_vdW_surface=use_vdW_surface,
            ref_charge=ref_charge, spin_multiplicities=spin_multiplicities, dupefinder=DupeFinderEG(),
            priority=priority, parent_fwid=1, additional_user_tags=user_tags, qm_method=qm_method)
    elif workflow_type == 'solvation energy':
        solvents = parameters.get('solvents', default_solvents)
        fws_tasks, connections = solvation_energy_fws(
            mol=snl.structure, name=molname, mission=mission, solvents=solvents, solvent_method=solvent_method,
            use_vdW_surface=use_vdW_surface,
            dupefinder=DupeFinderEG(), priority=priority, parent_fwid=1, additional_user_tags=user_tags,
            qm_method=qm_method)
    elif workflow_type == "single point energy":
        solvent = parameters.get('solvent', "water")
        fws_tasks, connections = single_point_energy_fws(
            mol=snl.structure, name=molname, mission=mission, solvent=solvent, solvent_method=solvent_method,
            use_vdW_surface=use_vdW_surface,
            qm_method=qm_method, pop_method=population_method, dupefinder=DupeFinderEG(),
            priority=priority, parent_fwid=1, additional_user_tags=user_tags)
    elif workflow_type == "md relax":
        high_temperature = parameters.get("high_temperature", 323.15)
        low_temperature = parameters.get("low_temperature", 273.15)
        md_steps = parameters.get("md_steps", 500)
        time_step = parameters.get("time_step", 1.0)
        md_runs = parameters.get("md_runs", 3)
        normal_basis = parameters.get("normal_basis", "6-31G*")
        diffuse_basis = parameters.get("diffuse_basis", "6-31+G*")
        charge_threshold = parameters.get("charge_threshold", -0.5)
        fws_tasks, connections = md_relax_fws(
            mol=snl.structure, name=molname, mission=mission, qm_method=qm_method,
            high_temperature=high_temperature, low_temperature=low_temperature, md_steps=md_steps,
            time_step=time_step, md_runs=md_runs, normal_basis=normal_basis, diffuse_basis=diffuse_basis,
            charge_threshold=charge_threshold, dupefinder=DupeFinderEG(), priority=priority,
            parent_fwid=1, additional_user_tags=user_tags)
    elif workflow_type == "bsse":
        charge = snl.structure.charge
        spin_multiplicity = snl.structure.spin_multiplicity
        fragments = parameters.get("fragments", None)
        fws_tasks, connections = counterpoise_correction_generation_fw(
            molname=molname, charge=charge, spin_multiplicity=spin_multiplicity, qm_method=qm_method,
            fragments=fragments, mission=mission, priority=priority, parent_fwid=1,
            additional_user_tags=user_tags)
    else:
        raise ValueError('Workflow "{}" is not supported yet'.
                         format(workflow_type))
    fws_all.extend(fws_tasks)

    wf_meta = get_meta_from_structure(snl.structure)
    wf_meta['run_version'] = 'Jan 27, 2014'

    if '_electrolytegenome' in snl.data and \
            'submission_id' in snl.data['_electrolytegenome']:
        wf_meta['submission_id'] = snl.data['_electrolytegenome'][
            'submission_id']
    return Workflow(fws_all, connections,
                    name=molname,
                    metadata=wf_meta)
