from fireworks.core.firework import FireWork, Workflow
from fireworks.utilities.fw_utilities import get_slug
from pymatgen import Composition
from rubicon.dupefinders.dupefinder_eg import DupeFinderEG
from rubicon.firetasks.egsnl_tasks import AddEGSNLTask
from rubicon.utils.snl.egsnl import EGStructureNL, get_meta_from_structure
from rubicon.workflows.multi_solvent_ipea_wf import multi_solvent_ipea_fws
from rubicon.workflows.multistep_ipea_wf import multistep_ipea_fws
from rubicon.workflows.solvation_energy_wf import solvation_energy_fws


def snl_to_eg_wf(snl, parameters=None):
    fws_all = []
    parameters = parameters if parameters else {}

    snl_priority = parameters.get('priority', 1)
    mission = parameters.get('mission', 'Electron Genome Production')
    priority = snl_priority * 2  # once we start a job, keep going!

    f = Composition.from_formula(snl.structure.composition.reduced_formula).\
        alphabetical_formula
    molname = parameters.get("nick_name", f)

    # add the SNL to the SNL DB and figure out duplicate group
    tasks = [AddEGSNLTask()]
    spec = {'task_type': 'Add to SNL database',
            'snl': snl.to_dict,
            '_priority': snl_priority}
    if 'snlgroup_id' in parameters and isinstance(snl, EGStructureNL):
        spec['force_egsnl'] = snl.to_dict
        spec['force_snlgroup_id'] = parameters['snlgroup_id']
        del spec['snl']
    fws_all.append(FireWork(tasks, spec,
                        name=get_slug(molname + ' -- Add to SNL database'),
                        fw_id=1))

    default_solvents = ['diglym', 'acetonitrile', 'dmso', 'thf',
                        'dimethylamine', 'dimethoxyethane',
                        'dimethylaniline', 'tetraglyme']

    workflow_type = parameters.get('workflow', 'ipea')
    if workflow_type == 'ipea':
        fws_tasks, connections = multistep_ipea_fws(
            snl.structure, molname, mission, DupeFinderEG(), priority, 1)
    elif workflow_type == 'multiple solvent ipea':
        solvents = parameters.get('solvents', default_solvents)
        fws_tasks, connections = multi_solvent_ipea_fws(
            snl.structure, molname, mission, solvents, DupeFinderEG(),
            priority, 1)
    elif workflow_type == 'solvation energy':
        solvents = parameters.get('solvents', default_solvents)
        fws_tasks, connections = solvation_energy_fws(
            snl.structure, molname, mission, DupeFinderEG(), priority, 1,
            solvents)
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
