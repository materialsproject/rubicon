from fireworks.core.firework import FireWork, Workflow
from fireworks.utilities.fw_utilities import get_slug
from pymatgen import Composition
from rubicon.dupefinders.dupefinder_eg import DupeFinderEG
from rubicon.firetasks.egsnl_tasks import AddEGSNLTask
from rubicon.utils.snl.egsnl import EGStructureNL, get_meta_from_structure
from rubicon.workflows.multistep_ipea_wf import multistep_ipea_fws


def snl_to_eg_wf(snl, parameters=None):
    fws = []
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
    fws.append(FireWork(tasks, spec,
                        name=get_slug(molname + ' -- Add to SNL database'),
                        fw_id=1))

    ipea_fws, connections = multistep_ipea_fws(snl.structure, molname, mission,
                                               DupeFinderEG(), priority, 1)
    fws.extend(ipea_fws)

    wf_meta = get_meta_from_structure(snl.structure)
    wf_meta['run_version'] = 'Oct 29, 2013'

    if '_electrolytegenome' in snl.data and \
            'submission_id' in snl.data['_electrolytegenome']:
        wf_meta['submission_id'] = snl.data['_electrolytegenome'][
            'submission_id']
    return Workflow(fws, connections,
                    name=molname,
                    metadata=wf_meta)
