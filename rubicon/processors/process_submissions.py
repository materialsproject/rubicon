import time
import traceback
from fireworks.core.launchpad import LaunchPad
from pymatgen.matproj.snl import StructureNL
from rubicon.submission.submission_mongo import SubmissionMongoAdapterJCESR
from rubicon.utils.snl.egsnl import EGStructureNL
from rubicon.workflows.snl_to_wf import snl_to_wf


class SubmissionProcessor():
    MAX_SITES = 200

    # This is run on the server end
    def __init__(self, sma, launchpad):
        self.sma = sma
        self.jobs = sma.jobs
        self.launchpad = launchpad

    def run(self, sleep_time=None, infinite=False):
        sleep_time = sleep_time if sleep_time else 30
        while True:
            self.submit_all_new_workflows()
            self.update_existing_workflows()
            if not infinite:
                break
            print 'sleeping', sleep_time
            time.sleep(sleep_time)

    def submit_all_new_workflows(self):
        last_id = -1
        while last_id:
            last_id = self.submit_new_workflow()

    def submit_new_workflow(self):
        # finds a submitted job, creates a workflow, and submits it to FireWorks
        job = self.jobs.find_and_modify({'state': 'SUBMITTED'},
                                        {'$set': {'state': 'WAITING'}})
        if job:
            submission_id = job['submission_id']
            # noinspection PyBroadException
            try:
                if 'snl_id' in job:
                    snl = EGStructureNL.from_dict(job)
                else:
                    snl = StructureNL.from_dict(job)
                if len(snl.structure.sites) > SubmissionProcessor.MAX_SITES:
                    self.sma.update_state(submission_id, 'REJECTED',
                                          'too many sites', {})
                    print 'REJECTED WORKFLOW FOR {} - too many sites ' \
                          '({})'.format(snl.structure.formula,
                                        len(snl.structure.sites))
                elif not job['is_valid']:
                    self.sma.update_state(submission_id, 'REJECTED',
                                          'invalid structure (atoms too close)',
                                          {})
                    print 'REJECTED WORKFLOW FOR {} - invalid ' \
                          'structure'.format(snl.structure.formula)
                else:
                    snl.data['_electrolytegenome'] = \
                        snl.data.get('_electrolytegenome', {})
                    snl.data['_electrolytegenome']['submission_id'] \
                        = submission_id

                    # create a workflow
                    wf = snl_to_wf(snl, job['parameters'])
                    self.launchpad.add_wf(wf)
                    print 'ADDED WORKFLOW FOR {}'.format(snl.structure.formula)
            except:
                self.jobs.find_and_modify({'submission_id': submission_id},
                                          {'$set': {'state': 'ERROR'}})
                traceback.print_exc()

            return submission_id

    def update_existing_workflows(self):
        # updates the state of existing workflows by querying the FireWorks
        # database
        for submission in self.jobs.find({'state': {'$nin': ['COMPLETED',
                                                             'ERROR',
                                                             'REJECTED']}},
                                         {'submission_id': 1}):
            submission_id = submission['submission_id']
            # noinspection PyBroadException
            try:
                # get a wf with this submission id
                fw_id = self.launchpad.get_wf_ids({'metadata.submission_id':
                                                   submission_id},
                                                  limit=1)[0]
                # get a workflow
                wf = self.launchpad.get_wf_by_fw_id(fw_id)
                # update workflow
                self.update_wf_state(wf, submission_id)
            except:
                print 'ERROR while processing s_id', submission_id
                traceback.print_exc()

    def update_wf_state(self, wf, submission_id):
        # state of the workflow

        details = '(none available)'
        for fw in wf.fws:
            if fw.state == 'READY':
                details = 'waiting to run: {}'.format(fw.spec['task_type'])
            elif fw.state in ['RESERVED', 'RUNNING', 'FIZZLED']:
                machine_name = 'unknown'
                for l in fw.launches:
                    if l.state == fw.state:
                        machine_name = 'unknown'
                        if 'hopper' in l.host or 'nid' in l.host:
                            machine_name = 'hopper'
                        elif 'c' in l.host:
                            machine_name = 'mendel/carver'
                        break
                if fw.state == 'RESERVED':
                    details = 'queued to run: {} on {}'.format(
                        fw.spec['task_type'], machine_name)
                if fw.state == 'RUNNING':
                    details = 'running: {} on {}'.format(
                        fw.spec['task_type'], machine_name)
                if fw.state == 'FIZZLED':
                    details = 'fizzled while running: {} on {}'.format(
                        fw.spec['task_type'], machine_name)

        m_taskdict = {}
        states = [fw.state for fw in wf.fws]
        if any([s == 'COMPLETED' for s in states]):
            for fw in wf.fws:
                if fw.state == 'COMPLETED':
                    for l in fw.launches:
                        # if task_id is not there, it means we went on DETOUR...
                        if l.state == 'COMPLETED' and \
                                      'task_id' in l.action.stored_data:
                            t_id = l.action.stored_data['task_id']
                            if 'task_type' in fw.spec:
                                m_taskdict[fw.spec['task_type']] = t_id
                            break

        self.sma.update_state(submission_id, wf.state, details, m_taskdict)
        return wf.state, details, m_taskdict

    @classmethod
    def auto_load(cls):
        sma = SubmissionMongoAdapterJCESR.auto_load()
        lp = LaunchPad.auto_load()

        return SubmissionProcessor(sma, lp)