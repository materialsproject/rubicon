from fireworks.features.dupefinder import DupeFinderBase

__author__ = 'xiaohuiqu'


class DupeFinderEG(DupeFinderBase):
    """
    General Duplication finder for Electrolyte Genome project
    """

    _fw_name = 'Dupe Finder EG'

    def verify(self, spec1, spec2):
        # assert: task_type and snlgroup_id have already been checked through query
        return spec1['user_tags'] == spec2['user_tags']

    def query(self, spec):
        return {'spec.task_type': spec['task_type'],
                'spec.snlgroup_id': spec['snlgroup_id']}