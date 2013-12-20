from fireworks.features.dupefinder import DupeFinderBase

__author__ = 'xiaohuiqu'


class DupeFinderEG(DupeFinderBase):
    """
    General Duplication finder for Electrolyte Genome project
    """

    _fw_name = 'Dupe Finder EG'

    def verify(self, spec1, spec2):
        # assert: task_type and snlgroup_id have already been checked through
        # query
        comp_tags1 = dict()
        comp_tags2 = dict()
        comp_tags1['run_tags'] = spec1['run_tags']
        comp_tags1['implicit_solvent'] = spec1['implicit_solvent']
        comp_tags2['run_tags'] = spec2['run_tags']
        comp_tags2['implicit_solvent'] = spec2['implicit_solvent']
        return comp_tags1 == comp_tags2

    def query(self, spec):
        return {'spec.task_type': spec['task_type'],
                'spec.snlgroup_id': spec['snlgroup_id'],
                'spec.charge': spec['charge'],
                'spec.spin_multiplicity': spec['spin_multiplicity']}