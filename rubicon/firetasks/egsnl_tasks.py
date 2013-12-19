from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen.matproj.snl import StructureNL
from rubicon.utils.snl.egsnl_mongo import EGSNLMongoAdapter


class AddEGSNLTask(FireTaskBase, FWSerializable):
    """
    Add a new SNL into the SNL database, and build duplicate groups
    """

    _fw_name = "Add EG SNL Task"

    def run_task(self, fw_spec):
        # pass-through option for when we start with an egsnl and don't
        # actually want to add
        if 'force_egsnl' in fw_spec and 'force_snlgroup_id' in fw_spec:
            print 'USING FORCED EGSNL'
            return FWAction(update_spec={'egsnl': fw_spec['force_egsnl'],
                                         'snlgroup_id':
                                         fw_spec['force_snlgroup_id']})

        sma = EGSNLMongoAdapter.auto_load()
        snl = StructureNL.from_dict(fw_spec['snl'])
        egsnl, snlgroup_id = sma.add_snl(snl)

        return FWAction(update_spec={'egsnl': egsnl.to_dict,
                                     'snlgroup_id': snlgroup_id})