"""
Build derived collections for Electrolyte Genome.
"""
import logging
import traceback
import time
import eg_molecules_builder as mols
from rubicon.builders import eg_shared

__author__ = 'Xiaohui Qu <xqu@lbl.gov>'
__date__ = '1/1/14'



_logname = 'eg.' + __name__

def run(colls, args):
    """Run with information from command-line args.
    """
    _log = logging.getLogger(_logname)
    ncores = args.num_cores
    # set default collection to be the real name of the tasks collection
    # (may not necessarily be 'tasks')
    colls.database.set_collection(colls.get_collection_name('tasks'))
    kw = dict(ncores=ncores)

    _log.info("create_builders.start")
    try:
        builders = [mols.MoleculesBuilder(colls, **kw)]
    except Exception, err:
        tb = traceback.format_exc()
        _log.error("create_builders.end.error msg={} traceback="
                   "{}".format(err, tb))
        raise
    _log.info("create_builders.end num={}".format(len(builders)))
    map(_run_one, builders)
    return 0


def _run_one(obj):
    _log = logging.getLogger(_logname)
    name = str(obj)
    _log.info("run.start builder={}".format(name))
    t0 = time.time()
    try:
        status = obj.run()
    except Exception, err:
        tb = traceback.format_exc()
        _log.error("run.end.error builder={} msg={}".format(name, tb))
        status = -1
    dur = time.time() - t0
    if status == 0:
        _log.info("run.end builder={} status=0 duration_sec="
                  "{:g}".format(name, dur))
    else:
        _log.error("run.end builder={} status={:d} duration_sec="
                   "{:g}".format(name, status, dur))
        raise eg_shared.BuildError(name, "unknown")   # TODO: get better error
        # message