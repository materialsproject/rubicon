"""
Run derived collection builder.

This should be run from the "egbuild" shell script.

See the -h option for details on the options.
Can run in two modes, 'merge' for building sandbox+core task collections,
and 'eg' for building the derived collections.
"""
from pymongo import MongoClient
from rubicon.builders import eg_builders_run, eg_shared

__author__ = 'Dan Gunter <dkgunter@lbl.gov>'
__date__ = '5/22/13'

# System
import argparse
import json
import logging
import os
import sys
import traceback


class ConfigurationError(Exception):
    def __init__(self, where, why):
        Exception.__init__(self, "Failed to load configuration {}: "
                                 "{}".format(where, why))


def tell_user(message):
    "Flexible way to control output for the user."
    print(message)


def get_settings(config_file):
    """Read settings from a configuration file.
    """
    try:
        if os.path.exists(config_file):
            return json.load(open(config_file))
        elif "DB_LOC" in os.environ and os.path.exists(os.environ["DB_LOC"]):
            with open(os.path.join(os.environ["DB_LOC"],
                                   config_file)) as f:
                return json.load(f)
        else:
            raise ValueError("configuration file not found")
    except Exception, err:
        raise ConfigurationError(config_file, err)

def get_db(**db_creds):
    conn = MongoClient(db_creds['host'], db_creds['port'])
    db = conn[db_creds['database']]
    if db_creds['admin_user']:
        db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
    return db

def main():
    # Setup args
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-c", "--config", dest="config_file", type=str,
                        metavar='FILE', default="tasks_db.json",
                        help="Configure database connection from FILE "
                             "(%(default)s)")
    common.add_argument("-p", "--prefix", dest="coll_prefix", type=str,
                        metavar='PREFIX', default=None,
                        help="Collection name prefix, ie namespace")
    common.add_argument('--verbose', '-v', dest='vb', action="count", default=0,
                        help="Print more verbose messages to standard error. "
                             "Repeatable. (default=ERROR)")
    common.add_argument('-P', '--ncores', dest="num_cores", type=int,
                        default=16,
                        help="Number of cores, thus number of threads to runin "
                             "parallel (%(default)d)")
    p = argparse.ArgumentParser(description="Build derived collections")
    subparsers = p.add_subparsers(description="Actions")
    # Action for VASP builder
    subp = subparsers.add_parser("eg",
                                 help="Build EG derived collections",
                                 parents=[common])
    subp.set_defaults(func=eg_builders_run.run, func_name='eg')
    subp.add_argument("-W", "--wipe", dest="wipe_target", action="store_true",
                      help="Wipe target collection, removing all data in it, "
                           "before merge")
    # Parse

    args = p.parse_args()

    # Configure logging

    _log = logging.getLogger("eg")  # parent
    _log.propagate = False
    hndlr = logging.StreamHandler()
    hndlr.setFormatter(logging.Formatter("[%(levelname)-6s] %(asctime)s "
                                         "%(name)s :: %(message)s"))
    _log.addHandler(hndlr)
    if args.vb > 1:
        lvl = logging.DEBUG
    elif args.vb > 0:
        lvl = logging.INFO
    else:
        lvl = logging.WARN
    _log.setLevel(lvl)
    _log = logging.getLogger("eg.build")
    # don't send logs up

    # Configure core database

    try:
        settings = get_settings(args.config_file)
    except ConfigurationError, err:
        p.error(str(err))
    core_db = get_db(**settings)
    suffix = None
    collections = eg_shared.Collections(core_db, task_suffix=suffix)

    # Run

    status = 0
    _log.info("run.start")
    if args.func_name == 'sandbox':
        pass
    else:
        # Build derived collection
        try:
            status = args.func(collections, args)
        except Exception, err:
            tb = traceback.format_exc()
            _log.error("Failed to run '{}': {}".format(args.func_name, tb))
            status = -1
    _log.info("run.end status={}".format(status))
    return status

if __name__ == '__main__':
    sys.exit(main())
