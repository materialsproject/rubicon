"""
Shared code for EG builders.

Adapted from Dan Gynter and Wei Chen's Materials Genome builders

For developers implementing a new builder,
you should inherit from `ParallelBuilder`. See the
documentation of this class for details.
"""
__author__ = 'Xiaohui Qu <xqu@lbl.gov>'
__date__ = '12/31/13'


import logging
import Queue
import threading
import multiprocessing

_log = logging.getLogger("eg.build")


class BuildError(Exception):
    def __init__(self, who, why):
        errmsg = "Builder {} failed: {}".format(who, why)
        Exception.__init__(self, errmsg)


class Collections:
    """Interface to normalized names for collections.
    
    After initialization with a MongoDB database and optional parameters,
    you can access collections in `known_collections` as attributes.
    """
    #: Collection names that are accessible as attributes
    #: of an instance of this class
    known_collections = ['tasks', 'molecules', 'reactions']

    MIN_VER = 1
    MAX_VER = 1

    def __init__(self, db, version=1, prefix=None, task_suffix=None):
        """Set collections from database.

        db:
            MongoDB database, but really anything that acts like a dict
            Type (dict-like object)
        version:
            Version of naming scheme for the collections
            Type (int)
        prefix:
            Prefix string to put before collection names,
            e.g. "dahn". Full
            collection name will be <prefix>.<name>; don't include
            '.' in the input.
            Type (str)
        task_suffix:
            Add this suffix to the tasks collection. Used for merged
            collections for sandboxes.
            Type (str)

        raise ValueError if `version` is not known
        """
        if not self.MIN_VER <= version <= self.MAX_VER:
            raise ValueError("Bad version ({v:d}) not in range {v0} .. {v1} ".
                             format(v=version, v0=self.MIN_VER,
                                    v1=self.MAX_VER))
        self._names, self._coll = {}, {}
        if version == 1:
            for name in self.known_collections:
                full_name = "{}.{}".format(prefix, name) if prefix else name
                if name == 'tasks' and task_suffix is not None:
                    full_name = "{}.{}".format(full_name, task_suffix)
                self._names[name] = full_name
                self._coll[full_name] = None
        self._db = db
        self._prefix = prefix

    def __getattr__(self, item):
        if item in self._names:
            coll_name = self._names[item]
            coll = self._coll[coll_name]
            if coll is None:
                self._coll[coll_name] = coll = self._db[coll_name]
            return coll
        return self.__dict__[item]

    def get_collection_name(self, alias):
        return self._names[alias]

    @property
    def database(self):
        """Return the current database object.
        """
        return self._db


class Builder(object):
    @staticmethod
    def combine_status(codes):
        """Combine integer status codes.

        Return: -1 if any is nonzero, else 0
        """
        return (-1, 0)[not filter(None, codes)]

    def __str__(self):
        return self.__class__.__name__


class ParallelBuilder(Builder):
    """Parallel builder base class.
        All the builder classes should inherit from this.

    Subclasses should define two methods:

    1) run() - Do any serial tasks, then put all objects that must be operated
        on in parallel into the queue, one at a time, by calling `add_item`.
        Then call `run_parallel`. The return value of this function is a list
        of status codes, one per thread/process. You can join the status codes
        with `combine_status`.

    2) process_item(item) - Do the work for one item. Return integer status
            code, 0 for OK and non-zero for failure.
    """
    def __init__(self, ncores=0, threads=False):
        """Create new builder for threaded or multiprocess execution.

        Args:
            ncores: Desired number of threads/processes to run
                Type (int)
            threads: Use threads (True) or processes (False)
                Type (bool)
        """
        self._ncores = ncores if ncores else 15
        self._threaded = threads
        if threads:
            self._queue = Queue.Queue()
            self._states_lock = threading.Lock()
            self._run_parallel_mode = self._run_parallel_threaded
        else:
            self._queue = multiprocessing.Queue()
            self._run_parallel_mode = self._run_parallel_multiprocess
        self._states = []

    def run(self):
        """This method should:
        (1) put work in the queue with add_item()
        (2) call run_parallel() to do the work

        Returns: Status code, 0 for OK
            Type (int)
        """
        raise NotImplementedError()

    def add_item(self, item):
        """Put an item of work in the queue.
        Subclasses do not need to modify this.
        """
        self._queue.put(item)

    def process_item(self, item):
        """Implement the analysis for each item of work here.

        Args:
            item: One item of work from the queue
                Type (object)

        Returns: Status code, 0 for OK
            Type (int)
        """
        raise NotImplementedError()

    def run_parallel(self):
        """Called to run threads, once queue is filled.

        Returns: Multiple integer status codes, 1 per thread
            Type (list of int)
        """
        return self._run_parallel_mode()

    def _run_parallel_threaded(self):
        """Called to run threads, once queue is filled.

        Returns: array of integer status codes, 1 per thread
        """
        _log.debug("run.parallel.threaded.start")
        threads = []
        for i in xrange(self._ncores):
            thr = threading.Thread(target=self._run)
            thr.start()
            threads.append(thr)
        for i in xrange(self._ncores):
            threads[i].join()
            if threads[i].isAlive():    # timed out
                _log.warn("run.parallel.threaded: timeout for thread="
                          "{:d}".format(i))
        _log.debug("run.parallel.threaded.end states="
                   "{}".format(self._set_status()))
        return self._set_status()

    def _run_parallel_multiprocess(self):
        """Called to run processes, once queue is filled.

        Returns: array of integer status codes, 1 per thread
        """
        _log.debug("run.parallel.multiprocess.start")
        processes = []
        ProcRunner.instance = self
        for i in xrange(self._ncores):
            proc = multiprocessing.Process(target=ProcRunner.run)
            proc.start()
            processes.append(proc)
        states = []
        for i in xrange(self._ncores):
            processes[i].join()
            states.append(processes[i].exitcode)
        _log.debug("run.parallel.multiprocess.end states="
                   "{}".format(','.join(map(str, states))))
        return states

    def _run(self):
        """Run method for one thread or process
        Just pull an item off the queue and process it,
        until the queue is empty.
        """
        while 1:
            try:
                item = self._queue.get(timeout=2)
                self.process_item(item)
            except Queue.Empty:
                break
            except Exception, err:
                _log.error("Processing exception: {}".format(err))
                self._set_status(-1)
                raise
        self._set_status(0)

    def _set_status(self, value=None):
        if self._threaded:
            self._states_lock.acquire()
            if value is not None:
                self._states.append(value)
            result = self._states[:]
            self._states_lock.release()
            return result
        else:
            return value


class ProcRunner:
    """This is a work-around to the limitation of multiprocessing that the
    function executed in the new module cannot be a method of a class.
    We simply set the instance (self) into the class before forking each
    process, and the class' method calls the instance method for us.
    """

    def __init__(self):
        pass

    instance = None

    @classmethod
    def run(cls):
        # noinspection PyProtectedMember
        cls.instance._run()
