import os
import subprocess
from custodian.custodian import Job, gzip_dir
import shutil

__author__ = 'xiaohuiqu'

class MopacJob(Job):
    """
    A basic MOPAC job.
    """

    def __init__(self, mopac_cmd, input_file="mol.mop", gzipped=False, backup=True):
        self.mopac_cmd = mopac_cmd
        self.input_file = input_file
        basename = os.path.splitext(self.input_file)[0]
        self.output_file = basename + ".out"
        self.arc_file = basename + ".arc"
        self.gzipped = gzipped
        self.backup = backup

    def setup(self):
        if self.backup:
            i = 0
            while os.path.exists("{}.{}.orig".format(self.input_file, i)):
                i += 1
            shutil.copy(self.input_file, "{}.{}.orig".format(self.input_file, i))
            if os.path.exists(self.output_file):
                shutil.copy(self.output_file, "{}.{}.orig".format(self.output_file, i))
            if os.path.exists(self.arc_file):
                shutil.copy(self.arc_file, "{}.{}.orig".format(self.arc_file, i))

    def run(self):
        cmd = [self.mopac_cmd, self.input_file]
        returncode = subprocess.call(cmd)
        return returncode

    def postprocess(self):
        if self.gzipped:
            gzip_dir(".")