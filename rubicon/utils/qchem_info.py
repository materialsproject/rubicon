import os

from pkg_resources import parse_version


def get_qchem_version():
    if "QC" in os.environ:
        version_file_path = os.path.abspath(
            os.path.join(
                os.environ["QC"], "version.txt"))
        if os.path.exists(version_file_path):
            with open(version_file_path) as f:
                text = f.readlines()
            version = parse_version(text[0].strip())
            return version
        else:
            raise Exception("Can't find QChem version "
                            "file \"{}\"".format(version_file_path))
    else:
        raise Exception("Environmental Variable $QC is not defined")

