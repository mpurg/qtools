#########################
# py.test test functions
#########################

from __future__ import absolute_import
from __future__ import print_function
import re
import pytest
import shutil
import os

from Qpyl.qgeninp import genfeps, genrelax

class TestGenFeps:
    # TODO: needs more tests (see qscripts-cli test2)

    def test_genfeps(self, tmpdir):
        # make inputs in temporary directory and compare to reference
        try:
            cwd = os.getcwd()
            shutil.copytree("data/qgeninp/genfep.1", str(tmpdir.join("genfep.1")))
            tmpdir.join("genfep.1").chdir()
            genfeps("genfeps.proc", "relax_003.inp", "relax", "random_name.list",
                    51, 1, 0.5, "test_", False, pdb_file="probr_cl_start.pdb",
                    runscript_file="run_feps_q.sh")

            reference_files = os.listdir("test_000_ref")
            generated_files = os.listdir("test_000")
            for rf in reference_files:
                if rf not in generated_files:
                    print("Missing file: ", rf)
                    assert 0
                if ".re" in rf or ".top" in rf or ".sh" in rf:
                    continue
                genstr = open(os.path.join("test_000", rf)).read()
                refstr = open(os.path.join("test_000_ref", rf)).read()
                genstr = re.sub("(\*|\!|#).*", "", genstr)

                if genstr != refstr:
                    print("File '{}' is bad.".format(rf))
                    assert genstr == refstr
        except:
            raise
        finally:
            os.chdir(cwd)


class TestGenRelax:
    # TODO: needs more tests (see qscripts-cli test2)

    def test_genrelax(self, tmpdir):
        # make inputs in temporary directory and compare to reference
        try:
            cwd = os.getcwd()
            shutil.copytree("data/qgeninp/genrelax.1", str(tmpdir.join("genrelax.1")))
            tmpdir.join("genrelax.1").chdir()
            genrelax("genrelax.proc", "test_relax", "top",
                     top_file="probr_cl.top", fep_file="probr_cl.fep",
                     runscript_file="run_relax_q.sh",
                     pdb_file="probr_cl_start.pdb")

            reference_files = os.listdir("test_relax_ref")
            generated_files = os.listdir("test_relax")
            for rf in reference_files:
                if rf not in generated_files:
                    print("Missing file: ", rf)
                    assert 0
                if ".re" in rf or ".top" in rf or ".sh" in rf:
                    continue
                genstr = open(os.path.join("test_relax", rf)).read()
                refstr = open(os.path.join("test_relax_ref", rf)).read()
                genstr = re.sub("(\*|\!|#).*", "", genstr)

                if genstr != refstr:
                    print("File '{}' is bad.".format(rf))
                    assert genstr == refstr
        except:
            raise
        finally:
            os.chdir(cwd)
