#########################
# py.test test functions
#########################


from __future__ import absolute_import
import pytest
import shutil
import os
from Qpyl.common import np, backup_file
from six.moves import range

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def test_backupfile(tmpdir):
    shutil.copy("data/qdyn.inp.1", tmpdir.dirname)
    fn = backup_file(os.path.join(tmpdir.dirname, "qdyn.inp.1"))
    inpstr = open("data/qdyn.inp.1").read()
    inpstr2 = open(os.path.join(tmpdir.dirname, fn)).read()
    assert inpstr == inpstr2


class TestStats:
    def test_mean(self):
        vals = [i**2 for i in range(1, 21)]
        assert is_close(np.mean(vals), 143.5)

    def test_median(self):
        vals = [i**2 for i in range(1, 21)]
        assert is_close(np.median(vals), 110.5)

    def test_std(self):
        vals = [i**2 for i in range(1, 21)]
        assert is_close(np.std(vals), 127.9023064686)

    def test_sem(self):
        vals = [i**2 for i in range(1, 21)]
        assert is_close(np.std_error(vals), 28.5998251743)
