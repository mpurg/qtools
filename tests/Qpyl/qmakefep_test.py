#########################
# py.test test functions
#########################


import re
import pytest

from Qpyl.qmakefep import make_fep, QMakeFepError

class TestMakeFep:
    def test_makefep(self):
        libs = ["data/qmakefep/3hp.lib", "data/qmakefep/3h2.lib"]
        fep_file_str = make_fep("data/qmakefep/3hp.qmap",
                                "data/qmakefep/3hp_start.pdb",
                                "oplsaa",
                                ["data/qmakefep/3hp.prm",],
                                libs)
        fep_file_str = re.sub("(\*|\!|#).*", "", fep_file_str)
        ref_fep = open("data/qmakefep/3hp.fep.tmplt.gen").read()
        assert fep_file_str == ref_fep

