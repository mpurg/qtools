#########################
# py.test test functions
#########################

import pytest
import os
import re

from Qpyl.qmapping import QMapper, QMapperError

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

class TestQMapper:
    def test_mapall(self):

        try:
            qfep = os.path.join(os.environ["QBIN_DIR"], "Qfep6")
        except KeyError:
            raise Exception("QBIN_DIR environment variable not defined.")

        dirs = ["data/qmapping/testrep1",
                "data/qmapping/testrep2",
                "data/qmapping/testrep3"]

        qmapper_parms = {"mapdirs": dirs,
                         "hij": 100.0,
                         "alpha": 10.0,
                         "nthreads": 1,
                         "temperature": 298.0,
                         "points_skip": 1,
                         "minpts_bin": 1,
                         "gap_bins": 20,
                         "qfep_exec": qfep,
                         "en_list_fn": "q_enfiles.list",
                         "gas_const": 0.0019872041}

        qmapper = QMapper(**qmapper_parms)
        try:
            qmapper.mapall()
        except KeyboardInterrupt:
            qmapper.kill_event.set()
            raise

        assert qmapper.input_parms_str == \
            "q_mapper.py 100.0 10.0 --bins 20 --skip 1 --min 1 --temp 298.0 "

        for md, (qfep_inp, qfep_out) in qmapper.mapped.iteritems():
            fn = "data/qmapping/qfep_" + os.path.basename(md) + ".inp"
            #open(fn, "w").write(qfep_inp)
            assert open(fn, "r").read() == qfep_inp

            fn = "data/qmapping/qfep_" + os.path.basename(md) + ".out"
            #open(fn, "w").write(qfep_out)
            assert open(fn, "r").read() == re.sub("Current .*", "", qfep_out)

    def test_fit_to_reference(self):

        try:
            qfep = os.path.join(os.environ["QBIN_DIR"], "Qfep6")
        except KeyError:
            raise Exception("QBIN_DIR environment variable not defined.")

        dirs = ["data/qmapping/testrep1",
                "data/qmapping/testrep2",
                "data/qmapping/testrep3"]

        qmapper_parms = {"mapdirs": dirs,
                         "hij": 0.0,
                         "alpha": 0.0,
                         "nthreads": 1,
                         "temperature": 298.0,
                         "points_skip": 1,
                         "minpts_bin": 1,
                         "gap_bins": 50,
                         "qfep_exec": qfep,
                         "en_list_fn": "q_enfiles.list",
                         "gas_const": 0.0019872041}

        qmapper = QMapper(**qmapper_parms)

        try:
            qmapper.fit_to_reference(12, -5, 10, 0.001, 10)
        except KeyboardInterrupt:
            qmapper.kill_event.set()
            raise

        assert qmapper.input_parms_str == "q_mapper.py 82.1201056503 " \
                "4.56708284345 --bins 50 --skip 1 --min 1 --temp 298.0 "

        assert is_close(qmapper.parms["hij"], 82.1201056503)
        assert is_close(qmapper.parms["alpha"], 4.56708284345)

