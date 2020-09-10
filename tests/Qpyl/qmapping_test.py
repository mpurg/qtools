#########################
# py.test test functions
#########################

from __future__ import absolute_import, unicode_literals
from __future__ import print_function, division
from io import open

import pytest
import os
import re

from Qpyl.qmapping import QMapper, QMapperError
import six

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

        for md, (qfep_inp, qfep_out) in six.iteritems(qmapper.mapped):
            fn = "data/qmapping/qfep_" + os.path.basename(md) + ".inp"
            #open(fn, "w").write(qfep_inp)
            assert open(fn, "r").read() == qfep_inp


            qfep_out = re.sub(".*Number of energy files",
                              "Number of energy files",
                              qfep_out, 0, re.DOTALL)

            fn = "data/qmapping/qfep_" + os.path.basename(md) + ".out"

            #open(fn, "w").write(qfep_out)
            qfep_out = re.sub(" +", " ", qfep_out)

            ref_out = open(fn, "r").read()
            ref_out = re.sub(" +", " ", ref_out)

            assert ref_out == qfep_out

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

        assert is_close(qmapper.parms["hij"], 82.120105)
        assert is_close(qmapper.parms["alpha"], 4.567082)
        assert is_close(qmapper.parms["gas_const"], 0.0019872041)
        assert is_close(qmapper.parms["temperature"], 298.0)
        assert qmapper.parms["gap_bins"] == 50
        assert qmapper.parms["points_skip"] == 1
        assert qmapper.parms["minpts_bin"] == 1

