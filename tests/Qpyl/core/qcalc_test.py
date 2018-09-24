#########################
# py.test test functions
#########################

from __future__ import absolute_import
import pytest

from Qpyl.core.qcalc import QCalcOutput

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestQFepOutput:

    @pytest.fixture
    def qco1(self):
        qco_str = open("data/qcalc.out.1", "r").read()
        return QCalcOutput(qco_str)

    @pytest.fixture
    def qco2(self):
        qco_str = open("data/qcalc.out.2", "r").read()
        return QCalcOutput(qco_str)

    def test_header(self, qco1):
        assert qco1.qcalc_version == "5.10.7"

    def test_gc1(self, qco1):
        lj, el = qco1.results["gc"].get_columns(["E_LJ", "E_EL"])
        assert is_close(sum(lj), -10.40)
        assert is_close(sum(el), -212.96)

    def test_gc2(self, qco1):
        els = dict(qco1.results["gc"].get_rows(columns=["Residue", "E_EL"]))
        assert is_close(els[35], 24.27)
        assert is_close(els[227], -13.65)
        assert is_close(els[286], 22.21)

    def test_calcs(self, qco2):
        assert is_close(sum(qco2.results["1"].get_columns()[1]), 14.189316)
        assert is_close(sum(qco2.results["2"].get_columns()[1]), 16.44)
        assert is_close(sum(qco2.results["3"].get_columns()[1]), 32.81)
        assert is_close(sum(qco2.results["4"].get_columns()[1]), 179.94)
        assert is_close(sum(qco2.results["5"].get_columns()[1]), 1183.01)

