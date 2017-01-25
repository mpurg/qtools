#########################
# py.test test functions
#########################

import pytest

from Qpyl.core.qcalc import QCalcOutput


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
        assert abs(sum(lj) - -10.40) < 1e-7
        assert abs(sum(el) - -212.96) < 1e-7

    def test_gc2(self, qco1):
        els = dict(qco1.results["gc"].get_rows(columns=["Residue", "E_EL"]))
        assert abs(els[35] - 24.27) < 1e-7
        assert abs(els[227] - -13.65) < 1e-7
        assert abs(els[286] - 22.21) < 1e-7

    def test_calcs(self, qco2):
        assert abs(sum(qco2.results["1"].get_columns()[1]) - 14.189316) < 1e-7
        assert abs(sum(qco2.results["2"].get_columns()[1]) - 16.44) < 1e-7
        assert abs(sum(qco2.results["3"].get_columns()[1]) - 32.81) < 1e-7
        assert abs(sum(qco2.results["4"].get_columns()[1]) - 179.94) < 1e-7
        assert abs(sum(qco2.results["5"].get_columns()[1]) - 1183.01) < 1e-7

