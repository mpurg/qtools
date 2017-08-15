#########################
# py.test test functions
#########################

import pytest

from Qpyl.core.qfep import QFepOutput, QFepOutputError

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestQFepOutput:

    @pytest.fixture(scope='session')
    def qfo1(self):
        qfo_str = open("data/qfep.out.1", "r").read()
        return QFepOutput(qfo_str)

    @pytest.fixture(scope='session')
    def qfo2(self):
        qfo_str = open("data/qfep.out.2", "r").read()
        return QFepOutput(qfo_str)

    def test_bad_output(self):
        with pytest.raises(QFepOutputError):
            QFepOutput(open("data/qcalc.out.1", "r").read())

    def test_header(self, qfo1):
        # no point in having multiple tests for this
        assert qfo1.header.qfep_version == "5.10.1"
        assert qfo1.header.qfep_moddate == "20151002"
        assert str(qfo1.header.Hij) == "76.54"
        assert str(qfo1.header.alpha) == "2.35"
        assert str(qfo1.header.kT) == "0.596"
        assert qfo1.header.bins == 51
        assert qfo1.header.nfiles == 51
        assert qfo1.header.nstates == 2
        assert qfo1.header.pts_skip == 50
        assert qfo1.header.min_pts_bin == 10


    def test_part0(self, qfo1):
        # EQtotal
        eq1_first = qfo1.part0.data_state[0].get_columns(["EQtot"])[0][0]
        eq2_last = qfo1.part0.data_state[1].get_columns(["EQtot"])[0][-1]
        assert str(eq1_first) == "-74.59"
        assert str(eq2_last) == "-87.63"

    def test_part1(self, qfo1):
        sum_dg = sum(qfo1.part1.data.get_columns(["dG"])[0])
        assert is_close(sum_dg, 1651.77)

    def test_dgfep1(self, qfo1):
        assert is_close(qfo1.part1.dg, -10.204)

    def test_lra1(self, qfo1):
        # LRA of EQbond (1.0 -> 0.0)
        # pen&paper
        # E1_10 = 1.0*-65.90 + 0.0*-8.21
        # E2_10 = 0.0*-65.90 + 1.0*-8.21
        # E1_01 = 1.0*-6.81  + 0.0*-77.97
        # E2_01 = 0.0*-6.81  + 1.0*-77.97
        # E = 0.5*( <E2-E1>_10 + <E2-E1>_01 )
        # E = 0.5*( (-8.21 --65.90) + (-77.97 --6.81) )
        # E = -6.7350
        lras = qfo1.part0.calc_lra(1.0, 0.0)
        print lras
        lra_eqbond = lras.get_columns(["LRA"])[0][1]
        assert is_close(lra_eqbond, -6.7350)

    def test_lra2(self, qfo2):
        lras = qfo2.part0.calc_lra(1.0, 0.0)
        print lras
        lras = " ".join(["{:.2f}".format(x) for x in lras.get_columns()[3]])
        expected = "206.07 139.64 -2.59 -6.68 -0.01 69.30 6.41 "\
                   "43.62 6.75 30.10 -0.61 -4.43 0.29 0.00"
        assert lras == expected

    def test_lra3(self, qfo1):
        # LRA of EQbond (0.9 -> 0.6)
        # pen&paper
        # E1_10 = 0.9*-65.52 + 0.1*-11.18 = -60.086
        # E2_10 = 0.6*-65.52 + 0.4*-11.18 = -43.784
        # E1_01 = 0.9*-56.26 + 0.1*-28.54 = -53.488
        # E2_01 = 0.6*-56.26 + 0.4*-28.54 = -45.172
        # E = 0.5*( <E2-E1>_10 + <E2-E1>_01 )
        # E = 0.5*( (-43.784 --60.086) + (-45.172 --53.488) )
        # E = 12.309
        lras = qfo1.part0.calc_lra(0.9, 0.6)
        print lras
        lra_eqbond = lras.get_columns(["LRA"])[0][1]
        assert is_close(lra_eqbond, 12.309)

    def test_lra4(self, qfo2):
        lras = qfo2.part0.calc_lra(0.9, 0.6)
        print lras
        lras = " ".join(["{:.2f}".format(x) for x in lras.get_columns()[3]])
        expected = "187.28 89.52 51.98 0.27 -0.00 30.78 14.74 19.51 14.87 "\
                   "10.84 -0.20 0.42 0.07 0.00"
        assert lras == expected

    def test_dgs1(self, qfo2):
        assert is_close(qfo2.part3.dga, 23.98)
        assert is_close(qfo2.part3.dg0, -21.98)

    def test_dgs2(self, qfo2):
        sub_calc = qfo2.sub_calcs["ex_full_49"]
        assert is_close(sub_calc.part3.dga, 23.28)
        assert is_close(sub_calc.part3.dg0, -23.45)

    def test_dgs3(self, qfo2):
        sub_calc = qfo2.sub_calcs["ex_full_200"]
        assert is_close(sub_calc.part3.dga, 25.82)
        assert is_close(sub_calc.part3.dg0, -15.95)

    def test_dgs4(self, qfo2):
        sub_calc = qfo2.sub_calcs["ex_full_221"]
        assert is_close(sub_calc.part3.dga, 9.77)
        assert is_close(sub_calc.part3.dg0, -52.08)

    def test_dgs5(self, qfo2):
        sub_calc = qfo2.sub_calcs["ex_full_49_200_221"]
        assert is_close(sub_calc.part3.dga, 10.69)
        assert is_close(sub_calc.part3.dg0, -47.67)

    def test_dgs6(self, qfo2):
        sub_calc = qfo2.sub_calcs["ex_el_49_200_221"]
        assert is_close(sub_calc.part3.dga, 10.64)
        assert is_close(sub_calc.part3.dg0, -47.71)

    def test_dgs7(self, qfo2):
        sub_calc = qfo2.sub_calcs["ex_vdw_49_200_221"]
        assert is_close(sub_calc.part3.dga, 24.04)
        assert is_close(sub_calc.part3.dg0, -21.94)

    def test_dgs8(self, qfo2):
        sub_calc = qfo2.sub_calcs["QCP"]
        assert is_close(sub_calc.part3.dga, 23.31)
        assert is_close(sub_calc.part3.dg0, -22.12)

    def test_dgs9(self, qfo2):
        sub_calc = qfo2.sub_calcs["QCP_mass"]
        assert is_close(sub_calc.part3.dga, 23.60)
        assert is_close(sub_calc.part3.dg0, -21.95)











