#########################
# py.test test functions
#########################

import re
import pytest

from Qpyl.qanalysis import QAnalyseFeps
from Qpyl import plotdata

class TestAnalyseFeps:
    def test_plotdata(self):
        # regression test, see if the data outputed is the same
        ref_values = open("data/qaf.PlotData.json").read().strip()

        qfep_outs = ["data/qfep.out.1",
                     "data/qfep.out.2"]
        qafs = QAnalyseFeps(qfep_outs, lra_lambdas=(1.0, 0.0))
        jsonenc = plotdata.PlotDataJSONEncoder(indent=2)
        assert jsonenc.encode(qafs.plotdata) == ref_values

    def test_lra_stats(self):
        # regression test, see if the data outputed is the same
        qfep_outs = ["data/qfep.out.1",
                     "data/qfep.out.2"]
        qafs = QAnalyseFeps(qfep_outs, lra_lambdas=(1.0, 0.0))
        assert  str(qafs.lra_stats).split() == """E_type      (E2-E1)_10_mean  (E2-E1)_10_std  (E2-E1)_01_mean  (E2-E1)_01_std  LRA_mean    LRA_std     REORG_mean  REORG_std
 EQtot                495.59          453.60          -307.57          136.64       94.01      158.48      401.58      295.12
 EQbond               223.15          233.99           -90.24           26.98       66.45      103.50      156.69      130.49
 EQang                164.00          111.67          -165.61          116.72       -0.81        2.52      164.80      114.20
 EQtor                  0.78            0.17            -7.38            9.72       -3.30        4.78        4.08        4.95
 EQimp                 -0.01            0.01             0.00            0.00       -0.01        0.01       -0.01        0.01
 EQel                  62.03           84.78             2.57           19.87       32.30       52.33       29.73       32.46
 EQvdW                 45.73           22.90           -47.00            2.98       -0.63        9.96       46.36       12.94
 Eel_qq                36.90           51.12             1.83           17.50       19.36       34.31       17.53       16.81
 EvdW_qq               44.93           25.08           -45.86            4.68       -0.47       10.20       45.39       14.88
 Eel_qp                20.77           27.51             9.51           14.79       15.14       21.15        5.63        6.36
 EvdW_qp                0.77            2.23            -1.38            1.36       -0.31        0.43        1.07        1.79
 Eel_qw                 4.36            6.16            -8.79           12.42       -2.21        3.13        6.57        9.29
 EvdW_qw                0.04            0.05             0.25            0.35        0.14        0.20       -0.11        0.15
 Eqrstr                -0.07            0.10             0.07            0.10
 0.00        0.00       -0.07        0.10""".split()

    def test_dg_all(self):
        # regression test, see if the data outputed is the same
        qfep_outs = ["data/qfep.out.1",
                     "data/qfep.out.2"]
        qafs = QAnalyseFeps(qfep_outs, lra_lambdas=(1.0, 0.0))
        assert str(qafs.dg_all).split() == """Qfep_output  dG*         dG0         dG_lambda   QCP_dG*     QCP_dG0     QCP_dG_lambda  QCP_mass_dG*  QCP_mass_dG0  QCP_mass_dG_lambda  ex_el_49_200_221_dG*  ex_el_49_200_221_dG0  ex_el_49_200_221_dG_lambda  ex_full_200_dG*  ex_full_200_dG0  ex_full_200_dG_lambda  ex_full_221_dG*  ex_full_221_dG0  ex_full_221_dG_lambda  ex_full_49_dG*  ex_full_49_dG0  ex_full_49_dG_lambda  ex_full_49_200_221_dG*  ex_full_49_200_221_dG0  ex_full_49_200_221_dG_lambda  ex_vdw_49_200_221_dG*  ex_vdw_49_200_221_dG0  ex_vdw_49_200_221_dG_lambda
 data/qfep.out.1       12.21       -6.44      -10.20  None        None        None           None          None          None                None                  None                  None                        None             None             None                   None             None             None                   None            None            None                  None                    None                    None                          None                   None                   None
 data/qfep.out.2       23.98      -21.98      -46.71       23.31      -22.12         -48.33         23.60        -21.95              -47.31                 10.64                -47.71                      -83.24            25.82           -15.95                 -38.98             9.77           -52.08                 -88.85           23.28          -23.45                -48.75                   10.69                  -47.67                        -83.16                  24.04                 -21.94                       -46.63""".split()
