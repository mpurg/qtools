#########################
# py.test test functions
#########################

import pytest
import os
import re

from Qpyl.qgroupcontrib import QGroupContrib, QGroupContribError

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

class TestQGroupContrib:
    def test_calcall(self):

        try:
            qcalc = os.path.join(os.environ["QBIN_DIR"], "Qcalc6")
        except KeyError:
            raise Exception("QBIN_DIR environment variable not defined.")

        dirs = ["data/qgroupcontrib/testrep1",
                "data/qgroupcontrib/testrep2",
                "data/qgroupcontrib/testrep3"]

        pdb = "data/qgroupcontrib/dfpase_dfp_start.pdb"

        qgc = QGroupContrib(qcalc, dirs, pdb, "q_enfiles.list",
                            [0.86, 0.14], [0.46, 0.54], 15, 20, 1, 1, None)
        try:
            qgc.calcall()
        except KeyboardInterrupt:
            qgc.kill_event.set()
            raise

        de1, de2, lras, reorgs = qgc.gcs_stats.get_columns([5, 9, 13, 17])
        # Glu19
        # qcalc manual calculation + pen&paper
        # REP1: 
        #    E1_1 = 79.51, E2_1 = 68.34, E1_2 = 79.42, E2_2 = 65.88
        #    de1 = -11.17, de2 = -13.54, lra = -12.355, reorg = 1.185
        # REP2:
        #    E1_1 = 76.02, E2_1 = 64.62, E1_2 = 79.88, E2_2 = 67.55
        #    de1 = -11.4, de2 = -12.33, lra = -11.865, reorg = 0.465
        # means:
        #    de1 = -11.285, de2 = -12.935, lra = -12.11, reorg = 0.825
        assert is_close(de1[4], -11.285)
        assert is_close(de2[4], -12.935)
        assert is_close(lras[4], -12.11)
        assert is_close(reorgs[4], 0.825)


    def test_calcall2(self):

        try:
            qcalc = os.path.join(os.environ["QBIN_DIR"], "Qcalc6")
        except KeyError:
            raise Exception("QBIN_DIR environment variable not defined.")

        dirs = ["data/qgroupcontrib/testrep1",
                "data/qgroupcontrib/testrep2",
                "data/qgroupcontrib/testrep3"]

        pdb = "data/qgroupcontrib/dfpase_dfp_start.pdb"

        # set the iscale to 4.0
        qgc = QGroupContrib(qcalc, dirs, pdb, "q_enfiles.list",
                            [0.86, 0.14], [0.46, 0.54], 15, 20, 4.0, 1, None)
        try:
            qgc.calcall()
        except KeyboardInterrupt:
            qgc.kill_event.set()
            raise

        de1, de2, lras, reorgs = qgc.gcs_stats.get_columns([5, 9, 13, 17])
        try:
            qgc.calcall()
        except KeyboardInterrupt:
            qgc.kill_event.set()
            raise

        de1, de2, lras, reorgs = qgc.gcs_stats.get_columns([5, 9, 13, 17])
        # Glu19 scaled down by 4 (see above)
        assert is_close(de1[4], -2.82125)
        assert is_close(de2[4], -3.23375)
        assert is_close(lras[4], -3.0275)
        assert is_close(reorgs[4], 0.20625)
        

    def test_calcall3(self):

        try:
            qcalc = os.path.join(os.environ["QBIN_DIR"], "Qcalc6")
        except KeyError:
            raise Exception("QBIN_DIR environment variable not defined.")

        dirs = ["data/qgroupcontrib/testrep1",
                "data/qgroupcontrib/testrep2",
                "data/qgroupcontrib/testrep3"]

        pdb = "data/qgroupcontrib/dfpase_dfp_start.pdb"

        qgc = QGroupContrib(qcalc, dirs, pdb, "q_enfiles.list",
                            [0.86, 0.14], [0.46, 0.54], 15, 20, 1, 1,
                            [4841, 4842, 4843, 4844])
        try:
            qgc.calcall()
        except KeyboardInterrupt:
            qgc.kill_event.set()
            raise

        de1, de2, lras, reorgs = qgc.gcs_stats.get_columns([5, 9, 13, 17])
        assert is_close(de1[4], -0.475)
        assert is_close(de2[4], -0.84)
        assert is_close(lras[4], -0.6575)
        assert is_close(reorgs[4], 0.1825)
        

