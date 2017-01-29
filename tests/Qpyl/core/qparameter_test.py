#########################
# py.test test functions
#########################


import re
import pytest
from Qpyl.core.qparameter import QPrm, QPrmError


class TestQ:
# TODO: tests for _PrmAtom/Bond/Angle/Torsion/Improper's functions
    def test_read_write_prm(self):
        qp_str = open("data/qamber14.prm", "r").read()
        qp_str = re.sub("(\*|\!|#).*", "", qp_str)
        qp_str = re.sub("\s+$", "", qp_str, 0, re.MULTILINE)
        qp_str = re.sub("^\n", "", qp_str, 0, re.MULTILINE)

        qprm = QPrm("amber")
        qprm.read_prm("data/qamber14.prm")
        qp_str2 = qprm.get_string()
        qp_str2 = re.sub("\s+$", "", qp_str2, 0, re.MULTILINE)
        qp_str2 = re.sub("^\n", "", qp_str2, 0, re.MULTILINE)

        assert qp_str == qp_str2

    def test_wrong_ff_fail(self):
        qprm = QPrm("amber")
        with pytest.raises(QPrmError):
            qprm.read_ffld("data/ace_ash_nma.ffld11")
        qprm = QPrm("oplsaa")
        with pytest.raises(QPrmError):
            qprm.read_amber_parm("data/ff-amber14/parm/parm10.dat")
        with pytest.raises(QPrmError):
            qprm.read_amber_frcmod("data/ff-amber14/parm/frcmod.ff14SB")

    def test_types(self):
        qprm = QPrm("amber")
        qprm.read_prm("data/qamber14.prm")
        assert qprm.atom_types["CK"].lj_R == 1.908
        assert qprm.atom_types["CK"].lj_eps == 0.086
        assert qprm.atom_types["CK"].mass == 12.01

        assert qprm.bonds["Br CA"].fc == 344.0
        assert qprm.bonds["Br CA"].r0 == 1.89

        assert qprm.angles["C4 C4 Nstar"].fc == 140.0
        assert qprm.angles["C4 C4 Nstar"].theta0 == 121.2

        assert qprm.torsions["Cstar CT CX N3"].fcs == [0.031, 0.234,
                                                       0.313, 0.079]


class TestAmber:
    def test_amber_conversion(self):
        qp_str = open("data/qamber14.prm", "r").read()
        qp_str = re.sub("(\*|\!|#).*", "", qp_str)
        qp_str = re.sub("\s+$", "", qp_str, 0, re.MULTILINE)
        qp_str = re.sub("^\n", "", qp_str, 0, re.MULTILINE)

        qprm = QPrm("amber", ignore_errors=True) # duplicates
        qprm.read_amber_parm("data/ff-amber14/parm/parm10.dat")
        ow = qprm.read_amber_frcmod("data/ff-amber14/parm/frcmod.ff14SB")

        # add options section manually and compare with official lib
        for line in """\
name                           Q-Amber14SB
type                           AMBER
vdw_rule                       arithmetic
scale_14                       0.8333
switch_atoms                   off
improper_potential             periodic
improper_definition            explicit\
""".splitlines():
            lf = line.split()
            qprm.options[lf[0]] = " ".join(lf[1:])

        qp_str2 = qprm.get_string()
        qp_str2 = re.sub("\s+$", "", qp_str2, 0, re.MULTILINE)
        qp_str2 = re.sub("^\n", "", qp_str2, 0, re.MULTILINE)

        assert qp_str == qp_str2

    def test_read_amber_parm(self):
        qprm = QPrm("amber")
        qprm.read_amber_parm("data/gaff.dat")
        assert len(qprm.atom_types) == 71
        assert len(qprm.bonds) == 832
        assert len(qprm.angles) == 4618
        assert len(qprm.generic_torsions) == 587
        assert len(qprm.torsions) == 66
        assert len(qprm.generic_impropers) == 8
        assert len(qprm.impropers) == 27

    def test_read_amber_parm_fail(self):
        qprm = QPrm("amber")
        with pytest.raises(QPrmError):
            qprm.read_amber_parm("data/ff-amber14/parm/frcmod.ff14SB")

    def test_read_amber_frcmod(self):
        qprm = QPrm("amber", ignore_errors=True) # duplicates
        qprm.read_amber_frcmod("data/ff-amber14/parm/frcmod.ff14SB")
        assert len(qprm.atom_types) == 4
        assert len(qprm.bonds) == 27
        assert len(qprm.angles) == 92
        assert len(qprm.generic_torsions) == 10
        assert len(qprm.torsions) == 128
        assert len(qprm.generic_impropers) == 1
        assert len(qprm.impropers) == 2

    def test_overwritten_prm(self):
        qprm = QPrm("amber", ignore_errors=True) # duplicates
        qprm.read_amber_parm("data/ff-amber14/parm/parm10.dat")
        ow = qprm.read_amber_frcmod("data/ff-amber14/parm/frcmod.ff14SB")
        # check overwritten parm
        assert ow[0].prm_id == "C N CX CT"
        assert ow[0].fcs == [0.0, 0.4, 2.0, 2.0]

    def test_overwritten_prm_fail(self):
        qprm = QPrm("amber") # no ignore_errors
        qprm.read_amber_parm("data/ff-amber14/parm/parm10.dat")
        with pytest.raises(QPrmError):
            qprm.read_amber_frcmod("data/ff-amber14/parm/frcmod.ff14SB")



class TestOplsaa:
    def test_read_ffld(self):
        qprm = QPrm("oplsaa")
        qprm.read_ffld("data/ace_ash_nma.ffld11")
        assert len(qprm.atom_types) == 11
        assert len(qprm.bonds) == 14
        assert len(qprm.angles) == 26
        assert len(qprm.torsions) == 35
        assert len(qprm.impropers) == 5

    def test_types2(self):
        qprm = QPrm("oplsaa")
        qprm.read_ffld("data/ace_ash_nma.ffld11")
        print qprm.torsions.keys()
        lj_A = 4*0.17*((3.25)**12)
        lj_B = 4*0.17*((3.25)**6)
        at = qprm.atom_types["N_N1_238"]
        assert abs(at.lj_A - lj_A) < 1e-7
        assert abs(at.lj_B - lj_B) < 1e-7

        bond = qprm.bonds["CT_C1_135 C_C2_235"]
        assert abs(bond.fc/2.0 - 317.0) < 1e-7
        assert abs(bond.r0 - 1.522) < 1e-7

        ang = qprm.angles["CT_C1_135 CT1_C1_224 C_C2_235"]
        assert abs(ang.fc/2.0 - 63.0) < 1e-7
        assert abs(ang.theta0 - 111.1) < 1e-7

        ang = qprm.angles["CT_C1_135 CT1_C1_224 C_C2_235"]
        assert abs(ang.fc/2.0 - 63.0) < 1e-7
        assert abs(ang.theta0 - 111.1) < 1e-7

        tors = qprm.torsions["C_C2_235 CT1_C1_224 N_N1_238 C_C2_235"]
        assert abs(tors.fcs[0]*2.0 - -2.365) < 1e-7
        assert abs(tors.fcs[1]*2.0 - 0.912) < 1e-7
        assert abs(tors.fcs[2]*2.0 - -0.850) < 1e-7
        assert abs(tors.periodicities[0] - 1.0) < 1e-7
        assert abs(tors.periodicities[1] - 2.0) < 1e-7
        assert abs(tors.periodicities[2] - 3.0) < 1e-7
        assert abs(tors.phases[0] - 0.0) < 1e-7
        assert abs(tors.phases[1] - 180.0) < 1e-7
        assert abs(tors.phases[2] - 0.0) < 1e-7

