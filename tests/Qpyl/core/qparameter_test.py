#########################
# py.test test functions
#########################


import re
import pytest
from Qpyl.core.qparameter import QPrm, QPrmError
from Qpyl.core.qstructure import QStruct

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


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

    def test_read_write_prm2(self):
        qp_str = open("data/ace_ash_nma.prm", "r").read()
        qp_str = re.sub("(\*|\!|#).*", "", qp_str)
        qp_str = re.sub("\s+$", "", qp_str, 0, re.MULTILINE)
        qp_str = re.sub("^\n", "", qp_str, 0, re.MULTILINE)

        qprm = QPrm("oplsaa")
        qprm.read_prm("data/ace_ash_nma.prm")
        qp_str2 = qprm.get_string()
        qp_str2 = re.sub("(\*|\!|#).*", "", qp_str2)
        qp_str2 = re.sub("\s+$", "", qp_str2, 0, re.MULTILINE)
        qp_str2 = re.sub("^\n", "", qp_str2, 0, re.MULTILINE)

        assert qp_str == qp_str2


    def test_wrong_ff_fail(self):
        qprm = QPrm("amber")
        with pytest.raises(QPrmError):
            qprm.read_ffld("data/ace_ash_nma.ffld11", None)
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

    def test_types_opls(self):
        qprm = QPrm("oplsaa")
        qprm.read_prm("data/ace_ash_nma.prm")
        assert qprm.atom_types["ace.C"].lj_A == 1802.2385
        assert qprm.atom_types["ace.C"].lj_B == 34.1758
        assert qprm.atom_types["ace.C"].mass == 12.011

        assert qprm.bonds["ace.CH3 ace.HH31"].fc == 680.0
        assert qprm.bonds["ace.CH3 ace.HH31"].r0 == 1.09

        assert qprm.angles["ash.OD1 ash.CG ash.OD2"].fc == 160.0
        assert qprm.angles["ash.OD1 ash.CG ash.OD2"].theta0 == 121.0

        t = qprm.torsions["ash.CA ash.CB ash.CG ash.OD2"]
        assert t.fcs == [0.225, 0.273, 0.5]
        assert t.multiplicities == [3.0, 2.0, 1.0]

        i = qprm.impropers["ash.CA ash.C ash.O nma.N"]
        assert i.fc == 10.5
        assert i.phi0 == 180.0



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
        qstruct = QStruct("data/ace_ash_nma.pdb", "pdb")
        qprm = QPrm("oplsaa")
        qprm.read_ffld("data/ace_ash_nma.ffld11", qstruct)
        assert len(qprm.atom_types) == 25
        assert len(qprm.bonds) == 24
        assert len(qprm.angles) == 40
        assert len(qprm.torsions) == 49
        assert len(qprm.impropers) == 5

    def test_types_ffld(self):
        qstruct = QStruct("data/ace_ash_nma.pdb", "pdb")
        qprm = QPrm("oplsaa")
        qprm.read_ffld("data/ace_ash_nma.ffld11", qstruct)
        print qprm.torsions.keys()
        lj_A_i= ( 4*0.17*((3.25)**12) )**0.5
        lj_B_i = ( 4*0.17*((3.25)**6) )**0.5
        at = qprm.atom_types["nma.N"]
        assert is_close(at.lj_A, lj_A_i)
        assert is_close(at.lj_B, lj_B_i)

        bond = qprm.bonds["ace.C ace.CH3"]
        assert is_close(bond.fc/2.0, 317.0)
        assert is_close(bond.r0, 1.522)

        ang = qprm.angles["ash.C ash.CA ash.CB"]
        assert is_close(ang.fc/2.0, 63.0)
        assert is_close(ang.theta0, 111.1)

        tors = qprm.torsions["ace.C ash.N ash.CA ash.C"]
        assert is_close(tors.fcs[0]*2.0, -2.365)
        assert is_close(tors.fcs[1]*2.0, 0.912)
        assert is_close(tors.fcs[2]*2.0, -0.850)
        assert is_close(tors.multiplicities[0], 1.0)
        assert is_close(tors.multiplicities[1], 2.0)
        assert is_close(tors.multiplicities[2], 3.0)
        assert is_close(tors.phases[0], 0.0)
        assert is_close(tors.phases[1], 180.0)
        assert is_close(tors.phases[2], 0.0)

