#########################
# py.test test functions
#########################


from __future__ import absolute_import
import pytest

from Qpyl.core.qlibrary import QLib, QLibError
from Qpyl.core.qstructure import QStruct

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

class TestQ:
    def test_read_write_lib(self):
        qlib = QLib("amber")
        qlib.read_lib("data/qamber14.lib")
        qlib.residue_dict.pop("HOH")  # different order
        ql_str2 = qlib.get_string()
        assert ql_str2 in open("data/qamber14.lib", "r").read()

    def test_wrong_ff_fail(self):
        qlib = QLib("amber", ignore_errors=True)
        with pytest.raises(QLibError):
            qlib.read_ffld("data/ace_ash_nma.ffld11", None)

        qlib = QLib("oplsaa", ignore_errors=True)
        with pytest.raises(QLibError):
            qlib.read_amber_lib("data/ff-amber14/lib/amino12.lib")
        with pytest.raises(QLibError):
            qlib.read_prepin_impropers("data/ff-amber14/lib/amino12.lib")
        with pytest.raises(QLibError):
            qlib.read_mol2("data/all_amino_acids.mol2")

    def test_atom(self):
        qlib = QLib("amber")
        qlib.read_lib("data/qamber14.lib")
        trp_CG = qlib.residue_dict["TRP"].atoms[7]
        assert trp_CG.name == "CG"
        assert trp_CG.atom_type == "Cstar"
        assert is_close(trp_CG.charge, -0.1415)
        assert trp_CG.residue.name == "TRP"

    def test_residue(self):
        qlib = QLib("amber")
        qlib.read_lib("data/qamber14.lib")
        asp = qlib.residue_dict["ASP"]
        assert asp.bonds[5] == ("CB", "HB2")
        assert asp.impropers[0] == ["-C", "N", "CA", "H"]
        assert asp.connections == ["head N", "tail C"]
        hoh = qlib.residue_dict["HOH"]
        assert int(hoh.info["solvent"]) == 1
        assert is_close(float(hoh.info["density"]), 0.0335)


    def test_rounding(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        old_charge = prc.atoms[0].charge
        prc.atoms[0].charge -= 1e-7
        prc.rescale(prc.charge_groups[0], 1)
        assert is_close(prc.atoms[0].charge, old_charge)

    def test_rescale(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc.atoms[0].charge -= 0.5
        prc.rescale(prc.charge_groups[0], 1)
        assert is_close(prc.atoms[0].charge, -0.72561699999)

    def test_rescale_badatom_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        with pytest.raises(QLibError):
            prc.rescale(prc.charge_groups[0] + ["BadAtom"], 0.5)

    def test_rescale_over_threshold_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc.atoms[0].charge -= 0.6
        with pytest.raises(QLibError):
            prc.rescale(prc.charge_groups[0], 0.3)

    def test_dupatom_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc.atoms.append(prc.atoms[0])
        with pytest.raises(QLibError):
            prc.check_valid()

    def test_bad_atom_bonds_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc.bonds.append(["C1", "BadAtom"])
        with pytest.raises(QLibError):
            prc.check_valid()

    def test_bad_atom_imps_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc._add_improper("C1", ["C4", "C7", "BadAtom"])
        with pytest.raises(QLibError):
            prc.check_valid()

    def test_bad_atom_cg_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc.charge_groups[1].append("BadAtom")
        with pytest.raises(QLibError):
            prc.check_valid()

    def test_net_charge_cg_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc.atoms[0].charge -= 1e-7
        with pytest.raises(QLibError):
            prc.check_valid()

    def test_net_charge_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc.charge_groups = []
        prc.atoms[0].charge -= 1e-7
        with pytest.raises(QLibError):
            prc.check_valid()

    def test_build_rules_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc.build_rules = ["torsion C1 C4 C7 H9 0"] # passes
        prc.check_valid()

        prc.build_rules = ["torsion C1 C4 C7 H9"]
        with pytest.raises(QLibError):
            prc.check_valid()
        prc.build_rules = ["bad_build_rule C1 C4 C7 H9 0"]
        with pytest.raises(QLibError):
            prc.check_valid()
        prc.build_rules = ["torsion C1 BadAtom C7 H9 0"]
        with pytest.raises(QLibError):
            prc.check_valid()

    def test_info_fail(self):
        qlib = QLib("oplsaa")
        qlib.read_lib("data/prc.lib")
        prc = qlib.residue_dict["PRC"]
        prc.info["random_number_generation"] = "true"
        with pytest.raises(QLibError):
            prc.check_valid()



class TestAmber:
    def read_amber_conversion(self):
        qlib = QLib("amber")
        qlib.read_amber_lib("data/ff-amber14/amber12_mod.lib")
        qlib.read_amber_lib("data/ff-amber14/arn.lib")
        qlib.read_prepin_impropers("data/ff-amber14/prep/amino12.in")
        qlib.read_prepin_impropers("data/ff-amber14/arn.prepi")

        # remove head from ACE and tail from NME
        cons = qlib.residue_dict["ACE"].connections
        cons = [con for con in cons if "head" not in con]
        qlib.residue_dict["ACE"].connections = cons

        cons = qlib.residue_dict["NME"].connections
        cons = [con for con in cons if "tail" not in con]
        qlib.residue_dict["NME"].connections = cons

        assert qlib.get_string() in open("data/qamber14.lib").read()

    def test_read_amber_lib_fail(self):
        # no residues found
        qlib = QLib("amber")
        with pytest.raises(QLibError):
            qlib.read_amber_lib("data/ff-amber14/parm/parm10.dat")

    def test_read_mol2_fail(self):
        # no residues found
        qlib = QLib("amber")
        with pytest.raises(QLibError):
            qlib.read_amber_lib("data/ff-amber14/parm/parm10.dat")

    def test_read_prepin_impropers_fail(self):
        # no residues in library
        qlib = QLib("amber")
        with pytest.raises(QLibError):
            qlib.read_amber_lib("data/ff-amber14/prep/amino12.in")


class TestOplsaa:
    def test_read_ffld(self):
        qlib = QLib("oplsaa")
        qstruct = QStruct("data/ace_ash_nma.pdb", "pdb")
        qlib.read_ffld("data/ace_ash_nma.ffld11", qstruct)
        assert len(qlib.residue_dict) == 3
        assert len(qlib.residue_dict["ACE"].atoms) == 6
        assert len(qlib.residue_dict["ASH"].atoms) == 13
        assert len(qlib.residue_dict["NMA"].atoms) == 6
        ash = qlib.residue_dict["ASH"]
        assert ash.atoms[1].atom_type == "ash.CA"
        assert is_close(ash.atoms[1].charge, 0.14)
        assert "tail C" in ash.connections
        assert "head N" in ash.connections

    def test_convert_oplsaa(self):
        qlib = QLib("oplsaa")
        qstruct = QStruct("data/ace_ash_nma.pdb", "pdb")
        qlib.read_ffld("data/ace_ash_nma.ffld11", qstruct)
        ql_str = qlib.get_string()
        assert ql_str == open("data/ace_ash_nma.lib").read()

    def test_read_ffld_fail(self):
        # no residues found
        qlib = QLib("oplsaa")
        qstruct = QStruct("data/ace_ash_nma.pdb", "pdb")
        with pytest.raises(QLibError):
            qlib.read_ffld("data/ace_ash_nma.pdb", qstruct)

    def test_read_ffld_wrong_order_fail(self):
        # see if it fails with PDB with wrong atom-order
        qlib = QLib("oplsaa")
        qstruct = QStruct("data/ace_ash_nma_bad.pdb", "pdb")
        with pytest.raises(QLibError):
            qlib.read_ffld("data/ace_ash_nma.ffld11", qstruct)


