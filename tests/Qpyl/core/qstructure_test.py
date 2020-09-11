#########################
# py.test test functions
#########################


from __future__ import absolute_import
import pytest
from Qpyl.core.qstructure import QStruct, QStructError

def test_notsupported():
    with pytest.raises(QStructError):
        QStruct("data/all_amino_acids.pdb", "txt")


class TestPDB: 
    # cannot possibly account for all types of garbage inputs...

    def test_read_pdb(self):
        # check the basics (number of atoms&residues, names...)
        qstruct = QStruct("data/all_amino_acids.pdb", "pdb")
        assert "HH31 CH3 HH32 HH33 C" == \
               " ".join([a.name for a in qstruct.atoms[0:5]])
        assert "ACE ALA ARG" == \
               " ".join([r.name for r in qstruct.residues[0:3]])
        assert len(qstruct.atoms) == 456
        assert len(qstruct.residues) == 30
        assert len(qstruct.molecules) == 1

    def test_read_pdb_fail(self):
        # bad atom indexes will fail parsing
        with pytest.raises(QStructError):
            QStruct("data/all_amino_acids_bad.pdb", "pdb")

    def test_read_pdb_fail2(self):
        # bad residue indexes will fail parsing
        with pytest.raises(QStructError):
            QStruct("data/all_amino_acids_bad2.pdb", "pdb")

    def test_read_pdb_fail3(self):
        # random input
        with pytest.raises(QStructError):
            QStruct("data/all_amino_acids_bad.mol2", "pdb")

    def test_read_pdb_unsafe(self):
        # ignore_errors will skip bad indexes but can produce crap
        qstruct_bad = QStruct("data/all_amino_acids_bad.pdb", "pdb",
                              ignore_errors=True)
        qstruct = QStruct("data/all_amino_acids.pdb", "pdb")
        assert len(qstruct.atoms) == len(qstruct_bad.atoms)
        assert len(qstruct.residues) == len(qstruct_bad.residues) - 2
        assert len(qstruct.molecules) == len(qstruct_bad.molecules) - 3

    def test_convert_placeholders(self):
        inp_string = """
$$ $AAA$ $.AA$ $AA.$ $A.A $ $ A.A$
$2.CB$ $3.CA$ $LAST.ID$
"""
        out_string = """
$$ $AAA$ $.AA$ $AA.$ $A.A $ $ A.A$
11     21     25
""".strip()
        qstruct = QStruct("data/ace_ash_nma.pdb", "pdb")
        assert qstruct.convert_placeholders(inp_string).strip() == out_string

    def test_convert_placeholders_fail(self):
        qstruct = QStruct("data/ace_ash_nma.pdb", "pdb")
        with pytest.raises(QStructError):
            qstruct.convert_placeholders("$1.CB$")




class TestMOL2:

    def test_read_mol2(self):
        # check the basics (number of atoms&residues, names...)
        qstruct = QStruct("data/all_amino_acids.mol2", "mol2")
        assert "HH31 CH3 HH32 HH33 C" == \
               " ".join([a.name for a in qstruct.atoms[0:5]])
        assert "ACE ALA ARG" == \
               " ".join([r.name for r in qstruct.residues[0:3]])
        assert len(qstruct.atoms) == 456
        assert len(qstruct.residues) == 30
        assert len(qstruct.molecules) == 1

    def test_read_mol2_fail(self):
        # bad atom indexes will fail parsing
        with pytest.raises(QStructError):
            QStruct("data/all_amino_acids_bad.mol2", "mol2")

    def test_read_mol2_fail2(self):
        # bad residue indexes will fail parsing
        with pytest.raises(QStructError):
            QStruct("data/all_amino_acids_bad2.mol2", "mol2")

    def test_read_mol2_unsafe(self):
        # ignore_errors will skip bad indexes but can produce crap
        qstruct_bad = QStruct("data/all_amino_acids_bad.mol2", "mol2",
                              ignore_errors=True)
        qstruct = QStruct("data/all_amino_acids.mol2", "mol2")
        assert len(qstruct.atoms) == len(qstruct_bad.atoms)
        assert len(qstruct.residues) == len(qstruct_bad.residues) - 2
        assert len(qstruct.molecules) == len(qstruct_bad.molecules)

    def test_convert_placeholders(self):
        inp_string = """
$$ $AAA$ $.AA$ $AA.$ $A.A $ $ A.A$
$2.CB$ $3.CA$ $LAST.ID$
"""
        out_string = """
$$ $AAA$ $.AA$ $AA.$ $A.A $ $ A.A$
11 19 456
""".split()
        qstruct = QStruct("data/all_amino_acids.mol2", "mol2")
        assert qstruct.convert_placeholders(inp_string).split() == out_string

        # FAIL
        with pytest.raises(QStructError):
            qstruct.convert_placeholders("$1.CB$")



