#########################
# py.test test functions
#########################


import pytest

from Qpyl.core.qpotential import *
from Qpyl.core.qstructure import PosVector as P

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestPotFunc:
    def test_bond_energy(self):
        e = bond_energy(3, 1000, 1)
        assert is_close(e, 2000)

    def test_angle_energy(self):
        e = angle_energy(100, 100, 120)
        assert is_close(e, 6.0923484)

    def test_torsion_energy(self):
        e = torsion_energy(0, 10, 3, 3, 0)
        assert is_close(e, 6.6666666666)

    def test_improper_energy(self):
        e = improper_energy_periodic(150, 10.5, 2, 180)
        assert is_close(e, 5.25)

class TestGeomFunc:
    def test_distance(self):
        a1, a2 = P(0, 0, 0), P(5, 6, 7)
        d = bond_distance(a1, a2)
        assert is_close(d, 10.488088482)

    def test_angle(self):
        a1, a2, a3 = P(1, 0, 0), P(0, 0, 0), P(0, 1, 0)
        d = angle_angle(a1, a2, a3)
        assert is_close(d, 90)
        d = angle_angle(a2, a1, a3)
        assert is_close(d, 45)

    def test_dihedral(self):
        a1, a2, a3, a4 = P(0, 0, 0), P(1, 0, 0), P(1, 1, 0), P(2, 1, 0)
        d = torsion_angle(a1, a2, a3, a4)
        assert is_close(d, 180)

        a1, a2, a3, a4 = P(0, 0, 0), P(1, 0, 0), P(1, 1, 0), P(2, 1, 1)
        d = torsion_angle(a1, a2, a3, a4)
        assert is_close(d, 135)





