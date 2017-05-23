#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# MIT License
#
# Copyright (c) 2017  Miha Purg <miha.purg@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
#

import math

from Qpyl.core import qstructure
from Qpyl.common import __version__, raise_or_log


def bond_distance(ac1, ac2):
    """
    Calculate the bond distance.

    Args:
      ac1, ac2 (qstructure.PosVector):   atom coordinates

    Returns:
      r (float):   distance in Angstrom
    """

    r21 = qstructure.PosVector(x=ac1.x - ac2.x,
                               y=ac1.y - ac2.y,
                               z=ac1.z - ac2.z)

    return math.sqrt(sum([a**2 for a in r21]))


def angle_angle(ac1, ac2, ac3):
    """
    Calculate the angle (of the angle).

    Args:
      ac1, ac2, ac3 (qstructure.PosVector):   atom coordinates

    Returns:
      theta (float):   angle in Degrees
    """

    r21 = qstructure.PosVector(x=ac1.x - ac2.x,
                               y=ac1.y - ac2.y,
                               z=ac1.z - ac2.z)

    r23 = qstructure.PosVector(x=ac3.x - ac2.x,
                               y=ac3.y - ac2.y,
                               z=ac3.z - ac2.z)

    # get the angle from the dot product equation
    # ( A*B = |A|*|B|*cos(theta) )
    # where A and B are vectors bewteen atoms (2->1 and 2->3)
    dr21 = math.sqrt(sum([a**2 for a in r21]))
    dr23 = math.sqrt(sum([a**2 for a in r23]))
    dot_product = sum([a*b for (a, b) in zip(r21, r23)])
    cos_theta = 1.0*dot_product/(dr21*dr23)
    theta = math.acos(cos_theta)   # in radians
    return 180.0 * theta/math.pi


def torsion_angle(ac1, ac2, ac3, ac4):
    """
    Calculate the dihedral angle.

    Args:
      ac1, ac2, ac3, ac4 (qstructure.PosVector):  atom coordinates

    Returns:
      theta (float):   angle in Degrees

    """

    r21 = qstructure.PosVector(x=ac1.x - ac2.x,
                               y=ac1.y - ac2.y,
                               z=ac1.z - ac2.z)

    r23 = qstructure.PosVector(x=ac3.x - ac2.x,
                               y=ac3.y - ac2.y,
                               z=ac3.z - ac2.z)

    r34 = qstructure.PosVector(x=ac4.x - ac3.x,
                               y=ac4.y - ac3.y,
                               z=ac4.z - ac3.z)

    # vector product n2 = r21 x r23
    n2 = qstructure.PosVector(x=r21.y*r23.z - r21.z*r23.y,
                              y=r21.z*r23.x - r21.x*r23.z,
                              z=r21.x*r23.y - r21.y*r23.x)
    # vector product n3 = r34 x r23
    n3 = qstructure.PosVector(x=r34.y*r23.z - r34.z*r23.y,
                              y=r34.z*r23.x - r34.x*r23.z,
                              z=r34.x*r23.y - r34.y*r23.x)

    # get the angle from the dot product equation ( A*B = |A|*|B|*cos(phi) )
    # where A and B are normal vectors n2 and n3
    dn2 = math.sqrt(sum([a**2 for a in n2]))
    dn3 = math.sqrt(sum([a**2 for a in n3]))
    dot_product = sum([a*b for (a, b) in zip(n2, n3)])

    cos_phi = 1.0*dot_product/(dn2*dn3)
    if cos_phi < -1: cos_phi = -1
    elif cos_phi > 1: cos_phi = 1
    phi = math.acos(cos_phi) # in radians
    return 180.0 * phi/math.pi


def improper_angle(ac1, ac2, ac3, ac4):
    """
    Calculate the dihedral angle.

    Args:
      ac1, ac2, ac3, ac4 (qstructure.PosVector):  atom coordinates

    Returns:
      theta (float):   angle in Degrees
    """
    return torsion_angle(ac1, ac2, ac3, ac4)



def bond_energy(r, fc, r0):
    """
    Calculate the bond energy using the harmonic potential.

    Args:
      r (float):  distance between atoms [angstrom]
      fc (float):  force constant [kcal/mol]
      r0 (float):  equilibrium distance [angstrom]

    Returns:
      e_bond (float):  energy of bond [kcal/mol]
    """

    return 0.5 * fc * (r - r0)**2


def angle_energy(theta, fc, theta0):
    """
    Calculate the angle energy using the harmonic potential.

    Args:
      theta (float):  angle between atoms [degrees]
      fc (float):  force constant [kcal/mol]
      theta0 (float):  equilibrium angle [degrees]

    Returns:
      e_angle (float):  energy of angle [kcal/mol]
    """

    theta = math.pi/180.0 * theta   # degrees to radians
    theta0 = math.pi/180.0 * theta0   # degrees to radians
    return 0.5 * fc * (theta - theta0)**2


def torsion_energy(phi, fc, periodicity, npaths, phi0):
    """
    Calculate the torsion energy using the periodic potential.

    One torsion term at a time, using this function:
    E = fc/npaths * [1 + cos(periodicity*phi - phi0)]

    Args:
      phi (float):  dihedral angle [degrees]
      fc (float):  force constant [kcal/mol]
      periodicity (float):  periodicity [/]
      npaths (float):  division factor [/]
      phi0 (float):  phase shift [degrees]

    Returns:
      e_torsion (float):  energy of torsion [kcal/mol]
    """
    phi = math.pi/180.0 * phi   # degrees to radians
    phi0 = math.pi/180.0 * phi0   # degrees to radians
    return 1.0 * fc/npaths * (1 + math.cos(periodicity*phi - phi0))

def improper_energy_periodic(phi, fc, periodicity, phi0):
    """
    Calculate the improper energy using the periodic potential.

    Using the same potential as torsion_energy(), with npaths==1.
    E = fc/npaths * [1 + cos(periodicity*phi - phi0)]

    Args:
      phi (float):  dihedral angle [degrees]
      fc (float):  force constant [kcal/mol]
      periodicity (float):  periodicity [/]
      phi0 (float):  phase shift [degrees]

    Returns:
      e_improper (float):  energy of improper [kcal/mol]
    """
    phi = math.pi/180.0 * phi   # degrees to radians
    phi0 = math.pi/180.0 * phi0   # degrees to radians
    return fc * (1 + math.cos(periodicity*phi - phi0))



