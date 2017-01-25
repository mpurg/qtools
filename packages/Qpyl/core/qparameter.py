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

import re
import math
import logging
from collections import OrderedDict

from Qpyl.common import __version__, raise_or_log

logger = logging.getLogger(__name__)


ATOM_MASSES = { "H": 1.0079,
                "C": 12.011,
                "N": 14.007,
                "O": 15.999,
                "F": 18.988,
                "P": 30.974,
                "S": 32.065,
                "CL": 35.453,
                "BR": 79.904,
                "I": 126.90,
                "DU": 0 }


class QPrmError(Exception):
    pass

class QPrm(object):
    """
    Class for reading and writing Q parameter files.

    Also supports parsing oplsaa_ffld, amber_parm and amber_frcmod files.

    Args:
       ff_type (string):  Either 'oplsaa' or 'amber'
       ignore_errors (boolean):  Optional, default is False.
                                 If set to True, some non-vital
                                 exceptions are logged instead.
    """
    def __init__(self, ff_type, ignore_errors=False):
        self.ignore_errors = ignore_errors
        supported_ff = ['oplsaa', 'amber']
        ff_type = ff_type.lower()
        if ff_type not in supported_ff:
            raise QPrmError("Force field type '{}' not supported. "
                            "Use {}".format(ff_type,
                                            " or ".join(supported_ff)))
        self.ff_type = ff_type
        self.options = OrderedDict()
        self.amber_masses = OrderedDict()
        self.atom_types = OrderedDict()
        self.bonds = OrderedDict()
        self.angles = OrderedDict()
        self.generic_torsions = OrderedDict()   # with wildcards
        self.torsions = OrderedDict()
        self.generic_impropers = OrderedDict()  # with wildcards
        self.impropers = OrderedDict()



    def read_prm(self, parm_fn):
        """
        Read and parse a Q parameter file.

        Args:
            parm_fn (string):   name of parameter file

        Returns:
            overwritten_prms (list):  list of overwritten parameters

        Same parameters with equal values (duplicates) are overwritten,
        same parameters with different values raise QPrmError, unless
        QPrm.ignore_errors is set to True.

        """

        section = None
        prms = {"atoms": [],
                "bonds": [],
                "angles": [],
                "torsions": [],
                "impropers": []}

        with open(parm_fn, 'r') as prmfile:
            for lnumber, line in enumerate(prmfile.readlines()):
                lnumber += 1

                line, comment = map(str.strip, re.split("#|\*|\!", line+"!", 1))
                comment = " ".join(comment.strip("!").split())

                if line == "":
                    continue
                if line[0] == "[":
                    # it is apparently allowed to write text after
                    # the section identifier, so a simple strip isn't enough
                    section = line.split("]")[0].strip(" [").lower()
                    continue
                if not section:
                    raise QPrmError("Line {} in PARM file '{}' is not a "
                                    "comment and is not inside any section "
                                    "([atom_types], [bonds], ...):\n{}"
                                    .format(lnumber, parm_fn, line))


                if section == "options":
                    key, value = line.split()
                    self.options[key] = value


                elif section == "atom_types":
                    parms = line.split()
                    try:
                        atom_type = parms[0]
                        lj_Ar, lj_Beps = float(parms[1]), float(parms[3])
                        mass = float(parms[6])
                    except Exception as e:
                        raise QPrmError("Could not parse line {} in "
                                        "[atom_types] section of parm file "
                                        "'{}':\n{}"
                                        .format(lnumber, parm_fn, line))

                    if self.ff_type == "oplsaa":
                        atom_type = _PrmAtom(atom_type, mass,
                                             lj_A=lj_Ar, lj_B=lj_Beps,
                                             comment=comment)
                    elif self.ff_type == "amber":
                        atom_type = _PrmAtom(atom_type, mass,
                                             lj_R=lj_Ar, lj_eps=lj_Beps,
                                             comment=comment)

                    prms["atoms"].append(atom_type)


                elif section == "bonds":
                    parms = line.split()
                    try:
                        atom_types = parms[0:2]
                        fc, r0 = float(parms[2]), float(parms[3])
                    except Exception as e:
                        raise QPrmError("Could not parse line {} in [bonds] "
                                        "section of parm file '{}':\n{}"
                                        .format(lnumber, parm_fn, line))

                    bond = _PrmBond(atom_types, fc, r0, comment=comment)
                    prms["bonds"].append(bond)

                elif section == "angles":
                    parms = line.split()
                    try:
                        atom_types = parms[0:3]
                        fc, theta0 = float(parms[3]), float(parms[4])
                    except Exception as e:
                        raise QPrmError("Could not parse line {} in [angles] "
                                        "section of parm file '{}':\n{}"
                                        .format(lnumber, parm_fn, line))

                    angle = _PrmAngle(atom_types, fc, theta0,
                                      comment=comment)
                    prms["angles"].append(angle)

                elif section == "torsions":
                    parms = line.split()
                    try:
                        atom_types = parms[0:4]
                        fc, periodicity, phase, npaths = map(float, parms[4:8])
                    except Exception as e:
                        raise QPrmError("Could not parse line {} in [torsions]"
                                        " section of parm file '{}':\n{}"
                                        .format(lnumber, parm_fn, line))

                    torsion = _PrmTorsion(atom_types, comment=comment)

                    # torsion parameters belonging to the same torsion
                    # must be sequential in the parameter file,
                    # otherwise there is no way of knowing which is correct
                    if prms["torsions"] and torsion.prm_id == \
                            prms["torsions"][-1].prm_id:
                        torsion = prms["torsions"][-1]
                    else:
                        prms["torsions"].append(torsion)

                    try:
                        torsion.add_prm(fc,
                                        periodicity,
                                        phase,
                                        npaths)
                    # in case of two parms sharing same periodicity
                    except ValueError:
                        raise QPrmError("Duplicate parameter found: {}"
                                        .format(torsion))



                elif section == "impropers":
                    parms = line.split()
                    try:
                        atom_types = parms[0:4]
                        center_atom = atom_types.pop(1)
                        fc, phi0 = map(float, parms[4:6])
                    except Exception as e:
                        raise QPrmError("Could not parse line {} in "
                                        "[impropers] section of parm file '{}'"
                                        ":\n{}".format(lnumber, parm_fn, line))


                    improper = _PrmImproper(center_atom, atom_types, fc, phi0,
                                            comment=comment)
                    prms["impropers"].append(improper)

                else:
                    raise QPrmError("Unknown section found in the parm file "
                                    "{}: {}".format(parm_fn, section))

        dups = []
        for type_, params in prms.iteritems():
            for prm in params:
                dup = self._add_prm(prm)
                if dup != None:
                    dups.append(dup)
        return dups



    def read_amber_parm(self, parm_fn):
        """
        Read and parse an Amber parameter file (.parm).

        Args:
            parm_fn (string):   name of parameter file

        Returns:
            overwritten_prms (list):  list of overwritten parameters

        Same parameters with equal values (duplicates) are overwritten,
        same parameters with different values raise QPrmError, unless
        QPrm.ignore_errors is set to True.

        At the moment, all comments are ignored.
        """

        if self.ff_type != "amber":
            raise QPrmError("Function not supported with "
                            "force field '{}'".format(self.ff_type))

        prms = {"atoms": [],
                "bonds": [],
                "angles": [],
                "torsions": [],
                "impropers": []}

        with open(parm_fn) as parmf:
            # line with FF description
            parmf.readline()

            # get the masses
            while True:
                line = parmf.readline().split()
                if not line: break
                try:
                    atom_name, mass = line[0].replace("*", "star"), float(line[1])
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))
                try:
                    old_mass = self.amber_masses[atom_name]
                except KeyError:
                    self.amber_masses[atom_name] = mass
                else:
                    if abs(mass - old_mass) > 0.00001:
                        raise_or_log("Different masses for same type ({}): "
                                     "{}, {}"
                                     .format(atom_name, mass, old_mass),
                                     QPrmError, logger, self.ignore_errors)

            # useless line (hydrophilic atoms?)
            parmf.readline()

            # get the bond types
            while True:
                line = parmf.readline().replace("-", " ", 1).split()
                if not line: break
                try:
                    a_types = [x.strip().replace("*", "star") for x in line[:2]]
                    fc, r0 = float(line[2]), float(line[3])
                    fc *= 2   # In Q the fc is divided by 2
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                bond = _PrmBond(a_types, fc, r0)
                prms["bonds"].append(bond)

            # get the angle types
            while True:
                line = parmf.readline().replace("-", " ", 2).split()
                if not line: break
                try:
                    a_types = [x.strip().replace("*", "star") for x in line[:3]]
                    fc, t0 = float(line[3]), float(line[4])
                    fc *= 2   # In Q the fc is divided by 2
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                angle = _PrmAngle(a_types, fc, t0)
                prms["angles"].append(angle)

            # get the torsion (dihedral) types
            while True:
                line = parmf.readline().replace("-", " ", 3).split()
                if not line: break
                try:
                    a_types = [x.strip().replace("*", "star") for x in line[:4]]
                    a_types = map(lambda x: "?" if x == "X" else x, a_types)
                    npaths, fc, phase, periodicity = map(float, line[4:8])
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                torsion = _PrmTorsion(a_types)

                # torsion parameters belonging to the same torsion
                # must be sequential in the parameter file,
                # otherwise there is no way of knowing which is correct
                if prms["torsions"] and torsion.prm_id == \
                        prms["torsions"][-1].prm_id:
                    torsion = prms["torsions"][-1]
                else:
                    prms["torsions"].append(torsion)

                try:
                    torsion.add_prm(fc,
                                    periodicity,
                                    phase,
                                    npaths)
                # in case of two parms sharing same periodicity
                except ValueError:
                    raise QPrmError("Duplicate parameter found: {}"
                                    .format(torsion))


            # get the improper types
            while True:
                line = parmf.readline().replace("-", " ", 3).split()
                if not line: break
                try:
                    a_types = [x.strip().replace("*", "star") for x in line[:4]]
                    a_types = map(lambda x: "?" if x == "X" else x, a_types)
                    center_atom = a_types.pop(2)
                    fc, phi0, periodicity = map(float, line[4:7])
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                improper = _PrmImproper(center_atom, a_types, fc,
                                        phi0, periodicity=periodicity)
                prms["impropers"].append(improper)

            # two useless lines
            parmf.readline()
            parmf.readline()

            # list of same parameters
            same_types = {}
            while True:
                line = parmf.readline().replace("*", "star").split()
                if not line: break
                same_types[line[0]] = line

            # another useless line  (MOD4 RE)
            parmf.readline()

            # get the nonbonding parameters
            while True:
                line = parmf.readline().split()
                if not line: break
                try:
                    a_type = line[0].replace("*", "star")
                    same_a_types = same_types.get(a_type, [a_type,])
                    lj_R, lj_eps = map(float, line[1:3])
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                try:
                    mass = self.amber_masses[a_type]
                except KeyError:
                    raise QPrmError("No mass for atom type: {}".format(a_type))

                for a_type in same_a_types:
                    atom_type = _PrmAtom(a_type, mass, lj_R=lj_R, lj_eps=lj_eps)
                    prms["atoms"].append(atom_type)


        dups = []
        for type_, params in prms.iteritems():
            for prm in params:
                dup = self._add_prm(prm)
                if dup != None:
                    dups.append(dup)
        return dups


    def read_amber_frcmod(self, frcmodfile):
        """
        Read and parse an Amber parameter file (.frcmod).

        Args:
            parm_fn (string):   name of parameter file

        Returns:
            overwritten_prms (list):  list of overwritten parameters

        Same parameters with equal values (duplicates) are overwritten,
        same parameters with different values raise QPrmError, unless
        QPrm.ignore_errors is set to True.

        At the moment, all comments are ignored.
        """

        if self.ff_type != "amber":
            raise QPrmError("Function not supported with "
                            "force field '{}'".format(self.ff_type))

        prms = {"atoms": [],
                "bonds": [],
                "angles": [],
                "torsions": [],
                "impropers": []}

        with open(frcmodfile) as parmf:
            # line with FF description
            parmf.readline()

            # get the masses
            parmf.readline()    # MASS
            while True:
                line = parmf.readline().split()
                if not line: break
                try:
                    atom_name, mass = line[0].replace("*", "star"), float(line[1])
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))
                try:
                    old_mass = self.amber_masses[atom_name]
                except KeyError:
                    self.amber_masses[atom_name] = mass
                else:
                    if abs(mass - old_mass) > 0.00001:
                        raise QPrmError("Different masses for same type ({}): "
                                        "{}, {}".format(atom_name, mass,
                                                        old_mass))

            # get the bond types
            parmf.readline()    # BOND
            while True:
                line = parmf.readline().replace("-", " ", 1).split()
                if not line: break
                try:
                    a_types = [x.strip().replace("*", "star")
                               for x in line[:2]]
                    fc, r0 = float(line[2]), float(line[3])
                    fc *= 2   # In Q the fc is divided by 2
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                bond = _PrmBond(a_types, fc, r0)
                prms["bonds"].append(bond)

            # get the angle types
            parmf.readline()    # ANGL
            while True:
                line = parmf.readline().replace("-", " ", 2).split()
                if not line: break
                try:
                    a_types = [x.strip().replace("*", "star") for x in line[:3]]
                    fc, t0 = float(line[3]), float(line[4])
                    fc *= 2   # In Q the fc is divided by 2
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                angle = _PrmAngle(a_types, fc, t0)
                prms["angles"].append(angle)

            # get the torsion (dihedral) types
            parmf.readline()    # DIHE
            while True:
                line = parmf.readline().replace("-", " ", 3).split()
                if not line: break
                try:
                    a_types = [x.strip().replace("*", "star") for x in line[:4]]
                    a_types = map(lambda x: "?" if x == "X" else x, a_types)
                    npaths, fc, phase, periodicity = map(float, line[4:8])
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                torsion = _PrmTorsion(a_types)
                # torsion parameters belonging to the same torsion
                # must be sequential in the parameter file,
                # otherwise there is no way of knowing which is correct
                if prms["torsions"] and torsion.prm_id == \
                        prms["torsions"][-1].prm_id:
                    torsion = prms["torsions"][-1]
                else:
                    prms["torsions"].append(torsion)

                try:
                    torsion.add_prm(fc,
                                    periodicity,
                                    phase,
                                    npaths)
                # in case of two parms sharing same periodicity
                except ValueError:
                    raise QPrmError("Duplicate parameter found: {}"
                                    .format(torsion))


            # get the improper types
            parmf.readline()    # IMPR
            while True:
                line = parmf.readline().replace("-", " ", 3).split()
                if not line: break
                try:
                    a_types = [x.strip().replace("*", "star") for x in line[:4]]
                    a_types = map(lambda x: "?" if x == "X" else x, a_types)
                    center_atom = a_types.pop(2)
                    fc, phi0, periodicity = map( float, line[4:7] )
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                improper = _PrmImproper(center_atom, a_types, fc, phi0,
                                        periodicity=periodicity)
                prms["impropers"].append(improper)


            # get the improper types
            parmf.readline()    # NONB
            while True:
                line = parmf.readline().split()
                if not line: break
                try:
                    a_type = line[0].replace("*", "star")
                    lj_R, lj_eps = map(float, line[1:3])
                except Exception as e:
                    raise QPrmError("Could not parse line '{}'".format(line))

                try:
                    mass = self.amber_masses[a_type]
                except KeyError:
                    raise QPrmError("No mass for atom type: {}".format(a_type))

                atom_type = _PrmAtom(a_type, mass, lj_R=lj_R, lj_eps=lj_eps)
                prms["atoms"].append(atom_type)

        dups = []
        for type_, params in prms.iteritems():
            for prm in params:
                dup = self._add_prm(prm)
                if dup != None:
                    dups.append(dup)
        return dups


    def read_ffld(self, ffld_file):
        """
        Read and parse a Macromodel's FFLD file for oplsaa parameters.

        Args:
            ffld_file (string):  path/name of ffld file

        Returns:
            overwritten_prms (list):  list of overwritten parameters

        Same parameters with equal values (duplicates) are overwritten,
        same parameters with different values raise QPrmError, unless
        QPrm.ignore_errors is set to True.

        Atom types are created by combining 'symbol', 'vdw' and 'type'.
        See qlibrary.QLib.read_ffld for more details.
        """

        if self.ff_type != "oplsaa":
            raise QPrmError("Function not supported with "
                            "force field '{}'".format(self.ff_type))

        # key is ffld atom name (C1), value is atom type (CA_145)
        lookup_aname = {}
        section = None

        prms = {"atoms": [],
                "bonds": [],
                "angles": [],
                "torsions": [],
                "impropers": []}

        for line in open(ffld_file, 'r').readlines():
            line = line.strip()
            if (line == "") or ("------" in line):
                continue
            elif line.startswith("atom   type  vdw  symbol"):
                section = "ATOMS"
                continue
            elif line.startswith("Stretch            k"):
                section = "BONDS"
                continue
            elif line.startswith("Bending                      k"):
                section = "ANGLES"
                continue
            elif line.startswith("proper Torsion"):
                section = "TORSIONS"
                continue
            elif line.startswith("improper"):
                section = "IMPROPERS"
                continue

#
#  C1      135  C1   CT      -0.0175   3.5000   0.0660 high   C: alkanes
#
            if section == "ATOMS":
                lf = line.split()
                try:
                    name, type_, vdw, symbol = lf[0:4]
                    charge, sigma, epsilon = map(float, lf[4:7])
                    comment = "FFLD: " + " ".join(lf[8:])
                except Exception as e:
                    raise QPrmError("Could not parse line: '{}'"
                                    .format(line))

                lj_A = 4*epsilon*((sigma)**12)
                lj_B = 4*epsilon*((sigma)**6)
                element = "".join(c for c in vdw if c.isalpha())
                try:
                    mass = ATOM_MASSES[element.upper()]
                except KeyError:
                    logger.warning("Mass for element '{}' (atom '{}') "
                                   "not found, set it manually."
                                   .format(element, name))
                    mass = "<FIX>"


                atype = "{}_{}_{}".format(symbol, vdw, type_)
                lookup_aname[name] = atype
                atom_type = _PrmAtom(atype, mass,
                                     lj_A=lj_A, lj_B=lj_B,
                                     comment=comment)
                prms["atoms"].append(atom_type)

#
#  C1      H2      340.00000    1.09000   high      140  0   CT  -HC    ==> CT  -HC
#
            if section == "BONDS":
                lf = line.split()
                try:
                    anames = lf[0:2]
                    fc, r0 = map(float, lf[2:4])
                    comment = "FFLD: " + " ".join(lf[4:])
                except Exception as e:
                    raise QPrmError("Could not parse line: '{}'"
                                    .format(line))

                atypes = [lookup_aname[name] for name in anames]
                bond = _PrmBond(atypes, fc*2.0, r0, comment=comment)
                prms["bonds"].append(bond)

#
#  H2      C1      H3        33.00000  107.80000   high      ...
#
            if section == "ANGLES":
                lf = line.split()
                try:
                    anames = lf[0:3]
                    fc, theta0 = map(float, lf[3:5])
                    comment = "FFLD: " + " ".join(lf[5:])
                except Exception as e:
                    raise QPrmError("Could not parse line: '{}'"
                                    .format(line))

                atypes = [lookup_aname[name] for name in anames]
                angle = _PrmAngle(atypes, fc*2.0, theta0, comment=comment)
                prms["angles"].append(angle)

#
#  O2      P1      O3      C4        0.000   0.000   0.562   0.000    high ...
#
            if section == "TORSIONS":
                lf = line.split()
                try:
                    anames = lf[0:4]
                    fcs = map(float, lf[4:8])
                    comment = "FFLD: " + " ".join(lf[8:])
                except Exception as e:
                    raise QPrmError("Could not parse line: '{}'"
                                    .format(line))


                periodicity = (1, 2, 3, 4)
                phi0s = (0.0, 180.0, 0.0, 180.0)
                paths = (1.0, 1.0, 1.0, 1.0)

                atypes = [lookup_aname[name] for name in anames]
                torsion = _PrmTorsion(atypes, comment=comment)
                for fc, nmin, phi0, path in zip(fcs, periodicity, phi0s, paths):
                    if abs(fc) > 0.000001:
                        torsion.add_prm(fc/2.0, nmin, phi0, path)
                if not torsion.get_prms():
                    torsion.add_prm(0.0, 1.0, 0.0, 1.0)
                prms["torsions"].append(torsion)

#
#  C21     C22     C20     O19       2.200   high   ...
#
            if section == "IMPROPERS":
                lf = line.split()
                try:
                    anames = lf[0:4]
                    fc = float(lf[4])
                    comment = "FFLD: " + " ".join(lf[5:])
                except Exception as e:
                    raise QPrmError("Could not parse line: '{}'"
                                    .format(line))

                atypes = [lookup_aname[name] for name in anames]
                center_atom = atypes.pop(2)
                improper = _PrmImproper(center_atom, atypes, fc/2.0,
                                        180.0, periodicity=2, comment=comment)
                prms["impropers"].append(improper)

        dups = []
        for type_, params in prms.iteritems():
            for prm in params:
                dup = self._add_prm(prm)
                if dup != None:
                    dups.append(dup)
        return dups



    def _add_prm(self, prm):
        # add the parameter 'prm' to the approprate dictionary
        # Depending on the value of QPrm.ignore_errors:
        #   True - overwrite parameter types with different values
        #   False - raise QPrmError
        if isinstance(prm, _PrmAtom):
            prm_dict = self.atom_types
        elif isinstance(prm, _PrmBond):
            prm_dict = self.bonds
        elif isinstance(prm, _PrmAngle):
            prm_dict = self.angles
        elif isinstance(prm, _PrmTorsion):
            if prm.is_generic:
                prm_dict = self.generic_torsions
            else:
                prm_dict = self.torsions
        elif isinstance(prm, _PrmImproper):
            if prm.is_generic:
                prm_dict = self.generic_impropers
            else:
                prm_dict = self.impropers
        else:
            raise TypeError("'prm' is not of type "
                            "_PrmAtom/Bond/Angle/Torsion/Improper")

        old_prm = None
        if prm.prm_id in prm_dict:
            sa1 = str(prm)
            sa2 = str(prm_dict[prm.prm_id])
            if sa1 != sa2:
                raise_or_log("Same parameter types with different "
                             "values:\n{}\n{}".format(sa1, sa2),
                             QPrmError, logger, self.ignore_errors)
                old_prm = prm_dict[prm.prm_id]
        prm_dict[prm.prm_id] = prm
        return old_prm






    def get_string(self, atom_types=None, bonds=None,
                   angles=None, torsions=None, impropers=None):
        """Return the Q parameter file in string format.

        Optional arguments are lists of _PrmAtom, _PrmBond, 
        _PrmAngle, _PrmTorsion and _PrmImproper objects.
        If unset, all parameters will be printed.

        Example:
        get_string(atom_types=["CA", "NA", "CB"], bonds=["CA_CB"],
                   angles=[], torsions=[], impropers=[])
        """
        
        prm_l = {"options": [],
                 "atom_types": [],
                 "bonds": [],
                 "angles": [],
                 "torsions": [],
                 "impropers": []}

        gentorsions = []
        genimpropers = []

        if atom_types == None:
            atom_types = self.atom_types.values()
        if bonds == None:
            bonds = self.bonds.values()
        if angles == None:
            angles = self.angles.values()
        if torsions == None:
            gentorsions = self.generic_torsions.values()
            torsions = self.torsions.values()
        if impropers == None:
            genimpropers = self.generic_impropers.values()
            impropers = self.impropers.values()

        for k, v in self.options.iteritems():
            prm_l["options"].append("{:<30s} {:<30s}".format(k,v))

        for v in sorted(set(atom_types), key=lambda x: x.prm_id):
            if self.ff_type == "amber":
                prm_l["atom_types"].append("{:<12} {:>10} {:>10} "
                                           "{:>10} {:>10} {:>10} "
                                           "{:>10} {}"
                                           .format(v.atom_type, v.lj_R, 0.0,
                                                   v.lj_eps, v.lj_R,
                                                   v.lj_eps/2, v.mass,
                                                   v.comment))
            elif self.ff_type == "oplsaa":
                lj_A = round(v.lj_A**0.5, 4)
                lj_B = round(v.lj_B**0.5, 4)
                prm_l["atom_types"].append("{:<12} {:>10} {:>10} "
                                           "{:>10} {:>10} {:>10} "
                                           "{:>10} {}"
                                           .format(v.atom_type, lj_A, lj_A,
                                                   lj_B, round(lj_A/2**0.5, 4),
                                                   round(lj_B/2**0.5, 4),
                                                   v.mass, v.comment))

        for v in sorted(set(bonds), key=lambda x: x.prm_id):
            a1, a2 = v.atom_types
            prm_l["bonds"].append("{:<12} {:<12} {:>10} {:>10} {}"
                                  .format(a1, a2, v.fc, v.r0, v.comment))

        for v in sorted(set(angles), key=lambda x: x.prm_id):
            a1, a2, a3 = v.atom_types
            prm_l["angles"].append("{:<12} {:<12} {:<12} {:>10} {:>10} {}"
                                   .format(a1, a2, a3, v.fc, v.theta0,
                                           v.comment))

        for v in sorted(set(gentorsions), key=lambda x: x.prm_id):
            a1, a2, a3, a4 = v.atom_types
            for fc, periodicity, phase, npaths in v.get_prms():
                prm_l["torsions"].append("{:<12} {:<12} {:<12} {:<12} "
                                         "{:>10} {:>5} {:>10} {:>5} {}"
                                         .format(a1, a2, a3, a4, fc,
                                                 periodicity, phase, npaths,
                                                 v.comment))

        for v in sorted(set(torsions), key=lambda x: x.prm_id):
            a1, a2, a3, a4 = v.atom_types
            for fc, periodicity, phase, npaths in v.get_prms():
                prm_l["torsions"].append("{:<12} {:<12} {:<12} {:<12} "
                                         "{:>10} {:>5} {:>10} {:>5} {}"
                                         .format(a1, a2, a3, a4, fc,
                                                 periodicity, phase, npaths,
                                                 v.comment))

        # sorting function includes empty spaces before the prm_id so
        # parameters with more wildcards come first
        for v in sorted(set(genimpropers), key=lambda x: \
                x.prm_id.count("?") * " " + x.prm_id):
            a1, a2, a3, a4 = v.atom_types
            prm_l["impropers"].append("{:<12} {:<12} {:<12} {:<12} "
                                      "{:>10} {:>10} {}"
                                      .format(a1, a2, a3, a4, v.fc,
                                              v.phi0, v.comment))

        for v in sorted(set(impropers), key=lambda x: x.prm_id):
            a1, a2, a3, a4 = v.atom_types
            prm_l["impropers"].append("{:<12} {:<12} {:<12} {:<12} "
                                      "{:>10} {:>10} {}"
                                      .format(a1, a2, a3, a4, v.fc,
                                              v.phi0, v.comment))

        for k in prm_l.keys():
            prm_l[k] = "\n".join(prm_l[k])

        return """[options]
{options}

[atom_types]
{atom_types}

[bonds]
{bonds}

[angles]
{angles}

[torsions]
{torsions}

[impropers]
{impropers}
""".format(**prm_l)


class _PrmAtom(object):
    def __init__(self, atom_type, mass,
                 lj_A=None, lj_B=None,
                 lj_R=None, lj_eps=None,
                 comment=None):
        """
        Stores atom parameters like type, mass, LJ params.

        Args:
            atom_type (string):  atom type (CT, CA, Cstar)
            mass (float):  mass of atom
            lj_A (float):  A_ii parameter for 12-6 pot. function
                           (geometric rules)  [ kcal/A^12 ]
            lj_B (float):  B_ii parameter
                           (geometric rules)   [ kcal/A^6 ]
            lj_R (float):  Rm_ii parameter
                           (arithmetic rules)  [ A ]
            lj_eps (float):  epsilon_ii parameter
                             (arithmetic rules)  [ kcal ]
            comment (string):  comment

        Note: Both A and B or both R and eps have to be given.
        """
        # atom_type is a string like "CA" or "OW"
        self.atom_type = atom_type
        self.prm_id = self.get_id(atom_type)
        self.lj_A = lj_A
        self.lj_B = lj_B
        self.lj_R = lj_R
        self.lj_eps = lj_eps
        self.mass = mass
        self._comment = comment

    @property
    def comment(self):
        if self._comment:
            return " # {}".format(self._comment) 
        else:
            return ""

    def __repr__(self):
        return "_PrmAtom({}, {})".format(self.prm_id, self.strval)

    @property
    def strval(self):
        """Return parameter values in string format."""
        if self.lj_A != None:
            prms = "lj_A={:.3f}, lj_B={:.3f}".format(self.lj_A, self.lj_B)
        else:
            prms = "lj_R={:.3f}, lj_eps={:.3f}".format(self.lj_R, self.lj_eps)
        return "{}, mass={:.3f}".format(prms, self.mass)


    @staticmethod
    def get_id(atom_type):
        """Return the unique identifier (atom type)"""
        return atom_type


class _PrmBond(object):
    def __init__(self, atom_types, fc, r0, comment=None):
        self.fc = fc
        self.r0 = r0
        self.prm_id = self.get_id(atom_types)
        self.atom_types = self.prm_id.split()   # list of atom_type strings
        self._comment = comment

    @property
    def comment(self):
        if self._comment:
            return " # {}".format(self._comment) 
        else:
            return ""

    def __repr__(self):
        return "_PrmBond({}, fc={:.3f}, r0={:.3f})".format(self.prm_id,
                                                           self.fc,
                                                           self.r0)
    @property
    def strval(self):
        """Return parameter values in string format."""
        return "fc={:.3f}, r0={:.3f}".format(self.fc, self.r0)

    @staticmethod
    def get_id(atom_types):
        """Return the unique identifier (sorted atom types)"""
        # prm_id = atom_types - sorted, "CA CB" instead of "CB CA", to prevent double entries
        return " ".join(sorted(atom_types))


class _PrmAngle(object):
    def __init__(self, atom_types, fc, theta0, comment=None):
        self.fc = fc
        self.theta0 = theta0
        self.prm_id = self.get_id(atom_types)
        self.atom_types = self.prm_id.split()
        self._comment = comment

    @property
    def comment(self):
        if self._comment:
            return " # {}".format(self._comment) 
        else:
            return ""

    def __repr__(self):
        return "_PrmAngle({}, {})".format(self.prm_id, self.strval)

    @property
    def strval(self):
        """Return parameter values in string format."""
        return "fc={:.3f}, th0={:.3f})".format(self.fc, self.theta0)

    @staticmethod
    def get_id(atom_types):
        """Return the unique identifier (sorted atom types)"""
        # smaller of forward or reverse lists, "CA OS CT" not "CT OS CA"
        return " ".join(min(atom_types, list(reversed(atom_types))))


class _PrmTorsion(object):
    def __init__(self, atom_types, comment=None):
        self.is_generic = True if "?" in atom_types else False
        self.periodicities = []
        self.fcs = []
        self.phases = []
        self.npaths = []
        self.prm_id = self.get_id(atom_types)
        self.atom_types = self.prm_id.split()
        self._comment = comment

    def __repr__(self):
        fcs = ", ".join("{:.4f}".format(fc) for fc in self.fcs)
        prds = ",".join("{:.1f}".format(prd) for prd in self.periodicities)
        phases = ",".join("{:.1f}".format(phase) for phase in self.phases)
        npaths = ",".join("{:.1f}".format(path) for path in self.npaths)

        return "_PrmTorsion({}, {})".format(self.prm_id, self.strval)
            
    @property
    def strval(self):
        """Return parameter values in string format."""
        fcs = ", ".join("{:.4f}".format(fc) for fc in self.fcs)
        prds = ",".join("{:.1f}".format(prd) for prd in self.periodicities)
        phases = ",".join("{:.1f}".format(phase) for phase in self.phases)
        npaths = ",".join("{:.1f}".format(path) for path in self.npaths)

        return "fcs=({}), periodicities=({}), phi0=({}), "\
               "npaths=({})".format(fcs, prds, phases, npaths)

    @property
    def comment(self):
        if self._comment:
            return " # {}".format(self._comment) 
        else:
            return ""

    @staticmethod
    def get_id(atom_types):
        """Return the unique identifier (sorted atom types)"""
        # smaller of forward or reverse lists, "CA CB CG OG1"
        # instead of "OG1 CG CB CA"
        return " ".join(min(atom_types, list(reversed(atom_types))))


    def add_prm(self, fc, periodicity, phase, npaths):
        # a torsion is comprised of several functions with different
        # periodicities
        # if periodicity already exists, raise exception
        prd = abs(periodicity)
        if prd in self.periodicities:
            raise ValueError("Duplicate parameter - periodicity")

        self.periodicities.append(prd)
        self.fcs.append(fc)
        self.phases.append(phase)
        self.npaths.append(npaths)

    def get_prms(self):
        """Return a list of the torsion parameters.
        
        The list is in reverse order of periodicity,
        all but last periodicities are set to negative value.

        [ (fc3, -periodicity3, phase3, npaths3),
          (fc2, -periodicity2, phase2, npaths2),
          (fc1, periodicity1, phase1, npaths1) ]
        """
        prms = sorted(zip(self.fcs, self.periodicities,
                          self.phases, self.npaths), 
                          key=lambda x: x[1], reverse=True)
        rl = []
        for i,prm in enumerate(prms):
            prm = list(prm)
            if i != len(prms)-1:
                prm[1] = -prm[1]
            rl.append(prm)
        return rl


class _PrmImproper(object):
    def __init__(self, center_atom_type, other_atom_types,
                 fc, phi0, periodicity=2.0, comment=None):
        self.is_generic = True if "?" in other_atom_types else False
        self.fc = fc
        self.phi0 = phi0
        if abs(periodicity) - 2.0 > 0.00001:
            raise QPrmError("Q does not support periodicity != 2.0 in impropers "
                            "({})".format(self))
        self.periodicity = periodicity
        self.prm_id = self.get_id(center_atom_type, other_atom_types)
        self.atom_types = self.prm_id.split()
        self._comment = comment

    @property
    def comment(self):
        if self._comment:
            return " # {}".format(self._comment) 
        else:
            return ""

    @staticmethod
    def get_id(center_atom_type, other_atom_types):
        """Return the unique identifier (sorted atom types)"""
        # sorted, with center (second) fixed, "CB CG OG1 OG2"
        # instead of "OG1 CG CB OG2" or "CB CG OG2 OG1" or ...

        # have to fix for general impropers (with wildcards)
        #
        # get the "not center" and "not wildcard" atoms
        ats = [x for x in sorted(other_atom_types) if x != "?"]
        if len(ats) == 3:
            # switch 2 and 3
            a_types = [ats[0], center_atom_type, ats[1], ats[2]]
        elif len(ats) == 2:
            a_types = ["?", center_atom_type, ats[0], ats[1]]
        elif len(ats) == 1:
            a_types = ["?", center_atom_type, ats[0], "?"]
        else:
            raise QPrmError("Improper with three wildcard atom types?")
        return " ".join(a_types)

    def __repr__(self):
        return "_PrmImproper({}, {})".format(self.prm_id, self.strval)

    @property
    def strval(self):
        """Return parameter values in string format."""
        return "fc={:.3f}, phi0={:.3f}".format(self.fc, self.phi0)


