#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# MIT License
#
# Copyright (c) 2018  Miha Purg <miha.purg@gmail.com>
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
"""
Module for wrapping Qcalc functionality.
Contains classes for running Qcalc (QCalc), generating inputs (QCalcInput),
parsing output (QCalcOutput).
"""

from __future__ import absolute_import, unicode_literals, division
import six
from six.moves import zip

import subprocess
import logging
import re
from collections import OrderedDict

from Qpyl.common import __version__, raise_or_log, DataContainer

logger = logging.getLogger(__name__)


class QCalcError(Exception):
    pass

class QCalc(object):
    """Class for running Qcalc.

    Args:
        qcalc_exec (string):  Qcalc executable filename
    """

    def __init__(self, qcalc_exec):
        self.qcalc_exec = qcalc_exec
        self.process = None

    def run(self, qcalc_input_str, workdir=None):
        """Run the calculation and parse the output.

        Args:
            qcalc_input_str (string):  qcalc input
            workdir (string, optional):  working directory
        Returns:
            qcalc_output_str (string):  qcalc stdout

        Raises QCalcError.
        """
        try:
            self.process = subprocess.Popen(self.qcalc_exec,
                                            stdin=subprocess.PIPE,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            cwd=workdir)

        except OSError as error_msg:
            raise QCalcError("Problem when running qcalc: {}"
                             "".format(error_msg))

        # "\n" is added to fix the blocking qcalc issue
        qcalc_input_str += "\n"
        # convert to bytes (Py3+)
        stdin = qcalc_input_str.encode("utf-8")
        stdout, stderr = self.process.communicate(stdin)

        # stderr is useless in qcalc
        return stdout.decode("utf-8")


class QCalcInput(object):
    """Class for creating qcalc inputs.

    Args:
        top (string): topology filename
        dcds (list of strings): trajectory filenames
        fep (string, optional): FEP file filename
        lambdas (list of floats, optional): lambda values

    """
    CALC_RMSD = 1
    CALC_DIST = 3
    CALC_ANGLE = 4
    CALC_TORSION = 5
    CALC_RES_NB = 8
    CALC_RMSF = 13

    def __init__(self, top, dcds, fep=None, lambdas=None):
        self.top = top
        self.dcds = dcds
        self.fep = fep
        self.lambdas = lambdas
        self.actions = []
        # [ (3, ["1231 123"]),
        #   (8, ["144", "144", "res 314", "1554", "1555", "."]) ]

    def add_dist(self, atom1, atom2):
        """Add a distance/bond_energy calculation.

        Args:
            atom1 (int): index of atom in topology
            atom2 (int): index of atom in topology

        """
        self.actions.append((self.CALC_DIST, ["{} {}".format(atom1, atom2)]))

    def add_angle(self, atom1, atom2, atom3):
        """Add a angle/angle_energy calculation.

        Args:
            atom1 (int): index of atom in topology
            atom2 (int): index of atom in topology
            atom3 (int): index of atom in topology

        """
        self.actions.append((self.CALC_ANGLE, ["{} {} {}".format(atom1,
                                                                 atom2,
                                                                 atom3)]))
    def add_torsion(self, atom1, atom2, atom3, atom4):
        """Add a torsion/torsion_energy calculation.

        Args:
            atom1 (int): index of atom in topology
            atom2 (int): index of atom in topology
            atom3 (int): index of atom in topology
            atom4 (int): index of atom in topology

        """
        self.actions.append((self.CALC_TORSION,
                             ["{} {} {} {}".format(atom1, atom2,
                                                   atom3, atom4)]))


    def add_rmsd(self, masks):
        """Add a RMSD calculation.

        Args:
            masks (list of strings): masks used for the calculation
                                     (eg. ["res 314", "1245", "1246", ...])
        """
        if isinstance(masks, six.string_types):
            masks = [masks,]
        # rmsd.out is due to the unfortunate feature added in Q c79ef672400
        self.actions.append((self.CALC_RMSD, masks + [".",] + ["rmsd.out",]))

    def add_residue_nb_mon(self, resid_first, resid_last, masks):
        """Add a residue nonbond monitor calculation.

        This is used mostly for calculating LRA group contributions.

        Args:
            resid_first (int):  index of first residue
            resid_las (int) : index of last residue
            masks (list of strings):  masks used to describe the Q region\
                                      (eg. ["res 314", "1245", "1246", ...])

        """
        if isinstance(masks, six.string_types):
            masks = [masks,]
        lines = [resid_first, resid_last] + masks + ["."]
        self.actions.append((self.CALC_RES_NB, lines))


    def get_string(self):
        """Return string representation of qcalc input.

        Raises QCalcError if no actions were defined.
        """

        if not self.actions:
            raise QCalcError("No actions defined in input.")

        inp_str = []
        inp_str.append(self.top)
        if self.fep and self.lambdas:
            inp_str.append(self.fep)
            inp_str.append(" ".join([str(x) for x in self.lambdas]))
        else:
            inp_str.append(".")

        for action_num, action_lines in self.actions:
            inp_str.append(action_num)
            for line in action_lines:
                inp_str.append(line)

        inp_str.append("go")
        inp_str.extend(self.dcds)
        inp_str.append(".")


        return "\n".join([str(x) for x in inp_str])



class QCalcOutput(object):
    """Class for parsing qcalc output.

    It stores all the data in DataContainer objects inside the
    dictionary 'results'.

    Args:
        qcalc_output (string):  qcalc output

    """
    _CALCLIST_RE = re.compile("List of defined calculations.*?"
                              "(^.*$).*?Calculation results",
                              re.DOTALL | re.MULTILINE)

    _RES_RE = re.compile("Calculation results.*?"
                         "(?=Average nonbonded|\Z)",
                         re.DOTALL)

    _RESNB_RE = re.compile("Average nonbonded.*", re.DOTALL)

    _VERSION_RE = re.compile("Build number\s*(\S*)")


    def __init__(self, qcalc_output):
        self.qcalc_output = qcalc_output
        self.qcalc_version = None
        self.results = OrderedDict()
        self._parse()

    def _parse(self):
        # find the version
        try:
            self.qcalc_version = self._VERSION_RE.findall(self.qcalc_output)[0]
        except IndexError:
            self.qcalc_version = "Unknown, likely ancient"
        # look for errors
        err = "\n".join(re.findall("ERROR.*", self.qcalc_output))
        if err:
            raise QCalcError("Errors in qcalc output: {}".format(err))

        # parse the list of calculations
        calc_list = self._CALCLIST_RE.findall(self.qcalc_output)
        if not calc_list:
            raise QCalcError("Failed to parse qcalc output")

        for line in calc_list[0].split("\n"):
            lf = line.split()
            calc_i = lf[0]
            if "Root Mean Square Deviation" in line:
                self.results[calc_i] = DataContainer(["Frame", "RMSD"])
            elif "distance between" in line:
                self.results[calc_i] = DataContainer(["Frame", "distance"])
            # TODO: extract the energy as well
            elif "distance, bond energy between" in line:
                self.results[calc_i] = DataContainer(["Frame", "distance"])
            # TODO: extract the energy as well
            elif "distance, qbond energy between" in line:
                self.results[calc_i] = DataContainer(["Frame", "distance"])
            elif "angle between" in line:
                self.results[calc_i] = DataContainer(["Frame", "angle"])
            elif "angle, angle energy between" in line:
                self.results[calc_i] = DataContainer(["Frame", "angle"])
            elif "angle, qangle energy between" in line:
                self.results[calc_i] = DataContainer(["Frame", "angle"])
            elif "torsion between" in line:
                self.results[calc_i] = DataContainer(["Frame", "torsion"])
            elif "torsion, torsion energy between" in line:
                self.results[calc_i] = DataContainer(["Frame", "torsion"])
            elif "torsion, qtorsion energy between" in line:
                self.results[calc_i] = DataContainer(["Frame", "torsion"])
            elif "nonbond monitor for residues" in line:
                pass
            else:
                logger.warning("Ignoring unknown QCalc5 results: {}"
                               "".format(line))


        # parse the normal results (distances, rmsds)
        res_list = self._RES_RE.findall(self.qcalc_output)
        if not res_list:
            raise QCalcError("Failed to parse qcalc output")

        # skip first row (--- Calculation results ---)
        res_list = res_list[0].split("\n")[1:]
        colheaders = res_list.pop(0)

        coltitles = []
        colheaders = colheaders.replace(": ", ":") #fix
        for colheader in colheaders.split():
            if ":" in colheader:
                colheader, calctype = colheader.split(":")
                if not calctype:
                    continue  # residue nonbond calc
            coltitles.append(colheader)

        if coltitles and res_list:
            tmpdata = DataContainer(coltitles)
            for line in res_list:
                lf = line.split()
                if not lf:
                    continue
                tmpdata.add_row(lf)

        for k, datac in self.results.items():
            for i, v in enumerate(zip(*tmpdata.get_columns(columns = [k,]))):
                datac.add_row((i, float(v[0])))

        # parse the average residue nonbond energies (if they exist)
        res_resnb = self._RESNB_RE.findall(self.qcalc_output)
        if res_resnb:
            self.results["gc"] = DataContainer(["Residue", "E_LJ", "E_EL"])
            # skip two lines 
            # TODO: extract qatoms indexes?
            res_resnb = res_resnb[0].split("\n")[2:]
            for line in res_resnb:
                lf = line.split()
                if lf:
                    resid, elj, eel = int(lf[0]), float(lf[1]), float(lf[2])
                    self.results["gc"].add_row((resid, elj, eel))

