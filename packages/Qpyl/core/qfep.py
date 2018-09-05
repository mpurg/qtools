#!/usr/bin/env python2
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
Module for wrapping Qfep functionality.
Contains classes for generating Qfep inputs (QFepInput),
performing system calls to Qfep (QFep), and parsing Qfep output (QFepOutput).
"""


import subprocess
import re
import logging

from Qpyl.common import __version__, DataContainer

logger = logging.getLogger(__name__)

class QFepError(Exception):
    pass

class QFepOutputError(Exception):
    pass

class _QFepHeader(object):
    """Class for parsing and storing data from the Qfep output header.

    Contains the Qfep version, modification date, mapping parameters.

    Usage:
    version = _QFepHeader.qfep_version
    nstates = _QFepHeader.nstates


    Args:
        header_string (string):  header string in qfep output

    """

    # TODO: works only for 2 states and Hij=const 
    _MAP_PRMS_RE = {"qfep_version": re.compile(r"version\s*(\S+)"),
                    "qfep_moddate": re.compile(r"modified on\s*(\d+)"),
                    "nfiles": re.compile(r"Number of files\s*=\s*(\d+)"),
                    "nstates": re.compile(r"Number of states\s*=\s*(\d+)"),
                    "kT": re.compile(r"kT\s*=\s*(\S+)"),
                    "pts_skip": re.compile(r"points to skip\s*=\s*(\d+)"),
                    "bins": re.compile(r"Number of gap-bins\s*=\s*(\d+)"),
                    "min_pts_bin": re.compile(r"points per bin=\s*(\d+)"),
                    "alpha": re.compile(r"Alpha for state  2\s*=\s*(\S+)"),
                    "Hij": re.compile(r"eta_ij, r_xy0: #\s*\d+\s*\d+\s*(\S+)")  }

    def __init__(self, header_string):
        self._header_string = header_string
        self.qfep_version = None
        self.qfep_moddate = None
        self.nfiles = None
        self.nstates = None
        self.kT = None
        self.pts_skip = None
        self.bins = None
        self.min_pts_bin = None
        self.alpha = None
        self.Hij = None

        self._parse()

    def _parse(self):
        for map_prm_key in self._MAP_PRMS_RE:
            c = self._MAP_PRMS_RE[map_prm_key].search(self._header_string)
            if c:
                mprm = c.group(1).strip()
                try:
                    if map_prm_key in ["nfiles", "nstates", "pts_skip",
                                       "bins", "min_pts_bin"]:
                        mprm = int(mprm)
                    elif map_prm_key in ["kT", "alpha", "Hij"]:
                        mprm = float(mprm)
                except ValueError:
                    pass
                finally:
                    setattr(self, map_prm_key, mprm)

        # nstates is a must-have
        if self.nstates == None:
            raise QFepOutputError("Failed to read number of states.")

        # set version to more informative than None
        if not self.qfep_version:
            self.qfep_version = "Unknown, likely ancient"

        


###############################################################################


class _QFepPart0(object):
    """Class for parsing and storing data from Part0 in Qfep output.

    Part0 contains average energies and number of points for every
    EVB state in every frame (energy file).

    If parsing is unsuccessful QFepOutputError is raised,
    else all the data is stored in DataContainer object 'data_state'.
    LRA and reorganization energies can be calculated with functions
    calc_lra.

    Args:
        part0_string (string):  string of Part0 in qfep output
        num_evb_states (int):  number of EVB states
        calc_index (int):  calculation index (0 for normal full calc,\
                           1,2... for whatever follows (exclusion, qcp,..)

    Examples:
        >>> state1_data = _QFepPart0.data_state[0].get_rows()
        >>> cols = ["Lambda", "EQtot"]
        >>> state2_EQtot_lambda = _QFepPart0.data_state[1].get_rows(columns=cols)

    """

    _PART0_HEADER = "# file             state   pts   lambda    EQtot   "\
                    "EQbond  EQang   EQtor   EQimp    EQel   EQvdW  Eel_qq  "\
                    "EvdW_qq Eel_qp  EvdW_qp Eel_qw EvdW_qw Eqrstr"

    _COLUMN_TITLES = ["file", "state", "pts", "Lambda", "EQtot", "EQbond",
                      "EQang", "EQtor", "EQimp", "EQel", "EQvdW", "Eel_qq",
                      "EvdW_qq", "Eel_qp", "EvdW_qp", "Eel_qw", "EvdW_qw",
                      "Eqrstr"]

    def __init__(self, part0_string, num_evb_states, calc_index):
        self._part0_string = part0_string
        self._num_evb_states = num_evb_states
        self._calc_index = calc_index
        self.data_state = [DataContainer(self._COLUMN_TITLES) for _ in
                                                        range(num_evb_states)]
        self._parse(calc_index)
        for e_dc in self.data_state:
            if not e_dc.get_rows():
                raise QFepOutputError("Part0 is empty (no rows).")

    def calc_lra(self, lambda_a, lambda_b):
        """Calculate LRA and reorganization energies between two states.

        LRA = 0.5*(<E2-E1>_10+<E2-E1>_01)
        REO = 0.5*(<E2-E1>_10-<E2-E1>_01)

        E1 == Potential energy of state A
        E2 == Potential energy of state B
        <>_10 == Configuration space A (lambda_a)
        <>_01 == Configuration space B (lambda_b)

        E2_10 == Potential energy of state B at lambda_a

        Args:
            lambda_a (float):  lambda value of first state, usually 1.0
            lambda_b (float):  lambda value of second state, usually 0.0

        Returns:
            lra (DataContainer):  LRA and reorganization energies,\
                                  as well as contributions from\
                                  individual states
        """

        if self._num_evb_states != 2:
            raise QFepOutputError("LRA works only with two states")

        lra = DataContainer(["E_type", "(E2-E1)_10", "(E2-E1)_01",
                             "LRA", "REORG"])

        e1_a, e1_b, e2_a, e2_b = None, None, None, None
        # get the appropriate rows of energies
        # note that these energies are not scaled by lambda
        # [4:] ignores 'file', 'state', 'points' and 'lambda'
        for row in self.data_state[0].get_rows():
            if abs(row[3] - lambda_a) < 1e-7:
                e1_a = row[4:]
            if abs(row[3] - lambda_b) < 1e-7:
                e1_b = row[4:]
        # lambda2 in data_state[1] is actually (1-lambda), correct for that
        for row in self.data_state[1].get_rows():
            if abs((1 - row[3]) - lambda_a) < 1e-7:
                e2_a = row[4:]
            if abs((1 - row[3]) - lambda_b) < 1e-7:
                e2_b = row[4:]

        if not e1_a:
            raise QFepOutputError("LRA: No energy values for lambda == '{}'"
                                  "".format(lambda_a))
        if not e1_b:
            raise QFepOutputError("LRA: No energy values for lambda == '{}'"
                                  "".format(lambda_b))

        la, lb = lambda_a, lambda_b
        # calculate total E=(l1*E1 + l2*E2) energies
        e1_state1 = [la*e1a + (1-la)*e2a for e1a, e2a in zip(e1_a, e2_a)]
        e1_state2 = [la*e1b + (1-la)*e2b for e1b, e2b in zip(e1_b, e2_b)]

        e2_state1 = [lb*e1a + (1-lb)*e2a for e1a, e2a in zip(e1_a, e2_a)]
        e2_state2 = [lb*e1b + (1-lb)*e2b for e1b, e2b in zip(e1_b, e2_b)]

        # (E2-E1)_10    (reactant state) = First row E2 - E1
        # (E2-E1)_01    (products state) = Last row E2 - E1
        des_st1 = [e2 - e1 for e1, e2 in zip(e1_state1, e2_state1)]
        des_st2 = [e2 - e1 for e1, e2 in zip(e1_state2, e2_state2)]

        # LRA=0.5*(<E2-E1>_10+<E2-E1>_01)
        # REO=0.5*(<E2-E1>_10-<E2-E1>_01)
        des_st1_st2 = zip(des_st1, des_st2)
        es_lra = [0.5 * (de_st1 + de_st2) for de_st1, de_st2 in des_st1_st2]
        es_reo = [0.5 * (de_st1 - de_st2) for de_st1, de_st2 in des_st1_st2]

        e_types = self.data_state[0].column_titles[4:]

        for row in zip(e_types, des_st1, des_st2, es_lra, es_reo):
            lra.add_row(row)

        return lra


    def _parse(self, calc_index):
        lines = self._part0_string.split('\n')
        # the first line is a comment
        lines.pop(0)
        # comment with column names
        header = lines.pop(0).strip()
        if header != self._PART0_HEADER:
            raise QFepOutputError("Part0 has a wrong header, did the qfep "
                                  "binary change?")
        n_lines_parsed = 0
        lines_to_read = range(calc_index*self._num_evb_states,
                              (calc_index+1)*self._num_evb_states)
        for line in lines:
            # fix for Q version df165865
            if "Could not read file header!" in line:
                continue
            # remove comments
            line = re.split("#|\!", line)[0].strip()
            # fix 'sticky' values in old versions of Q
            line = re.sub("(\d)-(\d)", "\g<1> -\g<2>", line)

            if not line:
                continue
            if "-->" in line:
                n_lines_parsed = 0
                continue
            

            if n_lines_parsed not in lines_to_read:
                n_lines_parsed += 1
                continue
            n_lines_parsed += 1

            cols = line.split()
            enfile = cols[0]
            state = int(cols[1])
            points = int(cols[2])
            lamb = float(cols[3])
            energies = [float(x) for x in cols[4:]]
            row = [enfile, state, points, lamb,]
            row.extend(energies)
            self.data_state[state-1].add_row(row)


###############################################################################


class _QFepPart1(object):
    """Class for parsing and storing data from Part1 in Qfep output.

    Part1 contains free energies vs. lambda (FEP).

    If parsing is unsuccessful QFepOutputError is raised,
    else all the data is stored in DataContainer object 'data'.

    Args:
        part1_string (string):  string of Part1 in qfep output

    Usage:
    >>> cols = ["Lambda", "dG"]
    >>> dG_lambda = _QFepPart1.data.get_rows(columns=cols)

    """

    _PART1_HEADER = "# lambda(1)      dGf sum(dGf)      dGr sum(dGr)     <dG>"
    _COLUMN_TITLES = ["Lambda", "dGf", "sum_dGf", "dGr", "sum_dGr", "dG"]

    def __init__(self, part1_string):
        self._part1_string = part1_string
        self.data = DataContainer(self._COLUMN_TITLES)

        self._parse()
        if not self.data.get_rows():
            raise QFepOutputError("Part1 is empty (no rows).")

    @property
    def dg(self):
        """Return final dG(lambda)   (FEP)"""
        return self.data.get_columns(["dG"])[0][-1]

    def _parse(self):

        lines = self._part1_string.split('\n')
        # the first line is a comment
        lines.pop(0)

        ## In newer versions of Q, two additional lines are printed
        # to distinguish between 'full', 'exclusions' and 'qcp'
        # check for the two extra lines and remove them
        if "Calculation" in lines[1]:
            lines = lines[2:]
        # comment with column names
        header = lines.pop(0).strip()
        if header != self._PART1_HEADER:
            raise QFepOutputError("Part1 has a wrong header, did the qfep "
                                  "binary change?")
        for line in lines:
            line = re.split("#|\!", line)[0].strip()
            if not line:
                continue
            row = [float(x) for x in line.split()]
            self.data.add_row(row)



###############################################################################



class _QFepPart2(object):
    """Class for parsing and storing data from Part2 in Qfep output.

    Part2 contains the results of the binning process - EVB ground state
    free energy vs. lambda and Egap, as well as the diabatic free energy
    profiles.

    If parsing is unsuccessful QFepOutputError is raised,
    else all the data is stored in DataContainer object 'data'.

    Args:
        part2_string (string):  string of Part2 in qfep output

    Usage:
    >>> cols = ["Lambda", "dGg"]
    >>> dGg_lambda = _QFepPart2.data.get_rows(columns=cols)


    """

    _PART2_HEADER = "# Lambda(1)  bin Energy gap      dGa     dGb     dGg    "\
                    "# pts    c1**2    c2**2"

    _COLUMN_TITLES = ["Lambda", "bin", "Egap", "dGa", "dGb", "dGg", "points",
                      "c1**2", "c2**2"]

    def __init__(self, part2_string):
        self._part2_string = part2_string
        self.data = DataContainer(self._COLUMN_TITLES)

        self._parse()
        if not self.data.get_rows():
            raise QFepOutputError("Part2 is empty (no rows).")


    def _parse(self):
        lines = self._part2_string.split('\n')
        # the first line is a comment
        lines.pop(0)
        # comment with column names
        header = lines.pop(0).strip()
        if header != self._PART2_HEADER:
            raise QFepOutputError("Part2 has a wrong header, did the qfep "
                                  "binary change?")
        for line in lines:
            line = re.split("#|\!", line)[0].strip()
            if not line:
                continue
            row = [float(x) for x in line.split()]
            self.data.add_row(row)


###############################################################################


class _QFepPart3(object):
    """Class for parsing and storing data from Part3 in Qfep output.

    Part3 contains the bin-averaged dGg values, points and squared
    eigenvectors from Part2.

    If parsing is unsuccessful QFepOutputError is raised,
    else all the data is stored in DataContainer object 'data'.

    Args:
        part3_string (string):  string of Part3 in qfep output

    Usage:
    >>> cols = ["Lambda", "dGg"]
    >>> dGg_lambda = _QFepPart3.data.get_rows(columns=cols)


    """

    _PART3_HEADER = "# bin  energy gap  <dGg> <dGg norm> pts  <c1**2> "\
                   "<c2**2> <r_xy>"

    _COLUMN_TITLES = ["bin", "Egap", "dGg", "dGg_norm", "points", "c1**2",
                      "c2**2", "r_xy"]

    def __init__(self, part3_string):
        self._part3_string = part3_string
        self.data = DataContainer(self._COLUMN_TITLES)
        self._dga = None
        self._dg0 = None
        self._maxima_bins = None
        self._minima_bins = None
        self.warning = None

        self._parse()
        if not self.data.get_rows():
            raise QFepOutputError("Part3 is empty (no rows).")


    def _parse(self):
        lines = self._part3_string.split('\n')
        # the first line is a comment
        lines.pop(0)
        # comment with column names
        header = lines.pop(0).strip()
        if header != self._PART3_HEADER:
            raise QFepOutputError("Part3 has a wrong header, did the qfep "
                                  "binary change?")
        for line in lines:
            line = re.split("#|\!", line)[0].strip()
            if not line:
                continue
            row = [float(x) for x in line.split()]
            self.data.add_row(row)

    @property
    def dga(self):
        if self._dga == None:
            self._get_dgs()
        return self._dga

    @property
    def dg0(self):
        if self._dg0 == None:
            self._get_dgs()
        return self._dg0

    @property
    def minima_bins(self):
        if self._minima_bins == None:
            self._get_dgs()
        return self._minima_bins

    @property
    def maxima_bins(self):
        if self._maxima_bins == None:
            self._get_dgs()
        return self._maxima_bins

    def _get_dgs(self):
        # Get minima and maxima without any smoothing.
        # If there is more than one maxima and less or more than 2 minima,
        # raise an exception search for maxima only between 0.2*nbins and
        # 0.8*nbins (bad sampling on the edges can raise an error)
        # Also, save the bins of the minima.

        bins, des, dgs = self.data.get_columns(["bin", "Egap", "dGg_norm"])
        minima, maxima = [], []
        nbins = len(bins)
        for i in range(1, nbins-1):     # from the second to the second last

            dg, dgnext, dgprev = dgs[i], dgs[i+1], dgs[i-1]
            if dgprev >= dg and dg < dgnext:
                minima.append(i)
            elif dgprev <= dg and dg > dgnext and \
                                 i > nbins*0.2 and i < nbins*0.8:
                maxima.append(i)

        if len(minima) > 2 or len(maxima) > 1:
            # Bad sampling, more minima and maxima than wanted.
            # Get the highest maxima from those found so far.
            # Get the absolute minima to the left and to the right of this
            # maxima. Save the warning.
            max1 = max(maxima, key=lambda i: dgs[i])
            react = [(dgs[i], i) for i in minima if i < max1]
            prod = [(dgs[i], i) for i in minima if i > max1]
            try:
                min1 = min(react)[1]   # min() will return tuple with lowest dg
                min2 = min(prod)[1]
            except ValueError:
                # multiple minima on one side, none on the other
                # (starts/ends at the lowest point)
                raise QFepOutputError("Bad reaction free energy profile - "
                                      "reactants minima: {}, products minima: "
                                      "{}".format(len(react), len(prod)))

            self.warning = "Rough Free energy profile ({} minima and {} "\
                           "maxima found), look at the graphs!"\
                           "".format(len(minima), len(maxima))
            maxima = [max1,]
            minima = [min1, min2]

        if len(minima) != 2:
            raise QFepOutputError("Bad reaction free energy profile - {} "
                                  "local minima (instead of 2)"
                                  "".format(len(minima)))
        elif len(maxima) != 1:
            raise QFepOutputError("Bad reaction free energy profile - {} "
                                  "local maxima (instead of 1)"
                                  "".format(len(maxima)))

        self._dga = dgs[maxima[0]] - dgs[minima[0]]
        self._dg0 = dgs[minima[1]] - dgs[minima[0]]
        self._minima_bins = [bins[mini] for mini in minima]
        self._maxima_bins = [bins[maxi] for maxi in maxima]

        # adjust the values in data so that the reactants are zero
        colindex = self.data.column_titles.index("dGg_norm")
        for row in self.data.get_rows():
            row[colindex] = row[colindex] - dgs[minima[0]]


###############################################################################


class QFepOutput(object):
    """Class for parsing and analysing Qfep output.

    All data is stored in DataContainer objects in separate
    _QFepPart0/1/2/3 objects:
    QFepOutput.part0.data[1].get_rows(columns=["EQtot"])
    QFepOutput.part2.data.get_rows(columns=["Egap", "dGg"])

    Some special data can be obtained directly:
    QFepOutput.dGa  # activation free energy
    QFepOutput.dG0  # reaction free energy

    Results of exclusion and QCP calculations are stored in dictionaries
    QFepOutput.exclusions and QanalyseMap.QCP as QFepOutput objects.

    Args:
        qfep_output (string):  qfep output
        _calc_index (int):  used internally for exclusions and QCP

    """

    _PART0_RE = re.compile(r"(# Part 0.*?)# Part 1", re.DOTALL)
    _PART1_RE = re.compile(r"(# Part 1.*?)# Part 2", re.DOTALL)
    _PART2_RE = re.compile(r"(# Part 2.*?)# Part 3", re.DOTALL)
    _PART3_RE = re.compile(r"(# Part 3.*)", re.DOTALL)

    # To extract the full, exclusions and QCP results separately
    # '# Part 1' is added manually to the end of the logfile
    _PART1to1_RE = re.compile(r"(# Part 1.*?)(?=# Part 1)", re.DOTALL)

    _EXCL_RE = re.compile(r"Calculation for system with (\w+) exclusion, "
                          "residues (.*)")

    _QCP_RE = re.compile(r"Calculation for QCP, number of atoms \s+(\d+)")

    _QCP_MASS_RE = re.compile(r"Calculation for QCP Mass Perturbation, number "
                              r"of atoms \s+(\d+)")

    # Header regular expression
    _HEADER_RE = re.compile(r"(# Qfep.*?|qfep version.*?)# Part 0", re.DOTALL)


    def __init__(self, qfep_output, _calc_index=0):
        self._qfep_output = qfep_output + "\n# Part 1"
        self._calc_index = _calc_index

        # parts
        self.header = None
        self.part0 = None
        self.part1 = None
        self.part2 = None
        self.part3 = None

        # exclusions and QCP
        # {"ex_full_386" : QFepOutput_instance,
        #  "QCP": QFepOutput_instance, ... }
        #  "QCP_mass": QFepOutput_instance, ... }
        self.sub_calcs = {}

        self._parse()


    def _parse(self):
        # Parse the qfep output file.
        # Assign _QFepPart0-3 objects to self.part0-3.
        # Store exclusions and qcp as QFepOutput objects

        # parse the header
        c = self._HEADER_RE.search(self._qfep_output)
        if not c:
            raise QFepOutputError("Missing header, ancient Qfep "
                                  "executable?")
        self.header = _QFepHeader(c.group(1))

        # Part 0
        # All of the data (full, excl, qcp) is in one string,
        # but is separated in QFepPart0 based on self._calc_index
        c = self._PART0_RE.search(self._qfep_output)
        if not c:
            raise QFepOutputError("Part 0 not found in qfep output.")
        self.part0 = _QFepPart0(c.group(1), self.header.nstates,
                                self._calc_index)

        # get the other parts as one string (for each calc)
        calcs = self._PART1to1_RE.findall(self._qfep_output)
        if not calcs:
            raise QFepOutputError("Couldn't parse qfep output.")

        # parse parts 1 to 3
        calc_this = calcs.pop(0)
        c1 = self._PART1_RE.search(calc_this)
        c2 = self._PART2_RE.search(calc_this)
        c3 = self._PART3_RE.search(calc_this)
        if not c1:
            raise QFepOutputError("Part 1 not found in qfep output.")
        if not c2:
            raise QFepOutputError("Part 2 not found in qfep output.")
        if not c3:
            raise QFepOutputError("Part 3 not found in qfep output.")
        self.part1 = _QFepPart1(c1.group(1))
        self.part2 = _QFepPart2(c2.group(1))
        self.part3 = _QFepPart3(c3.group(1))

        # parse sub-calcs (exclusions and qcp)
        calc_index = 1
        for calc in calcs:
            # add header and part 0
            full_s = "{}\n{}\n{}".format(self.header._header_string,
                                         self.part0._part0_string,
                                         calc)
            excl = self._EXCL_RE.search(calc)
            if excl:
                # "Calculation for system with full exclusion, "
                # "residues    88  177 " returns ("full", "    88 177")
                #
                # convert to "ex_full_88_177", "ex_el_100", "ex_vdw_221",...
                typ, residues = excl.groups()
                typ = typ.replace("electrostatic", "el").lower()
                residues = "_".join(residues.split())
                k = "ex_{}_{}".format(typ, residues)
                self.sub_calcs[k] = QFepOutput(full_s, _calc_index=calc_index)
            elif self._QCP_RE.search(calc):
                self.sub_calcs["QCP"] = QFepOutput(full_s,
                                                   _calc_index=calc_index)
            elif self._QCP_MASS_RE.search(calc):
                self.sub_calcs["QCP_mass"] = QFepOutput(full_s,
                                                        _calc_index=calc_index)
            else:
                raise QFepOutputError("Impossible error - debug")
            calc_index += 1


###############################################################################


class QFep(object):
    """Class for running Qfep.

    Args:
        qfep_exec (string):  Qfep executable path

    """

    def __init__(self, qfep_exec):
        self.qfep_exec = qfep_exec
        self.process = None

    def run(self, qfep_input_str, workdir=None):
        """Run qfep on the input string and return output.

        Args:
            qfep_input_str (string):  qfep input
            workdir (string, optional):  working directory

        Returns:
            qfep_output_str  (string):  qfep stdout

        Raises QFepError on OSError or stderr.
        """
        try:
            self.process = subprocess.Popen(self.qfep_exec,
                                            stdin=subprocess.PIPE,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            cwd=workdir)
        except OSError as error_msg:
            raise QFepError("Problem when running qfep: {}"
                            "".format(error_msg))

        stdout, stderr = self.process.communicate(qfep_input_str)

        # not sure if this ever happens, but will keep it anyway
        if stderr:
            raise QFepError("QFep wrote to STDERR: {}".format(stderr))

        return stdout

###############################################################################

class QFepInput(object):
    """Class for creating qfep inputs.

    At the moment, supports only:
    - 2 states
    - no predefined Hij
    - single, constant Hij (mu=eta=r0 = 0)

    Args:
        energy_files (list of strings): list of energy files
        hij (float): H12_A offdiagonal constant
        alpha (float): state2 alpha shift
        temperature (float): temperature in Kelvin
        gas_const (float): gas_constant in kcal/mol/K
        points_skip (int): number of energy points to skip in each frame
        gap_bins (int): number of gap-bins
        minpts_bin (int): minimum number of points per bin

    """

    def __init__(self, energy_files, hij, alpha, temperature,
                 gas_const, points_skip, gap_bins, minpts_bin):

        self.energy_files = energy_files
        self.hij = float(hij)
        self.alpha = float(alpha)
        self.temperature = float(temperature)
        self.gas_const = float(gas_const)
        self.points_skip = int(points_skip)
        self.gap_bins = int(gap_bins)
        self.minpts_bin = int(minpts_bin)

    def get_string(self):
        """Return input string for Qfep
        """

        return """\
{frames}         # number of files/frames
2  0                 # number of states and predefined off-diagonals
{RT} {points_skip}        # RT and number of points to skip
{gap_bins}              # number of bins
{minpts_bin}    # minimum points for bin
{alpha}         # state2 shift (alpha)
1                    # number of diagonal elements
1 2 {hij} 0 0 0     # states 1 and 2, A=const=Hij, mu=eta=r0= 0
1 -1                 # linear combination of states ( E = e1 - e2 )
{en_files}
stop""".format(frames=len(self.energy_files),
               RT=self.temperature*self.gas_const,
               points_skip=self.points_skip,
               gap_bins=self.gap_bins,
               minpts_bin=self.minpts_bin,
               alpha=self.alpha,
               hij=self.hij,
               en_files="\n".join(self.energy_files))


