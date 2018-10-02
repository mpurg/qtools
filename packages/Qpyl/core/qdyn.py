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
Module for wrapping Qdyn functionality.
Contains classes for generating and parsing Qdyn inputs (QDynInput)
and parsing Qdyn output (QDynOutput).
"""


from __future__ import absolute_import, unicode_literals, division
import six
from six.moves import range
import re
import copy
import logging
from collections import OrderedDict as ODict

from Qpyl.common import __version__, raise_or_log, DataContainer, gzopen

logger = logging.getLogger(__name__)

def _str_in(x, list_of_values):
    if str(x).lower() not in list_of_values:
        raise ValueError

_logical    = lambda x: _str_in(x, ["on", "off", "1", "0"])
_integrator = lambda x: _str_in(x, ["leap-frog", "velocity-verlet"])
_thermostat = lambda x: _str_in(x, ["berendsen", "langevin", "nose-hoover"])

Q_PARAMETERS = ODict( [ ("md", ODict([("steps", int), 
                                      ("random_seed", int),
                                      ("temperature", float),
                                      ("stepsize", float),
                                      ("bath_coupling", float),
                                      ("initial_temperature", float),
                                      ("separate_scaling", _logical),
                                      ("lrf", _logical),
                                      ("shake_solvent", _logical),
                                      ("shake_solute", _logical),
                                      ("shake_hydrogens", _logical),
                                      ("shake_heavy", _logical),
                                      ("shake_all_solvent", _logical),
                                      ("shake_all_solute", _logical),
                                      ("shake_all_hydrogens", _logical),
                                      ("shake_all_heavy", _logical),
                                      ("force_rms", _logical),
                                      ("integrator", _integrator),
                                      ("thermostat", _thermostat),
                                      ("langevin_random", _logical),
                                      ("langevin_friction", float),
                                      ("nhchains", int),
                                      ("nose-hoover_mass", int),
                                     ])),
                        ("cut-offs", ODict([("solute_solute", float),
                                            ("solvent_solvent", float),
                                            ("solute_solvent", float),
                                            ("q_atom", float),
                                            ("lrf", float)
                                           ])),
                        ("sphere", ODict([("centre", str),
                                          ("radius", float),
                                          ("shell_radius", float),
                                          ("shell_force", float),
                                          ("shell_force", float),
                                          ("excluded_force", float),
                                          ("excluded_freeze", _logical),
                                          ("exclude_bonded", _logical)
                                         ])),
                        ("pbc", ODict([("pressure_seed", int),
                                       ("rigid_box_centre", _logical),
                                       ("constant_pressure", _logical),
                                       ("max_volume_displ", float),
                                       ("pressure", float),
                                       ("atom_based_scaling", _logical),
                                       ("control_box", str),
                                       ("put_solvent_back_in_box", _logical),
                                       ("put_solute_back_in_box", _logical)
                                      ])),
                        ("solvent", ODict([("radius", float),
                                           ("centre", str),
                                           ("pack", float),
                                           ("radial_force", float),
                                           ("polarisation", _logical),
                                           ("charge_correction", _logical),
                                           ("polarisation_force", float),
                                           ("morse_depth", float),
                                           ("morse_width", float),
                                           ("model", str)
                                          ])),
                        ("qcp", ODict([("selection", str),
                                       ("qcp_size", str),
                                       ("qcp_kie", _logical),
                                       ("qcp_pdb", str),
                                       ("qcp_write", _logical),
                                       ("qcp_seed", int)
                                      ])),
                        ("intervals", ODict([("non_bond", int),
                                             ("output", int),
                                             ("temperature", int),
                                             ("energy", int),
                                             ("trajectory", int),
                                             ("volume_change", int)
                                            ])),
                        ("files", ODict([("topology", str),
                                         ("restart", str),
                                         ("final", str),
                                         ("trajectory", str),
                                         ("traj_input", str),
                                         ("energy", str),
                                         ("fep", str),
                                         ("restraint", str),
                                         ("water", str)
                                        ])),
                        ("lambdas", str),
                        ("group_contribution", list),
                        ("atom_restraints", list),
                        ("sequence_restraints", list),
                        ("distance_restraints", list),
                        ("angle_restraints", list),
                        ("wall_restraints", list)
                      ])


class QDynInputError(Exception):
    pass


class QDynInput(object):
    """Used for parsing, modifying and generating QDyn inputs.

    Args:
        input_string (str):   string of a qdyn input file
        parameters (dict):   { "MD": { "steps":10000, \
                                       "stepsize":1.00, ... }, \
                                       ... }
        ignore_errors (boolean, optional):   if set, write error messages to \
                                             logger.warning instead of raising\
                                             QDynInputError

    Examples:
        >>> try:
        # load and parse
        >>>     inp = QDynInput( input_file_string )
        # or inp=QDynInput( parameters={ "md": ... } )
        # update with another input and save the overridden paramaters
        >>>     overridden_parms = inp.update(input_file_2_string)
        # update with a dictionary
        >>>     new_parameters = { "md":  { "steps" : 100000 },\
                                          { "temperature" : 50 } } )
        >>>     inp.update(new_parameters)
        >>>     # check the input
        >>>     inp.check()
        >>>     # get the input string
        >>>     new_inp_str = inp.get_string()
        >>> except QDynInputError as e:
        >>>     print("Problem with input file:", str(e))

    """

    def __init__(self, input_string="", parameters={}, ignore_errors=False):
        self._ignore_errors = ignore_errors
        self.parameters = {}
        self.update(input_string=input_string, parameters=parameters)


    def _parse_inp(self, input_string):
        # just extracts the keyword:value pairs into self.parameters,
        # no checking is done at this point

        if not input_string:
            return {}
        parms = {}
        qsection = ""
        for line in input_string.split("\n"):
            # remove comments and strip whitespaces.
            line = line.split("#")[0].strip()
            line = line.split("!")[0].strip()
            # empty lines are useless
            if line == "":
                continue
            # found a qsection
            if line[0] == "[":
                qsection = line.strip("[").strip("]").lower()
                if qsection in parms:
                    raise QDynInputError("Section '{}' appears more than once"
                                         "".format(qsection))
                if "group_contribution" in qsection or "restraints" in qsection:
                    # make sure the restraints are cleared if
                    # an empty section is defined
                    parms[qsection] = []
                continue

            if not qsection:
                raise QDynInputError("Line '{}' not in any qsection."
                                     "".format(line))

            if "group_contribution" in qsection:
                parms[qsection].append(line)
            elif "restraints" in qsection:
                # prettify it
                rest = " ".join(["{:<6}".format(x) for x in line.split()])
                parms[qsection].append(rest)
            elif "lambdas" in qsection:
                parms[qsection] = line
            else:
                c = line.strip().split()
                key = c[0]
                try:
                    value = " ".join(c[1:])
                except IndexError:
                    value = None   # checking is done later in _check_parms
                if qsection not in parms:
                    parms[qsection] = {}
                parms[qsection][key] = value

        return parms


    def _check_parms(self, parms):
        # Checks if parameters are supported (typos and such)
        # and if they are of correct type.

        for qsection, qsec_parms in six.iteritems(parms):
            if qsection not in Q_PARAMETERS:
                raise_or_log("Unsupported section: '{}'".format(qsection),
                             QDynInputError, logger, self._ignore_errors)
            try:
                if isinstance(qsec_parms, dict):
                    for key, value in six.iteritems(qsec_parms):
                        exp_type = Q_PARAMETERS[qsection][key]
                        exp_type(value)
            except KeyError:
                raise_or_log("Unknown keyword '{}' in section '{}'"
                             "".format(key, qsection),
                             QDynInputError, logger, self._ignore_errors)
            except ValueError:
                raise_or_log("Bad value '{}' for parameter '{}' in section "
                             "'{}'".format(value, key, qsection),
                             QDynInputError, logger, self._ignore_errors)


    def _update_dict(self, d1, d2):
        # Updates values in dictionary d1 with values in dictionary d2

        # contains parameters that were overwritten as tuples (old,new)
        overridden = {}

        for section, prms in six.iteritems(d2):
            if "group_contribution" in section \
                    or "restraints" in section \
                    or "lambdas" in section:
                if section in d1:
                    overridden[section] = (d1[section], prms)
                d1[section] = prms
            else:
                if section not in d1:
                    d1[section] = {}
                for keyword, prm in six.iteritems(prms):
                    if keyword in d1[section]:
                        if d1[section][keyword] != prm:
                            tmpk = section + "/" + keyword
                            overridden[tmpk] = (d1[section][keyword], prm)
                    d1[section][keyword] = prm
        return overridden



    def update(self, input_string=None, parameters=None):
        """Update/modify the parameters.

        Updates the parameters with either (or both) an input string or
        a parameter dictionary. Either argument works, parameters
        overwrites, no argument fails.
        """

        parms = self._parse_inp(input_string)
        if parameters:
            self._update_dict(parms, copy.deepcopy(parameters))
        elif not input_string:
            raise ValueError("Function requires at least one argument")

        overwritten = self._update_dict(self.parameters, parms)
        self._check_parms(self.parameters)
        return overwritten



    def check(self):
        """Check for missing parameters.

        Raises QDynInputError if required parameters are missing.
        (is not all knowing, please don't rely to much on it)
        If 'ignore_errors' is set, it logs warnings instead.
        """
        # check for nonsense or missing mandatory parameters
        mdp = self.parameters.get("md", [])
        fp = self.parameters.get("files", [])
        ip = self.parameters.get("intervals", [])

        for keyword in ("temperature", "steps", "stepsize"):
            if keyword not in mdp:
                raise_or_log("Missing parameter '{}'".format(keyword),
                             QDynInputError, logger, self._ignore_errors)

        # fep file and lambdas require each other
        if ("fep" in fp and "lambdas" not in self.parameters) or \
           ("fep" not in fp and "lambdas" in self.parameters):
            raise_or_log("Parameter 'fep' requires the 'lambdas' section "
                         "and vice versa", QDynInputError,
                         logger, self._ignore_errors)

        # when generating new velocities, both parms need to be present
        if ("initial_temperature" in mdp and "random_seed" not in mdp) or \
             ("initial_temperature" not in mdp and "random_seed" in mdp):
            raise_or_log("Parameter 'initial_temperature' requires "
                         "'random_seed' and vice versa",
                         QDynInputError, logger, self._ignore_errors)

        # if a restart file is not defined, we have to generate new velocities
        if "restart" not in fp and "initial_temperature" not in mdp:
            raise_or_log("No restart file, please set 'initial_temperature' "
                         "and 'random_seed' to generate velocities",
                         QDynInputError, logger, self._ignore_errors)

        # since energies are important let's not rely on default values in Q...
        # if an energy file is defined, energy interval must be defined
        if ("energy" not in fp and "energy" in ip) or \
           ("energy" in fp and "energy" not in ip):
            raise_or_log("'energy' must be defined in both 'intervals' "
                         "and 'files' sections",
                         QDynInputError, logger, self._ignore_errors)



    def get_string(self, check=True, sort=True):
        """Returns the input as a string.

        Args:
            check (boolean, optional):  if True (default), call self.check()
            sort (boolean, optional):  if True (default), sort the sections \
                                       and keywords according to the order \
                                       in which they appear in Q_PARAMETERS

        """

        if check:
            self.check()

        qsections = list(self.parameters.keys())
        if sort:
            qsections = sorted(qsections,
                        key=lambda x: (list(Q_PARAMETERS.keys()) + [x]).index(x))

        # generate the string such that all the sections and keywords
        s = []
        for qsection in qsections:
            s.append("[{}]".format(qsection))
            if "group_contribution" in qsection or "restraints" in qsection:
                s.extend(self.parameters[qsection])
            elif "lambda" in qsection:
                s.append(self.parameters[qsection])
            else:
                keywords = list(self.parameters[qsection].keys())

                if sort:
                    qkeys = list(Q_PARAMETERS[qsection].keys()) 
                    keywords = sorted(keywords,
                                      key=lambda x: (qkeys + [x]).index(x))

                for key in keywords:
                    if key in self.parameters[qsection]:
                        val = self.parameters[qsection][key]
                        s.append("{:<20} {:>30}".format(key, val))

            s.append("")

        return "\n".join(s)



################################################################################


class QDynOutputError(Exception):
    pass

class _QDynHeader(object):
    """Class for parsing and storing data from the Qdyn output header.

    Private, used indirectly through QDynOutput.header.
    Contains the Qdyn version, modification date, various MD parameters.

    Args:
        header_string (string): Qdyn output header (to "Initializing dynamics")
        step_size (float): use in case the output reads 0.000

    """

    _BUILD_RE = re.compile(r"Build number\s*([\d\.]+)")
    _QV506_RE = re.compile(r"QDyn version 5.06")
    _TOP_RE = re.compile(r"Topology file      =\s*(\S+)")
    _STEPS_RE = re.compile(r"Number of MD steps =\s*(\d+)")
    _STEPSIZE_RE = re.compile(r"Stepsize \(fs\)    =\s*([\d\.]+)")
    _FEP_RE = re.compile(r"FEP input file     =\s*(\S+)")
    _NSTATES_RE = re.compile(r"No. of fep/evb states    =\s*(\d+)")
    _OFFDSEC_RE = re.compile(r"(No. of offdiagonal \(Hij\) functions =.*?^$)",
                             re.MULTILINE | re.DOTALL)
    _OFFD_RE = re.compile(r"\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+[\d\.]+\s+[\d\.]+")

    def __init__(self, header_string, step_size=None):
        m = self._BUILD_RE.search(header_string)
        if m:
            self.qdyn_version = m.group(1)
        elif self._QV506_RE.search(header_string):
            self.qdyn_version = "5.06"
        else:
            raise QDynOutputError("Not a valid Q output or Q version "
                                  "is very old...")

        m = self._TOP_RE.search(header_string)
        if m: self.top_file = m.group(1)
        else: raise QDynOutputError("Couldn't find the topology filenamet!?")

        m = self._STEPS_RE.search(header_string)
        if m: self.nsteps = int(m.group(1))
        else: raise QDynOutputError("Couldn't find number of steps!?")

        m = self._STEPSIZE_RE.search(header_string)
        if m: self.stepsize = float(m.group(1))
        else: raise QDynOutputError("Couldn't find the stepsize!?")

        if not step_size:
            if abs(self.stepsize - 0.0) < 1e-8:
                raise QDynOutputError("Can't convert steps to time, stepsize "
                                      "is 0.0 in the Qdyn output (Q sucks). Set "
                                      "the stepsize manually.")
        else:
            if self.stepsize:
                raise QDynOutputError("Will not override the non-zero "
                                      "stepsize in the Qdyn output...")
            else:
                self.stepsize = stepsize

        m = self._FEP_RE.search(header_string)
        if m: self.fep_file = m.group(1)
        else: self.fep_file = None

        if self.fep_file:
            m = self._NSTATES_RE.search(header_string)
            if m: self.nstates = int(m.group(1))
            else: raise QDynOutputError("Couldn't find the number of states!?")
        else:
            self.nstates = 0

        offdsection = self._OFFDSEC_RE.search(header_string).group(1)
        offdgs = self._OFFD_RE.findall(offdsection)
        self.offdiagonals = [(a1, a2) for a1, a2 in offdgs]


class QDynOutput(object):
    """Class for parsing Qdyn output and storing the data.

    Supports Qdyn versions 5.10 or higher.
    Typically used indirectly by wrapper QAnalyseDyns.

    Args:
        qdyn_output (string): Qdyn output filename
        time_unit (string): fs,ps,ns (optional, default is ps)
        step_size (float): use in case the output reads 0.000
        start_time (float): redefine the start time in given units\
                            in case of continuation simulation (default is 0)

    Examples:
        # Load a qdyn output
        >>> qdo = QDynOutput("qdyn.log")
        # list 
        >>> print qdo.data_EQ_Q[0].get_rows(["Time", "El"])
        # print out the Q-Q electrostatic energy
        >>> print qdo.data_EQ_Q[0].get_rows(["Time", "El"])
    """
    # TODO: write examples above



    def __init__(self, qdyn_output,
                 time_unit="ps", step_size=None, start_time=0):

        self._qdyn_output = qdyn_output

        _MAP_TIME = {"fs": 1.0, "ps": 1e-3, "ns": 1e-6}
        if time_unit.lower() not in _MAP_TIME:
            raise QDynOutputError("Timeunit has to be either 'fs',"
                                  "'ps' or 'ns'")
        self._timeconv = _MAP_TIME[time_unit.lower()]
        self._stepsize_user = step_size


        # parse the header
        self.time_begin = start_time
        self.time_unit = time_unit.lower()
        self._parse_header()


        ###  Datacontainer variables for storing all the data
        # temperature
        self.data_temp = DataContainer(["Time", "T_tot", "T_free",
                                        "T_free_solute", "T_free_solvent"])
        # energies
        columns1 = ["Time", "El", "VdW", "Bond",
                    "Angle", "Torsion", "Improper"]
        columns2 = ["Time", "Total", "Fix", "Solvent_rad",
                    "Solvent_pol", "Shell", "Solute"]
        columns3 = ["Time", "Total", "Potential", "Kinetic"]
        self.data_E_solute = DataContainer(columns1)
        self.data_E_solvent = DataContainer(columns1)
        self.data_E_solute_solvent = DataContainer(["Time", "El", "VdW"])
        self.data_E_LRF = DataContainer(["Time", "El"])
        self.data_E_Q_atom = DataContainer(columns1)
        self.data_E_restraints = DataContainer(columns2)
        self.data_E_SUM = DataContainer(columns3)
        # Q energies
        q_columns1 = ("Time", "Lambda", "El", "VdW")
        q_columns2 = ("Time", "Lambda", "El", "VdW", "Bond",
                      "Angle", "Torsion", "Improper")
        q_columns3 = ("Time", "Lambda", "Total", "Restraint")

        self.data_EQ_Q, self.data_EQ_prot = [], []
        self.data_EQ_wat, self.data_EQ_surr = [], []
        self.data_EQ_any, self.data_EQ_SUM = [], []
        for i in range(self.header.nstates):
            self.data_EQ_Q.append(DataContainer(q_columns1))
            self.data_EQ_prot.append(DataContainer(q_columns1))
            self.data_EQ_wat.append(DataContainer(q_columns1))
            self.data_EQ_surr.append(DataContainer(q_columns1))
            self.data_EQ_any.append(DataContainer(q_columns2))
            self.data_EQ_SUM.append(DataContainer(q_columns3))

        # mapping of energy types (label in the output) with containers
        self.map_en_section = ODict([("solute", self.data_E_solute),
                                     ("solvent", self.data_E_solvent),
                                     ("solute-solvent",
                                         self.data_E_solute_solvent),
                                     ("LRF", self.data_E_LRF),
                                     ("Q-atom", self.data_E_Q_atom),
                                     ("restraints", self.data_E_restraints),
                                     ("SUM", self.data_E_SUM)])

        self.map_qen_section = ODict([("Q-Q", self.data_EQ_Q),
                                      ("Q-prot", self.data_EQ_prot),
                                      ("Q-wat", self.data_EQ_wat),
                                      ("Q-surr.", self.data_EQ_surr),
                                      ("Q-any", self.data_EQ_any),
                                      ("Q-SUM", self.data_EQ_SUM)])


        # parse the rest
        self._parse_dyn()
        self.time_end = self.header.nsteps \
                      * self.header.stepsize \
                      * self._timeconv \
                      + self.time_begin



    def _parse_header(self):
        """Parses the header of the Qdyn output (called by init)
        """
        header_string = ""
        try:
            with gzopen(self._qdyn_output) as qdo:
                for line in qdo:
                    header_string += line
                    if "Initialising dynamics" in line:
                        break
        except IOError as e:
            raise QDynOutputError("Could not read the Qdyn output: {}"
                                  "".format(e))

        self.header = _QDynHeader(header_string, step_size=self._stepsize_user)
        self._header_length = len(header_string)




    def _parse_dyn(self):
        """Parses the dynamics part of the Qdyn output (called by init)

        Extracts all the temperatures, energies, Q energies and off-diagonals.
        """

        # tmp temperature vars
        t_free, t_tot = None, None
        temps_q6 = {"Total": [], "Free": [], "Solute": [],
                    "Solvent": [], "time": []}
        # tmp offdiagonal vars
        tmp_offdiags = ODict()
        for atom1, atom2 in self.header.offdiagonals:
            k = "{}_{}".format(atom1, atom2)
            tmp_offdiags[k] = []

        time = self.time_begin
        insection = False
        step = 0
        with gzopen(self._qdyn_output) as qdyn_output:
            qdyn_output.seek(self._header_length)
            for line in qdyn_output:
                lf = line.split()
                if not lf:
                    continue
                if "Initialising dynamics" in line:
                    raise QDynOutputError("Found more than one qdyn_output...",
                                           "Please don't concatenate...")

                # Temperature
                if self.header.qdyn_version > "6":
                    if "temperature at step" in line:
                        # fix for large step numbers
                        lf = line.replace("step", "step ")
                        lf = lf.replace("System", "").split()
                        t_type, t, step = lf[0], float(lf[6]), int(lf[4])
                        temps_q6[t_type].append(t)
                        if t_type == "Total":
                            time = step * self.header.stepsize \
                                 * self._timeconv  + self.time_begin
                            temps_q6["time"].append(time)
                else:
                    # second line with temps (pre Q6)
                    if t_free != None: # second line with temps
                        try:
                            tf_solute = float(lf[1])
                        except: # gas phase
                            tf_solute = 0
                        try:
                            tf_solvent = float(lf[3])
                        except: # gas phase
                            tf_solvent = 0
                        self.data_temp.add_row((time, t_tot, t_free,
                                                tf_solute, tf_solvent))
                        t_free, t_tot = None, None
                    # first line with temps (pre Q6)
                    elif "Temperature at step" in line:
                        # fix for large step numbers
                        lf = line.replace("step", "step ").split()
                        step = int(lf[3].strip(":"))
                        time = step * self.header.stepsize \
                             * self._timeconv  + self.time_begin
                        t_tot, t_free = float(lf[5]), float(lf[7])

                if "Energy summary at step" in line or \
                        "Q-atom energies at step" in line:
                    insection = True
                    step = int(lf[5])
                    time = step * self.header.stepsize \
                         * self._timeconv  + self.time_begin

                elif "FINAL  Energy summary" in line or \
                        "FINAL Q-atom energies" in line:
                    insection = True
                    time = self.header.nsteps * self.header.stepsize \
                         * self._timeconv  + self.time_begin

                elif "===================================================="\
                     "======================" in line:
                    insection = False

                # skip the 0th step
                if step == 0:
                    continue
                elif insection:
                    key = lf[0]
                    if key in self.map_en_section:
                        row = [time,] + [float(x) for x in lf[1:]]
                        self.map_en_section[key].add_row(row)
                    elif key in self.map_qen_section:
                        evb_index = int(lf[1]) - 1
                        row = [time,] + [float(x) for x in lf[2:]]
                        self.map_qen_section[key][evb_index].add_row(row)
                    elif "dist. between" in line:
                        atom1, atom2, dist = lf[8], lf[9], float(lf[11])
                        k = "{}_{}".format(atom1, atom2)
                        tmp_offdiags[k].append([time, dist])

        # join temperatures to one DataContainer (Q6+)
        for i, time in enumerate(temps_q6["time"]):
            try:
                t_solv = temps_q6["Solvent"][i]
            except IndexError:
                # gas phase
                t_solv = 0

            self.data_temp.add_row((time,
                                    temps_q6["Total"][i],
                                    temps_q6["Free"][i],
                                    temps_q6["Solute"][i],
                                    t_solv))

        # join Offdiagonal distances to single DataContainer
        offd_keys = list(tmp_offdiags.keys())
        cts = ["Time",] + offd_keys
        self.data_offdiags = DataContainer(cts)
        for i, (time, _) in enumerate(list(tmp_offdiags.values())[0]):
            row = [time,] + [tmp_offdiags[k][i][1] for k in offd_keys]
            self.data_offdiags.add_row(row)



