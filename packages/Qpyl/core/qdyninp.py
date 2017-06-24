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


import copy
import logging
from collections import OrderedDict as ODict

from Qpyl.common import __version__, raise_or_log

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
        parameters (dict):   { "MD": { "steps":10000,
                                       "stepsize":1.00, ... },
                                       ... }
        ignore_errors (boolean, optional):   if set, write error messages to
                                             logger.warning instead of raising
                                             QDynInputError

    Usage:

    try:
        # load and parse
        inp = QDynInput( input_file_string )
        # or inp=QDynInput( parameters={ "md": ... } )

        # update with another input and save the overridden paramaters
        overridden_parms = inp.update( input_file_2_string )

        # update with a dictionary
        new_parameters = { "md":  { "steps" : 100000 },
                                  { "temperature" : 50 } } )
        inp.update( new_parameters )

        # check the input
        inp.check()

        # get the input string
        new_inp_str = inp.get_string()
    except QDynInputError as e:
        print "Problem with input file: " + str(e)


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
                raise QDynInputError("Line '%s' not in any qsection" % line)

            if "group_contribution" in qsection:
                parms[qsection].append(line)
            elif "restraints" in qsection:
                # prettify it
                rest = " ".join(["%-6s" % x for x in line.split()])
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
                if qsection not in parms.keys():
                    parms[qsection] = {}
                parms[qsection][key] = value

        return parms


    def _check_parms(self, parms):
        # Checks if parameters are supported (typos and such)
        # and if they are of correct type.

        for qsection, qsec_parms in parms.iteritems():
            if qsection not in Q_PARAMETERS:
                raise_or_log("Unsupported section: '{}'".format(qsection),
                             QDynInputError, logger, self._ignore_errors)
            try:
                if isinstance(qsec_parms, dict):
                    for key, value in qsec_parms.iteritems():
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

        for section, prms in d2.iteritems():
            if "group_contribution" in section \
                    or "restraints" in section \
                    or "lambdas" in section:
                if section in d1:
                    overridden[section] = (d1[section], prms)
                d1[section] = prms
            else:
                if section not in d1:
                    d1[section] = {}
                for keyword, prm in prms.iteritems():
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

        Arguments:
            check (boolean, optional):  if True (default), call self.check()
            sort (boolean, optional):  if True (default), sort the sections and
                                       keywords according to the order in which
                                       they appear in Q_PARAMETERS
        """

        if check:
            self.check()

        qsections = self.parameters.keys()
        if sort:
            qsections = sorted(qsections,
                        key=lambda x: (Q_PARAMETERS.keys() + [x]).index(x))

        # generate the string such that all the sections and keywords
        s = []
        for qsection in qsections:
            s.append("[%s]" % qsection)
            if "group_contribution" in qsection or "restraints" in qsection:
                s.extend(self.parameters[qsection])
            elif "lambda" in qsection:
                s.append(self.parameters[qsection])
            else:
                keywords = self.parameters[qsection].keys()

                if sort:
                    qkeys = Q_PARAMETERS[qsection].keys() 
                    keywords = sorted(keywords,
                                      key=lambda x: (qkeys + [x]).index(x))

                for key in keywords:
                    if key in self.parameters[qsection]:
                        val = self.parameters[qsection][key]
                        s.append("{:<20} {:>30}".format(key, val))

            s.append("")

        return "\n".join(s)
