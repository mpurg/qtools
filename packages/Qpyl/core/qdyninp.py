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
from collections import OrderedDict as ODict

from Qpyl.common import __version__

def logical(x):
    if str(x).lower() not in [ "on", "off", "1", "0" ]: raise ValueError

Q_PARAMETERS = ODict( [ ("md", ODict( [ ("steps", int), 
                                        ("random_seed", int),
                                        ("temperature", float), 
                                        ("stepsize", float),
                                        ("bath_coupling", float),
                                        ("initial_temperature", float),
                                        ("separate_scaling", logical ),
                                        ("lrf", logical ),
                                        ("shake_solvent", logical ),  
                                        ("shake_solute", logical ),
                                        ("shake_hydrogens", logical ),
                                        ("shake_heavy", logical ),
                                        ("shake_all_solvent", logical ),
                                        ("shake_all_solute", logical ),
                                        ("shake_all_hydrogens", logical ),
                                        ("shake_all_heavy", logical ),
                                        ("force_rms", logical ) 
                                      ])),
                        ("cut-offs", ODict( [ ("solute_solute", float),
                                              ("solvent_solvent", float),
                                              ("solute_solvent", float),
                                              ("q_atom", float),
                                              ("lrf", float) 
                                            ])),
                        ("sphere", ODict( [ ("centre", str),
                                            ("radius", float),
                                            ("shell_radius", float),
                                            ("shell_force",float),
                                            ("shell_force",float),
                                            ("excluded_force",float),
                                            ("excluded_freeze",logical),
                                            ("exclude_bonded",logical) 
                                          ])),
                        ("pbc", ODict( [ ("pressure_seed", int),
                                         ("rigid_box_centre", logical),
                                         ("constant_pressure", logical),
                                         ("max_volume_displ", float),
                                         ("pressure", float),
                                         ("atom_based_scaling", logical),
                                         ("control_box", str),
                                         ("put_solvent_back_in_box", logical),
                                         ("put_solute_back_in_box", logical)
                                       ])),
                        ("solvent", ODict( [ ("radius", float),
                                             ("centre", str),
                                             ("pack", float),
                                             ("radial_force", float),
                                             ("polarisation", logical),
                                             ("charge_correction", logical),
                                             ("polarisation_force", float),
                                             ("morse_depth", float),
                                             ("morse_width", float),
                                             ("model", str)
                                           ])),
                        ("intervals", ODict( [ ("non_bond", int),
                                               ("output", int),
                                               ("temperature", int),
                                               ("energy", int),
                                               ("trajectory", int),
                                               ("volume_change", int)
                                             ])),
                        ("files", ODict( [ ("topology", str),
                                           ("restart", str),
                                           ("final", str),
                                           ("trajectory", str),
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
                        ("angle_restraints", list)
                      ]) 

class QDynInputError(Exception):
    pass

class QDynInput(object):
    """ 
    Base class for parsing, modifying and generating QDyn inputs.
    The api should look like this:

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

    def __init__(self, input_string="", parameters={}):
        """
        Arguments:
        input_string (str):   string of the qdyn input file
        parameters (dict):   { "MD": { "steps":10000, "stepsize":1.00, ... }, ... }
        """

        self.parameters = {}

        self.update(input_string=input_string, parameters=parameters)


    def _parse_inp(self, input_string):
        """ just extracts the keyword:value pairs into self.parameters, no checking is done at this point """

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
                    raise QDynInputError("Section '%s' appears more than once" % qsection)
                if "group_contribution" in qsection or "restraints" in qsection:
                    parms[ qsection ] = []    # make sure the restraints are cleared if an empty section is defined
                continue
          
            if not qsection:
                raise QDynInputError("Line '%s' not in any qsection" % line)

            if "group_contribution" in qsection:
                parms[ qsection ].append( line )
            elif "restraints" in qsection:
                rest = " ".join( [ "%-6s" % x for x in line.split() ] )  # prettify it
                parms[ qsection ].append( rest )
            elif "lambdas" in qsection:
                parms[ qsection ] = line
            else:
                c = line.strip().split()
                key = c[0]
                try:
                    value = " ".join( c[1:] )
                except IndexError:
                    value = None   # checking is done later in _check_parms
                if qsection not in parms.keys():
                    parms[ qsection ] = {}
                parms[ qsection ][ key ] = value 

        return parms


    def _check_parms(self, parms):
        """ 
        Checks if parameters are supported (typos and such) and if they are of correct type.
        """
        for qsection, qsec_parms in parms.iteritems():
            if qsection not in Q_PARAMETERS:
                raise QDynInputError("Unsupported section: '%s'" % qsection) 
            try:
                if isinstance(qsec_parms, dict):
                    for key,value in qsec_parms.iteritems():
                        exp_type = Q_PARAMETERS[qsection][key]
                        exp_type(value)
            except KeyError:
                raise QDynInputError("Keyword '%s' in section '%s' unsupported" % (key,qsection))
            except ValueError:
                raise QDynInputError("Bad value '%s' for parameter '%s' in Q-section '%s'" % (value, key, qsection) )
           

    def _update_dict(self, d1, d2):
        """
        Updates values in dictionary d1 with values in dictionary d2
        """
        overridden = {}   # contains parameters that were overwritten as tuples (old,new)

        for section, prms in d2.iteritems():
            if "group_contribution" in section or "restraints" in section or "lambdas" in section:
                if section in d1:
                    overridden[section] = (d1[section], prms)
                d1[section] = prms
            else:
                if section not in d1:
                    d1[section] = {}
                for keyword,prm in prms.iteritems():
                    if keyword in d1[section]:
                        if d1[section][keyword] != prm:
                            overridden[section + "/" + keyword] = (d1[section][keyword], prm)
                    d1[section][keyword] = prm
        return overridden



    def update(self, input_string=None, parameters=None):
        """ 
        Updates the parameters with either (or both) an input string or a parameter dictionary.
        Either argument works, parameters overwrites, no argument fails 
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
        """ 
        Raises QDynInputError if required parameters are missing.
        (is not all knowing, please don't rely to much on it)
        """
        # check for nonsense or missing mandatory parameters
        mdp = self.parameters.get( "md", [] )
        fp = self.parameters.get( "files", [] )
        ip = self.parameters.get( "intervals", [] )

        for keyword in ("temperature", "steps", "stepsize"):
            if keyword not in mdp:
                raise QDynInputError("Missing parameter '%s'" % keyword)

        # fep file and lambdas require each other
        if ("fep" in fp and "lambdas" not in self.parameters) or \
           ("fep" not in fp and "lambdas" in self.parameters):
            raise QDynInputError("Parameter 'fep' requires the 'lambdas' section and vice versa")

        # when generating new velocities, both parms need to be present
        if ("initial_temperature" in mdp and "random_seed" not in mdp) or \
             ("initial_temperature" not in mdp and "random_seed" in mdp):
            raise QDynInputError("Parameter 'initial_temperature' requires 'random_seed' and vice versa")

        # if a restart file is not defined, we have to generate new velocities
        if "restart" not in fp and "initial_temperature" not in mdp:
            raise QDynInputError("No restart file, please set 'initial_temperature' and 'random_seed' to generate velocities")

        # since energies are important let's not rely on default values in Q...
        # if an energy file is defined, energy interval must be defined
        # (there is no room for libertarian politics in stupidville)
        if ("energy" not in fp and "energy" in ip) or \
           ("energy" in fp and "energy" not in ip):
            raise QDynInputError("'energy' must be defined in both 'intervals' and 'files' sections")



    def get_string(self, check=True):
        """
        Returns the input as a string.
        Arguments:
            check (boolean):  if True (default),call self.check()
        """

        if check:
            self.check()

        # generate the string
        s = []
        for qsection, qsec_parms in Q_PARAMETERS.iteritems():
            if not qsection in self.parameters:
                continue
            s.append("[%s]" % qsection)
            if "group_contribution" in qsection or "restraints" in qsection:
                s.extend(self.parameters[qsection])
            elif "lambda" in qsection:
                s.append(self.parameters[qsection])
            else:
                for key,value in qsec_parms.iteritems():
                    if key in self.parameters[qsection]:
                        s.append("%-20s %30s" % (key,self.parameters[qsection][key]))

            s.append("")
        return "\n".join(s)

