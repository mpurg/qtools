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
#


from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import argparse
import logging

from Qpyl.qgeninp import genfeps, QGenfepsError
from Qpyl.common import init_logger, get_version_full

logger = init_logger('Qpyl')

parser = argparse.ArgumentParser(description="""
Script for generating FEP inputs from the 'procedure' file (example can be
found in 'template_examples' folder). It's not pretty but it gets the jobs
done.  
""", add_help=False)

reqarg = parser.add_argument_group("Required")
reqarg.add_argument("fep_proc",
                    help="Input file containing the FEP procedure.")

reqarg.add_argument("relax_input",
                    help="Path to the input file from the last "
                         "relaxation step (to extract all the relevant "
                         "filenames - re,top,fep and lambda values).")

reqarg.add_argument("restraint",
                    help="Sequence restraints applied to topology ('top') "
                         "or relaxed structure ('relax') or "
                         "whatever is used in the relaxation input ('inp').")

optarg = parser.add_argument_group("Optional")
optarg.add_argument("--rs", dest="runscript", default=None,
                    help="shell runscript for Q")

optarg.add_argument("--frames", dest="frames", type=int,
                    help="Number of frames (31,51,101,...). Default={}."
                         "".format(QScfg.get("inputs", "fep_frames")),
                    default=QScfg.get("inputs", "fep_frames"))

optarg.add_argument("--repeats", dest="repeats", type=int,
                    help="number of repeats/replicas. Default={}."
                         "".format(QScfg.get("inputs", "num_repeats")),
                    default=QScfg.get("inputs", "num_repeats"))

optarg.add_argument("--fep", dest="fep", default=None,
                    help="FEP file (default is the one in the input file)."
                         "It can contain atom placeholders.")

optarg.add_argument("--fromlambda", dest="fromlambda",
                    type=float, default=None,
                    help="Starting lambda for state 1. Example: "
                         "--fromlambda 0.45  will go from 0.45,0.55 in "
                         "both directions, to 1.0,0.0 and 0.0,1.0. "
                         "Example2: --fromlambda 0.0 will drive the "
                         "reaction in reverse direction. Default is the "
                         "one in the relaxation input (usually 1.0 - "
                         "starting from the reactants state).")

optarg.add_argument("--pdb", dest="pdb", default=None,
                    help="PDB file created with qprep. Used to replace "
                         "$RESID.ATOMNAME$ placeholders with atom indices "
                         "(eg. $512.N1$ -> 8356).")

optarg.add_argument("--prefix", dest="prefix",
                    help="Prefix for repeat/replica folder names."
                         "(default='{}')".format(QScfg.get("inputs",
                                                           "prefix_rep")),
                    default=QScfg.get("inputs", "prefix_rep"))

optarg.add_argument("--first_frame_eq", dest="first_frame_eq",
                    action="store_true", default=False,
                    help="If set, the first FEP frame will be replaced "
                         "by the last equilibration step (CADEE stuff).")

optarg.add_argument("--ignore_errors", action="store_true", default=False,
                    help="Keyword/parameter checks will no longer be fatal."
                         "Use with care.")

optarg.add_argument("-v", "--version", action="version",
                    version=get_version_full())

optarg.add_argument("-h", "--help", action="help", help="show this "
                    "help message and exit")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

kwargs = {"fep_proc_file": args.fep_proc,
          "relax_input_file": args.relax_input,
          "restraint": args.restraint,
          "energy_list_fn": QScfg.get("files", "en_list_fn"),
          "pdb_file": args.pdb,
          "fep_file": args.fep,
          "runscript_file": args.runscript,
          "frames": args.frames,
          "repeats": args.repeats,
          "fromlambda": args.fromlambda,
          "prefix": args.prefix,
          "first_frame_eq": args.first_frame_eq,
          "ignore_errors": args.ignore_errors}

try:
    gen_dirs = genfeps(**kwargs)
    #print gen_dirs
except QGenfepsError as e:
    print "ERROR: {}".format(e)
    sys.exit(1)

