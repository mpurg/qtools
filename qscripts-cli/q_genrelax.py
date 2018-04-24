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

from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import argparse
import logging

from Qpyl.qgeninp import genrelax, QGenrelaxError
from Qpyl.common import init_logger, get_version_full

logger = init_logger('Qpyl')

parser = argparse.ArgumentParser(description="""
Script for generating QDyn MD inputs from the 'procedure' file (example can be
found in 'template_examples' folder). It's not pretty but it gets the jobs
done.  
""", add_help=False)

reqarg = parser.add_argument_group("Required")
reqarg.add_argument("relax_proc",
                    help="Input file containing the relaxation procedure")

optarg = parser.add_argument_group("Optional")
optarg.add_argument("--top", dest="top", default=None,
                    help="Path to topology file")

optarg.add_argument("--fep", dest="fep", default=None,
                    help="Path to fep file (if there is one). "
                         "It can contain atom placeholders.")

optarg.add_argument("--rs", dest="runscript", default=None,
                    help="Submission script (for Slurm,SGE,Torque,...) ")

optarg.add_argument("--cont", dest="cont", default=None,
                    help="Continue a previous relaxation, argument is the "
                         "name of the last input file (e.g. 'relax_012.inp')")

optarg.add_argument("--rest", dest="restraint", default=None,
                    help="Sequence restraints applied to: 'top' - topology, "
                         "'cont_inp' - same as in the --cont input, "
                         "'cont_final' - final restart in --cont input. "
                         "Required if --cont is specified, otherwise, "
                         "defaults to 'top'.")

optarg.add_argument("--pdb", dest="pdb", default=None,
                    help="PDB file created with qprep. Used to replace "
                         "$RESID.ATOMNAME$ placeholders with atom indices "
                         "(eg. $512.N1$ -> 8356).")

optarg.add_argument("--outdir", dest="outdir",
                    help="Output directory name. Default='{}'"
                         "".format(QScfg.get("inputs", "relax_dir")),
                    default=QScfg.get("inputs", "relax_dir"))

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

if args.cont == None and args.restraint == None:
    args.restraint = 'top'
elif args.cont != None and args.restraint == None:
    print "FATAL! --rest is required with --cont"
    sys.exit(1)

kwargs = {"relax_proc_file": args.relax_proc,
          "top_file" : args.top,
          "fep_file" : args.fep,
          "cont_file" : args.cont,
          "runscript_file": args.runscript,
          "restraint": args.restraint,
          "pdb_file": args.pdb,
          "outdir": args.outdir,
          "ignore_errors": args.ignore_errors}

try:
    gen_inps = genrelax(**kwargs)
    #print gen_inps
except QGenrelaxError as err:
    print "ERROR: {}".format(err)
    sys.exit(1)


