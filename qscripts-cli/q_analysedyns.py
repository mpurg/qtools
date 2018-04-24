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
# Standalone CLI script for Qpyl.qanalysis.QAnalyseDyns
# Extracts dynamics information from the logfile and saves to json
#

from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import argparse
import os
from collections import OrderedDict as ODict
import logging

from Qpyl.qanalysis import QAnalyseDyns, QAnalyseDynsError
from Qpyl import plotdata
from Qpyl.common import backup_file, init_logger, get_version_full

logger = init_logger('Qpyl')

parser = argparse.ArgumentParser(description="""
Script for extracting temperatures and energies from QDyn outputs.
Mostly used for debugging.
""", add_help=False)
reqarg = parser.add_argument_group("Required")
reqarg.add_argument("outputs", nargs="+", help="Qdyn output files")
optarg = parser.add_argument_group("Optional")
optarg.add_argument("--plots_out", dest="plots_out",
                    help="Output filename for plot data (default='{}')"
                         "".format(QScfg.get("files", "analysedyns_plots")),
                    default=QScfg.get("files", "analysedyns_plots"))

optarg.add_argument("--stepsize", dest="stepsize", default=None, type=float,
                    help="If the stepsize in your Qdyn output is 0.000, "
                         "define it with this flag.")

optarg.add_argument("--timeunit", dest="timeunit", default="ps",
                    help="Time unit. Options are 'fs', 'ps', 'ns'. "
                         "Default is ps.")

optarg.add_argument("--stride", dest="stride", type=int, default=1,
                    help="Read only every Nth point. Default=1")

optarg.add_argument("-v", "--version", action="version",
                    version=get_version_full())
optarg.add_argument("-h", "--help", action="help", help="show this help "
                    "  message and exit")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

for qdynout in args.outputs:
    if not os.path.lexists(qdynout):
        print "FATAL! File '{}' doesn't exist".format(qdynout)
        sys.exit(1)

try:
    qads = QAnalyseDyns(args.outputs,
                        time_unit=args.timeunit,
                        step_size=args.stepsize)
except QAnalyseDynsError as e:
    print "Error: {}".format(e)
    sys.exit(1)

print qads.get_temp_stats()

plots = qads.get_plotdata(stride=args.stride)

jsonenc = plotdata.PlotDataJSONEncoder(indent=2)
backup = backup_file(args.plots_out)
if backup:
    print "Backed up '{}' to '{}'".format(args.plots_out, backup)
open(args.plots_out, 'w').write(jsonenc.encode(plots))
print "\nWrote '{}'. Use q_plot.py to visualize the plots.".format(args.plots_out)

