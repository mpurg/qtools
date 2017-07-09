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
# Standalone CLI script for Qpyl.qanalysis.QAnalyseDyns
# Extracts dynamics information from the logfile and saves to json
#

from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import argparse
import os
from collections import OrderedDict as ODict
import logging

from Qpyl.qanalysis import QAnalyseDyns, QAnalyseDynError
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
                    help="If the stepsize in your log is 0.000, define it "
                         "with this flag.")

optarg.add_argument("--timeunit", dest="timeunit", default="ps",
                    help="Which units of time should the results be in? fs, "
                         "ps or ns? Default is ps.")

optarg.add_argument("--stride", dest="stride", type=int, default=1,
                    help="Read only every Nth point. Default=1")

optarg.add_argument("--skip", dest="skip", type=int, default=0,
                    help="Skip percentage of data points in each log. "
                         "Default=0") 
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
                        timeunit=args.timeunit,
                        stepsize=args.stepsize)
except QAnalyseDynError as e:
    print "Error: " + str(e)
    sys.exit(1)


print qads.get_temp_stats()

time_label = "Time [{}]".format(args.timeunit)
plots = ODict()
plots["temp"] = plotdata.PlotData("Temperature",
                                  xlabel=time_label,
                                  ylabel="T [K]")

plots["offdiags"] = plotdata.PlotData("Offdiagonal distances",
                                      xlabel=time_label,
                                      ylabel="Distance [A]")

t_dc = qads.get_temps(percent_skip=args.skip, stride=args.stride)
t_cs, t_cts = t_dc.get_columns(), t_dc.get_column_titles()
for i, t_ct in enumerate(t_cts[1:]):
    plots["temp"].add_subplot(t_ct, t_cs[0], t_cs[i+1])   # 0==Time


d_dc = qads.get_offdiags(percent_skip=args.skip, stride=args.stride)
d_cs, d_cts = d_dc.get_columns(), d_dc.get_column_titles()
for i, d_ct in enumerate(d_cts[1:]):
    plots["offdiags"].add_subplot(d_ct, d_cs[0], d_cs[i+1])   # 0==Time


for k in qads.en_section_keys:
    key = "E_{}".format(k)
    plots[key] = plotdata.PlotData("Energy: " + k,
                                   xlabel=time_label,
                                   ylabel="Energy [kcal/mol]")
    e_dc = qads.get_energies(k, percent_skip=args.skip,
                             stride=args.stride)
    e_cs, e_cts = e_dc.get_columns(), e_dc.get_column_titles()
    if e_cs:
        for i, e_ct in enumerate(e_cts[1:]):
            plots[key].add_subplot(e_ct, e_cs[0], e_cs[i+1])   # 0==Time


for k in qads.qen_section_keys:
    for evb_state in range(1, qads.n_evb_states + 1):
        key = "EQ{}_{}".format(evb_state, k)
        plots[key] = plotdata.PlotData("Q Energy: {} (state {})"
                                       "".format(k, evb_state),
                                       xlabel=time_label,
                                       ylabel="Energy [kcal/mol]")
        qe_dc = qads.get_q_energies(k, evb_state,
                                    percent_skip=args.skip,
                                    stride=args.stride)
        qe_cs, qe_cts = qe_dc.get_columns(), qe_dc.get_column_titles()
        if qe_cs:
            for i, qe_ct in enumerate(qe_cts[1:]):
                plots[key].add_subplot(qe_ct, qe_cs[0], qe_cs[i+1]) # 0==Time



jsonenc = plotdata.PlotDataJSONEncoder(indent=2)
backup = backup_file(args.plots_out)
if backup:
    print "Backed up '{}' to '{}'".format(args.plots_out, backup)
open(args.plots_out, 'w').write(jsonenc.encode(plots))
print "\nWrote '{}'. Use q_plot.py to visualize the plots.".format(args.plots_out)

