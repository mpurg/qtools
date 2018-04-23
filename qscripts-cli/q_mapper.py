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

from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import os
import argparse
import logging

from Qpyl.qanalysis import QAnalyseFeps
from Qpyl.qmapping import QMapper
from Qpyl.common import backup_file, init_logger, get_version_full

def main():
    logger = init_logger('Qpyl')

    parser = argparse.ArgumentParser(description="""
    Command-line interface for mapping EVB (or just plain old FEP) simulations
    with QFep. At the moment, supports only single Hij (constant) and alpha.
    For FEP, use Hij=alpha=0.

    By default, all subdirectories in current dir will be used
    for mapping, or current dir if no subdirs are found.
    This can be changed with --dirs.
    """, add_help=False)

    reqarg = parser.add_argument_group("Required")
    reqarg.add_argument("hij", type=float, help="Hij coupling constant")

    reqarg.add_argument("alpha", type=float, help="state 2 shift (alpha)")

    optarg = parser.add_argument_group("Optional")
    optarg.add_argument("--nt", dest='nthreads', type=int,
                        default=QScfg.get("mapping", "nthreads"),
                        help="Number of threads (default = {})"
                             "".format(QScfg.get("mapping", "nthreads")))

    optarg.add_argument("--bins", dest="gap_bins", type=int,
                        default=QScfg.get("mapping", "gap_bins"),
                        help="Number of gap-bins (default={})."
                             "".format(QScfg.get("mapping", "gap_bins")))

    optarg.add_argument("--skip", dest="points_skip", type=int,
                        default=QScfg.get("mapping", "points_skip"),
                        help="Number of points to skip in each frame "
                             "(default={})."
                             "".format(QScfg.get("mapping", "points_skip")))

    optarg.add_argument("--min", dest="minpts_bin", type=int,
                        default=QScfg.get("mapping", "minpts_bin"),
                        help="Minimum points for gap-bin (default={})."
                             "".format(QScfg.get("mapping", "minpts_bin")))

    optarg.add_argument("--temp", dest="temperature", type=float,
                        default=QScfg.get("mapping", "temperature"),
                        help="Temperature (default={})."
                             "".format(QScfg.get("mapping", "temperature")))

    optarg.add_argument("--dirs", nargs="+", dest="mapdirs", default=[],
                        help="Directories to map (default=all subdirs "
                             "in cwd or current dir)")

    optarg.add_argument("--out", dest="outfile",
                        default=QScfg.get("files", "mapper_log"),
                        help="Logfile name (default={})."
                             "".format(QScfg.get("files", "mapper_log")))

    optarg.add_argument("--qfep_exec", dest="qfep_exec",
                        default=QScfg.get("qexec", "qfep"),
                        help="qfep executable path (default={})."
                             "".format(QScfg.get("qexec", "qfep")))

    optarg.add_argument("-v", "--version", action="version",
                        version=get_version_full())
    optarg.add_argument("-h", "--help", action="help", help="show this "
                        "help message and exit")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    mapdirs = args.mapdirs

    # if mapping directories were not passed in as an argument,
    # store all directories in the current folder
    if not mapdirs:
        lsdir = os.listdir(os.getcwd())
        mapdirs = [md for md in lsdir if os.path.isdir(md)]

    mapdirs.sort()

    # if there are no folders in the current working directory,
    # map just the current one
    if not mapdirs:
        mapdirs = [os.getcwd(),]
        print "No subdirectories. Mapping files in current directory only."
    else:
        print "Will use these directories for mapping (use --dirs to "\
              "change this): {}".format(", ".join(mapdirs))

    qmapper_parms = {"mapdirs": mapdirs,
                     "hij": args.hij,
                     "alpha": args.alpha,
                     "nthreads": args.nthreads,
                     "temperature": args.temperature,
                     "points_skip": args.points_skip,
                     "minpts_bin": args.minpts_bin,
                     "gap_bins": args.gap_bins,
                     "qfep_exec": args.qfep_exec,
                     "en_list_fn": QScfg.get("files", "en_list_fn"),
                     "gas_const": QScfg.get("mapping", "gas_const")}

    qmapper = QMapper(**qmapper_parms)
    try:
        qmapper.mapall()
    except KeyboardInterrupt:
        qmapper.kill_event.set()
        raise

    qfep_inp_fn = QScfg.get("files", "qfep_inp")
    qfep_out_fn = QScfg.get("files", "qfep_out")

    # write out the inputs and outputs
    for mapdir, (qfep_inp_str, qfep_out_str) in qmapper.mapped.iteritems():
        qfep_inp = os.path.join(mapdir, qfep_inp_fn)
        qfep_out = os.path.join(mapdir, qfep_out_fn)
        open(qfep_inp, "w").write(qfep_inp_str)
        open(qfep_out, "w").write(qfep_out_str)

    # analyse the outputs
    output_files = [os.path.join(md, qfep_out_fn) for md in qmapper.mapped]
    qafs = QAnalyseFeps(output_files)

    outstr = """
{mapper_details}
Analysis Stats:
{analysis_stats}
Analysis Fails:
FEP: {fails}, EVB: {fails_dg}

Run q_analysefeps.py for more info...
""".format(mapper_details=qmapper.details, analysis_stats=qafs.stats_str,
           fails=len(qafs.failed) or "None",
           fails_dg=len(qafs.failed_dg) or "None")

    print outstr
    backup = backup_file(args.outfile)
    if backup:
        print "# Backed up '{}' to '{}'".format(args.outfile, backup)
    open(args.outfile, "w").write(outstr)
    print "Wrote '{}'...".format(args.outfile)



if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print "\nCtrl-C detected. Quitting..."
        sys.exit(1)

