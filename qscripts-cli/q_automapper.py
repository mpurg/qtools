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
# Automapper for Q
# Run the script with no arguments first.
#
# Varies EVB parameters Hij and Gas_shift until desired (experiment or qm)
# activation and reaction free energies are obtained.
# This script should be run in the directory that contains
# the reference reaction's replicas (gas or water), each in its own directory
# (by default all subdirectories in the current dir are used for mapping
# or if there are no subdirectories, the current dir is mapped.
# This can be changed with --dirs dir1 dir2 dir3 ...
#
# Initial guess values for Hij and gas_shift should be relatively close to
# their correct values (+-50) otherwise qfep crashes. If it doesn't converge,
# change the step size (--step), number of iterations (--iter) or the threshold
# (--threshold). For information about other parameters (bins, nt, skip,
# temperature...) see q_mapper.py
#

from qscripts_config import QScriptsConfig as QScfg

import sys
import os
import argparse
import logging
import inspect

from Qpyl.common import backup_file, SpecialFormatter
from Qpyl.qanalysis import QAnalyseFeps
from Qpyl.qmapping import QMapper, QMapperError


def main():
    logger = logging.getLogger('Qpyl')
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(SpecialFormatter())
    logger.addHandler(handler)

    parser = argparse.ArgumentParser(description="""
    Bored of changing your Hij and alpha manually when calibrating 
    your EVB potential? Look no further, this script will do it for you.
    Give it the two reference values (activation and reaction free 
    energy), and initial guesses for hij and alpha and enjoy reddit.
    """, add_help=False)
    reqarg = parser.add_argument_group("Required")
    reqarg.add_argument("ref_dga", type=float,
                        help="Reference activation free energy.")

    reqarg.add_argument("ref_dg0", type=float,
                        help="Reference reaction free energy.")

    reqarg.add_argument("init_hij", type=float,
                        help="Initial guess for Hij (offdiagonal)")

    reqarg.add_argument("init_alpha", type=float,
                        help="Initial guess for alpha (state 2 shift)")

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
                             "in cwd that contain the energy-files list {})."
                             "".format(QScfg.get("files", "en_list_fn")))

    optarg.add_argument("--out", dest="outfile",
                        default=QScfg.get("files", "automapper_log"),
                        help="Logfile name (default={})."
                             "".format(QScfg.get("files", "automapper_log")))

    _args, _, _, _defaults = inspect.getargspec(QMapper.fit_to_reference)
    defs = dict(zip(_args[-len(_defaults):], _defaults))

    optarg.add_argument("--step", dest="step_size", type=float,
                        help="Step size (default={})."
                        "".format(defs["step_size"]),
                        default=defs["step_size"])

    optarg.add_argument("--threshold", dest="threshold", type=float,
                        help="Convergence threshold for dG# and dG0 "
                             "(default={}).".format(defs["threshold"]),
                        default=defs["threshold"])

    optarg.add_argument("--iter", dest="max_iterations", type=int,
                        help="Max number of iterations (default={})."
                        "".format(defs["max_iterations"]),
                        default=defs["max_iterations"])

    optarg.add_argument("--nosingle", dest="nosingle", action="store_true",
                        help="Do not run the first iteration on only 1 dir.")

    optarg.add_argument("--qfep_exec", dest="qfep_exec",
                        default=QScfg.get("qexec", "qfep"),
                        help="qfep5 executable path (default={})."
                             "".format(QScfg.get("qexec", "qfep")))

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    print """\
Attempting to fit to dG# = {} and dG0 = {}
(stepsize = {}, threshold = {}, max iterations = {})
""".format(args.ref_dga, args.ref_dg0, args.step_size,
           args.threshold, args.max_iterations)

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

    qmapper_parms = {"hij": args.init_hij,
                     "alpha": args.init_alpha,
                     "nthreads": args.nthreads,
                     "temperature": args.temperature,
                     "points_skip": args.points_skip,
                     "minpts_bin": args.minpts_bin,
                     "gap_bins": args.gap_bins,
                     "qfep_exec": args.qfep_exec,
                     "en_list_fn": QScfg.get("files", "en_list_fn"),
                     "gas_const": QScfg.get("mapping", "gas_const")}

    # automap with only the first replica (when we have 3 or more)
    # to get a better init guess quickly
    if not args.nosingle and len(mapdirs) > 2:
        print "\nInitial fit, using only the first folder (disable this "\
                "with --nosingle)."
        # create QMapper instance with all arguments
        qmapper_parms["mapdirs"] = mapdirs[:1]
        qmapper_single = QMapper(**qmapper_parms)
        try:
            qmapper_single.fit_to_reference(args.ref_dga, args.ref_dg0,
                                            step_size=args.step_size,
                                            threshold=args.threshold,
                                            max_iterations=1)
                                            
        except QMapperError:
            print "...failed, will try with all dirs anyhow..."
        except KeyboardInterrupt:
            qmapper_single.kill_event.set()
            raise
        else:
            qmapper_parms.update({"hij": qmapper_single.parms["hij"],
                                  "alpha": qmapper_single.parms["alpha"]})

        print "\nSwitching to all directories..."

    qmapper_parms["mapdirs"] = mapdirs
    qmapper = QMapper(**qmapper_parms)

    try:
        rcode = qmapper.fit_to_reference(args.ref_dga, args.ref_dg0,
                                         step_size=args.step_size,
                                         threshold=args.threshold,
                                         max_iterations=args.max_iterations)
    except QMapperError as error_msg:
        print "\nMassive fail:\n{}\n".format(error_msg)
        sys.exit(1)
    except KeyboardInterrupt:
        qmapper.kill_event.set()
        raise

    if not rcode:
        print "Did not converge. Try changing the step (--step), increasing "\
              "number of iterations (--iter) or lowering the treshold "\
              "(--threshold)\n"


    else:
        print """

Well done! Use this on your non-reference simulations:
{}

""".format(qmapper.input_parms_str)

        # write out the inputs and outputs from the last step
        qfep_inp_fn = QScfg.get("files", "qfep_inp")
        qfep_out_fn = QScfg.get("files", "qfep_out")
        for mapdir, (qfep_inp_str, qfep_out_str) in qmapper.mapped.iteritems():
            qfep_inp = os.path.join(mapdir, qfep_inp_fn)
            qfep_out = os.path.join(mapdir, qfep_out_fn)
            open(qfep_inp, "w").write(qfep_inp_str)
            open(qfep_out, "w").write(qfep_out_str)

        # analyse the outputs
        output_files = [os.path.join(md, qfep_out_fn) for md in qmapper.mapped]
        qafs = QAnalyseFeps(output_files)
        fails = "\n".join(["{}: {}".format(qfo, err) for qfo, err in
                                            qafs.failed.iteritems()])

        outstr = """
{mapper_details}
Analysis Stats:
{analysis_stats}
Analysis Fails:
{analysis_fails}
""".format(mapper_details=qmapper.details, analysis_stats=qafs.stats_str,
           analysis_fails=fails or "None")

        if fails or qmapper.failed:
            print """
WARNING! Some dirs failed to map/analyse! Look at the log!

"""

        print "Writting out the logfile..."
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
