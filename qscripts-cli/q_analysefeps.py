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

from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import time
import os
import argparse
import logging

from Qpyl.common import backup_file, SpecialFormatter
from Qpyl.qanalysis import QAnalyseFeps
from Qpyl import plotdata


def main():
    logger = logging.getLogger('Qpyl')
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(SpecialFormatter())
    logger.addHandler(handler)

    parser = argparse.ArgumentParser(description="""
Tool for analysing QFep outputs - extracting FEP results, activation and reaction free
energies, calculating LRA contributions, calculating statistics over all
outputs, and exporting all the data into JSON format. Should be used after
every mapping.
    """, add_help=False)
    reqarg = parser.add_argument_group("Required")
    reqarg.add_argument("fepdirs", nargs="+",
                        help="Directories to scan for qfep output.")

    optarg = parser.add_argument_group("Optional")
    def_lra_l = QScfg.get("analysis", "lambdas_state1").split(",")

    optarg.add_argument("--lra_l", dest="lra_l", nargs=2,
                        metavar=("l1", "l2"),
                        default=def_lra_l,
                        help="Specify lambdas (state 1) at which LRA and "
                             "REORG are calculated. Default is '{}'."
                             "".format(str(def_lra_l)))


    optarg.add_argument("--qfep_out", dest="qfep_out",
                        help="Qfep output filename (default='{}')"
                             "".format(QScfg.get("files", "qfep_out")),
                        default=QScfg.get("files", "qfep_out"))

    optarg.add_argument("--out", dest="output_fn",
                        help="Output filename (default='{}')"
                             "".format(QScfg.get("files", "analysefeps_log")),
                        default=QScfg.get("files", "analysefeps_log"))

    optarg.add_argument("--plots_out", dest="plots_out",
                        help="Output filename for plot data (default='{}')"
                             "".format(QScfg.get("files", "analysefeps_plots")),
                        default=QScfg.get("files", "analysefeps_plots"))

    optarg.add_argument("--po_qcp", dest="plots_out_qcp", default=None,
                        help="Output filename for QCP plot data. By default "
                             "this data is not written out.")

    optarg.add_argument("--po_qcpm", dest="plots_out_qcp_mass", default=None,
                        help="Output filename for QCP_mass plot data. "
                             "By default this data is not written out.")


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    lra_l = []
    for lamb in args.lra_l:
        try:
            lamb = float(lamb)
            if lamb < 0 or lamb > 1:
                raise ValueError
        except ValueError:
            print "FATAL! LRA lambdas make no sense. 0<lambda<1 please."
            sys.exit(1)
        lra_l.append(lamb)

    qos = [os.path.join(md, args.qfep_out) for md in sorted(args.fepdirs)]
    qaf = QAnalyseFeps(qos, lra_lambdas=lra_l)

    stats, fails = [], []

    stats.append(qaf.stats_str)
    for sub_calc_key, sub_calc in sorted(qaf.sub_calcs.iteritems()):
        stats.append(sub_calc.stats_str)

    if qaf.failed:
        fails.append("Failed to parse:")
    for failed_path, failed_msg in sorted(qaf.failed.iteritems()):
        relp = os.path.relpath(failed_path)
        fails.append("-> {}: {}".format(relp, failed_msg))

    if qaf.failed_dg:
        fails.append("Failed to produce dGa/dG0:")
    for failed_path, failed_msg in sorted(qaf.failed_dg.iteritems()):
        relp = os.path.relpath(failed_path)
        fails.append("-> {}: {}".format(relp, failed_msg))

    stats = "\n".join(stats)
    fails = "\n".join(fails) or None

    summary = """
----------------------------------- SUMMARY -----------------------------------
# Analysed with: QTools/q_analysefeps.py ({version})
# Work dir: {cwd}
# Date: {date}
# CMDline: {cmdline}


----- Statistics -----

{stats}

------- Fails --------

{fails}
-------------------------------------------------------------------------------
""".format(version=__version__, date=time.ctime(), cwd=os.getcwd(),
           stats=stats, fails=fails, cmdline=" ".join(sys.argv))

    print summary

    if not qaf.qfos:
        print "\nFATAL! None of the outputs could be parsed!"
        print "Are you running an ancient Q version? Then don't..."
        print "If not, report a bug."
        sys.exit(1)


    # save some useful data
    output_string = """
------------------------------- Free energies ---------------------------------
{}
{}
""".format(qaf.dg_all, summary)

    fn_out = args.output_fn
    backup = backup_file(fn_out)
    if backup:
        print "# Backed up '{}' to '{}'".format(fn_out, backup)
    open(fn_out, "w").write(output_string)
    print "Wrote '{}'...".format(fn_out)

    # convert plots to json and write them out
    fn_out = args.plots_out
    plots = qaf.plotdata
    jsonenc = plotdata.PlotDataJSONEncoder(indent=2)
    backup = backup_file(fn_out)
    if backup:
        print "# Backed up '{}' to '{}'".format(fn_out, backup)
    open(fn_out, 'w').write(jsonenc.encode(plots))
    print "Wrote '{}'... (q_plot.py is your "\
          "friend)".format(fn_out)

    if args.plots_out_qcp:
        try:
            plots = qaf.sub_calcs["QCP"].plotdata
        except KeyError:
            print "No QCP data found, can't write it out..."
        else:
            fn_out = args.plots_out_qcp
            jsonenc = plotdata.PlotDataJSONEncoder(indent=2)
            backup = backup_file(fn_out)
            if backup:
                print "# Backed up '{}' to '{}'".format(fn_out, backup)
            open(fn_out, 'w').write(jsonenc.encode(plots))
            print "Wrote '{}'... (q_plot.py is your "\
                "friend)".format(fn_out)

    if args.plots_out_qcp_mass:
        try:
            plots = qaf.sub_calcs["QCP_mass"].plotdata
        except KeyError:
            print "No QCP_mass data found, can't write it out..."
        else:
            fn_out = args.plots_out_qcp_mass
            jsonenc = plotdata.PlotDataJSONEncoder(indent=2)
            backup = backup_file(fn_out)
            if backup:
                print "# Backed up '{}' to '{}'".format(fn_out, backup)
            open(fn_out, 'w').write(jsonenc.encode(plots))
            print "Wrote '{}'... (q_plot.py is your "\
                "friend)".format(fn_out)



if __name__ == "__main__":
    main()


