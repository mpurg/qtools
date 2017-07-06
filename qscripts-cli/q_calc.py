#!/usr/bin/env python2

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
import logging
import argparse

from Qpyl.qgroupcontrib import QGroupContrib, QGroupContribError
from Qpyl import plotdata
from Qpyl.common import backup_file, init_logger

def gc(args):

    if not os.path.lexists(args.pdb):
        print "This file went missing: {}".format(args.pdb)
        sys.exit(1)

    if args.qmaskfile:
        try:
            qmask = open(args.qmaskfile, "r").read().split()
            if not qmask: raise IOError
        except IOError:
            print "Can't read '{}' or file empty".format(args.qmaskfile)
            sys.exit(1)
    else:
        qmask = None

    lambdas = []
    for lamb in args.lra_l:
        try:
            lamb = float(lamb)
            if lamb < 0 or lamb > 1:
                raise ValueError
        except ValueError:
            print "FATAL! Lambda values make no sense. 0<lambda<1 please."
            sys.exit(1)
        lambdas.append((lamb, 1-lamb))

    calcdirs = args.dirs
    if not calcdirs:
        lsdir = os.listdir(os.getcwd())
        calcdirs = [f for f in lsdir if os.path.isdir(f)]
    if not calcdirs:
        calcdirs = [os.getcwd(),]
        print "No subdirectories. Calculating in current directory only.\n"
    else:
        print "Will use these directories for calculating GCs (use --dirs to "\
              "change this): {}\n".format(", ".join(calcdirs))

    qgc = QGroupContrib(args.qcalc_exec, calcdirs, args.pdb,
                        QScfg.get("files", "en_list_fn"),
                        lambdas[0], lambdas[1],
                        args.resid_first, args.resid_last,
                        args.scale_ionized,
                        args.nthreads,
                        qmask)

    try:
        qgc.calcall()
    except QGroupContribError as error_msg:
        print "\nMassive fail:\n{}\n".format(error_msg)
        sys.exit(1)
    except KeyboardInterrupt:
        qgc.kill_event.set()
        raise

    # writeout QCalc inputs and outputs
    if args.writeout:
        for calcdir, (qcinps, qcouts) in qgc._qcalc_io.iteritems():
            for i, qci in enumerate(qcinps):
                fn = os.path.join(calcdir, "q_calc.gc.{}.inp".format(i+1))
                try:
                    open(fn, 'w').write(qci)
                except OSError as err:
                    print "Error when writing to {}: {}".format(fn, err)
                else:
                    print "Wrote {}".format(fn)
            for i, qco in enumerate(qcouts):
                fn = os.path.join(calcdir, "q_calc.gc.{}.out".format(i+1))
                try:
                    open(fn, 'w').write(qco)
                except OSError as err:
                    print "Error when writing to {}: {}".format(fn, err)
                else:
                    print "Wrote {}".format(fn)



    # write out details and top 20 El GCs to stdout and outfile
    if not qgc.gcs:
        top_gcs = "None, all directories failed..."
    else:
        top_rows = sorted(qgc.gcs_stats.get_rows(),
                          key=lambda x: -abs(x[5]))[:10]

        out_l = ["{:<15} {:>15} {:>15}".format("Residue",
                                               "GC_El_mean",
                                               "GC_El_std")]

        for rid, rn, _, _, _, el, els in top_rows:
            tmp = "{}_{}".format(rn.capitalize(), rid)
            tmp2 = "{:<15} {:15.2f} {:15.2f}"\
                   "".format(tmp, el, els)
            out_l.append(tmp2)
        top_gcs = "\n".join(out_l)

    outstr = """
{gc_details}
Top group contributions:
{top_gcs}
""".format(gc_details=qgc.details, top_gcs=top_gcs)

    print outstr
    fn_out = args.output_fn
    backup = backup_file(fn_out)
    if backup:
        print "# Backed up '{}' to '{}'".format(fn_out, backup)
    open(fn_out, "w").write(outstr)
    print "Wrote '{}'...".format(fn_out)

    # convert plots to json and write them out
    fn_out = args.plots_out
    plots = qgc.plotdata
    jsonenc = plotdata.PlotDataJSONEncoder(indent=2)
    backup = backup_file(fn_out)
    if backup:
        print "# Backed up '{}' to '{}'".format(fn_out, backup)
    open(fn_out, 'w').write(jsonenc.encode(plots))
    print "Wrote '{}'... (q_plot.py is your "\
          "friend)".format(fn_out)

    # writeout the pdbgc if requested
    if args.pdbgc_out:
        backup = backup_file(args.pdbgc_out)
        if backup:
            print "# Backed up '{}' to '{}'".format(args.pdbgc_out, backup)
        open(args.pdbgc_out, 'w').write(qgc.get_pdbgc())
        print "Wrote '{}'... (use Pymol/Chimera/VMD and color by B-factor)"\
            "".format(args.pdbgc_out)


def main():
    logger = init_logger('Qpyl')

    parser = argparse.ArgumentParser(description="""
    A friendly command-line interface for calculating distances, angles, rmsds,
    group contributions, etc. with Qcalc. To get help on a sub-command, just run
    it without additional arguments.""")

    subp = parser.add_subparsers(title="subcommands", dest="command")

    subps = {}
    subps["gc"] = subp.add_parser("gc", help="Group contribution calculation.",
                                  description="""
                                  Calculate group contributions - nonbonded
                                  Linear Response Approximation energies
                                  between (protein) residues and the 
                                  reactive Q region.""", add_help=False)

    gc_reqarg = subps["gc"].add_argument_group("Required")
    gc_reqarg.add_argument("pdb", help="PDB structure created with Qprep.")

    gc_reqarg.add_argument("resid_first", type=int,
                           help="Indexes of first residue to be "
                                "included in the calculation.")

    gc_reqarg.add_argument("resid_last", type=int,
                             help="Index of last residue to be "
                                  "included in the calculation. ")

    gc_optarg = subps["gc"].add_argument_group("Optional")
    gc_optarg.add_argument("--iscale", dest="scale_ionized",
                           type=float, default=QScfg.get("gc",
                                                         "scale_ionized"),
                           help="Scale down electrostatic interactions of "
                                "ionized residues (ASP, GLU, HIP, LYS, ARG) "
                                "(see doi:10.1021/jp962478o). Default is {}"
                                "".format(QScfg.get("gc", "scale_ionized")))

    gc_optarg.add_argument("--qmask", dest="qmaskfile", default=None,
                           help="File containing the Q-atom mask (line "
                                "or space separated atom indexes). By default, "
                                "this is extracted from the [atoms] section "
                                "in the FEP file.")

    def_lra_l = QScfg.get("gc", "lambdas_state1").split(",")
    gc_optarg.add_argument("--lra_l", dest="lra_l", nargs=2,
                           metavar=("l1", "l2"),
                           default=def_lra_l,
                           help="Specify lambdas (state 1) for GC LRA."
                                "Default is '{}'."
                                "".format(str(def_lra_l)))

    gc_optarg.add_argument("--nt", dest='nthreads', type=int,
                           default=QScfg.get("gc", "nthreads"),
                           help="Number of threads (default = {})"
                           "".format(QScfg.get("gc", "nthreads")))

    gc_optarg.add_argument("--plots_out", dest="plots_out",
                           help="Output filename for plot data (default="
                                "'{}').".format(QScfg.get("files",
                                                          "gc_plots")),
                           default=QScfg.get("files", "gc_plots"))

    gc_optarg.add_argument("--dirs", dest="dirs", nargs="+",
                           help="Directories to use (default is "
                                "all subdirs or current working dir).",
                           default=[])

    gc_optarg.add_argument("--out", dest="output_fn",
                           help="Output filename (default='{}')."
                                "".format(QScfg.get("files", "calcs_log")),
                           default=QScfg.get("files", "calcs_log"))

    gc_optarg.add_argument("--pdbgc", dest="pdbgc_out",
                           help="Output filename of PDB structure file "
                                "with group contributions in place of the "
                                "B-factor. Default=Don't output.",
                           default=None)

    gc_optarg.add_argument("--writeout", action="store_true", default=False,
                           help="Write out QCalc inputs and outputs."
                                "Default=Don't")

    gc_optarg.add_argument("--qcalc_exec", dest="qcalc_exec",
                           default=QScfg.get("qexec", "qcalc"),
                           help="qcalc5 executable path (default={})."
                                "".format(QScfg.get("qexec", "qcalc")))

    gc_optarg.add_argument("-h", "--help", action="help", help="show this "
                           "help message and exit")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif len(sys.argv) == 2 and sys.argv[1] not in ['-h', '--help']:
        try:
            subps[sys.argv[1]].print_help()
        except KeyError:
            print "Subcommand '{}' not available...".format(sys.argv[1])
            print
        sys.exit(1)

    args = parser.parse_args()

    if args.command == "gc":
        gc(args)



if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print "\nCtrl-C detected. Quitting..."
        sys.exit(1)
