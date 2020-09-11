#!/usr/bin/env python
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
#
# Converts FFLD output to library and parameter files for Q.
# Two arguments are required:
# - the FFLD file created with 'ffld_server' (SCHRODINGER's Maestro)
# - PDB structure from which the FFLD_OUTPUT was created (QM optimized)
#
# Creates three files: .lib, .prm and .prm.chk
#
# ffld_server usage example:
#   $SCHRODINGER/utilities/ffld_server -ipdb paraoxon.pdb -print_parameters -version 14 > paraoxon.ffld11
#
#

from __future__ import absolute_import, print_function
from __future__ import division, unicode_literals
import six
from six.moves import range

from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import os
import time
import argparse
import logging

from Qpyl.core.qlibrary import QLib, QLibError
from Qpyl.core.qparameter import QPrm, QPrmError
from Qpyl.core.qpotential import torsion_energy
from Qpyl.core.qstructure import QStruct, QStructError
from Qpyl.core.qtopology import QTopology, QTopologyError
from Qpyl.common import backup_file, init_logger, get_version_full

logger = init_logger('Qpyl')

parser = argparse.ArgumentParser(description="""
Command-line tool for converting OPLS-AA force-field parameters from FFLD
(Schrodinger MacroModel) to Q format.

Additionally, it checks the quality of the parameters by calculating all the bonding
energies in the given structure.

Also, it complains about stuff like duplicate or overwritten parameters
and non-integer residue charges.
""", add_help=False)
reqarg = parser.add_argument_group("Required")
reqarg.add_argument("ffld_output", help = "ffld_server output")
reqarg.add_argument("pdb", help="PDB structure file WHICH WAS USED TO CREATE "
                                "THE FFLD_OUTPUT (used to copy "
                                "atom names and check parameter energies)")
optarg = parser.add_argument_group("Optional")
optarg.add_argument("-o", dest="output_basename",
                    help="Basename for output files (.lib, .prm and .prm.chk)."
                         " Default is 'XXX'.", default="XXX")
optarg.add_argument("--ignore_errors", action="store_true", default=False,
                    help="Use in case you have double parameter definitions, "
                         "non-integer residue charge (from MCPB.py for "
                         "instance), or other weird stuff, but PLEASE don't "
                         "ignore the output message and triple check your "
                         "outputs.")
optarg.add_argument("-v", "--version", action="version",
                    version=get_version_full())
optarg.add_argument("-h", "--help", action="help", help="show this "
                    "help message and exit")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

for k, v in six.iteritems(vars(args)):
    if k in ['ffld_output', 'pdb'] and not os.path.lexists(v):
        print("FATAL! File '{}' doesn't exist.".format(v))
        sys.exit(1)


#
# create QLib, QPrm, QStruct and QTopology objects
#

qlib = QLib("oplsaa", ignore_errors=args.ignore_errors)
qprm = QPrm("oplsaa", ignore_errors=args.ignore_errors)

try:
    qstruct = QStruct(args.pdb, "pdb", ignore_errors=args.ignore_errors)
except QStructError as err:
    print("FATAL! Problem with pdb: {}".format(err))
    sys.exit(1)

try:
    qprm.read_ffld(args.ffld_output, qstruct)
except QPrmError as err:
    print("FATAL! Problem with ffld file: {}".format(err))
    sys.exit(1)


try:
    qlib.read_ffld(args.ffld_output, qstruct)
except QLibError as err:
    print("FATAL! Problem with ffld file: {}".format(err))
    sys.exit(1)

try:
    qtop = QTopology(qlib, qprm, qstruct)
except QTopologyError as err:
    print("FATAL! Problem building the topology: {}".format(err))
    sys.exit(1)

#
# get total and max energies and number of parameters
#

data, total_e, max_e, nprm = {}, {}, {}, {}

for bati_name in ["bonds", "angles", "torsions", "impropers"]:
    total_e[bati_name] = 0
    max_e[bati_name] = 0
    nprm[bati_name] = 0
    data[bati_name] = []
    for bati in getattr(qtop, bati_name):
        if bati_name == "bonds":
            e0 = 0
            v0 = bati.prm.r0
        elif bati_name == "angles":
            e0 = 0
            v0 = bati.prm.theta0
        elif bati_name == "torsions":
            # get the torsion minimum
            energy_profile = []
            for phi in range(0, 181, 30):
                energy = 0
                for fc, multiplicity, phase, npaths in bati.prm.get_prms():
                    energy += torsion_energy(phi, fc, multiplicity,
                                             npaths, phase)
                energy_profile.append((energy, phi))
            energy_profile.sort()
            e0, v0 = min(energy_profile)
        elif bati_name == "improper":
            e0 = 0
            v0 = bati.prm.phi0

        e, v = bati.calc()
        de = e - e0

        data[bati_name].append((bati, de, v, v0))
        total_e[bati_name] += de
        if de > max_e[bati_name]:
            max_e[bati_name] = de
        nprm[bati_name] += 1

nprm["diff_torsions"] = sum([len(list(tor.prm.get_prms())) for tor in qtop.torsions])

# sort the data according to energies
all_prms = data["bonds"] + data["angles"] + \
           data["torsions"] + data["impropers"]
all_prms.sort(key=lambda x: x[1]) # sort by energy
pchk = []
for p, de, v, v0 in all_prms:
    prmtype = "-".join(p.prm.prm_id.split())
    pchk.append("{!r:<40}  value:{:>8.2f}   eq.value:{:>8.2f}    "
                "dE:{:>6.2f}    # Prm: {}".format(p, v, v0, de, prmtype))


#
# output some information
#
print("""
Details about the system:

Bonds: {n[bonds]}
Angles: {n[angles]}
Torsions: {n[torsions]}
Q Torsions (diff parms): {n[diff_torsions]}
Impropers: {n[impropers]}
Bond energy: total {t[bonds]:.2f}, max {m[bonds]:.2f}
Angle energy: total {t[angles]:.2f}, max {m[angles]:.2f}
Torsion energy: total {t[torsions]:.2f}, max {m[torsions]:.2f}
Improper energy: total {t[impropers]:.2f}, max {m[impropers]:.2f}

""".format(n=nprm, t=total_e, m=max_e))

#
# write the files
#
libfn = args.output_basename + ".lib"
prmfn = args.output_basename + ".prm"
prmchkfn = args.output_basename + ".prm.chk"

print("Writing the library file: {}".format(libfn))
backup = backup_file(libfn)
if backup:
    print("# Backed up '{}' to '{}'".format(libfn, backup))
outstring = """# Generated with {}, version {}
# Date: {}
#
{}
""".format(os.path.basename(__file__), __version__,
           time.ctime(), qlib.get_string())
open(libfn, "w").write(outstring)

print("Writing the parameter file: {}".format(prmfn))
backup = backup_file(prmfn)
if backup:
    print("# Backed up '{}' to '{}'".format(prmfn, backup))
outstring = """# Generated with {}, version {}
# Date: {}
#
{}
""".format(os.path.basename(__file__), __version__,
           time.ctime(), qprm.get_string())
open(prmfn, "w").write(outstring)

print("Writing the parameter check file: {}".format(prmchkfn))
backup = backup_file(prmchkfn)
if backup:
    print("# Backed up '{}' to '{}'".format(prmchkfn, backup))
open(prmchkfn, "w").write("\n".join(pchk))
