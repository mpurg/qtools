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

from __future__ import absolute_import, print_function
from __future__ import division, unicode_literals
import six

from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import os
import argparse
import logging

from Qpyl.qmakefep import make_fep, QMakeFepError
from Qpyl.common import backup_file, init_logger, get_version_full


if __name__ == "__main__":
    logger = init_logger('Qpyl')

    parser = argparse.ArgumentParser(description="""
    Generates a FEP file for EVB simulations in Q. The changes in atom types,
    charges, bonding terms are determined automatically from the 'qmap' file,
    structure, libraries and parameter files. Parameters that cannot be
    determined (soft_pair C, Morse D,a,r, Hij) are replaced with <FIX>.
    """, add_help=False)

    reqarg = parser.add_argument_group("Required")
    reqarg.add_argument("-m", nargs=1, dest="qmap", required=True,
                        help="QMAP file. Used to map Q atoms between states. "
                             "Examples can be found in the test folder.")

    reqarg.add_argument("-s", nargs=1, dest="pdb", required=True,
                        help="PDB structure file (best to use the one "
                             "created by qprep when making the topology).")

    reqarg.add_argument("-f", dest="forcefield", required=True,
                        help="Forcefield type. Supported options are 'amber'"
                             " and 'oplsaa'.")

    reqarg.add_argument("-p", nargs="+", dest="prms",
                        help="Force field parameter files.", required=True)

    reqarg.add_argument("-l", nargs="+", dest="libs",
                        help="Force field library files.", required=True)

    optarg = parser.add_argument_group("Optional")
    optarg.add_argument("-o", nargs=1, dest="outfile",
                        default=["generated.fep"],
                        help="Output file (default=generated.fep).")

    optarg.add_argument("--ignore_errors", action="store_true", default=False,
                        help="Use in case you have double parameter "
                             "definitions, non-integer residue charge "
                             "(from MCPB.py for instance), or other weird "
                             "stuff, but PLEASE don't ignore the output "
                             "messages and PLEASE triple check your outputs.")
    optarg.add_argument("-v", "--version", action="version",
                        version=get_version_full())
    optarg.add_argument("-h", "--help", action="help", help="show this "
                        "help message and exit")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    for k, v in six.iteritems(vars(args)):
        if k in ["pdb", "parms", "libs", "qmap"]:
            for fn in v:
                if not os.path.lexists(fn):
                    print("File '{}' doesn't exist.".format(fn))
                    sys.exit(1)

    try:
        fepstr = make_fep(args.qmap[0], args.pdb[0], args.forcefield,
                          args.prms, args.libs,
                          ignore_errors=args.ignore_errors)
    except QMakeFepError as err:
        print("\nFATAL!\n{}\n".format(err))
        sys.exit(1)

    output_file = args.outfile[0]
    backup = backup_file(output_file)
    if backup:
        print("\n# Backed up '{}' to '{}'".format(output_file, backup))
    open(output_file, 'w').write(fepstr)
    print("Wrote '{}'... ".format(output_file))
    
