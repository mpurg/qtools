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

from __future__ import absolute_import, print_function
from __future__ import division, unicode_literals

from qscripts_config import __version__, QScriptsConfig as QScfg

import os
import sys
import time
import argparse
import logging

from Qpyl.core.qlibrary import QLib, QLibError
from Qpyl.common import backup_file, init_logger, get_version_full

logger = init_logger('Qpyl')

parser = argparse.ArgumentParser(description="""
Script for rescaling charges in OPLS-AA. The given library must have correctly
defined charge groups, the charges will then be rescaled such that the charge
groups have integer charge.
""", add_help=False)
reqarg = parser.add_argument_group("Required")
reqarg.add_argument("lib_file",
                    help="Q oplsaa library file (single residue only)")
optarg = parser.add_argument_group("Optional")
optarg.add_argument("-t", dest="threshold", type=float, default=0.3,
                    help="Threshold of difference between net charge and "
                         "nearest integer charge, above which the script "
                         "will fail. Default=0.3")
optarg.add_argument("--ignore_errors", action="store_true", default=False,
                    help="Use if nothing else works")

optarg.add_argument("-v", "--version", action="version",
                    version=get_version_full())
optarg.add_argument("-h", "--help", action="help", help="show this "
                    "help message and exit")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

if not os.path.lexists(args.lib_file):
    print("FATAL! File %s doesn't exist." % args.lib_file)
    sys.exit(1)


# load the library
try:
    qlib = QLib("oplsaa", ignore_errors=args.ignore_errors)
    qlib.read_lib(args.lib_file)
except QLibError as e:
    print("FATAL! Problem with library: {}".format(str(e)))
    sys.exit(1)

if len(qlib.residue_dict) > 1:
    print("FATAL! Please supply a library with just one residue entry.")
    sys.exit(1)

residue_lib = list(qlib.residue_dict.values())[0]

if not residue_lib.charge_groups:
    print("No charge groups found, using all atoms in residue '{}'"\
          "".format(residue_lib.name))
    ch_groups = [ [a.name for a in residue_lib.atoms] ]
         
else:
    print("Rescaling charge groups in '{}'...".format(residue_lib.name))
    ch_groups = residue_lib.charge_groups

# iterate over the charge groups and rescale the charges
for ch_group in ch_groups:
    try:
        residue_lib.rescale(ch_group, args.threshold)
    except QLibError as e:
        print("Library error: {}".format(str(e)))
        sys.exit(1)
        
print("Checking the library...")
try:
    qlib.check_valid()
except QLibError as e:
    print("Library error: {}".format(str(e)))
    sys.exit(1)
print("Done")
print() 

backup = backup_file(args.lib_file)
if backup:
    print("Backed up '{}' to '{}'".format(args.lib_file, backup))
outstring = """# Generated with {}, version {}
# Date: {}
#
{}
""".format(os.path.basename(__file__), __version__,
           time.ctime(), qlib.get_string())
open(args.lib_file, 'w').write(outstring)
print("Successfully modified the library.")


