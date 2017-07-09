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
# changes the placeholders inside template files ("$593.CA$") to pdb indexes  ("1123")
# takes in three arguments: PDB (after qprep), file containing the placeholders (q_makefep.py generated FEP file, input templates), output filename
# extra keyword that can be used instead of the placeholder is 'LAST.ID' (no explanation needed)


from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import os
import argparse

from Qpyl.core.qstructure import QStruct, QStructError
from Qpyl.common import backup_file, init_logger, get_version_full

logger = init_logger('Qpyl')

parser = argparse.ArgumentParser(description="""
Command-line tool for converting atom placeholders to indexes. The
placeholders have the following format: $RESID.ATOM_NAME$
""", add_help=False)
reqarg = parser.add_argument_group("Required")
reqarg.add_argument("pdb",
                    help="pdb structure file (created with qprep)")
reqarg.add_argument("inp",
                    help="input/fep file containing the placeholders")
reqarg.add_argument("out",
                    help="output filename")
optarg = parser.add_argument_group("Optional")
optarg.add_argument("--ignore_errors", action="store_true", default=False,
                    help="Don't use this unless it's an emergency")
optarg.add_argument("-v", "--version", action="version",
                    version=get_version_full())
optarg.add_argument("-h", "--help", action="help", help="show this "
                    "help message and exit")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

for k, v in vars(args).iteritems():
    if k in ["inp", "pdb"] and not os.path.lexists(v):
        print "FATAL! File '{}' doesn't exist.".format(v)
        sys.exit(1)

inpstr = open(args.inp, "r").read()

try:
    qstruct = QStruct(args.pdb, "pdb", ignore_errors=args.ignore_errors)
    outstring = qstruct.convert_placeholders(inpstr)
except QStructError as e:
    print "FATAL! Exception raised: {}".format(str(e))
    sys.exit(1)

backup = backup_file(args.out)
if backup:
    print "Backed up '{}' to '{}'".format(args.out, backup)
open(args.out, "w").write(outstring)
print "Created file '{}'...".format(args.out)

