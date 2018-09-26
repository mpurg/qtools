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
# Takes as arguments:
#    - pdb structure
#    - file containing residue indexes that should be charged
#    (each in new line or space separated)
#
# Renames the ARG,GLU,ASP,LYS that should be neutral to ARN,GLH,ASH and LYN
# and vice versa
#

from __future__ import absolute_import
from __future__ import print_function
from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import os
import argparse
from Qpyl.common import backup_file, get_version_full

parser = argparse.ArgumentParser(description="""
Simple script for ionizing residues. It basically changes all the specified
residues to ARG,GLU,ASP,LYS, while the rest are renamed to ARN,GLH,ASH,LYN.
""", add_help=False)
reqarg = parser.add_argument_group("Required")
reqarg.add_argument("pdb", help="PDB structure file")
reqarg.add_argument("resids", help="Text file with space or newline separated "
                                   "indexes of residues that should be "
                                   "IONIZED. All others will be set to "
                                   "their neutral form.")
reqarg.add_argument("outfn", help="Output filename")

optarg = parser.add_argument_group("Optional")
optarg.add_argument("-v", "--version", action="version",
                    version=get_version_full())
optarg.add_argument("-h", "--help", action="help", help="show this "
                    "help message and exit")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

if not os.path.lexists(args.pdb):
    print("FATAL! File '{}' doesn't exist.".format(args.pdb))
    sys.exit(1)

if not os.path.lexists(args.resids):
    print("FATAL! File '{}' doesn't exist.".format(args.resids))
    sys.exit(1)

chres = ("ARG", "GLU", "ASP", "LYS") # charged names
nres = ("ARN", "GLH", "ASH", "LYN") # neutral names

pdb_lines = open(args.pdb, 'r').readlines()
charged_residues = open(args.resids, 'r').read().split()

new_pdb = ""

for line in pdb_lines:
    ri = line[22:26].strip()
    rn = line[17:20]
    new_line=""

    if ri in charged_residues:
        if rn in nres:
            new_rn = chres[ nres.index(rn) ]
            print("Changing: %s%s to %s%s" % (rn,ri, new_rn, ri))
            new_line=line.replace(rn, new_rn)
        else:
            new_line=line
    else:
        if rn in chres:
            new_rn = nres[ chres.index(rn) ]
            new_line=line.replace(rn, new_rn)
        else:
            new_line = line
    new_pdb += new_line


backup = backup_file(args.outfn)
if backup:
    print("Backed up '{}' to '{}'".format(args.outfn, backup))
open(args.outfn, 'w').write(new_pdb)
print("Wrote " + args.outfn)



