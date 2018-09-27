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
import six.moves.configparser
from six.moves import input
#from io import open

import argparse
import sys
import os

# check for python version
if sys.version_info < (2, 7):
    print("Python 2.7+ required. Detected version {0}".format(sys.version))
    sys.exit(1)

try:
    from Qpyl.common import __version__
except ImportError:
    print("Can't import Qpyl module, your installation is messed up...")
    sys.exit(1)


# this should change when the config file format changes
# to indicate an incompatible config file
CFG_VERSION = 10


# absolute config path
try:
    _CFG_FILE = os.path.abspath(os.path.join(os.environ["QTOOLS_HOME"],
                                             "qscripts.cfg"))
except KeyError:
    print("Your QTools instalation is messed up, did you source "\
          "'qscripts_init.sh'?")
    sys.exit(1)

# absolute path of default config file
_CFG_FILE_DEFAULT = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 "qscripts.cfg.default"))



class _QScriptsConfig(object):
    def __init__(self, cfgfile):

        self.cfgfile = cfgfile
        if not os.path.lexists(self.cfgfile):
            print("Configuration file '{}' not found. Please run "\
                  "'qscripts_config.py'.".format(self.cfgfile))
            sys.exit(1)

        args = {"inline_comment_prefixes": (';',)}
        try:
            self.config = six.moves.configparser.ConfigParser(**args)
        except TypeError: # py27
            self.config = six.moves.configparser.ConfigParser()


        try:
            self.config.read(self.cfgfile)
        except six.moves.configparser.ParsingError:
            print("Configuration file '{}' could not be read. Fix/remove it "\
                  "and run 'qscripts_config.py'.".format(self.cfgfile))
            sys.exit(1)

        # check version of config file
        version = self.config.getint("other", "cfgversion")
        if version != CFG_VERSION:
            print("Your configuration file '{}' is outdated. Please remove "\
                  "it and run 'qscripts_config.py'.".format(self.cfgfile))
            sys.exit(1)

    def get(self, section, option):
        try:
            return self.config.get(section, option)
        except (six.moves.configparser.NoSectionError, six.moves.configparser.NoOptionError) as e:
            print("Somehow your configuration file '{}' got messed up. "\
                  "Please fix/remove it and run 'qscripts_config.py'.\n"\
                  "Details: {}".format(self.cfgfile, e))
            sys.exit(1)

    def set(self, section, option, value):
        self.config.set(section, option, value)



def get_exec_path(name):
    paths = []
    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        if not os.path.lexists(path):
            continue
        for ex in os.listdir(path):
            if name in ex:
                ex = os.path.join(path, ex)
                if not ex in paths:
                    paths.append(ex)

    if paths:
        print("These '{}' executables were found in your PATH. "\
              "Choose the correct one (select number).".format(name))
        for i, path in enumerate(paths):
            print("      [{}] {}".format(i, path))
        path = ""
        while not path:
            try:
                inp = input("? ")
                i = int(inp)
                path = paths[i]
            except (ValueError, IndexError):
                pass
        return path
    return ""



def main():
    parser = argparse.ArgumentParser(description="""
    Creates a config file for the qtools command-line interface (CLI).
    ($QTOOLS_HOME/qtools.cfg)
    It searches for Q binaries in the PATH variable and prompts the
    user to choose the correct ones.
    """, add_help=False)
    optarg = parser.add_argument_group("Optional")
    optarg.add_argument("--no_bin", dest="no_bin", action="store_true",
                        help="Don't look for Q binaries.", default=False)
    optarg.add_argument("-h", "--help", action="help", help="show this "
                        "help message and exit")

    args = parser.parse_args()

    if os.path.lexists(_CFG_FILE):
        print("Configuration file '{}' exists. Please remove it before "\
              "running this script.".format(_CFG_FILE))
        sys.exit(1)

    QScriptsConfig = _QScriptsConfig(_CFG_FILE_DEFAULT)
    print("Creating a new configuration file...")

    if not args.no_bin:
        print("\nLooking for Qfep6 executables in PATH...")
        qfep_path = get_exec_path("Qfep6")
        if not qfep_path:
            print("Not found, looking for qfep (v5) instead...")
            qfep_path = get_exec_path("qfep")
        if not qfep_path:
            print("No Qfep executable was found in your PATH. "\
                  "Please set the path manually in the config file.")

        QScriptsConfig.set("qexec", "qfep", qfep_path)

        print("Looking for Qcalc6 executables in PATH...")
        qcalc_path = get_exec_path("Qcalc6")
        if not qcalc_path:
            print("Not found, looking for qcalc (v5) instead...")
            qcalc_path = get_exec_path("qcalc")
        if not qcalc_path:
            print("No Qcalc executable was found in your PATH. "\
                "Please set the path manually in the config file.")

        QScriptsConfig.set("qexec", "qcalc", qcalc_path)

    QScriptsConfig.config.write(open(_CFG_FILE, "w+"))
    print("\nThe following config file was created with some default "\
          "values:\n{}".format(_CFG_FILE))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nCtrl-C detected. Quitting...")
        sys.exit(1)
else:  # imported
    # This variable is imported from other modules.
    QScriptsConfig = _QScriptsConfig(_CFG_FILE)

