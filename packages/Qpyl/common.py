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

"""
This module contains some common classes and functions,
including simple statistical methods and data structures.
"""

from __future__ import absolute_import, division, unicode_literals
from six.moves import zip
from io import open
import math
import sys
import os
import shutil
import logging
import gzip

try:
    import statistics
except ImportError:
    pass

logger = logging.getLogger(__name__)

from Qpyl import __version__

def get_version_full():
    try:
        gitdir = os.path.join(os.environ["QTOOLS_HOME"], ".git")
        head = open(os.path.join(gitdir, "HEAD")).read().split()[1]
        branch = head.split("/")[-1]
        ref = open(os.path.join(gitdir, head)).read().strip()[:8]
    except:
        ref, branch = "Unknown", "Unknown"
    return "Qtools version: {}, git id: {} ({})"\
           "".format(__version__, ref, branch)


class SpecialFormatter(logging.Formatter):
    FORMATS = {logging.DEBUG :"DBG: %(module)s: %(lineno)d: %(message)s",
               logging.WARNING : "\nWARNING: %(message)s\n",
               logging.CRITICAL : "\nCRITICAL: %(message)s\n",
               logging.INFO : "# %(message)s",
               'DEFAULT' : "%(message)s"}

    def __init__(self, *args, **kwargs):
        super(SpecialFormatter, self).__init__(*args, **kwargs)

    def format(self, record):
# a bit diff in py2 vs py3
# https://stackoverflow.com/questions/14844970/modifying-logging-message-format-based-on-message-logging-level-in-python3
        try:
            self._style
            self._style._fmt = self.FORMATS.get(record.levelno,
                                                self.FORMATS['DEFAULT'])
        except AttributeError:
            self._fmt = self.FORMATS.get(record.levelno,
                                         self.FORMATS['DEFAULT'])
        return logging.Formatter.format(self, record)


def init_logger(name,
                level=None,
                handler=None,
                formatter=None):
    """Helper function for initializing the logger.

    Args:
        name (string):  module name, usually root: 'Qpyl'
        level (int, optional):  logging level (DEBUG, INFO, WARNING...), \
                                default is INFO
        handler: (logging.Handler, optional):  default is \
                                               StreamHandler(sys.stdout) \
        formatter: (logging.Formatter, optional): default is \
                                                  SpecialFormatter

    Returns:
        lg (logging.Logger)

    """
    lg = logging.getLogger(name)
    if level == None:
        level = logging.INFO
    lg.setLevel(level)

    if handler == None:
        handler = logging.StreamHandler(sys.stdout)
    if formatter == None:
        handler.setFormatter(SpecialFormatter())

    lg.addHandler(handler)
    return lg


def raise_or_log(message, exception_class, logger, ignore_errors):
    """Method used for raising exceptions or writting them to logger instead

    This way one can bypass certain exceptions like non-integer charge groups,
    that may occur in weird scenarios (Amber bonding metal model).

    Critical level is always used for logging.
    """
    if ignore_errors:
        logger.critical(message)
    else:
        raise exception_class(message)
    

def backup_file(filename):
    """Check if a file exists, make a backup (#filename.1#, #filename.2#...).
    
    Args:
        filename (string):  name of file to backup

    Returns:
        backup_filename (string):  basename of the new filename or empty
                                   string if the file was not found.

    """
    if os.path.lexists(filename):
        di = os.path.dirname(filename)
        fn = os.path.basename(filename)
        backup_filename = fn
        i = 1
        while os.path.lexists(os.path.join(di, backup_filename)):
            backup_filename = "#%s.%d#" % (fn, i)
            i += 1
            if i > 20:
                logger.warning("You have more than 20 backed up files... "
                               "Cleaning time...")
        shutil.copy2(filename, os.path.join(di, backup_filename))
        return backup_filename
    return ""
    

# no need for numpy to do these basic stats
class stats(object):

    @staticmethod
    def mean(vals):
        """Calculate mean.

        Args:
            vals (list of float):  sample values

        Wraps statistics.mean() in Py3+.

        Returns float('nan') on empty array.
        """
        if len(vals) == 0:
            return float('nan')
        try:
            return statistics.mean(vals)
        except NameError:
            return 1.0 / len(vals) * sum(vals)


    @staticmethod
    def stdev(vals):
        """Calculate sample standard deviation.
        
        Args:
            vals (list of float):  sample values

        Wraps statistics.stdev() in Py3+.

        Returns float('nan') when fewer than two values.
        """
        if len(vals) < 2:
            return float('nan')

        try:
            return statistics.stdev(vals)
        except NameError:
            mean = stats.mean(vals)
            variance = [(x - mean)**2 for x in vals]
            return math.sqrt(sum(variance)*1.0/(len(vals)-1))

    @staticmethod
    def sem(vals):
        """Calculates standard error of mean.
        
        Args:
            vals (list of float):  sample values

        Returns float('nan') when fewer than two values.
        """
        if len(vals) < 2:
            return float('nan')

        return stats.stdev(vals) / math.sqrt(len(vals))

    @staticmethod
    def median(vals):
        """Calculate median

        Args:
            vals (list of float):  sample values

        Wraps statistics.median() in Py3+.

        Returns float('nan') on empty array.
        """
        N = len(vals)
        if N == 0:
            return float('nan')
        try:
            return statistics.median(vals)
        except NameError:
            vals = sorted(vals)
            if N % 2 == 0: #even
                return stats.mean((vals[N//2-1], vals[N//2]))
            else: #odd
                return vals[N//2]




class DataContainer(object):
    """
    Contains a two dimensional array of values:

    [ [ row1_column1, row1_column2, row1_column3, ...],
      [ row2_column1, row2_column2, row2_column3, ...],
      ...                                               ]

    and column titles.

    Args:
        coltitles (list): column titles

    Examples:
        >>> dg = DataContainer(['Energy_gap', 'dG', 'points'], comment="asd"))
        >>> dg.add_row([-300.0, 10.0, 2000])
        >>> dg.add_row([-200.0, 5.0, 1000])
        >>> dg
        DataContainer(['Energy_gap', 'dG', 'points'], comment='asd', Nrows=2)
        # get all rows
        >>> dg.get_rows()
        [[-300.0, 10.0, 2000], [-200.0, 5.0, 1000]]
        # get rows from specific columns
        >>> dg.get_rows(columns=('Energy_gap', 'points'))
        [[-300.0, 2000], [-200.0, 1000]]
        # get all columns
        >>> dg.get_columns()
        [(-300.0, -200.0), (10.0, 5.0), (2000, 1000)]
        # get specific columns
        >>> dg.get_columns(columns=(0, 1))
        [(-300.0, -200.0), (10.0, 5.0)]
        >>> dg.get_columns(columns=('Energy_gap', 'dG')
        [(-300.0, -200.0), (10.0, 5.0)]
        # clean up
        >>> dg.delete_rows()
    """

    def __init__(self, coltitles, comment=""):
        if not isinstance(coltitles, (list, tuple)):
            coltitles = [coltitles,]

        self.column_titles = list(coltitles)
        # a list containing rows of values
        # (each row is a list with length = len(coltitles))
        self._rows = []
        self.comment = comment

    def __repr__(self):
        return "DataContainer({}, comment='{}', Nrows={})" \
               "".format(self.column_titles, self.comment, len(self._rows))


    def get_columns(self, columns=None):
        """
        Transposes the array and returns the columns instead of rows.

        Args:
            columns (list), optional: return only columns with
                                      these indices and/or titles

        Returns:
            list of columns (list of lists)
        """
        if not columns:
            columns = []
        col_inds = []
        for col in columns:
            if type(col) == int:
                col_inds.append(col)
            else:
                col_inds.append(self.column_titles.index(str(col)))
        cols = list(zip(*self._rows))   # transpose
        if col_inds:
            return [cols[i] for i in col_inds]
        else:
            return cols


    def get_rows(self, columns=None):
        """Return the rows.

        Args:
            columns (list), optional: return only columns with
                                      these indices and/or titles

        Returns:
            list of rows (list of lists)
        """
        if columns:
            cols = self.get_columns(columns)
            return list(zip(*cols))
        else:
            return self._rows


    def add_row(self, row):
        """Add a row.

        Args:
            row (list): a list of values

        Raises:
            ValueError: if number of elements in row is not equal to
                        number of column titles
        """

        if len(row) != len(self.column_titles):
            raise ValueError("Number of elements is not equal to number "
                             "of columns, in row:\n{}".format(row))
        self._rows.append(list(row))


    def delete_rows(self):
        """
        Removes the rows.
        """
        self._rows = []

    def __str__(self):
        if self.comment:
            outs = "#" + self.comment + "\n"
        else:
            outs = ""
        for name in self.column_titles:
            width = len(name)
            if width < 10:
                width = 10
            outs += " {name:{width}} ".format(name=name, width=width)
        for row in self._rows:
            outs += "\n"
            for i, val in enumerate(row):
                try:
                    width = len(self.column_titles[i])
                    if width < 10:
                        width = 10
                except IndexError:
                    width = 20
                if type(val) == float:
                    outs += " {val:{width}.2f} ".format(val=val, width=width)
                else:
                    outs += " {val:{width}} ".format(val=str(val), width=width)
        return outs


# Shamelessly borrowed from
# http://www.genomearchitecture.com/2014/01/how-to-gunzip-on-the-fly-with-python
class gzopen(object):
    """Generic opener that decompresses gzipped files
    if needed. Encapsulates an open file or a GzipFile.
    Use the same way you would use 'open()'.
    """
    def __init__(self, fname):
        f = open(fname)
        # Read magic number (the first 2 bytes) and rewind.
        magic_number = f.read(2)
        f.seek(0)
        # Encapsulated 'self.f' is a file or a GzipFile.
        if magic_number == '\x1f\x8b':
            self.f = gzip.GzipFile(fileobj=f)
        else:
            self.f = f

    # Define '__enter__' and '__exit__' to use in
    # 'with' blocks. Always close the file and the
    # GzipFile if applicable.
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        try:
            self.f.fileobj.close()
        except AttributeError:
            pass
        finally:
            self.f.close()

    # Reproduce the interface of an open file
    # by encapsulation.
    def __getattr__(self, name):
        return getattr(self.f, name)
    def __iter__(self):
        return iter(self.f)
    def next(self):
        return next(self.f)



