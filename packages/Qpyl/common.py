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
# Some common classes and functions

import math
import sys
import os
import shutil
import logging

logger = logging.getLogger(__name__)

__version__ = "0.5.11"


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

    def format(self, record):
        self._fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        return logging.Formatter.format(self, record)


def init_logger(name,
                level=None,
                handler=None,
                formatter=None):
    """Helper function for initializing the logger.

    Args:
        name (string):  module name, usually root: 'Qpyl'
        level (int, optional):  logging level (DEBUG, INFO, WARNING...),
                                default is INFO
        handler: (logging.Handler, optional):  default is
                                               StreamHandler(sys.stdout)
        formatter: (logging.Formatter, optional): default is
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
class np():
    @staticmethod
    def kendall_tau(list_a, list_b):
        """Calculate Kendall's tau-a correlation coefficient.

        Args:
            list_a (list):  list of ints/floats
            list_b (list):  list of ints/floats
        """
        lists = zip(list_a, list_b)
        N, conc, disc = len(lists), 0.0, 0.0
        for i, (x, y) in enumerate(lists[:-1]):
            for (x2, y2) in lists[i+1:]:
                dx, dy = (x - x2), (y - y2)
                if dx * dy < 0:
                    disc += 1
                elif dx * dy > 0:
                    conc += 1
                else:
                    pass # ignore ties
        return (conc - disc)/(N * (N - 1)/2.0)

    @staticmethod
    def mean(vals):
        """Calculate mean.

        Args:
            vals (list of float):  sample values
        """
        N = len(vals)
        if N == 0: 
            return float('nan')
        return sum(vals) * 1.0 / N

    @staticmethod
    def std(vals, ddof=1):
        """Calculate standard deviation.
        
        Args:
            vals (list of float):  sample values
            ddof (int, optional):  delta degrees of freedom (default is 1,
                                   producing *sample* standard deviation)
        """
        N = len(vals)
        if N == 0 or N-ddof == 0: 
            return float('nan')
        mean = np.mean(vals)
        variance = map(lambda x: (x - mean)**2, vals)
        return math.sqrt(sum(variance)/(N - ddof))

    @staticmethod
    def std_error(vals, ddof=1):
        """Calculates standard error of mean.
        
        Args:
            vals (list of float):  sample values
            ddof (int, optional):  see np.std()
        """
        N = len(vals)
        if N == 0 or N-ddof == 0:
            return float('nan')
        return np.std(vals, ddof=ddof) / math.sqrt(N)

    @staticmethod
    def median(vals):
        """Calculate median

        Args:
            vals (list of float):  sample values
        """
        N = len(vals)
        if N == 0:
            return float('nan')
        vals = sorted(vals)
        if N % 2 == 0: #even
            return np.mean((vals[N/2-1], vals[N/2]))
        else: #odd
            return vals[N/2]





class DataContainer(object):
    """
    Contains a two dimensional array of values:

    [ [ row1_column1, row1_column2, row1_column3, ...],
      [ row2_column1, row2_column2, row2_column3, ...],
      ...                                               ]

    and column titles.

    Args:
        coltitles (list): column titles

    Example of usage:
    >>> dg_de = DataContainer( ['Energy_gap', 'dG'] )
    >>> dg_de.add_row( [-300.0, 10.0 ]
    # reversed rows
    >>> rows = dg_de.get_rows( reversed(dg_de.get_column_titles()) )
    >>> cols = dg_de.get_columns( columns=[0, 1] )
    """

    def __init__(self, coltitles):
        if not isinstance(coltitles, (list, tuple)):
            coltitles = [coltitles,]

        self._column_titles = list(coltitles)
        # a list containing rows of values
        # (each row is a list with length = len(coltitles))
        self._rows = []
        self.comment = None


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
                col_inds.append(self._column_titles.index(str(col)))
        cols = zip(*self._rows)   # transpose
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
            return zip(*cols)
        else:
            return self._rows


    def get_column_titles(self):
        """Return the list of column titles."""

        return self._column_titles


    def add_row(self, row):
        """Add a row.

        Args:
            row (list): a list of values

        Raises:
            ValueError: if number of elements in row is not equal to
                        number of column titles
        """

        if len(row) != len(self._column_titles):
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
        for name in self._column_titles:
            width = len(name)
            if width < 10:
                width = 10
            outs += " {name:{width}} ".format(name=name, width=width)
        for row in self._rows:
            outs += "\n"
            for i, val in enumerate(row):
                try:
                    width = len(self._column_titles[i])
                    if width < 10:
                        width = 10
                except IndexError:
                    width = 20
                if type(val) == float:
                    outs += " {val:{width}.2f} ".format(val=val, width=width)
                else:
                    outs += " {val:{width}} ".format(val=str(val), width=width)
        return outs

