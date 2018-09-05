#!/usr/bin/env python2
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

"""
This module contains the QMapper class for automating the calculation and
calibration of EVB reaction free profiles (via Qfep). 
"""

import sys
import os
import time
import logging
import threading
from collections import OrderedDict as ODict

from Qpyl.core.qfep import QFep, QFepError, QFepInput
from Qpyl.core.qfep import QFepOutput, QFepOutputError
from Qpyl.common import __version__, np

logger = logging.getLogger(__name__)

class _QMapperThread(threading.Thread):
    # threading class used for calling QMapper._mapsingle
    def __init__(self, qmapperinstance, semaphore, mapdir, supress_info):
        threading.Thread.__init__(self)
        self.qm = qmapperinstance
        self.semaphore = semaphore
        self.mapdir = mapdir
        self.qfep = QFep(qmapperinstance._qfep.qfep_exec)
        self.supress_info = supress_info
        self.qinp = None
        self.qout = None
        self.error = None

    def run(self):
        with self.semaphore:

            if self.qm.kill_event.is_set():
                return

            try:
                self.qinp, self.qout = self.qm.mapsingle(self.mapdir,
                                                         _qfep=self.qfep)
                if self.supress_info:
                    logger.debug("{} mapped.".format(self.mapdir))
                else:
                    logger.info("{} mapped.".format(self.mapdir))

            except QMapperError as error_msg:
                self.error = error_msg
                if self.supress_info:
                    logger.debug("{} failed.".format(self.mapdir))
                else:
                    logger.info("{} failed.".format(self.mapdir))

            except Exception as error_msg:
                self.error = error_msg
                logger.critical("{}: Uncaught exception in mapping thread: "
                                "{}".format(self.mapdir, error_msg))


class QMapperError(Exception):
    pass

class QMapper(object):
    """'Maps' EVB simulations in specified directories.

    Basically a wrapper around Qpyl.qfep.QFep

    At the moment, supports only:
    - 2 states
    - no predefined Hij
    - single, constant Hij (mu=eta=r0 = 0)

    Args:
        qfep_exec (string): qfep executable path
        mapdirs (list of strings): list of directories
        en_list_fn (string): energy-files-list filename
        hij (float): H12_A offdiagonal constant
        alpha (float): state2 alpha shift
        temperature (float): temperature in Kelvin
        gas_const (float): gas_constant in kcal/mol/K
        points_skip (int): number of energy points to skip in each frame
        gap_bins (int): number of gap-bins
        minpts_bin (int): minimum number of points per bin
        nthreads (int): number of threads

    QMapper.parms contains mapping parameters, any direct change will be
    applied in the mapping.
    QMapper.mapped contains the inputs and outputs of succesfully mapped dirs
    # { 'directory_path1': (qfep_input_string, qfep_output_string), ... }
    QMapper.failed contains the error messages of failed dirs
    # { 'directory_path2': 'error_message', ... }

    QMapper.mapall() maps the directories and updates .mapped and .failed
    QMapper.fit_to_reference() fits and updates 'hij' and 'alpha' parameters

    """
    def __init__(self, qfep_exec, mapdirs, en_list_fn, hij, alpha, temperature,
                 gas_const, points_skip, gap_bins, minpts_bin, nthreads):

        self._en_list_fn = en_list_fn
        self._qfep = QFep(qfep_exec)
        self._mapdirs = [os.path.relpath(m) for m in mapdirs]
        self._nthreads = nthreads

        self.parms = {"hij": hij,
                      "alpha": alpha,
                      "temperature": temperature,
                      "gas_const": gas_const,
                      "points_skip": points_skip,
                      "gap_bins": gap_bins,
                      "minpts_bin": minpts_bin}

        self.mapped = ODict()
        self.failed = ODict()

        self.kill_event = threading.Event()


    def mapall(self, _supress_info=False):
        """'Map' all directories and update QMapper.mapped and QMapper.failed
        """
        semaphore = threading.BoundedSemaphore(self._nthreads)

        self.mapped.clear()
        self.failed.clear()
        threads = []
        for mapdir in self._mapdirs:
            threads.append(_QMapperThread(self, semaphore,
                                          mapdir, _supress_info))
            threads[-1].start()

        for t in threads:
            while t.isAlive():
                t.join(1.0)
                if self.kill_event.is_set():
                    try:
                        t.qfep.process.terminate()
                    except:
                        pass
                    return

            if t.error:
                self.failed[t.mapdir] = t.error
            else:
                self.mapped[t.mapdir] = (t.qinp, t.qout)


    def mapsingle(self, mapdir, _qfep=None):
        """'Map' a single directory.

        Args:
            mapdir (string):  directory path
            _qfep (QFep): internal, for threading

        Returns:
            qfep_inp_str (string):  QFep input string used
            qfep_out_str (string):  QFep output string returned

        Raises QMapperError
        """

        qfep = _qfep or self._qfep

        # get the list of energy-files
        en_list_fn = os.path.join(mapdir, self._en_list_fn)
        try:
            en_fn_str = open(en_list_fn, 'r').read()
        except IOError:
            raise QMapperError("No energy-files list '{}'."
                               "".format(self._en_list_fn))

        en_list = [enf for enf in en_fn_str.split("\n") if enf.strip() != ""]
        if not en_list:
            raise QMapperError("No energy files in '{}'."
                               "".format(self._en_list_fn))

        for enf in en_list:
            if not os.path.lexists(os.path.join(mapdir, enf)):
                raise QMapperError("Energy file '{}' missing.".format(enf))

        # make qfep input
        qfep_inp_str = QFepInput(en_list, **self.parms).get_string()

        # run qfep
        try:
            qfep_out_str = qfep.run(qfep_inp_str, workdir=mapdir)
        except QFepError as error_msg:
            raise QMapperError(error_msg)

        return (qfep_inp_str, qfep_out_str)


    def fit_to_reference(self, dga_ref, dg0_ref, step_size=10.0,
                         threshold=0.005, max_iterations=10):
        """Fit Hij and alpha to obtain desired dG* and dG0 values.

        Updates 'hij' and 'alpha' parameters.

        Args:
            dga_ref (float): reference activation free energy
            dg0_ref (float): reference reaction free energy
            step_size (float, optional): step size for calculating the\
                                         gradient, default=10.0
            threshold (float, optional): convergence threshold, default=0.005
            max_iterations (int, optional): max number of iterations,\
                                            default=10

        Returns:
            (boolean):  True if succeeded, False if max_iterations was\
                        reached.

        Raises QMapperError if all of the directories fail to map.
        """

        iteration = 1
        while iteration <= max_iterations:
            logger.info("Iteration #{}".format(iteration))
            logger.info("{:>10} {:>10} {:>10} {:>10}"
                        "".format("Hij", "Alpha", "dG#", "dG0"))
            if iteration == 1:
                try:
                    dga_mean, dg0_mean = self._getmeans()
                except QMapperError:
                    raise
                logger.info("{:10.2f} {:10.2f} {:10.2f} {:10.2f}"
                            "".format(self.parms["hij"], self.parms["alpha"],
                                      dga_mean, dg0_mean))

            try:
                dga_mean, dg0_mean = self._do_iteration(dga_mean, dg0_mean,
                                                        dga_ref, dg0_ref,
                                                        step_size)
            except QMapperError:
                raise

            if abs(dga_mean - dga_ref) < threshold and \
                   abs(dg0_mean - dg0_ref) < threshold:
                return True
            iteration += 1
        return False

    def _getmeans(self):
        # called by fit_to_reference and _do_iteration
        # calls map_all()
        # analyses the outputs and returns dG* and dG0 means
        # raises QMapperError on total failure

        self.mapall(_supress_info=True)
        for (qfo, err) in self.failed.iteritems():
            logger.info("Failed to map '{}': {}".format(qfo, err))

        if not self.mapped:
            raise QMapperError("All directories failed to map! Try changing "
                               "the initial-guess values (Hij and alpha) "
                               "or step_size... Also, check the mapping "
                               "parameters (skip, bins, ...).")

        dga, dg0 = [], []
        for mapdir, (_, qfo_str) in self.mapped.iteritems():
            try:
                qfo = QFepOutput(qfo_str)
                dga.append(qfo.part3.dga)
                dg0.append(qfo.part3.dg0)
            except QFepOutputError as error_msg:
                logger.info("Failed to analyse '{}': {}"
                            "".format(mapdir, error_msg))
            except Exception as error_msg:
                logger.warning("Uncaught exception when analysing '{}': {}"
                               "".format(mapdir, error_msg))

        if not dga or not dg0:
            raise QMapperError("All directories failed to analyse! Try "
                               "changing the initial-guess values (Hij and "
                               "alpha) or step_size...")
        return np.mean(dga), np.mean(dg0)



    def _do_iteration(self, dga_init, dg0_init, dga_ref, dg0_ref, step_size):
        # called by fit_to_reference
        # changes Hij and alpha and calls _getmeans
        # calculates the gradients of the changes and updates the variables

        hij_init = self.parms["hij"]
        alpha_init = self.parms["alpha"]

        # Step1, changing alpha=alpha+step_size,
        self.parms["hij"] = hij_init
        self.parms["alpha"] = alpha_init + step_size
        means1 = self._getmeans()
        logger.info("{:10.2f} {:10.2f} {:10.2f} {:10.2f}"
                    "".format(self.parms["hij"], self.parms["alpha"],
                              means1[0], means1[1]))

        # Step2, changing hij=hij+step_size,
        self.parms["hij"] = hij_init + step_size
        self.parms["alpha"] = alpha_init
        means2 = self._getmeans()
        logger.info("{:10.2f} {:10.2f} {:10.2f} {:10.2f}"
                    "".format(self.parms["hij"], self.parms["alpha"],
                              means2[0], means2[1]))

        # calculate gradients for hij and gs
        ka1 = (means1[0] - dga_init) / step_size
        ka2 = (means2[0] - dga_init) / step_size
        k01 = (means1[1] - dg0_init) / step_size
        k02 = (means2[1] - dg0_init) / step_size

        # get the shifts
        delta_dga = dga_ref - dga_init
        delta_dg0 = dg0_ref - dg0_init
        try:
            dalpha = ((delta_dga/ka1)-(delta_dg0*ka2/ka1/k02)) \
                   / (1-k01*ka2/ka1/k02)
            dhij = (delta_dg0 - dalpha*k01)/k02
        except ZeroDivisionError:
            dalpha = +1
            dhij = +1

        # Step 3, get results
        self.parms["hij"] = hij_init + dhij
        self.parms["alpha"] = alpha_init + dalpha
        means3 = self._getmeans()
        logger.info("{:10.2f} {:10.2f} {:10.2f} {:10.2f}"
                    "".format(self.parms["hij"], self.parms["alpha"],
                              means3[0], means3[1]))
        return means3




    @property
    def input_parms_str(self):
        # TODO: CLI script doesn't belong in here
        inp_parms = "q_mapper.py {hij} {alpha} "\
                    "--bins {gap_bins} --skip {points_skip} "\
                    "--min {minpts_bin} --temp {temperature} "\
                    "".format(**self.parms)
        return inp_parms


    @property
    def details(self):

        fails = "\n".join(["{}: {}".format(md, e) for md, e in \
                                                 self.failed.iteritems()])

        qfep_version = "Unknown, likely ancient"
        for _, qfep_out in self.mapped.values():
            try:
                qfo = QFepOutput(qfep_out)
            except Exception:
                pass
            else:
                qfep_version = qfo.header.qfep_version
                break

        mapdirs = ", ".join(self._mapdirs)
        outstr = """
------------------------------- Mapping details -------------------------------
# Mapped with: Qtools ({version}), Qfep ({qfep_version})
# Qfep path: {qfep_exec}
# Work dir: {cwd}
# Date: {date}
# CMDline: {cmdline}

Directories:
{dirs}

Input:
{inp_parms}

Fails:
{fails}
-------------------------------------------------------------------------------
""".format(version=__version__, cwd=os.getcwd(), date=time.ctime(),
           inp_parms=self.input_parms_str, cmdline=" ".join(sys.argv),
           fails=fails or "None", dirs=mapdirs, qfep_version=qfep_version,
           qfep_exec=os.path.abspath(self._qfep.qfep_exec))

        return outstr

