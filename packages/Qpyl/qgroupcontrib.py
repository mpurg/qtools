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

import sys
import os
import time
import logging
import threading
from collections import OrderedDict as ODict

from Qpyl.core.qcalc import QCalc, QCalcError, QCalcInput, QCalcOutput
from Qpyl.core.qdyninp import QDynInput, QDynInputError
from Qpyl.core.qstructure import QStruct, QStructError
from Qpyl.common import __version__, np, DataContainer
from Qpyl.plotdata import PlotData

logger = logging.getLogger(__name__)

class _QGroupContribThread(threading.Thread):
    # threading class used for calling QGroupContrib._calcsingle
    def __init__(self, qgcinstance, semaphore, calcdir):
        threading.Thread.__init__(self)
        self.qgc = qgcinstance
        self.semaphore = semaphore
        self.calcdir = calcdir
        self.qcalc = QCalc(qgcinstance._qcalc_exec)
        self.qinps = None
        self.qouts = None
        self.error = None

    def run(self):
        with self.semaphore:
            if self.qgc.kill_event.is_set():
                return
            relp = os.path.relpath(self.calcdir)
            try:
                self.qinps, self.qouts = self.qgc._calcsingle(self.calcdir,
                                                              self.qcalc)
                logger.info("'{}' done.".format(relp))
            except QGroupContribError as error_msg:
                self.error = error_msg
                logger.info("'{}' failed.".format(relp))
            except Exception as error_msg:
                self.error = error_msg
                logger.critical("{}: Uncaught exception in calc thread: "
                                "{}".format(relp, error_msg))


class QGroupContribError(Exception):
    pass


class QGroupContrib(object):
    """Class for calculating LRA group contributions from EVB trajectories.

    Args:
        qcalc_exec (string): path of qcalc executable
        calcdirs (list of strings): list of directories
        pdb_file (string): PDB created with qprep
        en_list_fn (string): energy-files-list filename
        lambdas_A (tuple of floats): lambdas of state A (1.0, 0.0)
        lambdas_B (tuple of floats): lambdas of state B (0.5, 0.5)
        resid_first (int): index of first residue used for calcs
        resid_last (int): index of last residue used for calcs
        scale_ionized (float): scale down ionized residues (ARG, LYS,
                               HIP, GLU, ASP) by this factor
                               (see doi:10.1021/jp962478o)
        nthreads (int): number of threads
        qmask (list of ints): list of atom indexes to be used as the
                              Q mask for the GC calculations. By default,
                              this is obtained from the FEP file.

    """
    def __init__(self, qcalc_exec, calcdirs, pdb_file, en_list_fn,
                 lambdas_A, lambdas_B, resid_first, resid_last,
                 scale_ionized, nthreads, qmask=None):

        self._en_list_fn = en_list_fn
        self._qcalc_exec = qcalc_exec
        try:
            self._pdb_qstruct = QStruct(pdb_file, "pdb")
        except QStructError as error_msg:
            raise QGroupContribError("Can't parse PDB file '{}': {}"
                                     "".format(pdb_file, error_msg))

        self._calcdirs = calcdirs
        self._nthreads = nthreads
        self._lambdas_A = lambdas_A
        self._lambdas_B = lambdas_B
        self._resid_first = resid_first
        self._resid_last = resid_last
        self._scale_ionized = scale_ionized
        self._qmask = qmask


        self._qcalc_io = ODict()
        self.gcs = ODict()
        self.failed = ODict()
        self.qcalc_version = None

        self.kill_event = threading.Event()

        lambda1_st1, lambda2_st1 = lambdas_A[0], lambdas_B[0]
        sci = self._scale_ionized
        colnames = ["Residue id",
                    "Residue name",
                    "N",
                    "VdW(l={:5.4f}->l={:5.4f})_mean"
                    "".format(lambda1_st1, lambda2_st1),
                    "VdW(l={:5.4f}->l={:5.4f})_stdev"
                    "".format(lambda1_st1, lambda2_st1),
                    "El(l={:5.4f}->l={:5.4f})_(scale={})_mean"
                    "".format(lambda1_st1, lambda2_st1, sci),
                    "El(l={:5.4f}->l={:5.4f})_(scale={})_stdev"
                    "".format(lambda1_st1, lambda2_st1, sci)]
        self.gcs_stats = DataContainer(colnames)


    def calcall(self):
        """Run the GC calcs, update .gcs, .failed and .gcs_stats.
        """
        semaphore = threading.BoundedSemaphore(self._nthreads)

        self._qcalc_io.clear()
        self.gcs.clear()
        self.gcs_stats.delete_rows()
        self.failed.clear()
        threads = []
        for calcdir in self._calcdirs:
            threads.append(_QGroupContribThread(self, semaphore, calcdir))
            threads[-1].start()

        for t in threads:
            while t.isAlive():
                t.join(1.0)
                if self.kill_event.is_set():
                    try:
                        t.qcalc.process.terminate()
                    except Exception as e:
                        pass
                    return

            if t.error:
                self.failed[t.calcdir] = t.error
            else:
                self._qcalc_io[t.calcdir] = (t.qinps, t.qouts)

        # parse the output for results and
        # calculate LRAs for each dir
        for _dir, (_, qouts) in self._qcalc_io.iteritems():
            gcs = []
            failed_flag = False
            for qout in qouts:
                try:
                    qco = QCalcOutput(qout)
                    res = qco.results["gc"]
                    if not self.qcalc_version:
                        self.qcalc_version = qco.qcalc_version
                except (QCalcError, KeyError) as error_msg:
                    self.failed[_dir] = error_msg
                    failed_flag = True
                    break
                gc = {}
                for row in res.get_rows():
                    resid, vdw, el = int(row[0]), float(row[1]), float(row[2])
                    gc[resid] = {"vdw": vdw, "el": el}
                gcs.append(gc)

            if failed_flag:
                continue

            resids = sorted(gcs[0].keys())
            resnames = [self._pdb_qstruct.residues[ri-1].name for ri in resids]

            # do the LRA thingy
            # LRA=0.5*(<E2-E1>_conf1+<E2-E1>_conf2)
            e2e1_st1_vdw = [gcs[1][key]["vdw"] - gcs[0][key]["vdw"] for key in resids]
            e2e1_st1_el = [gcs[1][key]["el"] - gcs[0][key]["el"] for key in resids]
            e2e1_st2_vdw = [gcs[3][key]["vdw"] - gcs[2][key]["vdw"] for key in resids] 
            e2e1_st2_el = [gcs[3][key]["el"] - gcs[2][key]["el"] for key in resids]

            vdw_lra = [0.5*(a + b) for a, b in zip(e2e1_st1_vdw, e2e1_st2_vdw)]
            el_lra = [0.5*(a + b) for a, b in zip(e2e1_st1_el, e2e1_st2_el)]

            # scale the ionized residues
            for i, resname in enumerate(resnames):
                if resname in ("ARG", "LYS", "HIP", "ASP", "GLU"):
                    el_lra[i] = el_lra[i] / self._scale_ionized

            # write the DataContainer
            lambda1_st1 = self._lambdas_A[0]
            lambda2_st1 = self._lambdas_B[0]
            gc_lra = DataContainer(["Residue_id", "Residue name",
                                    "VdW(l={:5.4f}->l={:5.4f})"
                                    "".format(lambda1_st1, lambda2_st1),
                                    "El(l={:5.4f}->l={:5.4f})_(scale={})"
                                    "".format(lambda1_st1, lambda2_st1,
                                              self._scale_ionized)])

            for row in zip(resids, resnames, vdw_lra, el_lra):
                gc_lra.add_row(row)

            self.gcs[_dir] = gc_lra

        # get GC stats over all directories
        self.gcs_stats.delete_rows()
        gcs = {}
        for _, gc in self.gcs.iteritems():
            for row in gc.get_rows():
                resid, resname = row[0:2]
                res_key = "{}.{}".format(resid, resname)
                values = [[val,] for val in row[2:]]
                if not gcs.has_key(res_key):
                    gcs[res_key] = values
                else:
                    for i, val in enumerate(gcs[res_key]):
                        val.extend(values[i])

        # iterate through each residue and calculate
        # means and stdevs
        for res_key in sorted(gcs.keys()):
            rc = gcs[res_key]
            resid, resname = res_key.split(".")
            # get mean and stdev
            rc_stats = [int(resid), resname, len(rc[0]),
                        np.mean(rc[0]), np.std(rc[0]), # vdw
                        np.mean(rc[1]), np.std(rc[1])] # el

            self.gcs_stats.add_row(rc_stats)



    def _calcsingle(self, calcdir, qcalc):
        # find input files with given lambdas
        # (and correct energy files)
        # extract information and run qcalc for each combination
        #   fep_000_1.000.dcd, "1.00 0.00"
        #   fep_000_1.000.dcd, "0.50 0.50"
        #   fep_025_0.500.dcd, "1.00 0.00"
        #   fep_025_0.500.dcd, "0.50 0.50"
        # return input output strings as a tuple of lists of strings
        # ( [inp1, inp2, inp3, inp4], [out1, out2, out3, out4] )
        # or raise QGroupContribError on failure

        # get the list of energy-files
        try:
            en_list_fn = os.path.join(calcdir, self._en_list_fn)
            en_list_fn_str = open(en_list_fn, 'r').read()
        except IOError:
            raise QGroupContribError("No energy-files list '{}'."
                                     "".format(self._en_list_fn))

        en_list = [enf for enf in en_list_fn_str.split("\n") \
                                                if enf.strip() != ""]

        if not en_list:
            raise QGroupContribError("No energy files in '{}'."
                                     "".format(self._en_list_fn))

        # parse all input files in calcdir for
        # a valid energy file and lambda values
        inp_fns = [inp for inp in os.listdir(calcdir) if inp.endswith(".inp")]
        lambda_inp_map = {}
        for inp in inp_fns:
            try:
                inp_file = os.path.join(calcdir, inp)
                qdi = QDynInput(input_string=open(inp_file, "r").read())
            except (IOError, QDynInputError) as error_msg:
                logger.debug("Error reading Q input '{}': {}"
                             "".format(inp, error_msg))
                continue

            try:
                lambda_st1 = float(qdi.parameters["lambdas"].split()[0])
                en_file = qdi.parameters["files"]["energy"]
            except KeyError:
                logger.debug("Input '{}' missing lambda or energy file"
                             "".format(inp))
                continue

            if en_file not in en_list:
                continue

            lambda_key = "{:.6f}".format(lambda_st1)
            try:
                inp2 = lambda_inp_map[lambda_key][0]
            except KeyError:
                lambda_inp_map[lambda_key] = (inp, qdi)
            else:
                raise QGroupContribError("Same lambda values in Qdyn "
                                         "inputs: '{}', '{}' ??"
                                         "".format(inp, inp2))

        # get inputs that match specified state1 lambda values
        lambdas_st1 = (self._lambdas_A[0], self._lambdas_B[0])
        try:
            inputs = []
            for lamb_st1 in lambdas_st1:
                lamb_key = "{:.6f}".format(lamb_st1)
                inputs.append(lambda_inp_map[lamb_key])
        except KeyError:
            raise QGroupContribError("QDyn input with lambda=='{}' "
                                     "(and energy file in '{}') not found."
                                     "".format(lamb_st1, en_list_fn))

        # get topology, fep and trajectory filenames from the inputs
        top_fn, fep_fn, dcd_fns = None, None, []
        for inp, qdi in inputs:
            try:
                tmp_top_fn = qdi.parameters["files"]["topology"]
            except KeyError:
                raise QGroupContribError("Topology not found in Qdyn "
                                         "input '{}'.".format(inp))
            if top_fn and top_fn != tmp_top_fn:
                raise QGroupContribError("Qdyn inputs with different "
                                         "topologies: '{}', '{}' ??"
                                         "".format(top_fn, tmp_top_fn))

            try:
                tmp_fep_fn = qdi.parameters["files"]["fep"]
            except KeyError:
                raise QGroupContribError("Fep file not found in Qdyn "
                                         "input '{}'.".format(inp))
            if fep_fn and fep_fn != tmp_fep_fn:
                raise QGroupContribError("Qdyn inputs with different "
                                         "fep files: '{}', '{}' ??"
                                         "".format(fep_fn, tmp_fep_fn))

            try:
                tmp_dcd_fn = qdi.parameters["files"]["trajectory"]
            except KeyError:
                raise QGroupContribError("Trajectory file not found in Qdyn "
                                         "input '{}'.".format(inp))

            top_fn = tmp_top_fn
            fep_fn = tmp_fep_fn
            dcd_fns.append(tmp_dcd_fn)

        # check if files are missing
        for fn in [top_fn, fep_fn] + dcd_fns:
            if not os.path.lexists(os.path.join(calcdir, fn)):
                raise QGroupContribError("Missing file: {}".format(fn))

        if not self._qmask:
            # parse fep for q atom numbers
            with open(os.path.join(calcdir, fep_fn), "r") as fep:
                section = ""
                q_atoms = []
                for line in fep.readlines():
                    line = line.split("#")[0].split("!")[0].strip()
                    if line == "":
                        continue
                    elif line[0] == "[":
                        section = line
                    elif section == "[atoms]":
                        q_atoms.append(line.split()[1])
        else:
            q_atoms = self._qmask

        masks = ["{} {}".format(ai, ai) for ai in q_atoms]

        # make qcalc inputs for every combination of
        # configuration (dcd) and potential (lambda),
        # run them and return the inputs and outputs
        combs = ((dcd_fns[0], self._lambdas_A),   # E1_conf1
                 (dcd_fns[0], self._lambdas_B),   # E2_conf1
                 (dcd_fns[1], self._lambdas_A),   # E1_conf2
                 (dcd_fns[1], self._lambdas_B))  # E2_conf2
        # example with lambdas "1.00 0.00" and "0.50 0.50":
        #
        # fep_000_1.000.dcd, (1.00, 0.00)
        # fep_000_1.000.dcd, (0.50, 0.50)
        # fep_025_0.500.dcd, (1.00, 0.00)
        # fep_025_0.500.dcd, (0.50, 0.50)

        input_strings = []
        output_strings = []
        for dcdfile, lambdas in combs:
            qci = QCalcInput(top_fn, [dcdfile,], fep_fn, lambdas)

            qci.add_residue_nb_mon(self._resid_first,
                                   self._resid_last,
                                   masks)

            qcalc_inp_str = qci.get_string()

            try:
                qcalc_out_str = qcalc.run(qcalc_inp_str, workdir=calcdir)
            except QCalcError as error_msg:
                raise QGroupContribError(error_msg)

            input_strings.append(qcalc_inp_str)
            output_strings.append(qcalc_out_str)

        return (input_strings, output_strings)



    @property
    def details(self):

        fails = "\n".join(["{}: {}".format(cd, e) for cd, e in \
                                                 self.failed.iteritems()])

        calcdirs = ", ".join(self._calcdirs)
        outstr = """
---------------------------------- GC details ---------------------------------
# Calculated with: QTools ({version}), QCalc ({qcalc_version})
# Work dir: {cwd}
# Date: {date}
# CMDline: {cmdline}

Directories:
{dirs}

Fails:
{fails}
-------------------------------------------------------------------------------
""".format(version=__version__, cwd=os.getcwd(), date=time.ctime(),
           cmdline=" ".join(sys.argv), qcalc_version=self.qcalc_version,
           fails=fails or "None", dirs=calcdirs)

        return outstr


    @property
    def plotdata(self):
        """Return GC data as a dictionary of PlotData objects.

        Example keys in returned dictionary:
            'gc_el': PlotData of electrostatic GCs, one subplot with means vs
                  residues index

            'gc_el_top': 'el', sorted by absolute contribution, only first 20
                      means vs "resid.resname"

            'gc_vdw': PlotData of electrostatic GCs, one subplot with means vs
                   residue indexes

            'gc_vdw_top': 'vdw', sorted by absolute contribution, only first 20
                       means vs "resid.resname"

        """

        plots = ODict()

        # all failed
        if not self.gcs:
            return plots

        lamb1, lamb2 = self._lambdas_A[0], self._lambdas_B[0]
        # make PlotData objects

        cols = self.gcs_stats.get_columns()
        resids, _, _, vdws, vdwss, els, elss = cols

        N = len(self.gcs)
        plots["gc_el"] = PlotData("LRA GC (electrostatic): "
                                  "dG( l={} -> l={} ) (iscale={})"
                                  "".format(lamb1, lamb2,
                                            self._scale_ionized),
                                  xlabel="Residue index",
                                  ylabel="Free energy  [kcal/mol]",
                                  plot_type="bar")
        plots["gc_el"].add_subplot("mean_N={}".format(N),
                                   resids, els, yerror=elss)

        plots["gc_vdw"] = PlotData("LRA GC (VdW): "
                                   "dG( l={} -> l={} )".format(lamb1, lamb2),
                                   xlabel="Residue index",
                                   ylabel="Free energy  [kcal/mol]",
                                   plot_type="bar")
        plots["gc_vdw"].add_subplot("mean_N={}".format(N),
                                    resids, vdws, yerror=vdwss)


        sorted_rows = sorted(self.gcs_stats.get_rows(), \
                                            key=lambda x: -abs(x[5]))

        resids, resnames, _, _, _, els, elss = zip(*sorted_rows[:20])

        keys = ["{}_{}".format(rn.capitalize(), ri) \
                                for ri, rn in zip(resids, resnames)]

        plots["gc_el_top"] = PlotData("LRA GC (electrostatic): "
                                      "dG( l={} -> l={} ) (iscale={}), top 20"
                                      "".format(lamb1, lamb2,
                                                self._scale_ionized),
                                      xlabel="Residue",
                                      ylabel="Free energy  [kcal/mol]",
                                      plot_type="bar")
        plots["gc_el_top"].add_subplot("mean_N={}".format(N),
                                       keys, els,
                                       yerror=elss)

        return plots




