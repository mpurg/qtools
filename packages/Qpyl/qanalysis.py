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
# Module for analysing Q fep and Q dynamics output (logfiles).
# Two main classes are QAnalyseFeps and QAnalyseDyns.
# See q_analysefeps.py and q_analysedyns.py for usage examples.
#

import os
import re
import logging
from collections import OrderedDict as ODict

from Qpyl.core.qfep import QFepOutput, QFepOutputError
from Qpyl.common import DataContainer, np
from Qpyl.plotdata import PlotData

logger = logging.getLogger(__name__)

class QAnalyseDynError(Exception):
    pass

class QAnalyseFeps(object):
    """Wrapper class for QFepOutput for analysing multiple outputs.

    Calculates statistics (mean, stdev,...), makes pretty summaries,
    extracts all data to PlotData objects via neat functions.

    Args:
        qfep_outputs (list of strings or list of tuples(string, string)): 
                    pathnames of qfep outputs or (name, output_string) tuples

        lra_lambdas (list of floats, optional):  lambdas at which LRA energies
                                                 should be calculated at.
                                                 Usual case is (1.0, 0.0).
                                                 Default = None (no LRA calcs)
                                                                   instead of 
        _parent (QAnalyseFeps): used internally for QCP and exclusions
        _subcalcs (string): used internally for QCP and exclusions


    Usage (printing statistics and failures, see q_analysefeps.py):

    qos = [os.path.join(md, "qfep.log") for md in sorted(mapdirs)]
    qaf = QAnalyseFeps(qos, lra_lambdas=(1.0, 0.0))

    print "Statistics"
    print qaf.stats_str
    for sub_calc_key, sub_calc in sorted(qaf.sub_calcs.iteritems()):
        print sub_calc.stats_str

    if qaf.failed:
        print "Failed to parse:"
    for failed_path, failed_msg in sorted(qaf.failed.iteritems()):
        relp = os.path.relpath(failed_path)
        print "-> {}: {}".format(relp, failed_msg)

    """

    def __init__(self, qfep_outputs, lra_lambdas=None,
                 _parent=None, _subcalc_key=None):

        self._qfep_outputs = sorted(set(qfep_outputs))
        self._lra_lambdas = lra_lambdas
        self._parent = _parent
        self._subcalc_key = _subcalc_key

        self.qfos = ODict()
        self.dgas = ODict()
        self.dg0s = ODict()
        self.dgs_fep = ODict()
        self.lras = ODict()
        self.failed = ODict()
        self.failed_dg = ODict()
        self.sub_calcs = ODict()

        # main calculation (not subcalc)
        if self._parent == None:
            # check if the files exist and if they parse
            for qfep_output in self._qfep_outputs:
                if isinstance(qfep_output, basestring):
                    try:
                        qfep_out_string = open(qfep_output, "r").read()
                    except IOError as error_msg:
                        self.failed[qfep_output] = error_msg
                        continue
                    else:
                        qfep_out_name = qfep_output
                else:
                    qfep_out_name, qfep_out_string = qfep_output

                try:
                    qfo = QFepOutput(qfep_out_string)
                except QFepOutputError as error_msg:
                    self.failed[qfep_out_name] = error_msg
                except Exception as error_msg:
                    self.failed[qfep_out_name] = "UNCAUGHT EXCEPTION: {}"\
                                                 "".format(error_msg)
                else:
                    self.qfos[qfep_out_name] = qfo
        # exctract data (exclusion or QCP) from parent's sub_calcs
        else:
            self.qfos_failed = self._parent.failed
            for p_qfo_key, p_qfo in self._parent.qfos.iteritems():
                if self._subcalc_key in p_qfo.sub_calcs:
                    self.qfos[p_qfo_key] = p_qfo.sub_calcs[self._subcalc_key]

        # get dG_FEP, dGa, dG0, LRAs
        for qfep_output, qfo in self.qfos.iteritems():
            self.dgs_fep[qfep_output] = qfo.part1.dg
            try:
                dga = qfo.part3.dga
                dg0 = qfo.part3.dg0
                self.dgas[qfep_output] = dga
                self.dg0s[qfep_output] = dg0
                if qfo.part3.warning:
                    logger.warning("{}: {}".format(qfep_output,
                                                   qfo.part3.warning))
            except QFepOutputError as error_msg:
                self.failed_dg[qfep_output] = error_msg
            except Exception as error_msg:
                self.failed_dg[qfep_output] = "UNCAUGHT EXCEPTION: {}"\
                                              "".format(error_msg)
            # get LRA energies or log if fails (3 states, weird data)
            if self._lra_lambdas != None:
                try:
                    lra = qfo.part0.calc_lra(self._lra_lambdas[0],
                                             self._lra_lambdas[1])
                except Exception as e:
                    logger.warning("LRA failed on '{}': {}"
                                   "".format(qfep_output, e))
                else:
                    self.lras[qfep_output] = lra

            # add sub_calcs to self.sub_calcs, if they exist
            for subcalc_key in qfo.sub_calcs:
                if subcalc_key not in self.sub_calcs:
                    subqaf = QAnalyseFeps(self._qfep_outputs,
                                          self._lra_lambdas,
                                          _parent=self,
                                          _subcalc_key=subcalc_key)
                    self.sub_calcs[subcalc_key] = subqaf

        # TODO: check that parameters match in different outputs,
        # otherwise, log a warning


    @property
    def dg_all(self):
        """DataContainer with all main and subcalc free energies."""
        subcalcs = sorted(self.sub_calcs.keys())
        coltitles = ["Qfep_output", "dG*", "dG0", "dG_lambda"]
        for sc in subcalcs:
            coltitles.extend(["{}_dG*".format(sc),
                              "{}_dG0".format(sc),
                              "{}_dG_lambda".format(sc)])

        dc = DataContainer(coltitles)
        for qfep_output in sorted(self.qfos):
            relp = os.path.relpath(qfep_output)
            dga = self.dgas.get(qfep_output, None)
            dg0 = self.dg0s.get(qfep_output, None)
            dg_fep = self.dgs_fep.get(qfep_output, None)
            row = [relp, dga, dg0, dg_fep]

            for subcalc in subcalcs:
                dga = self.sub_calcs[subcalc].dgas.get(qfep_output, None)
                dg0 = self.sub_calcs[subcalc].dg0s.get(qfep_output, None)
                dg_fep = self.sub_calcs[subcalc].dgs_fep.get(qfep_output, None)
                row.extend([dga, dg0, dg_fep])

            dc.add_row(row)
        return dc


    @property
    def stats_str(self):
        """Free energy stats in string format."""
        dgas = self.dgas.values()
        dg0s = self.dg0s.values()
        dgs_fep = self.dgs_fep.values()

        allres = {}
        allres["calc_type"] = self._subcalc_key or ""
        allres["dg_n"] = len(dgas)
        allres["dga"] = (np.mean(dgas), np.std(dgas),
                         np.median(dgas), np.std_error(dgas))

        allres["dg0"] = (np.mean(dg0s), np.std(dg0s),
                         np.median(dg0s), np.std_error(dg0s))

        allres["dg_fep_n"] = len(dgs_fep)
        allres["dg_fep"] = (np.mean(dgs_fep), np.std(dgs_fep),
                            np.median(dgs_fep), np.std_error(dgs_fep))
        return """\
# {calc_type:<15} Mean      Std.dev    Median    Std.error       N
dG*         {dga[0]:10.2f} {dga[1]:10.2f} {dga[2]:10.2f} {dga[3]:10.2f} {dg_n:10}
dG0         {dg0[0]:10.2f} {dg0[1]:10.2f} {dg0[2]:10.2f} {dg0[3]:10.2f} {dg_n:10}
dG_lambda   {dg_fep[0]:10.2f} {dg_fep[1]:10.2f} {dg_fep[2]:10.2f} \
{dg_fep[3]:10.2f} {dg_fep_n:10}
""".format(**allres)



    @property
    def plotdata(self):
        """Return 'useful data' as a dictionary of PlotData objects.

        Each qfep_output will be a subplot in one PlotData, except in the case
        of LRA where there is only one subplot: the average and stdev over all
        outputs.

        Useful data:
        - All energies from part 0
        - FEP back, forward and average dG profiles vs lambda
        - Sampling profiles
        - LRA contributions (statistics)
        - Free energy profiles vs Egap (bin-averaged)
        - Coefficients vs Egap (part3)
        """

        plots = ODict()

        # no QFepOutput objects (all failed to parse)
        if not self.qfos:
            return plots

        # make PlotData objects
        plots["dgde"] = PlotData("Free energy profile",
                                 xlabel="E1-E2  [kcal/mol]",
                                 ylabel="Free energy  [kcal/mol]")
        if self._lra_lambdas:
            l1, l2 = self._lra_lambdas 
            lra_de_st1 = "lra_de_st1_{}".format(l1)
            lra_de_st2 = "lra_de_st2_{}".format(l2)
            lra_lra = "lra_lra_{}{}".format(l1, l2)
            lra_reo = "lra_reo_{}{}".format(l1, l2)
            plots[lra_de_st1] = PlotData("E2-E1 (lambda={})".format(l1),
                                         xlabel="Energy type",
                                         ylabel="Potential energy  [kcal/mol]",
                                         plot_type="bar")
            plots[lra_de_st2] = PlotData("E2-E1 (lambda={})".format(l2),
                                         xlabel="Energy type",
                                         ylabel="Potential energy  [kcal/mol]",
                                         plot_type="bar")
            plots[lra_lra] = PlotData("LRA (l={} -> l={})".format(l1, l2),
                                      xlabel="Energy type",
                                      ylabel="Potential energy  [kcal/mol]",
                                      plot_type="bar")
            plots[lra_reo] = PlotData("Reorganization energy (l={} -> "
                                      "l={})".format(l1, l2),
                                      xlabel="Energy type",
                                      ylabel="Potential energy  [kcal/mol]",
                                      plot_type="bar")

        plots["lambda_egap"] = PlotData("Sampling (binning): "
                                        "Check the overlap between lambda "
                                        "frames in each bin",
                                        xlabel="Egap [kcal/mol]",
                                        ylabel="Lambda",
                                        plot_type="scatter")
        plots["pts_egap"] = PlotData("Sampling (total counts): "
                                     "Check for breaks.",
                                     xlabel="Egap [kcal/mol]",
                                     ylabel="Number of points")
        plots["pts_egap_hists"] = PlotData("Sampling (histograms, 1st output "
                                           "only): Check overlap ",
                                           xlabel="Egap",
                                           ylabel="Number of points,")
        plots["pts_egap_l"] = PlotData("Sampling3D (1st output only)",
                                       xlabel="Egap",
                                       ylabel="Lambda",
                                       zlabel="Number of points",
                                       plot_type="wireframe")
        plots["dgl"] = PlotData("dG vs Lambda",
                                xlabel="Lambda",
                                ylabel="Free energy  [kcal/mol]")
        plots["dgl_forw"] = PlotData("dG vs Lambda (forward)",
                                     xlabel="Lambda",
                                     ylabel="Free energy  [kcal/mol]")
        plots["dgl_rev"] = PlotData("dG vs Lambda (reverse)",
                                    xlabel="Lambda",
                                    ylabel="Free energy  [kcal/mol]")
        plots["rxy"] = PlotData("Reactive distance",
                                xlabel="E1-E2  [kcal/mol]",
                                ylabel=u"Rxy  [Å]")

        # get the column names from the first output (0th is lambda)
        qfo0 = self.qfos.values()[0]
        evb_states = qfo0.header.nstates
        part0_coltitles = qfo0.part0.data_state[0].get_column_titles()

        for col in part0_coltitles[4:]:
            for evb_state in range(evb_states):
                est = evb_state + 1
                key = "e{}l_{}".format(est, col)
                plots[key] = PlotData("E{} vs Lambda ({})".format(est, col),
                                      xlabel="Lambda (state {})".format(est),
                                      ylabel="E{} ({})  [kcal/mol]"
                                             "".format(est, col))
                key = "e{}l_{}".format(est, col)
                plots[key] = PlotData("E{} vs Lambda ({})".format(est, col),
                                      xlabel="Lambda (state {})".format(est),
                                      ylabel="E{} ({})  [kcal/mol]"
                                             "".format(est, col))

        # populate PlotData subplots (each output is a subplot)
        for qfo_path, qfo in self.qfos.iteritems():

            relp = os.path.relpath(qfo_path)

            # Part 0 energies
            for evb_state in range(evb_states):
                est = evb_state + 1
                data = qfo.part0.data_state[evb_state].get_columns()
                for i, colname in enumerate(part0_coltitles[4:]):
                    key = "e{}l_{}".format(est, colname)
                    # 3rd column is lambda, 4,5,6,7.. are energies
                    plots[key].add_subplot(relp, data[3], data[i+4])


            # Part 1 FEP
            data = qfo.part1.data.get_columns(["Lambda", "sum_dGf",
                                               "sum_dGr", "dG"])
            plots["dgl_forw"].add_subplot(relp, data[0], data[1])

            plots["dgl_rev"].add_subplot(relp, data[0], data[2])
            plots["dgl"].add_subplot(relp, data[0], data[3])


            # Part 2 (sampling/binning)
            data = qfo.part2.data.get_columns(["Lambda", "Egap", "points"])
            plots["lambda_egap"].add_subplot(relp, data[1], data[0])

            ## use only the first one, too much data otherwise
            if not plots["pts_egap_hists"].subplots:
                rows = zip(*data) #transpose columns to rows
                for l in sorted(set(data[0])):
                    rows_f = [(eg, pts) for lam, eg, pts in rows if lam == l]
                    eg, pts = zip(*rows_f) #transpose rows to columns

                    plots["pts_egap_hists"].add_subplot("{}_{}".format(relp, l),
                                                        eg, pts)
            ## use only the first one, too much data otherwise
            if not plots["pts_egap_l"].subplots:
                plots["pts_egap_l"].add_subplot(relp, data[1], data[0], data[2])


            # Part 3
            data = qfo.part3.data.get_columns(["Egap", "dGg_norm",
                                               "r_xy", "points"])
            plots["dgde"].add_subplot(relp, data[0], data[1])
            plots["rxy"].add_subplot(relp, data[0], data[2])
            plots["pts_egap"].add_subplot(relp, data[0], data[3])

        if self.lras:
            data = self.lra_stats.get_columns()
            plots[lra_de_st1].add_subplot("average", data[0],
                                          data[1], yerror=data[2])
            plots[lra_de_st2].add_subplot("average", data[0],
                                          data[3], yerror=data[4])
            plots[lra_lra].add_subplot("average", data[0],
                                       data[5], yerror=data[6])
            plots[lra_reo].add_subplot("average", data[0],
                                       data[7], yerror=data[8])

        return plots



    @property
    def lra_stats(self):
        """Calculate average and st.dev of LRA and reorg energies."""

        average_lras = DataContainer(["E_type", "(E2-E1)_10_mean",
                                      "(E2-E1)_10_std", "(E2-E1)_01_mean",
                                      "(E2-E1)_01_std", "LRA_mean", "LRA_std",
                                      "REORG_mean", "REORG_std"])

        allvals = []
        for lra in self.lras.values():
            rows = lra.get_rows()
            for irow, row in enumerate(rows):
                try:
                    allvals[irow].append(row)
                except IndexError:
                    allvals.append([row,])

    # allvals now looks like this:
    # [
    #   [
    #     ["EQtot", EQtot_de_st1_1, EQtot_de_st2_1, EQtot_lra_1, EQtot_reorg_1],
    #     ["EQtot", EQtot_de_st1_2, EQtot_de_st2_2, ...], ...
    #   ],
    #   [
    #     ["EQbond", EQbond_de_st1_1, EQbond_de_st2_1, EQbond_lra_1, EQbond_reorg_1],
    #     ["EQbond", EQbond_de_st1_2, EQbond_de_st2_2, ...], ...
    #   ]
    # ]
    #
        for values in allvals:
            # transpose to get [ ["EQtot","EQtot"...],
            #                    [ EQtot_de_st1_1, EQtot_de_st1_2,...],
            #                    [ EQtot_de_st2_1, EQtot_de_st2_2,...], ...]

            values = zip(*values)  
            # now they can be easily averaged and std-ed
            e_type = values[0][0]
            de_st1_mean = np.mean(values[1])
            de_st2_mean = np.mean(values[2])
            lra_mean = np.mean(values[3])
            reo_mean = np.mean(values[4])
            de_st1_std = np.std(values[1])
            de_st2_std = np.std(values[2])
            lra_std = np.std(values[3])
            reo_std = np.std(values[4])

            average_lras.add_row([e_type, de_st1_mean, de_st1_std,
                                  de_st2_mean, de_st2_std, lra_mean,
                                  lra_std, reo_mean, reo_std])

        return average_lras






# TODO: refactor most of the code below to qdyn.QDynOutput
# (same principle as qfep.QFepOutput)
class QAnalyseDyns(object):
    def __init__(self, logfiles, timeunit="ps", stepsize=None):
        """
        Wrapper class for QanalyseDyn for analysing a sequence of log files.
        Args:
           logfile (list):  paths/filenames of Q logfiles
           timeunit (string):  fs,ps,ns (optional, default is ps)
           stepsize (float):  in case the on in Q is 0.000 (Q printout is a work of art)

        Usage:

        qads = QAnalyseDyns(["fep_000.log", "fep_001.log", ...], timeunit="ps")

        # get average temperatures by combining all logs and skipping 10% in each one
        temps_dc = qads.get_temps(percent_skip=10)   # returns Datacontainer with all logs combined
        temps = temps_dc.get_columns()   # returns the columns
        coltitles = temps_dc.get_column_titles()       # "Time", "T_tot", "T_free", "T_free_solute", "T_free_solvent"

        for i,colt in coltitles[1:]:
            print colt, np.mean( [ x for j,x in enumerate(temps[i]) if temps[0][j] >= midpoint ] )

        """
        if not logfiles:
            raise QAnalyseDynError("No logfiles given")
        self.analysed = []
        starttime = 0
        for i,logfile in enumerate(logfiles):
            try:
                qad = _QAnalyseDyn(logfile, timeunit, stepsize=stepsize, starttime=starttime)
                self.analysed.append(qad)
            except QAnalyseDynError as e:
                raise QAnalyseDynError("%s: %s" % (logfile, str(e)))
            starttime = qad.get_endtime()
        self.n_evb_states = self.analysed[0]._evb_states
        self.en_section_keys = self.analysed[0].map_en_section.keys()
        self.qen_section_keys = self.analysed[0].map_qen_section.keys()


    def get_temps(self, percent_skip=0, stride=1):
        """
        Get temperatures from all logfiles combined.
        Args:
           percent_skip (int, optional):  percent of datapoints in each
                                          logfile to skip, default=0
           stride (int, optional):  use only every Nth point, default=1
        Returns:
           temperatures (DataContainer)
        """

        # "Time", "T_tot", "T_free", "T_free_solute", "T_free_solvent"
        cts = list(self.analysed[0].data_temp.get_column_titles())
        temps = DataContainer(cts)

        for qad in self.analysed:
            rows = qad.data_temp.get_rows()
            skip = int(round(len(rows)*percent_skip/100.0))
            for row in rows[skip::stride]:
                temps.add_row(row)
        return temps



    def get_temp_stats(self, percent_skip=0, stride=1):
        """
        Returns temperature stats in string format (used for cmdline printout)
        for all logfiles combined
        """
        temps = self.get_temps(percent_skip=percent_skip, stride=stride)
        tt, tf, tf_solu, tf_solv = temps.get_columns(("T_tot", "T_free",
                                                      "T_free_solute",
                                                      "T_free_solvent"))
        tt_mean, tt_std = np.mean(tt), np.std(tt)
        tf_mean, tf_std = np.mean(tf), np.std(tf)
        tf_solu_mean, tf_solu_std = np.mean(tf_solu), np.std(tf_solu)
        tf_solv_mean, tf_solv_std = np.mean(tf_solv), np.std(tf_solv)
        tt_max_dev = max(map(lambda x: abs(x - tt_mean), tt))
        tf_max_dev = max(map(lambda x: abs(x - tf_mean), tf))
        tf_solu_max_dev = max(map(lambda x: abs(x - tf_solu_mean), tf_solu))
        tf_solv_max_dev = max(map(lambda x: abs(x - tf_solv_mean), tf_solv))

        outstr = """\
Temperature stats:
{0:20s}{1:>20s}{2:>20s}{3:>20s}
{4:20s}{5:>20.2f}{6:>20.2f}{7:>20.2f}
{8:20s}{9:>20.2f}{10:>20.2f}{11:>20.2f}
{12:20s}{13:>20.2f}{14:>20.2f}{15:>20.2f}
{16:20s}{17:>20.2f}{18:>20.2f}{19:>20.2f}
""".format("", "Mean", "Stdev", "Max.Abs.Dev.",
           "T_total", tt_mean, tt_std, tt_max_dev,
           "T_free", tf_mean, tf_std, tf_max_dev,
           "T_free_solute", tf_solu_mean, tf_solu_std, tf_solu_max_dev,
           "T_free_solvent", tf_solv_mean, tf_solv_std, tf_solv_max_dev)

        return outstr



    def get_offdiags(self, percent_skip=0, stride=1):
        """
        Get distances from all logfiles combined.
        Args:
           percent_skip (int, optional):  percent of datapoints in each
                                          logfile to skip, default=0
           stride (int, optional):  use only every Nth point, default=1
        Returns:
           distances (dict):  e.g: { "13_31": DataContainer,
                                     "13_18": DataContainer }
        """

        coltitles = list(self.analysed[0].data_offdiags.get_column_titles())
        dists = DataContainer(coltitles)

        for qad in self.analysed:
            rows = qad.data_offdiags.get_rows()
            skip = int(round(len(rows)*percent_skip/100.0))
            for row in rows[skip::stride]:
                dists.add_row(row)
        return dists


    def get_energies(self, e_type, percent_skip=0, stride=1):
        """
        Get energies from all logfiles combined.
        Args:
           e_type (string):  keys in QAnalyseDyn.map_en_section dictionary
           percent_skip (int, optional):  percent of datapoints in each
                                          logfile to skip, default=0
           stride (int, optional):  use only every Nth point, default=1
        Returns:
           energies (DataContainer)
        """

        cts = list(self.analysed[0].map_en_section[e_type].get_column_titles())
        energies = DataContainer(cts)

        for qad in self.analysed:
            rows = qad.map_en_section[ e_type ].get_rows()
            skip = int(round(len(rows)*percent_skip/100.0))
            for row in rows[skip::stride]:
                energies.add_row(row)
        return energies


    def get_q_energies(self, qe_type, evb_state, percent_skip=0, stride=1):
        """
        Get Q energies from all logfiles combined.
        Args:
           qe_type (string):  keys in QAnalyseDyn.map_qen_section dictionary
           evb_state (int):  1 or 2 or 3...
           percent_skip (int, optional):  percent of datapoints in each
                                          logfile to skip, default=0
           stride (int, optional):  use only every Nth point, default=1
        Returns:
           energies (DataContainer)
        """

        cts = list(self.analysed[0].map_qen_section[qe_type][evb_state-1]\
                                                        .get_column_titles())
        energies = DataContainer(cts)

        for qad in self.analysed:
            rows = qad.map_qen_section[qe_type][evb_state-1].get_rows()
            skip = int(round(len(rows)*percent_skip/100.0))
            for row in rows[skip::stride]:
                energies.add_row(row)
        return energies




class _QAnalyseDyn(object):
    def __init__(self, logfile, timeunit="ps", stepsize=None, starttime=0):
        """
        Parses a Q dynamics logfile and extracts data (temperature, energies...)
        For interfacing, use QAnalyseDyns.

        Args:
           logfile (string):  path/filename of Q logfile
           timeunit (string):  fs,ps,ns (optional, default is ps)
           stepsize (float):  in case the one in Q is 0.000 (Q printout is a work of art)


        Usage looks like this:

        # parse
        qad = QAnalyseDyns(.....).analysed[0]

        # print out nicely formatted temperature stats
        print qad.get_temp_stats()

        # get averages for seconds half (step >= 50% of steps) of all the temperatures
        temps = qad.data_temp.get_columns()

        coltitles = qad.data_temp.get_column_titles()
        # [ "Time", "T_tot", "T_free", "T_free_solute", "T_free_solvent" ]

        midpoint = int(temps[0][-1])/2        # 0 == "Time", -1 == last frame
        for i,colt in coltitles[1:]:
            print colt, np.mean( [ x for j,x in enumerate(temps[i]) if temps[0][j] >= midpoint ] )

        # get the potential energy data and just print it out
        Epot = qad.data_E_SUM.get_columns( ["Time", "Potential"] )
        print Epot

        """

        # parse the logfile:
        # first the header using RE
        # then dynamics (_parse_dyn()) line by line using the lazy generator in 'open' (less memory consumption and faster than regular expressions)

        self._logfile = logfile
        self._starttime = starttime

        self.MAP_TIME = {"fs": 1.0, "ps": 1e-3, "ns": 1e-6}
        if timeunit not in self.MAP_TIME:
            raise QAnalyseDynError("Timeunit has to be either 'fs', 'ps' or 'ns'")
        self._timeconv = self.MAP_TIME[timeunit]

        self._header=""
        try:
            with open(self._logfile,'r') as lf:
                for line in lf:
                    self._header += line
                    if "Initialising dynamics" in line:
                        break
        except IOError as e:
            raise QAnalyseDynError("Could not read the logfile: " + str(e))

        # use RE to get some info about the simulations
        m = re.search("Build number\s*([\d\.]+)", self._header)
        if m: self._qversion = m.group(1)
        else:
            m = re.search('QDyn version 5.06', self._header)
            if m:
                self._qversion = '5.06'
            else:
                raise QAnalyseDynError("Not a valid Q log file or Q version "
                                       "is very old...")

        m = re.search("Topology file      =\s*(\S+)", self._header)
        if m: self._topfile = m.group(1)
        else: raise QAnalyseDynError("Couldn't find the topology filename!?")

        m = re.search("Number of MD steps =\s*(\d+)", self._header)
        if m: self._md_steps = int(m.group(1))
        else: raise QAnalyseDynError("Couldn't find number of steps!?")

        m = re.search("Stepsize \(fs\)    =\s*([\d\.]+)", self._header)
        if m: self._stepsize = float(m.group(1))
        else: raise QAnalyseDynError("Couldn't find the stepsize!?")

        if not stepsize:
            if abs(self._stepsize - 0.0) < 1e-8:
                raise QAnalyseDynError("Can't convert steps to time, stepsize "
                                       "is 0.0 in the logfile (Q sucks). Set "
                                       "the stepsize please.")
        else:
            if self._stepsize:
                raise QAnalyseDynError("Will not override the non-zero "
                                       "stepsize in the logfile...")
            else:
                self._stepsize = stepsize

        m = re.search("FEP input file     =\s*(\S+)", self._header)
        if m: self._fepfile = m.group(1)
        else: self._fepfile = None

        if self._fepfile:
            m = re.search("No. of fep/evb states    =\s*(\d+)", self._header)
            if m: self._evb_states = int(m.group(1))
            else: raise QAnalyseDynError("Couldn't find the number of states!?")

        offdsection = re.search("(No. of offdiagonal \(Hij\) functions =.*?^$)",
                                self._header,
                                re.MULTILINE | re.DOTALL).group(1)
        offdgs = re.findall("\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+[\d\.]+\s+[\d\.]+",
                            offdsection)

        #
        # make datacontainer variables for storing all the data

        # offdiags
        offdiags = ["{}_{}".format(a1, a2) for a1, a2 in offdgs]
        self._tmp_offdiags = {}
        for k in offdiags:
            self._tmp_offdiags[k] = DataContainer(["Time", "Distance"])

        # temperature
        self.data_temp = DataContainer(["Time", "T_tot", "T_free",
                                        "T_free_solute", "T_free_solvent"])

        # energies
        self.data_E_solute = DataContainer(["Time", "El", "VdW", "Bond",
                                            "Angle", "Torsion", "Improper"])
        self.data_E_solvent = DataContainer(["Time", "El", "VdW", "Bond",
                                             "Angle", "Torsion", "Improper"])
        self.data_E_solute_solvent = DataContainer(["Time", "El", "VdW"])
        self.data_E_LRF = DataContainer(["Time", "El"])
        self.data_E_Q_atom = DataContainer(["Time", "El", "VdW", "Bond",
                                            "Angle", "Torsion", "Improper"])
        self.data_E_restraints = DataContainer(["Time", "Total", "Fix",
                                                "Solvent_rad", "Solvent_pol",
                                                "Shell", "Solute"])
        self.data_E_SUM = DataContainer(["Time", "Total",
                                         "Potential", "Kinetic"])

        # Q energies
        q_columns1 = ("Time", "Lambda", "El", "VdW")
        q_columns2 = ("Time", "Lambda", "El", "VdW", "Bond",
                      "Angle", "Torsion", "Improper")
        q_columns3 = ("Time", "Lambda", "Total", "Restraint")

        self.data_EQ_Q, self.data_EQ_prot = [], []
        self.data_EQ_wat, self.data_EQ_surr = [], []
        self.data_EQ_any, self.data_EQ_SUM = [], []
        for i in range(self._evb_states):
            self.data_EQ_Q.append(DataContainer(q_columns1))
            self.data_EQ_prot.append(DataContainer(q_columns1))
            self.data_EQ_wat.append(DataContainer(q_columns1))
            self.data_EQ_surr.append(DataContainer(q_columns1))
            self.data_EQ_any.append(DataContainer(q_columns2))
            self.data_EQ_SUM.append(DataContainer(q_columns3))

        # mapping of energy types (label in the output) with containers
        self.map_en_section = {"solute": self.data_E_solute,
                               "solvent": self.data_E_solvent,
                               "solute-solvent": self.data_E_solute_solvent,
                               "LRF": self.data_E_LRF,
                               "Q-atom": self.data_E_Q_atom,
                               "SUM": self.data_E_SUM}

        self.map_qen_section = {"Q-Q": self.data_EQ_Q,
                                "Q-prot": self.data_EQ_prot,
                                "Q-wat": self.data_EQ_wat,
                                "Q-surr.": self.data_EQ_surr,
                                "Q-any": self.data_EQ_any,
                                "Q-SUM": self.data_EQ_SUM}

        self._parse_dyn()
        d_dcs = self._tmp_offdiags.values()
        cts = ["Time",] + self._tmp_offdiags.keys()
        self.data_offdiags = DataContainer(cts)
        # TODO: clean up this magic below
        for d_row in zip(* [d_dcs[0].get_columns([0,])[0],]  + [d_dc.get_columns([1,])[0] for d_dc in d_dcs ]):
            self.data_offdiags.add_row(d_row)

        self._endtime = self._md_steps * self._stepsize * self._timeconv \
                      + self._starttime


    def get_endtime(self):
        return self._endtime


    def _parse_dyn(self):
        """
        Parses the dynamics part of the logfile (used by init)
        """
        time = self._starttime
        t_free, t_tot = None, None
        insection = False
        with open(self._logfile, 'r') as logfile:
            logfile.seek(len(self._header))
            for line in logfile:
                lf = line.split()
                if not lf:
                    continue
                if "Initialising dynamics" in line:
                    raise QAnalyseDynError("Found more than one logfile...",
                                           "Don't concatenate man...")

                if t_free != None: # second line with temps
                    try:
                        tf_solute = float(lf[1])
                    except: # gas phase
                        tf_solute = 0
                    try:
                        tf_solvent = float(lf[3])
                    except: # gas phase
                        tf_solvent = 0

                    self.data_temp.add_row((time, t_tot, t_free,
                                            tf_solute, tf_solvent))
                    t_free, t_tot = None, None

                # first line with temps
                elif "Temperature at step" in line:
                    # fix for large step numbers
                    line = line.replace("step", "step ")
                    lf = line.split()
                    step = int(lf[3].strip(":"))
                    time = step * self._stepsize * self._timeconv \
                         + self._starttime
                    t_tot, t_free = float(lf[5]), float(lf[7])

                elif "Energy summary at step" in line or \
                        "Q-atom energies at step" in line:
                    insection = True
                    step = int(lf[5])
                    time = step * self._stepsize * self._timeconv \
                         + self._starttime

                elif "FINAL  Energy summary" in line or \
                        "FINAL Q-atom energies" in line:
                    insection = True
                    time = self._md_steps * self._stepsize * self._timeconv \
                         + self._starttime

                elif "===================================================="\
                     "======================" in line:
                    insection = False

                if not insection:
                    continue

                # always skip the 0th step
                if step == 0:
                    continue

                key = lf[0]
                if key in self.map_en_section:
                    row = [time,] + [float(x) for x in lf[1:]]
                    self.map_en_section[key].add_row(row)
                elif key in self.map_qen_section:
                    evb_index = int(lf[1]) - 1
                    row = [time,] + [float(x) for x in lf[2:]]
                    self.map_qen_section[key][evb_index].add_row(row)
                elif "dist. between" in line:
                    atom1, atom2, dist = lf[8], lf[9], float(lf[11])
                    k = "{}_{}".format(atom1, atom2)
                    self._tmp_offdiags[k].add_row([time, dist])




