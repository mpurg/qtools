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
"""
Module for analysing Qfep and Qdyn output (logfiles).
Two main classes are QAnalyseFeps and QAnalyseDyns.
See q_analysefeps.py and q_analysedyns.py for usage.
"""


from __future__ import absolute_import, unicode_literals
from io import open

import os
import re
import logging
from collections import OrderedDict as ODict

from Qpyl.core.qfep import QFepOutput, QFepOutputError
from Qpyl.core.qdyn import QDynOutput, QDynOutputError
from Qpyl.common import DataContainer, np
from Qpyl.plotdata import PlotData
import six
from six.moves import range
from six.moves import zip

logger = logging.getLogger(__name__)

class QAnalyseFeps(object):
    """Wrapper class for QFepOutput for analysing multiple outputs.

    Calculates statistics (mean, stdev,...), makes pretty summaries,
    extracts all data to PlotData objects via neat functions.

    Args:
        qfep_outputs (list of strings or list of tuples(string, string)): \
                    pathnames of qfep outputs or (name, output_string) tuples

        lra_lambdas (list of floats, optional):  lambdas at which LRA energies\
                                                 should be calculated at.\
                                                 Usual case is (1.0, 0.0).\
                                                 Default = None (no LRA calcs)
        _parent (QAnalyseFeps): used internally for QCP and exclusions
        _subcalcs (string): used internally for QCP and exclusions 

    Examples:
        >>> qos = [os.path.join(md, "qfep.log") for md in sorted(mapdirs)]
        >>> qaf = QAnalyseFeps(qos, lra_lambdas=(1.0, 0.0))
        >>> print qaf.stats_str
        # prints statistics
        >>> for sub_calc_key, sub_calc in sorted(qaf.sub_calcs.iteritems()):
        >>>     print sub_calc.stats_str
        # prints subcalculation (QCP, GE) statistics
        >>> for failed_path, failed_msg in sorted(qaf.failed.iteritems()):
        >>>     relp = os.path.relpath(failed_path)
        >>>     print "-> {}: {}".format(relp, failed_msg)
        # prints failed

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
                if isinstance(qfep_output, six.string_types):
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
                    raise
                    self.failed[qfep_out_name] = "UNCAUGHT EXCEPTION: {}"\
                                                 "".format(error_msg)
                else:
                    self.qfos[qfep_out_name] = qfo
        # exctract data (exclusion or QCP) from parent's sub_calcs
        else:
            self.qfos_failed = self._parent.failed
            for p_qfo_key, p_qfo in six.iteritems(self._parent.qfos):
                if self._subcalc_key in p_qfo.sub_calcs:
                    self.qfos[p_qfo_key] = p_qfo.sub_calcs[self._subcalc_key]

        # get dG_FEP, dGa, dG0, LRAs
        for qfep_output, qfo in six.iteritems(self.qfos):
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
        dgas = list(self.dgas.values())
        dg0s = list(self.dg0s.values())
        dgs_fep = list(self.dgs_fep.values())

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
        - FEP delta (forward - reverse) vs lambda
        - Binning overlap/histograms
        - Full dGg profile (not bin-averaged) for first output
        - LRA contributions (statistics)
        - Free energy profiles vs Egap (bin-averaged)
        """

        plots = ODict()

        # no QFepOutput objects (all failed to parse)
        if not self.qfos:
            return plots

        # make PlotData objects

        # Part 3 (bin-averaged profile)
        plots["dgde"] = PlotData("Free-energy profile (bin-averaged, norm.)",
                                 xlabel="E1-E2  [kcal/mol]",
                                 ylabel="Free energy  [kcal/mol]")
        plots["dgde_unnorm"] = PlotData("Free-energy profile (bin-averaged)",
                                        xlabel="E1-E2  [kcal/mol]",
                                        ylabel="Free energy  [kcal/mol]")
        plots["rxy"] = PlotData("Off-diagonal distance",
                                xlabel="E1-E2  [kcal/mol]",
                                ylabel=u"Rxy  [Å]")
        plots["pts_egap"] = PlotData("Binning (total counts)",
                                     xlabel="E1-E2 [kcal/mol]",
                                     ylabel="Number of points")

        # Part 2 (dGg/binning)
        plots["egap_lambda"] = PlotData("Binning (overlap)",
                                        xlabel="Lambda (state 1)",
                                        ylabel="E1-E2 [kcal/mol]")
        plots["pts_egap_hists"] = PlotData("Binning (histograms, "
                                           "1st output only)",
                                           xlabel="E1-E2",
                                           ylabel="Number of points,")
        plots["pts_egap_l"] = PlotData("Binning3D (1st output only)",
                                       xlabel="E1-E2",
                                       ylabel="Lambda (state 1)",
                                       zlabel="Number of points",
                                       plot_type="wireframe")
        plots["dgg_egap"] = PlotData("Full dGg profile (1st output only)",
                                     xlabel="E1-E2 [kcal/mol]",
                                     ylabel="Lambda (state 1)",
                                     zlabel="Free energy [kcal/mol]",
                                     plot_type="wireframe")

        # Part 1 (FEP)
        plots["dgl"] = PlotData("FEP (dG_sum)",
                                xlabel="Lambda (state 1)",
                                ylabel="Free energy  [kcal/mol]")
        plots["dgl_delta"] = PlotData("FEP (dGf-dGr)",
                                      xlabel="Lambda (state 1)",
                                      ylabel="Free energy  [kcal/mol]")
        plots["dgl_forw"] = PlotData("FEP (dG_forward)",
                                     xlabel="Lambda (state 1)",
                                     ylabel="Free energy  [kcal/mol]")
        plots["dgl_rev"] = PlotData("FEP (dG_reverse)",
                                    xlabel="Lambda (state 1)",
                                    ylabel="Free energy  [kcal/mol]")

        # Part 0 (Avg. energies / LRAs)
        if self._lra_lambdas:
            l1, l2 = self._lra_lambdas 
            lra_de_st1 = "lra_de_st1_{}".format(l1)
            lra_de_st2 = "lra_de_st2_{}".format(l2)
            lra_lra = "lra_lra_{}{}".format(l1, l2)
            lra_reo = "lra_reo_{}{}".format(l1, l2)
            plots[lra_de_st1] = PlotData("E2-E1 (lambda={})".format(l1),
                                         xlabel="Energy type",
                                         ylabel="Free energy  [kcal/mol]",
                                         plot_type="bar")
            plots[lra_de_st2] = PlotData("E2-E1 (lambda={})".format(l2),
                                         xlabel="Energy type",
                                         ylabel="Free energy  [kcal/mol]",
                                         plot_type="bar")
            plots[lra_lra] = PlotData("LRA (l={} -> l={})".format(l1, l2),
                                      xlabel="Energy type",
                                      ylabel="Free energy  [kcal/mol]",
                                      plot_type="bar")
            plots[lra_reo] = PlotData("Reorganization energy (l={} -> "
                                      "l={})".format(l1, l2),
                                      xlabel="Energy type",
                                      ylabel="Free energy  [kcal/mol]",
                                      plot_type="bar")

        # get the column names from the first output (0th is lambda)
        qfo0 = list(self.qfos.values())[0]
        evb_states = qfo0.header.nstates
        part0_coltitles = qfo0.part0.data_state[0].column_titles

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
        for qfo_path, qfo in six.iteritems(self.qfos):

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
            data = qfo.part1.data.get_columns(["Lambda", "dGf",
                                               "dGr", "dG"])

            delta = [0,]
            for dgf, dgb in zip(data[1][1:], data[2][:-1]):
                dg = abs(dgf)-abs(dgb)
                delta.append(dg)

            plots["dgl_delta"].add_subplot(relp, data[0], delta)
            plots["dgl_forw"].add_subplot(relp, data[0], data[1])
            plots["dgl_rev"].add_subplot(relp, data[0], data[2])
            plots["dgl"].add_subplot(relp, data[0], data[3])


            # Part 2 (sampling/binning)
            data = qfo.part2.data.get_columns(["Lambda", "Egap",
                                               "dGg", "points"])
            plots["egap_lambda"].add_subplot(relp, data[0], data[1])

            ## use only the first one, too much data otherwise
            if not plots["pts_egap_hists"].subplots:
                rows = list(zip(*data)) #transpose columns to rows
                for l in sorted(set(data[0])):
                    rows_f = [(eg, pts) for lam, eg, dgg, pts in rows if lam == l]
                    eg, pts = list(zip(*rows_f)) #transpose rows to columns

                    plots["pts_egap_hists"].add_subplot("{}_{}".format(relp, l),
                                                        eg, pts)

                plots["dgg_egap"].add_subplot(relp, data[1], data[0], data[2])
                plots["pts_egap_l"].add_subplot(relp, data[1], data[0], data[3])


            # Part 3
            data = qfo.part3.data.get_columns(["Egap", "dGg", "dGg_norm",
                                               "r_xy", "points"])
            plots["dgde"].add_subplot(relp, data[0], data[2])
            plots["dgde_unnorm"].add_subplot(relp, data[0], data[1])
            plots["rxy"].add_subplot(relp, data[0], data[3])
            plots["pts_egap"].add_subplot(relp, data[0], data[4])

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

            values = list(zip(*values))  
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



class QAnalyseDynsError(Exception):
    pass

class QAnalyseDyns(object):
    """Wrapper class for QanalyseDyn for parsing a sequence of Qdyn outputs.

    Args:
        qdyn_outputs (list): Qdyn output pathnames
        time_unit (string): fs,ps,ns (optional, default is ps)
        step_size (float): use in case the output reads 0.000

    Examples:
        >>> qads = QAnalyseDyns(["fep_000.log", "fep_001.log"], timeunit="ps")
        # get average temperatures by combining all outputs
        >>> print qads.get_temp_stats()

    Raises QAnalyseDynsError  (if it fails to parse an output)

    """ 
    def __init__(self, qdyn_outputs, time_unit="ps", step_size=None):
        self.analysed = []
        start_time = 0
        self.time_unit = time_unit

        for i, output in enumerate(qdyn_outputs):
            try:
                qdo = QDynOutput(output, time_unit,
                                 step_size=step_size,
                                 start_time=start_time)
                self.analysed.append(qdo)
            except QDynOutputError as e:
                raise QAnalyseDynsError("{}: {}".format(output, e))

            start_time = qdo.time_end

        self.n_evb_states = self.analysed[0].header.nstates
        self.en_section_keys = list(self.analysed[0].map_en_section.keys())
        self.qen_section_keys = list(self.analysed[0].map_qen_section.keys())


    def get_temps(self, stride=1):
        """Get temperatures from all logfiles combined.

        Args:
            stride (int, optional):  use only every Nth point, default=1

        Returns:
            temperatures (DataContainer)
        """

        # "Time", "T_tot", "T_free", "T_free_solute", "T_free_solvent"
        cts = self.analysed[0].data_temp.column_titles
        temps = DataContainer(cts)

        for qdo in self.analysed:
            rows = qdo.data_temp.get_rows()
            for row in rows[::stride]:
                temps.add_row(row)
        return temps



    def get_temp_stats(self, stride=1):
        """Returns temperature stats in string format (used for cmdline printout)
        for all logfiles combined
        """
        temps = self.get_temps(stride=stride)
        tt, tf, tf_solu, tf_solv = temps.get_columns(("T_tot", "T_free",
                                                      "T_free_solute",
                                                      "T_free_solvent"))
        tt_mean, tt_std = np.mean(tt), np.std(tt)
        tf_mean, tf_std = np.mean(tf), np.std(tf)
        tf_solu_mean, tf_solu_std = np.mean(tf_solu), np.std(tf_solu)
        tf_solv_mean, tf_solv_std = np.mean(tf_solv), np.std(tf_solv)
        tt_max_dev = max([abs(x - tt_mean) for x in tt])
        tf_max_dev = max([abs(x - tf_mean) for x in tf])
        tf_solu_max_dev = max([abs(x - tf_solu_mean) for x in tf_solu])
        tf_solv_max_dev = max([abs(x - tf_solv_mean) for x in tf_solv])

        outstr = """\
Temperature stats:
{0:20s}{1:>15s}{2:>15s}{3:>15s}
{4:20s}{5:>15.2f}{6:>15.2f}{7:>15.2f}
{8:20s}{9:>15.2f}{10:>15.2f}{11:>15.2f}
{12:20s}{13:>15.2f}{14:>15.2f}{15:>15.2f}
{16:20s}{17:>15.2f}{18:>15.2f}{19:>15.2f}
""".format("", "Mean", "Stdev", "Max.Abs.Dev.",
           "T_total", tt_mean, tt_std, tt_max_dev,
           "T_free", tf_mean, tf_std, tf_max_dev,
           "T_free_solute", tf_solu_mean, tf_solu_std, tf_solu_max_dev,
           "T_free_solvent", tf_solv_mean, tf_solv_std, tf_solv_max_dev)

        return outstr



    def get_offdiags(self, stride=1):
        """Get distances from all logfiles combined.

        Args:
           stride (int, optional):  use only every Nth point, default=1

        Returns:
           distances (dict):  e.g: { "13_31": DataContainer,
                                     "13_18": DataContainer }
        """

        coltitles = self.analysed[0].data_offdiags.column_titles
        dists = DataContainer(coltitles)

        for qdo in self.analysed:
            rows = qdo.data_offdiags.get_rows()
            for row in rows[::stride]:
                dists.add_row(row)
        return dists


    def get_energies(self, e_type, stride=1):
        """Get energies from all logfiles combined.

        Args:
           e_type (string):  keys in QDynOutput.map_en_section dictionary
           stride (int, optional):  use only every Nth point, default=1

        Returns:
           energies (DataContainer)
        """

        cts = self.analysed[0].map_en_section[e_type].column_titles
        energies = DataContainer(cts)

        for qdo in self.analysed:
            rows = qdo.map_en_section[e_type].get_rows()
            for row in rows[::stride]:
                energies.add_row(row)
        return energies


    def get_q_energies(self, qe_type, evb_state, stride=1):
        """Get Q energies from all logfiles combined.

        Args:
           qe_type (string):  keys in QDynOutput.map_qen_section dictionary
           evb_state (int):  1 or 2 or 3...
           stride (int, optional):  use only every Nth point, default=1

        Returns:
           energies (DataContainer)
        """

        cts = self.analysed[0].map_qen_section[qe_type][evb_state-1]\
                                                        .column_titles
        energies = DataContainer(cts)

        for qdo in self.analysed:
            rows = qdo.map_qen_section[qe_type][evb_state-1].get_rows()
            for row in rows[::stride]:
                energies.add_row(row)
        return energies


    def get_plotdata(self, stride=1):
        """Return 'useful data' as a dictionary of PlotData objects.

        Useful data:
        - Temperatures
        - Offdiagonal distances
        - Energies (Q and non-Q)

        Args:
           stride (int, optional):  use only every Nth point, default=1
        """

        plots = ODict()

        # make PlotData objects

        time_label = "Time [{}]".format(self.time_unit)
        plots = ODict()
        plots["temp"] = PlotData("Temperature",
                                 xlabel=time_label,
                                 ylabel="T [K]")

        plots["offdiags"] = PlotData("Offdiagonal distances",
                                     xlabel=time_label,
                                     ylabel="Distance [A]")

        t_dc = self.get_temps(stride=stride)
        t_cs, t_cts = t_dc.get_columns(), t_dc.column_titles
        for i, t_ct in enumerate(t_cts[1:]):
            plots["temp"].add_subplot(t_ct, t_cs[0], t_cs[i+1]) # 0==Time


        d_dc = self.get_offdiags(stride=stride)
        d_cs, d_cts = d_dc.get_columns(), d_dc.column_titles
        for i, d_ct in enumerate(d_cts[1:]):
            plots["offdiags"].add_subplot(d_ct, d_cs[0], d_cs[i+1]) # 0==Time


        for k in self.en_section_keys:
            key = "E_{}".format(k)
            plots[key] = PlotData("Energy: " + k,
                                  xlabel=time_label,
                                  ylabel="Energy [kcal/mol]")
            e_dc = self.get_energies(k, stride=stride)
            e_cs, e_cts = e_dc.get_columns(), e_dc.column_titles
            if e_cs:
                for i, e_ct in enumerate(e_cts[1:]):
                    plots[key].add_subplot(e_ct, e_cs[0], e_cs[i+1]) # 0==Time


        for k in self.qen_section_keys:
            for evb_state in range(1, self.n_evb_states + 1):
                key = "EQ{}_{}".format(evb_state, k)
                plots[key] = PlotData("Q Energy: {} (state {})"
                                      "".format(k, evb_state),
                                      xlabel=time_label,
                                      ylabel="Energy [kcal/mol]")
                qe_dc = self.get_q_energies(k, evb_state, stride=stride)
                qe_cs, qe_cts = qe_dc.get_columns(), qe_dc.column_titles
                if qe_cs:
                    for i, qe_ct in enumerate(qe_cts[1:]):
                        plots[key].add_subplot(qe_ct, qe_cs[0], qe_cs[i+1])

        return plots



