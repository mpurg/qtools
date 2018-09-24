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

"""
This module contains functions for generating FEP and
equilibration inputs for Qdyn.
Example of procedure files can be found in
qtools/template_examples/
"""

from __future__ import absolute_import
import sys
import os
import copy
import re
import time
import locale
import shutil
import tempfile
import random
import logging
from collections import OrderedDict as ODict

from Qpyl.core.qdyn import QDynInput, QDynInputError
from Qpyl.core.qstructure import QStruct, QStructError, find_placeholders
from Qpyl.common import __version__, raise_or_log
import six
from six.moves import range

logger = logging.getLogger(__name__)


class QGenrelaxError(Exception):
    pass

class QGenfepsError(Exception):
    pass


# TODO: break these two functions into something digestable

def genrelax(relax_proc_file, outdir, restraint,
             top_file=None, fep_file=None, runscript_file=None,
             pdb_file=None, cont_file=None, ignore_errors=False):

    """Generates inputs for an MD simulation with Q (Qdyn).

    Arguments:
        relax_proc_file (string):  genrelax procedure file pathname
        outdir (string):  output directory
        restraint (string):  restraint coordinate (a)

    Optional arguments (b)
        top_file (string):  Q topology pathname
        fep_file (string):  fep file pathname
        runscript_file (string):  slurm/sge run script
        pdb_file (string):  pdb pathname (used to convert placeholders)
        cont_file (string):  pathname of previous Qdyn input (continuation)
        ignore_errors (boolean):  passed to QStruct and QDynInp - write to\
                                  logger instead of raising exceptions on\
                                  non-critical things

    (a) Restraint coordinate can be set to:
    'top' - topology
    'cont_inp' - whatever is defined in cont_file
    'cont_final' - endpoint of previous simulation\
                   (final restart of cont_file)

    (b) top_file and cont_file are mutually exclusive, one of them has to\
    be provided

    """

    # check if files exist
    for k, v in six.iteritems(locals()):
        if k in ["pdb_file", "cont_file", "relax_proc_file",
                 "fep_file", "top_file", "runscript_file", "relax_input"]:
            if v and not os.path.lexists(v):
                raise QGenrelaxError("File '{}' doesn't exist.".format(v))

    if restraint not in ["top", "cont_inp", "cont_final"]:
        raise QGenrelaxError("Argument 'restraint' has to be either "
                             "'cont_inp', 'top' or 'cont_final'")

    # constants
    PREFIX = "relax_"
    DIR = os.path.join(os.getcwd(), outdir)
    if os.path.lexists(DIR):
        raise QGenrelaxError("Directory '{}' exists. Please (re)move it "
                             "or set 'outdir'.".format(DIR))

    TMPDIR = tempfile.mkdtemp()

    header_comment = """\
# Generated with QTools, version {}
# Date: {}
# CWD: {}
# Cmdline: {}
""".format(__version__, time.ctime(), os.getcwd(), " ".join(sys.argv))

    # find and replace placeholders. if not PDB was given to replace them, exit
    relax_proc_str = open(relax_proc_file, 'r').read()
    c = find_placeholders(relax_proc_str)
    if c and not pdb_file:
        raise QGenrelaxError("Found placeholders in proc.file, but no PDB "
                             "was given: {}".format(", ".join(c)))
    elif c:
        logger.info("These placeholders will be replaced with atom indices: {}"
                    "".format(", ".join(c)))
        try:
            qstruct = QStruct(pdb_file, "pdb", ignore_errors=ignore_errors)
            relax_proc_str = qstruct.convert_placeholders(relax_proc_str)
        except QStructError as err_msg:
            raise QGenrelaxError("Failed to replace placeholders: "
                                 "{}".format(err_msg))

    # get topology and fep and others from previous input if given (--cont)
    if cont_file:
        if top_file:
            raise QGenrelaxError("'top_file' and 'cont_file' don't like each "
                                 "other. Difficult to continue with a "
                                 "different topology...")
        try:
            c = QDynInput(open(cont_file, 'r').read(),
                          ignore_errors=ignore_errors)
        except QDynInputError as err_msg:
            raise QGenrelaxError("There is something wrong with the given "
                                 "input file ({}): {}".format(cont_file,
                                                              err_msg))

        cont_files = c.parameters["files"]
        di = os.path.dirname(cont_file)
        top_fn = cont_files["topology"]
        cont_re_fn = cont_files["final"]
        re_fn = "cont_{}".format(cont_re_fn)
        shutil.copy2(os.path.join(di, top_fn), TMPDIR)
        shutil.copy2(os.path.join(di, cont_re_fn),
                     os.path.join(TMPDIR, re_fn))

        if restraint == "cont_inp" and "restraint" in cont_files:
            cont_rest_fn = cont_files["restraint"]
            rest_fn = "cont_{}".format(cont_rest_fn)
        elif restraint == "cont_final":
            cont_rest_fn = cont_re_fn
            rest_fn = "cont_{}.rest".format(cont_rest_fn)
        else:
            rest_fn = None

        if rest_fn:
            shutil.copy2(os.path.join(di, cont_rest_fn),
                         os.path.join(TMPDIR, rest_fn))

        if fep_file:
            logger.warning("Using the fep file '{}', instead of the one "
                           "found in the input".format(fep_file))

            fep_fn = os.path.basename(fep_file)
            shutil.copy2(fep_file, TMPDIR)
        else:
            try:
                fep_fn = cont_files["fep"]
                shutil.copy2(os.path.join(di, fep_fn), TMPDIR)
            except KeyError:
                logger.info("No FEP file found in the input")

    # or take the arguments
    else:
        if not top_file:
            raise QGenrelaxError("Please specify the topology file or "
                                 "a previous input for a continuation run.")

        cont_files = None
        top_fn = os.path.basename(top_file)
        shutil.copy2(top_file, TMPDIR)
        try:
            fep_fn = os.path.basename(fep_file)
            shutil.copy2(fep_file, TMPDIR)
        except AttributeError:
            logger.info("NOTE: No FEP file!")

        if restraint in ["cont_inp", "cont_final"]:
            raise QGenrelaxError("Can't restrain to '{}'. Specify 'cont_file'."
                                 "".format(restraint))
        else:
            rest_fn = None

    logger.info("Restraining to: '{}'".format(rest_fn or 'topology'))

    try:
        shutil.copy2(runscript_file, TMPDIR)
    except AttributeError:
        logger.info("No submission script was given.")


    general_inp = []
    steps_inps = [[],]
    script_vars = {}

    section = ""
    for line in relax_proc_str.split("\n"):
        # remove comments and strip whitespaces.
        line = re.split("#|\!", line)[0].strip()
        # empty lines are useless
        if line == "":
            continue
        # found a section
        if line[0] == "{":
            section = line.strip("{}").lower()
            continue

        if not section:
            raise QGenrelaxError("Failed to parse '{}'... this line - '{}' "
                                 "is not inside any section:"
                                 "".format(relax_proc_file, line))

        if section == "script_vars":
            c = line.split()
            var = c[0]
            value = " ".join(c[1:])
            script_vars[var] = value
        elif section == "general":
            general_inp.append(line)
        elif section == "steps":
            if "__________" in line:
                steps_inps.append([])
            else:
                steps_inps[-1].append(line)

    if "fep_fn" in locals():
        # find and replace atom placeholders in FEP file
        # if no PDB was given to replace them, exit
        fep_tmp = os.path.join(TMPDIR, fep_fn)
        fep_file_str = open(fep_tmp, 'r').read()
        c = find_placeholders(fep_file_str)
        if c and not pdb_file:
            raise QGenfepsError("Found placeholders in FEP file, but no "
                                "PDB was given: {}".format(", ".join(c)))
        elif c:
            logger.info("Replacing FEP file placeholders...")
            try:
                qstruct = QStruct(pdb_file, "pdb", ignore_errors=ignore_errors)
                fep_file_str = qstruct.convert_placeholders(fep_file_str)
            except QStructError as err_msg:
                raise QGenfepsError("Failed to replace placeholders: {}"
                                    "".format(err_msg))
            else:
                open(fep_tmp, 'w').write(fep_file_str)


    # check for steps with no parameters
    # (too many _________ lines)and remove them
    for i in range(len(steps_inps)-1, -1, -1):
        if not steps_inps[i]:
            steps_inps.pop(i)

    # join lists of lines to strings and replace the placeholders
    gen_inp_s = "\n".join(general_inp)
    for placeholder, value in six.iteritems(script_vars):
        gen_inp_s = gen_inp_s.replace(placeholder, value)

    step_inps_s = []
    for i, step_inp in enumerate(steps_inps):
        s = "\n".join(step_inp)
        for placeholder, value in six.iteritems(script_vars):
            s = s.replace(placeholder, value)
        step_inps_s.append(s)

    # make and save the inputs
    steps = []
    overridden_prms_all = []
    step_n = 1
    inp_fns = []  # to store the filenames and use the return value
    for step_inp_s in step_inps_s:
        # create the files section
        final = "{}{:03d}.re".format(PREFIX, step_n)
        dcd = "{}{:03d}.dcd".format(PREFIX, step_n)
        files = {"final"      : final,
                 "trajectory" : dcd,
                 "topology"   : top_fn}
        try:
            files["fep"] = fep_fn
        except NameError:
            pass
        if step_n != 1:
            prev_step = step_n - 1
            files["restart"] = "{}{:03d}.re".format(PREFIX, prev_step)
        elif cont_files:
            files["restart"] = re_fn

        if rest_fn != None:
            files["restraint"] = rest_fn

        try:
            # parse the general input
            inp = QDynInput(gen_inp_s, ignore_errors=ignore_errors)
            # update the general parameters with step input, printout the
            # overriden parms, update the files section
            overridden_prms = inp.update(step_inp_s)
            if overridden_prms:
                overridden_prms_all.append((step_n, ", ".join(
                    ["{}:{}->{}".format(key, value_old, value_new) \
                            for key, (value_old, value_new) in \
                            six.iteritems(overridden_prms)])))

            if "energy" in inp.parameters["intervals"]:
                files["energy"] = "{}{:03d}.en".format(PREFIX, step_n)

            inp.update(parameters={"files": files})

        except QDynInputError as err_msg:
            raise QGenrelaxError("Problem with step no. {}: {}"
                                 "".format(step_n, err_msg))

        # set the random seed
        mdp = inp.parameters["md"]
        if "random_seed" in mdp and int(mdp["random_seed"]) < 1:
            rs = random.randint(1, 1000000)
            inp.update(parameters={"md": {"random_seed": rs}})
            logger.info("Generated random seed in step {}: {}"
                        "".format(step_n, rs))

        # get the input string
        try:
            inpstr = inp.get_string()
        except QDynInputError as err_msg:
            raise QGenrelaxError("Error in step {}: {}"
                                 "".format(step_n, err_msg))

        inpfn = "{}{:03d}.inp".format(PREFIX, step_n)
        inp_fns.append(os.path.join(DIR, inpfn))
        s = header_comment + inpstr
        open(os.path.join(TMPDIR, inpfn), 'w').write(s)

        steps.append(inp)
        step_n += 1

    try:
        shutil.copytree(TMPDIR, DIR)
    except OSError:
        raise QGenrelaxError("Cannot create directory '{}'.".format(DIR))

    # remove temporary directory
    shutil.rmtree(TMPDIR)
    logger.info("Created inputs {}{:03d}.inp - {}{:03d}.inp"
                "".format(PREFIX, 1, PREFIX, len(steps)))

    # print some useful information
    if overridden_prms_all:
        logger.info("Overridden parameters:")
        for step_n, op in overridden_prms_all:
            logger.info("{}: {}".format(step_n, op))

    summary = """
Quick summary
{0:<10} {1:>5} {2:>10} {3:>10} {4:^10} {5:^10} {6:^10} {7:^30} {8:^10} {9:>10}
""".format("Step", "T", "Stepsize", "Steps", "Seq.rest", "Dist.rest",
           "Ang.rest", "Shake", "Rand.Seed", "Data (MB)")
    locale.setlocale(locale.LC_ALL, '')
    restraints = []
    total_time = 0

    # print out how much data this run will produce
    # for this we need the atom count from the topology
    for line in open(os.path.join(DIR, top_fn), 'r').readlines(1024):
        if "no. of atoms, no. of solute atoms" in line:
            num_atoms_all = int(line.strip().split()[0])
            break

    REST_B_PER_ATOM = 48.0
    TRJ_B_PER_ATOM = 12.0
    # very rough estimate, depends on Q version
    # it can double if group_contributions are calculated
    EN_B_PER_STEP = 370.0
    CONV_MB = 2**20
    # very rough estimate
    OUT_B_PER_STEP = 2000
    TEMP_B_PER_STEP = 160
    NB_B_PER_STEP = 80

    # intervals mapping: q_parameter_key, q_default_value, approx_bytes_per_frame
    qintervals = {"trj"  : ["trajectory",   100, num_atoms_all*TRJ_B_PER_ATOM],
                  "log"  : ["output",        10, OUT_B_PER_STEP],
                  "temp" : ["temperature",   10, TEMP_B_PER_STEP],
                  "en"   : ["energy",        10, EN_B_PER_STEP],
                  "nb"   : ["non_bond",      10, NB_B_PER_STEP]}

    total_data = {"trj": 0, "log": 0, "en": 0, "rest": 0}

    for i, step in enumerate(steps):
        nstep = i+1
        try:

            # get md parameters
            mdparms = step.parameters["md"]
            total_time += float(mdparms["stepsize"])*int(mdparms["steps"])
            random_seed = mdparms.get("random_seed", "")

            # get restraints
            step_rests = {"sequence_restraints" : [],
                          "distance_restraints" : [],
                          "angle_restraints" : []}
            for rest_type in step_rests.keys():
                for seqrest in step.parameters.get(rest_type, []):
                    if seqrest in restraints:
                        rest_num = restraints.index(seqrest)+1
                        step_rests[rest_type].append(str(rest_num))
                    else:
                        restraints.append(seqrest)
                        step_rests[rest_type].append(str(len(restraints)))
            seq = ",".join(step_rests["sequence_restraints"])
            dist = ",".join(step_rests["distance_restraints"])
            angle = ",".join(step_rests["angle_restraints"])

            # get shake parameters
            shake = []
            # this is a Q default, hopefully it will not change
            if mdparms.get("shake_solvent", "on") == "on":
                shake.append("solvent")
            if mdparms.get("shake_hydrogens", "off") == "on":
                shake.append("hydrogens")
            if mdparms.get("shake_solute", "off") == "on":
                shake.append("solute")
            shake = ",".join(shake)

            # calculate approx amount of data
            data = {}
            mdsteps = int(mdparms["steps"])
            for k, v in six.iteritems(qintervals):
                interval_key = v[0]
                default_interval = v[1]
                bytes_per_step = v[2]
                try:
                    interval = int(step.parameters["intervals"][interval_key])
                    data[k] = mdsteps / interval * bytes_per_step
                except KeyError:
                    # default
                    data[k] = mdsteps / default_interval * bytes_per_step
                except ZeroDivisionError:
                    data[k] = 0   # no printout
                finally:
                    # if energy or trajectory, check that files for output are
                    # defined, otherwise set the printout to 0
                    if interval_key in ("energy", "trajectory") and not \
                            interval_key in list(step.parameters["files"].keys()):
                        data[k] = 0


            trj_data = data["trj"]
            en_data = data["en"]
            log_data = (data["log"] + data["temp"] + data["nb"])
            rest_data = num_atoms_all * REST_B_PER_ATOM
            total_data["trj"] += trj_data
            total_data["log"] += log_data
            total_data["en"] += en_data
            total_data["rest"] += rest_data

            data = (trj_data + log_data + rest_data + en_data)/CONV_MB

            summary += "{:<10} {:>5} {:>10} {:>10} {:^10} {:^10} {:^10} "\
                       "{:^30} {:^10} {:>10.2f}\n" \
                       "".format(nstep, mdparms["temperature"],
                                 mdparms["stepsize"],
                                 locale.format('%d', mdsteps, 1),
                                 seq, dist, angle, shake, random_seed, data)

        except KeyError as err_msg:
            raise QGenrelaxError("You are missing either 'steps', "
                                 "'temperature' or 'stepsize' in one of your "
                                 "relaxation steps. These parameters are "
                                 "quite important you know...")

    summary += "Restraints:\n"
    for i, rest in enumerate(restraints):
        summary += "{}: {}\n".format(i+1, rest)

    summary += """
Total time: {} ps
Total wasted storage (wild approximation): \
{:.2f} MB (trj: {:.1f}, log: {:.1f}, en: {:.1f}, rest: {:.1f})
""".format(total_time/1000.0, sum(total_data.values())/CONV_MB,
           total_data["trj"]/CONV_MB, total_data["log"]/CONV_MB,
           total_data["en"]/CONV_MB, total_data["rest"]/CONV_MB)

    for l in summary.split("\n"):
        logger.info(l)

    return inp_fns



def genfeps(fep_proc_file, relax_input_file, restraint, energy_list_fn,
            frames, repeats, fromlambda, prefix, first_frame_eq,
            pdb_file=None, fep_file=None, runscript_file=None,
            ignore_errors=False):

    """Generates inputs for a FEP/MD simulation with Q (Qdyn).

    Args:
        fep_proc_file (string):  genfeps procedure file pathname
        relax_input_file (string):  pathname of last relaxation step input
        restraint (string):  restraint coordinate (a)
        energy_list_fn (string):  name of file that will contain the en.f.list
        frames (int):  number of FEP frames
        repeats (int):  number of repeats/replicas
        fromlambda (float):  starting lambda (0.0 - 1.0)
        prefix (string):  prefix for repeat directories
        first_frame_eq (boolean):  use equil instead of first frame (cadee)

    Optional args:
        pdb_file (string):  pdb pathname (used to convert placeholders)
        fep_file (string):  alternate fep file pathname (ignoring input's fep)
        runscript_file (string):  slurm/sge run script
        ignore_errors (boolean):  passed to QStruct and QDynInp - write to\
                                  logger instead of raising exceptions on\
                                  non-critical things

    Returns:
        rep_dirs (list):  list of created replica folders

    (a) Restraint coordinate can be set to:
    'inp' - whatever is defined in relax_input_file
    'top' - topology
    'relax' - endpoint of relaxation

    """

    frames = int(frames)
    repeats = int(repeats)
    if fromlambda != None:
        fromlambda = float(fromlambda)

    # constants
    PREFIX_EQ = "equil_"
    PREFIX_FEP = "fep_"


    # check if files exist
    for k, v in six.iteritems(locals()):
        if k in ["pdb_file", "fep_proc_file", "fep_file",
                 "runscript_file", "relax_input_file"]:
            if v and not os.path.lexists(v):
                raise QGenfepsError("File '{}' doesn't exist.".format(v))

    if restraint not in ["top", "relax", "inp"]:
        raise QGenfepsError("Argument 'restraint' has to be either "
                            "'inp', 'top' or 'relax'")

    # find and replace atom placeholders.
    # if no PDB was given to replace them, exit
    fep_proc_str = open(fep_proc_file, 'r').read()
    c = find_placeholders(fep_proc_str)
    if c and not pdb_file:
        raise QGenfepsError("Found placeholders in proc. file, but no PDB "
                            "was given: {}".format(", ".join(c)))
    elif c:
        logger.info("These placeholders will be replaced with atom indices: "
                    + ", ".join(c))
        try:
            qstruct = QStruct(pdb_file, "pdb", ignore_errors=ignore_errors)
            fep_proc_str = qstruct.convert_placeholders(fep_proc_str)
        except QStructError as err_msg:
            raise QGenfepsError("Failed to replace placeholders: {}"
                                "".format(err_msg))


    # make a nice header comment in each input file with the
    header_comment = """\
# Generated with QTools, version {}
# Date: {}
# CWD: {}
# Cmdline: {}
""".format(__version__, time.ctime(), os.getcwd(), " ".join(sys.argv))


    # get topology and fep and others from the last relaxation input
    top_file_abs, fep_file_abs, re_file_abs, rest_file = None, None, None, None
    lambda_initial = None
    try:
        c = QDynInput(open(relax_input_file, 'r').read(),
                      ignore_errors=ignore_errors)
    except QDynInputError as err_msg:
        raise QGenfepsError("There is something wrong with the given input "
                            "file ({}): {}".format(relax_input_file, err_msg))


    di = os.path.dirname(relax_input_file)
    try:
        files = c.parameters["files"]
        lambda_initial = float(c.parameters["lambdas"].split()[0])
        top_file_abs = os.path.join(di, files["topology"])
        re_file_abs = os.path.join(di, files["final"])
        fep_file_abs = os.path.join(di, files["fep"])
        if "restraint" in files:
            rest_file = os.path.join(di, files["restraint"])
    except KeyError as err_msg:
        raise QGenfepsError("Parsing the relaxation input file failed, "
                            "keyword missing... {}".format(err_msg))

    # check if the files actually exist
    for fn, descr in [(top_file_abs, "topology"),
                      (fep_file_abs, "fep"),
                      (re_file_abs, "final"),
                      (rest_file, "restraint")]:
        if fn and not os.path.lexists(fn):
            raise QGenfepsError("When parsing the input, found this filename "
                                "'{}' next to the '{}' command. Unfortunately,"
                                " the file doesnt exist...".format(fn, descr))

    # change the FEP (when debugging your system, you might want to
    # use an old relax and not waste 100M core hours when changing
    # a soft core value in the fep)
    if fep_file:
        fep_file_abs = os.path.abspath(fep_file)


    # find and replace atom placeholders in FEP file
    # if no PDB was given to replace them, exit
    fep_file_str = open(fep_file_abs, 'r').read()
    c = find_placeholders(fep_file_str)
    if c and not pdb_file:
        raise QGenfepsError("Found placeholders in FEP file, but no PDB was "
                            "given: {}".format(", ".join(c)))
    elif c:
        logger.info("Replacing FEP file placeholders...")
        try:
            qstruct = QStruct(pdb_file, "pdb", ignore_errors=ignore_errors)
            fep_file_str = qstruct.convert_placeholders(fep_file_str)
        except QStructError as err_msg:
            raise QGenfepsError("Failed to replace placeholders: {}"
                                "".format(err_msg))


    # change the inital lambda (this is not recommended, the system should
    # be properly relaxed at a particual lambda before doing FEP)
    if fromlambda != None:
        lambda_initial = float(fromlambda)
        if lambda_initial > 1.0 or lambda_initial < 0.0:
            raise QGenfepsError("Lambda value is bogus, are you on drugs?")

    # create lambda values, find the closest to the starting one and
    # rearrange accordingly: [0.0, 0.02, 0.04, ... 0.98, 1.0]  for frames==51
    lambdas = [float(num) / (frames - 1) for num in range(0, frames)]

    # [2,]   for lambda_initial == 0.04 (or close to 0.04) and frames==51
    l_i = [i for i in range(0, frames) if \
            abs(lambdas[i] - lambda_initial) <= (1.0 / frames)]
    # there should be only one
    l_i = l_i[0]
    lambda_initial = lambdas[l_i]

    # [0.02, 0.0,] for the case of lambda_initial == 0.04 and frames == 51
    forward_lambdas = list(reversed(lambdas[0:l_i]))
    # [0.06, 0.08, ..., 1.0] for the case of lambda_initial == 0.04, fr. == 51
    backward_lambdas = lambdas[l_i+1:]

    lambdas = [lambda_initial,] + forward_lambdas + backward_lambdas

    # print out some useful information
    logger.info("Using restart file: {}"
                "".format(os.path.relpath(re_file_abs)))
    logger.info("Using topology file: {}"
                "".format(os.path.relpath(top_file_abs)))
    logger.info("Using FEP file: {}"
                "".format(os.path.relpath(fep_file_abs)))
    logger.info("Starting from lambda value (state 1): {}"
                "".format(lambda_initial))
    logger.info("Number of FEP frames: {} ".format(frames))


    # create a temporary directory to store the files that are identical
    # in all replicas - top, fep, runscript, relax restart, restraint file
    # (if any) and copy the common files
    TMPDIR = tempfile.mkdtemp()
    top_fn = os.path.basename(top_file_abs)
    fep_fn = os.path.basename(fep_file_abs)
    relax_re_fn = "cont_" + os.path.basename(re_file_abs)
    shutil.copy2(top_file_abs, TMPDIR)
    shutil.copy2(re_file_abs, os.path.join(TMPDIR, relax_re_fn))

    open(os.path.join(TMPDIR, fep_fn), "w").write(fep_file_str)
    if runscript_file:
        shutil.copy2(runscript_file, TMPDIR)
    else:
        logger.info("No Q runscript given.")

    # handle the whole restraint coordinates crap...
    # rest_file is the file from the relaxation input (if any)
    # rest_fn is the basename of the restraints file (either from input
    # or relaxed.re.rest), or None if rest. to topology
    if restraint == "relax":
        logger.info("Restraining to: relaxation")
        rest_fn = "cont_" + os.path.basename(re_file_abs) + ".rest"
        shutil.copy2(re_file_abs, os.path.join(TMPDIR, rest_fn))
    elif restraint == "top":
        logger.info("Restraining to: topology")
        rest_fn = None
    else: # default, from input
        if rest_file:
            logger.info("Restraining to: {} (from input)"
                        "".format(os.path.relpath(rest_file)))
            rest_fn = "cont_" + os.path.basename(rest_file)
            shutil.copy2(rest_file, os.path.join(TMPDIR, rest_fn))
        else:
            logger.info("Restraining to: topology (from input)")
            rest_fn = None


    # parse the proc file
    general_inp = []
    eq_steps_inps = [[],]
    fep_inp = []
    script_vars = {}

    section = ""
    for line in fep_proc_str.split("\n"):
        # remove comments and strip whitespaces.
        line = re.split("#|\!", line)[0].strip()
        # empty lines are useless
        if line == "":
            continue
        # found a section
        if line[0] == "{":
            section = line.strip("{}").lower()
            continue

        if not section:
            raise QGenfepsError("Parsing the procedure file failed... This "
                                "line: '{}' is not inside any section:"
                                "".format(line))

        if section == "script_vars":
            c = line.split()
            var, value = c[0], " ".join(c[1:])
            script_vars[var] = value
        elif section == "general":
            general_inp.append(line)
        elif section == "steps_equil":
            if "__________" in line:
                eq_steps_inps.append([])
            else:
                eq_steps_inps[-1].append(line)
        elif section == "fep":
            fep_inp.append(line)
        else:
            raise QGenfepsError("Parsing the procedure file failed: "
                                "Unsupported section: '{}'".format(section))


    # check for steps with no parameters (too many _________ lines)
    # and remove them
    for i in range(len(eq_steps_inps)-1, -1, -1):
        if not eq_steps_inps[i]:
            eq_steps_inps.pop(i)

    # check for missing sections
    for l, n in ((general_inp, "GENERAL"),
                 (eq_steps_inps, "STEPS_EQUIL"),
                 (fep_inp, "FEP")):
        if not l:
            raise QGenfepsError("Parsing the procedure file failed: "
                                "Section '{}' is missing".format(n))


    # join lists of lines to strings and replace the placeholders
    script_variables = sorted(list(script_vars.items()), reverse=True)
    gen_inp_s = "\n".join(general_inp)
    fep_inp_s = "\n".join(fep_inp)
    eq_steps_inps_s = ["\n".join(eq_s_inp) for eq_s_inp in eq_steps_inps]

    for placeholder, value in script_variables:
        gen_inp_s = gen_inp_s.replace(placeholder, value)
        fep_inp_s = fep_inp_s.replace(placeholder, value)
        for step_i, eq_step_inp_s in enumerate(eq_steps_inps_s):
            eq_steps_inps_s[step_i] = eq_step_inp_s.replace(placeholder, value)


    ####################
    # make equil. inputs
    eq_steps = []
    for step_n, eq_step_inp_s in enumerate(eq_steps_inps_s):
        # create the files section
        final = "{}{:03d}_{:4.3f}.re".format(PREFIX_EQ, step_n, lambda_initial)
        dcd = "{}{:03d}_{:4.3f}.dcd".format(PREFIX_EQ, step_n, lambda_initial)
        files = {"final"      : final,
                 "trajectory" : dcd,
                 "topology"   : top_fn,
                 "fep"        : fep_fn}

        if first_frame_eq:
            files["energy"] = "{}{:03d}_{:4.3f}.en".format(PREFIX_EQ,
                                                           step_n,
                                                           lambda_initial)

        if rest_fn:
            files["restraint"] = rest_fn

        if step_n != 0:
            files["restart"] = "{}{:03d}_{:4.3f}.re".format(PREFIX_EQ,
                                                            step_n-1,
                                                            lambda_initial)
        else:
            files["restart"] = relax_re_fn

        # parse the general input and update with step input and files section
        try:
            inp = QDynInput(gen_inp_s, ignore_errors=ignore_errors)
            inp.update(eq_step_inp_s)
            if "energy" in inp.parameters["intervals"]:
                files["energy"] = "{}{:03d}_{:4.3f}.en".format(PREFIX_EQ,
                                                               step_n,
                                                               lambda_initial)
            elif first_frame_eq:
                raise QGenfepsError("Argument 'first_frame_eq' requires the "
                                    "energy printout defined in the intervals "
                                    "section of the equilibration "
                                    "(e.g. 'energy   10')")

            inp.update(parameters={"files": files})
            inp.update(parameters={"lambdas": "{:9.7f} {:9.7f}"
                                              "".format(lambda_initial,
                                                        1-lambda_initial)})
        except QDynInputError as err_msg:
            raise QGenfepsError("Problem with equil. step no. {}: {}"
                                "".format(step_n, err_msg))

        # test the input string
        try:
            _ = inp.get_string()
        except QDynInputError as err_msg:
            raise QGenfepsError("Error in equil. step {}: {}"
                                "".format(step_n, err_msg))

        # check if random seed is not defined or is fixed in the first step
        if step_n == 0:
            if repeats > 1:
                if ("random_seed" not in inp.parameters["md"]) or \
                        (int(inp.parameters["md"]["random_seed"]) > 0):
                    raise QGenfepsError("Fixed random seed (or restart "
                                        "velocities) works only with one "
                                        "repeat (others will be identical).\n"
                                        "Please use 'random_seed   -1' in "
                                        "your first equilibration step to "
                                        "generate random random seeds.")

            elif "random_seed" not in inp.parameters["md"]:
                logger.info("No random seed in first step of equilibration,"
                            "using restart velocities.")

                if (not rest_file and rest_fn) or (not rest_fn and rest_file) \
                  or (rest_file and (os.path.basename(rest_file) != rest_fn)):
                    logger.warning("This will not be a true continuation run! "
                                   "The relaxation restraint does not match "
                                   "yours. Use 'inp' instead of 'top' or "
                                   "'relax' for the restraint.")

        # append the input
        eq_steps.append(inp)



    #################
    # make FEP inputs
    en_filenames = []
    feps = []

    for step_n, lam in enumerate(lambdas):
        # create the files section
        final = "{}{:03d}_{:4.3f}.re".format(PREFIX_FEP, step_n, lam)
        dcd = "{}{:03d}_{:4.3f}.dcd".format(PREFIX_FEP, step_n, lam)
        en = "{}{:03d}_{:4.3f}.en".format(PREFIX_FEP, step_n, lam)
        files = {"final"      : final,
                 "trajectory" : dcd,
                 "topology"   : top_fn,
                 "energy"     : en,
                 "fep"        : fep_fn}

        # if this step is in new direction (backwards) then
        # set the previous lambda and step to initial
        if backward_lambdas and lam == backward_lambdas[0]:
            prev_fep = feps[0]
        elif step_n == 0:
            prev_fep = eq_steps[-1]
        else:
            prev_fep = feps[-1]

        # if this flag is set, all steps that point to the first step
        # should point to the last eq step
        if first_frame_eq:
            if step_n == 1 or (backward_lambdas and lam == backward_lambdas[0]):
                prev_fep = eq_steps[-1]

        if rest_fn:
            files["restraint"] = rest_fn

        files["restart"] = prev_fep.parameters["files"]["final"]


        # update the parameters and check the input
        try:
            inp = QDynInput(gen_inp_s, ignore_errors=ignore_errors)
            inp.update(fep_inp_s)
            if "energy" not in inp.parameters["intervals"]:
                raise QGenfepsError("FEP stage requires the energy printout "
                                    "defined in the intervals section "
                                    "(e.g. 'energy   10')")

            inp.update(parameters={"files": files})
            inp.update(parameters={"lambdas": "{:9.7f} {:9.7f}"
                                              "".format(lam, 1-lam)})
            inp.check()
        except QDynInputError as err_msg:
            raise QGenfepsError("Error in FEP step {}: {}"
                                "".format(step_n, err_msg))

        # append the input
        feps.append(inp)

        # add the energy filename to the list
        en_filenames.append(inp.parameters["files"]["energy"])


    # if first_frame_eq is set add the energy file and remove the
    # first fep frame
    if first_frame_eq:
        logger.info("Replacing the first FEP frame with the last "
                    "equilibration step")
        en_filenames[0] = eq_steps[-1].parameters["files"]["energy"]
        feps.pop(0)


    # check random seed in fep
    if "random_seed" in feps[0].parameters["md"] and \
            int(feps[0].parameters["md"]["random_seed"]) < 1:
        logger.warning("Generating random seeds in FEP inputs. "
                       "Are you sure this is ok?")

    # write a file that contains the names of all energy files in proper order
    # this file is used later by q_mapper.py
    # sort the enfiles according to lambda (1.0 -> 0.0) so that the mapping
    # will always go from reactants to products
    enfiles_lambdas = sorted([(enf.split("_")[-1], i) for i, enf in \
            enumerate(en_filenames)], reverse=True)
    en_filenames_sorted = [en_filenames[i] for l, i in enfiles_lambdas]
    enf = os.path.join(TMPDIR, energy_list_fn)
    open(enf, 'w').write("\n".join(en_filenames_sorted))

    # create directories for repeats/replicas (rep_000,rep_001,rep_002...)
    # copy everything from TMPDIR (topology, fep file, relax restart and
    # restraint file (if any)); create the eq and fep inputs
    #
    # first check for existing directories
    for num in range(0, repeats):
        rep = "{}{:03d}".format(prefix, num)
        if os.path.lexists(rep):
            raise QGenfepsError("Directory '{}' exists. Please (re)move it or "
                                "change the prefix with --prefix.".format(rep))

    lsdir = os.listdir(TMPDIR)
    rep_dirs = []
    for num in range(0, repeats):
        rep = "{}{:03d}".format(prefix, num)
        os.mkdir(rep)
        # copy stuff from TMPDIR
        for f in lsdir:
            shutil.copy2(os.path.join(TMPDIR, f), rep)

        # create eq inputs
        for step_n, eq_step in enumerate(eq_steps):

    # check if random seed is a fixed value or not (generate random or fail)
            eqs = copy.deepcopy(eq_step)  # a copy
            if "random_seed" in eqs.parameters["md"] and \
                    int(eqs.parameters["md"]["random_seed"]) < 1:
                rs = random.randint(1, 1e6)
                eqs.update(parameters={"md": {"random_seed" : rs}})

            try:
                s = eqs.get_string()
            except QDynInputError as err_msg:
                raise QGenfepsError("Error in step {}: {}"
                                    "".format(step_n, err_msg))
            fn = os.path.join(rep, "{}{:03d}_{:4.3f}.inp"
                                   "".format(PREFIX_EQ,
                                             step_n,
                                             lambda_initial))
            s = header_comment + s
            open(fn, 'w').write(s)

        last_eq_fn = fn
        # create FEP inputs
        for step_n, fep in enumerate(feps):
            if first_frame_eq:
                step_n += 1

            fs = copy.deepcopy(fep)  # a copy
            if "random_seed" in fs.parameters["md"] and \
                    int(fs.parameters["md"]["random_seed"]) < 1:
                rs = random.randint(1, 1e6)
                fs.update(parameters = {"md": {"random_seed" : rs}})

            try:
                s = fs.get_string()
            except QDynInputError as err_msg:
                raise QGenfepsError("Error in step {}: {}"
                                    "".format(step_n, err_msg))
            lam = lambdas[step_n]  # feps was created in lambdas iteration
            fn = os.path.join(rep, "{}{:03d}_{:4.3f}.inp"
                                   "".format(PREFIX_FEP, step_n, lam))
            s = header_comment + s
            open(fn, 'w').write(s)
      
        logger.info("Created inputs for repeat/replica '{}'.".format(rep))
        rep_dirs.append(rep)



    # get the amount of storage that will be wasted
    # for this we need the atom count from the topology
    for line in open(os.path.join(TMPDIR,top_fn), 'r').readlines(1024):
        if "no. of atoms, no. of solute atoms" in line:
            num_atoms_all = int(line.split()[0])
            break

    REST_B_PER_ATOM = 48.0
    TRJ_B_PER_ATOM = 12.0
    # very rough estimate, depends on Q version
    # it can double if group_contributions are calculated
    EN_B_PER_STEP = 370.0
    CONV_MB = 2**20
    # very rough estimate
    OUT_B_PER_STEP = 2000
    TEMP_B_PER_STEP = 160
    NB_B_PER_STEP = 80

    # intervals maps: q_parameter_key, q_default_value, approx_bytes_per_frame
    qintervals = {"trj"  : ["trajectory",   100, num_atoms_all*TRJ_B_PER_ATOM],
                  "log"  : ["output",        10, OUT_B_PER_STEP],
                  "temp" : ["temperature",   10, TEMP_B_PER_STEP],
                  "en"   : ["energy",        10, EN_B_PER_STEP],
                  "nb"   : ["non_bond",      10, NB_B_PER_STEP]}
    total_data = {"trj": 0, "log": 0, "en": 0, "rest": 0}
            # calculate approx amount of data
    for i, step in enumerate(eq_steps+feps):
        data = {}
        mdsteps = int(step.parameters["md"]["steps"])
        for k, v in six.iteritems(qintervals):
            interval_key = v[0]
            default_interval = v[1]
            bytes_per_step = v[2]
            try:
                interval = int(step.parameters["intervals"][interval_key])
                data[k] = mdsteps / interval * bytes_per_step
            except KeyError:
                # default
                data[k] = mdsteps / default_interval * bytes_per_step
            except ZeroDivisionError:
                data[k] = 0   # no printout
            finally:
                # if energy or trajectory, check that files for output are
                # defined, otherwise set the printout to 0
                if interval_key in ("energy", "trajectory") and not \
                   interval_key in list(step.parameters["files"].keys()):
                    data[k] = 0

        trj_data = data["trj"]
        en_data = data["en"]
        log_data = (data["log"] + data["temp"] + data["nb"])
        rest_data = num_atoms_all * REST_B_PER_ATOM
        total_data["trj"] += trj_data
        total_data["log"] += log_data
        total_data["en"] += en_data
        total_data["rest"] += rest_data

        data = (trj_data + log_data + rest_data + en_data)/CONV_MB

    logger.info("Your runs will waste approx. {:.2f} MB of storage. "
                "Per replica: {:.2f} MB (trj: {:.1f}, log: {:.1f}, "
                "en: {:.1f}, rest: {:.1f})"
                .format(sum(total_data.values())/CONV_MB*repeats,
                        sum(total_data.values())/CONV_MB,
                        total_data["trj"]/CONV_MB,
                        total_data["log"]/CONV_MB,
                        total_data["en"]/CONV_MB,
                        total_data["rest"]/CONV_MB))

    # remove temporary directory
    shutil.rmtree(TMPDIR)

    return rep_dirs


