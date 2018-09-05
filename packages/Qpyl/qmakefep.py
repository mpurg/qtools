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
This module contains the make_fep() function for generating Q FEP-files,
and implements a custom exception class QMakeFepError.
"""

import sys
import os
import re
import time
import tempfile
import logging
from collections import OrderedDict as ODict

from Qpyl.core.qlibrary import QLib, QLibError
from Qpyl.core.qparameter import QPrm, QPrmError
from Qpyl.core.qstructure import QStruct, QStructError
from Qpyl.core.qtopology import QTopology, QTopologyError
from Qpyl.common import raise_or_log, __version__

logger = logging.getLogger(__name__)

SUPPORTED_FF = ["amber", "oplsaa"]


class QMakeFepError(Exception):
    pass

class _FepPrmMorse(object):
    """Class for Morse parameters."""
    def __init__(self, harmonic_prm):
        self.harmonic_prm = harmonic_prm

def make_fep(qmap_file, pdb_file, forcefield,
             parm_files, lib_files, ignore_errors=False):
    """Generate a template FEP file for EVB simulations in Q.

    Parses a QMAP file (see below), state 1 structure file (PDB) and
    all states libraries and parameters, and determines the changes
    in connectivity/charges/parameters that occur between the states.

    QMAP is a text file that defines mappings of library ids (for each state)
    to state 1 structure/topology ids, best explained on an example:

    q   315.O     OHH.O     OHH.O
    q   315.H1    OHH.H1    HIP.HE2
    q   315.H2    OHH.H2    OHH.H2
    q   155.NE2   HID.NE2   HIP.NE2
    ...
    n   155.CA    HID.CA    HIP.CA

    The first column defines the atom as being a 'Q' atom or a 'neighbouring'
    atom. The latter will not be included in the 'Q-region' but will be
    included in the 'change_bonds/change_angles...' sections in case there is a
    change in bonding/parameters outside the Q region. Additionally, you can
    define a 'q' atom with 'q_qcp', which will create an additional section for
    isotopically clean masses used in QCP calculations.
    The second column is the PDB ID, comprised of residue index and atom name,
    separated by a dot.
    The third column is the library ID of this atom in state 1, comprised of
    residue name and atom name (should be the same as in the structure).
    The fourth column is the library ID of this atom in state 2.
    Additional columns can be added for other states.

    The returned template string contains several missing parts, denoted with
    <FIX>, which have to be manually replaced with appropriate values. These
    include the softpair C parameters, Morse parameters, Hij parameters.

    Args:
        qmap_file (string): QMAP file path
        pdb_file (string): state 1 PDB file path (the one built with qprep)
        forcefield (string): forcefield type (see SUPPORTED_FF)
        prms_files (list): Q parameter-file paths
        libs_files (list): Q library-file paths
        ignore_errors (boolean, optional): don't fail on certain non critical\
                                           errors

    Returns:
        fepstr (string):  fepfile template

    Raises:
        QMakeFepError

    """
    
    if forcefield not in SUPPORTED_FF:
        raise QMakeFepError("Force field '{}' not supported. Use {}"
                            "".format(forcefield, " or ".join(SUPPORTED_FF)))


    fep_types = {"atoms": [], "bonds": [], "angles": [],
                 "torsions": [], "impropers": []}
    fep_changes = {"atoms": [], "charges": [],
                   "bonds": ODict(), "angles": ODict(),
                   "torsions": ODict(), "impropers": ODict()}
    fep_qcp_atoms = []
    fep_morse_prms = {}
    fep_reacting_atoms = set()
    num_evb_states = None

    # parse the MAP file
    # pdb_ids_map = [ ('q', [pdbid1_state1,]),
    #                 ('q', [pdbid2_state1,]),
    #                 ...
    #                 ('n', [pdbid11_state1,]),
    #                 ...
    #               ]
    # lib_ids_map = [ [lib_id1_state1, lib_id2_state1...],
    #                 [lib_id1_state2, lib_id2_state2...],
    #                 ...
    #               ]
    #
    lib_ids_map = []
    pdb_ids_map = []
    with open(qmap_file, 'r') as qatom_map:
        for i, line in enumerate(qatom_map.readlines()):
            line = re.split("#|\*|\!", line, 1)[0].strip()   # remove comments
            if line == "":
                continue
            c = line.split()
            atom_type = c[0].lower()
            pdb_id = c[1]
            lib_ids = c[2:]

            if atom_type not in ["q", "n", "q_qcp"]:
                raise QMakeFepError("Lines in the QMAP file should begin "
                                    "with either 'q' (qatom) or 'n' "
                                    "(neighboring atom) or 'q_qcp' "
                                    "(QPI q atom)")
            try:
                resid, name = pdb_id.split(".")
                if not name or not int(resid): raise ValueError

            except ValueError:
                raise QMakeFepError("Invalid PDB ID '{}'. Should be "
                                    "RESID.ATOMNAME".format(pdb_id))

            tmp = (atom_type, [pdb_id,])
            if tmp in pdb_ids_map:
                raise QMakeFepError("Duplicate PDB ID: '{}'".format(pdb_id))
            pdb_ids_map.append(tmp)

            if num_evb_states == None:
                num_evb_states = len(lib_ids)
            elif len(lib_ids) != num_evb_states:
                raise QMakeFepError("Number of states in line '{}' not equal "
                                    "to number of PDB files".format(line))

            for state, lib_id in enumerate(lib_ids):
                try:
                    resname, name = lib_id.split(".")
                    if not resname or not name: raise ValueError
                except ValueError:
                    raise QMakeFepError("Invalid library ID '{}'. Should be "
                                        "RESNAME.ATOMNAME".format(lib_id))
                try:
                    if lib_id in lib_ids_map[state]:
                        raise QMakeFepError("The library IDs in one EVB state "
                                            "should be unique (double '{}'), "
                                            "otherwise proper bonding can't "
                                            "be determined.".format(lib_id))
                except IndexError:
                    lib_ids_map.append([])

                lib_ids_map[state].append(lib_id)


    # load libraries
    qlib = QLib(forcefield, ignore_errors=ignore_errors)
    for lib in lib_files:
        try:
            qlib.read_lib(lib)
        except QLibError as e:
            raise QMakeFepError("Problem parsing lib ({}): {}"
                                "".format(lib, e))


    # make dummy structures for other states
    structures = [None for _ in range(num_evb_states)]
    structures[0] = pdb_file
    libid_pdbid_map = [{} for _ in range(num_evb_states)]

    for state in range(1, num_evb_states):
        state_structure = []
        atom_index = 1
        processed_residues = []
        for i, (q_or_n, pdb_ids_all_states) in enumerate(pdb_ids_map):

            lib_id = lib_ids_map[state][i]
            resname, aname = lib_id.split(".")

            # add all atoms of current residue to the dummy structure
            # at the same time, storing the mapping lib_id:pdb_id
            # in libid_pdbid_map for later
            if resname not in processed_residues:
                try:
                    residue_lib = qlib.residue_dict[resname]
                except KeyError:
                    raise QMakeFepError("Residue '{}' not found in library."
                                        "".format(resname))
                processed_residues.append(resname)
                res_index = len(processed_residues)

                for atom in residue_lib.atoms:
                    lib_id2 = "{}.{}".format(resname, atom.name)
                    pdb_id2 = "{}.{}".format(res_index, atom.name)
                    state_structure.append("{:<6}{:5} {:4} {:3} {:5}    "
                                           "{:8.3f}{:8.3f}{:8.3f}"
                                           "".format("ATOM", atom_index,
                                                     atom.name, resname,
                                                     res_index, 0, 0, 0))
                    atom_index += 1

                    # map the newly created dummy atom's pdb_id to lib_id
                    libid_pdbid_map[state][lib_id2] = pdb_id2

            # add pdb_id of current atom in current (dummy structure)
            # state to pdb_ids_map (using its lib_id)
            try:
                pdb_id_this_state = libid_pdbid_map[state][lib_id]
            except KeyError:
                raise QMakeFepError("Library ID '{}' not valid.".format(lib_id))
            pdb_ids_all_states.append(pdb_id_this_state)

        _, structures[state] = tempfile.mkstemp()
        open(structures[state], "w").write("\n".join(state_structure))
        # DEBUG
        # print "Dummy PDB for st.{}: {}".format(state + 1, structures[state])



    # load parameters
    qprm = QPrm(forcefield, ignore_errors=ignore_errors)
    for parm in parm_files:
        try:
            qprm.read_prm(parm)
        except QPrmError as e:
            raise QMakeFepError("Problem with parm ({}): {}"
                                "".format(parm, e))

    # load structures and make topologies
    topologies = []
    for state in range(num_evb_states):
        try:
            qstruct = QStruct(structures[state], "pdb",
                              ignore_errors=ignore_errors)
        except QStructError as e:
            raise QMakeFepError("Problem parsing PDB file ({}): {} "
                                "".format(structures[state], e))

        try:
            topologies.append(QTopology(qlib, qprm, qstruct))
        except QTopologyError as e:
            raise QMakeFepError("Problem building the topology: {}"
                                "".format(e))



    # Make _TopoAtom (atoms in QTopology) maps out of qmap's lists
    # and extract types, type changes and charge changes
    #
    # atom_map = [ [_TopoAtom1_state1, _TopoAtom1_state2, ... ],
    #              [_TopoAtom2_state1, _TopoAtom2_state2, ... ], 
    #              [_TopoAtom3_state1, _TopoAtom3_state2, ... ], 
    #            ...
    #            ]
    atom_map = []
    for i, (q_or_n, pdb_id_all_states) in enumerate(pdb_ids_map):
        atom_all_states = []
        for state, pdb_id in enumerate(pdb_id_all_states):
            residue, aname = pdb_id.split(".")
            try:
                residue = topologies[state].residues[int(residue)-1]
                atom = [a for a in residue.atoms if a.name == aname][0]
            except (KeyError, IndexError) as e:
                raise QMakeFepError("Atom '{}' doesn't exist in PDB '{}'"
                                    "".format(pdb_id, structures[state]))
            atom_all_states.append(atom)
        atom_map.append(atom_all_states)

        # check for stupidity - lib_id in QMAP state 1 
        # not matching the structure/topology
        lib_id_qmap = lib_ids_map[0][i]
        lib_id = "{}.{}".format(atom_all_states[0].residue.name,
                                atom_all_states[0].name)
        if lib_id != lib_id_qmap:
            pdb_id = pdb_ids_map[i][1][0]
            raise QMakeFepError("QMAP state 1 library ID ({}) of atom '{}' "
                                "doesn't match topology library ID ({})."
                                "".format(lib_id_qmap, pdb_id, lib_id))


        # For Q atoms (and not the neighbor atoms):
        # get FEP atom types, type changes and charge changes
        if q_or_n in ["q", "q_qcp"]:
            for atom in atom_all_states:
                if atom.prm not in fep_types["atoms"]:
                    fep_types["atoms"].append(atom.prm)

            fep_changes["charges"].append([a.charge for a in atom_all_states])
            fep_changes["atoms"].append(atom_all_states)
            if q_or_n == "q_qcp":
                fep_qcp_atoms.append(atom_all_states)

    charge_sums = []
    for state in range(num_evb_states):
        charge_sum = sum([c[state] for c in fep_changes["charges"]])
        if abs(round(charge_sum) - charge_sum) > 1e-6:
            raise_or_log("Net charge in state {} not integer: {}"
                         "".format(state + 1, charge_sum),
                         QMakeFepError, logger, ignore_errors=ignore_errors)
        charge_sums.append(charge_sum)

    if any([abs(c-charge_sums[0]) > 1e-6 for c in charge_sums]):
        logger.warning("Net charge changes between states: {}"
                       "".format(" -> ".join([str(c) for c in charge_sums])))



    # get all Bonds, Angles, Torsions and Impropers which include
    # at least one atom defined in qmap
    batis = {"bonds": [], "angles": [], "torsions": [], "impropers": []}
    batis["bonds"] = [set() for _ in range(num_evb_states)]
    batis["angles"] = [set() for _ in range(num_evb_states)]
    batis["torsions"] = [set() for _ in range(num_evb_states)]
    batis["impropers"] = [set() for _ in range(num_evb_states)]
    for atom_all_states in atom_map:
        for state, atom in enumerate(atom_all_states):
            _ = [batis["bonds"][state].add(b) for b in atom.bonds]
            _ = [batis["angles"][state].add(a) for a in atom.angles]
            _ = [batis["torsions"][state].add(t) for t in atom.torsions]
            _ = [batis["impropers"][state].add(i) for i in atom.impropers]


    # map the bonds,angles,torsions,impropers (bati) in different states
    # to same key (ordered list of state1 PDB_IDs)
    #
    # bati_map =
    # { "bonds": {state1_bond1_key: [bond1_state1, bond1_state2,...],
    #             state1_bond2_key: [bond2_state1, bond2_state2,...], ...},
    #   "angles": {state1_angle1_key: [angle1_state1, angle1_state2,...],...}
    #  ... }
    #
    # also, include only batis which have all atoms defined in qmap
    # also, detect inter-residue batis and raies QMakeFepError

    bati_map = {"bonds": {}, "angles": {}, "torsions": {}, "impropers": {}}
    for state in range(num_evb_states):
        atoms_in_state = [a_all_st[state] for a_all_st in atom_map]
        for bati_type in bati_map:
            for bati in batis[bati_type][state]:

                # find the corresponding atoms in state1
                try:
                    atoms_st1 = [atom_map[atoms_in_state.index(a)][0] for a in
                                 bati.atoms]
                except ValueError:
                    # one of the Atoms is not defined in QMAP
                    continue

                pdbid_index = []
                for atom in atoms_st1:
                    pdbid_index.append((atom.index,
                                        "{}.{}".format(atom.residue.index,
                                                       atom.name)))
                # order the pdbids to prevent double entries
                if bati_type == "bonds":
                    pids = sorted(pdbid_index)
                elif bati_type == "angles":
                    pids = min(pdbid_index, list(reversed(pdbid_index)))
                elif bati_type == "torsions":
                    pids = min(pdbid_index, list(reversed(pdbid_index)))
                elif bati_type == "impropers":
                    # topology order == library order == correct order
                    pids = pdbid_index
                key = " ".join([p[1] for p in pids])

                # check for bonds/angles/torsions/impropers that are
                # shared between residues
                residue_ids = set(atom.residue.index for atom in bati.atoms)
                if len(residue_ids) > 1:
                    raise QMakeFepError("Inter-residue bond/angle/torsion '{}'"
                                        " not supported. Combine the residues "
                                        "into a single library entry if you "
                                        "want to make changes over the "
                                        "'head-tail' bond.".format(key))

                # add bati to bati_map
                try:
                    bati_map[bati_type][key][state] = bati
                except KeyError:
                    bati_map[bati_type][key] = [None for _ in
                                                range(num_evb_states)]
                    bati_map[bati_type][key][state] = bati


    # DEBUG
    # for k,v in bati_map.iteritems():
        # print k
        # for k2, v2 in v.iteritems():
            # print k2, v2[0], v2[1]

    def _bati_sort(key, bati_all_states):
        # to sort bonds/angles.. based on the key
        # also, breaking and forming bonds have priority
        try:
            return (-1 * bati_all_states.index(None), key)
        except:
            return (1, key)

    # find changes between states (add to fep_changes dict)
    for bati_type, batis in bati_map.iteritems():
        for bati_key, bati_all_states in sorted(batis.items(),
                                                key=lambda (key, val): \
                                                        _bati_sort(key, val)):

            # bond/angle/.. breaking or forming
            if None in bati_all_states:
                fep_changes[bati_type][bati_key] = bati_all_states

                # add bond atoms to "reactive atoms" set
                # and replace the bond parameter with a Morse type
                if bati_type == "bonds":
                    for bati in bati_all_states:
                        if bati != None:
                            fep_reacting_atoms |= set(bati.atoms)
                            # the bond parameter is replaced with a Morse
                            # parameter (_FepPrmMorse)
                            prm_id = bati.prm.prm_id
                            try:
                                bati.prm = fep_morse_prms[prm_id]
                            except KeyError:
                                bati.prm = _FepPrmMorse(bati.prm)
                                fep_morse_prms[prm_id] = bati.prm

            # the actual values of the parameters are not exactly the same
            else:
                tmp = [bati_all_states[0].prm.strval == bati.prm.strval
                                                for bati in bati_all_states]
                if not all(tmp):
                    fep_changes[bati_type][bati_key] = bati_all_states

    # DEBUG
    # for k,v in fep_changes.iteritems():
        # print k
        # try:
            # for k2,(v1,v2) in v.iteritems():
                # print k2,v1,v2
        # except:
            # for (v1,v2) in v:
                # print v1,v2


    # add parameters of changing batis to fep_types
    for bati_type in bati_map:
        for bati_all_states in fep_changes[bati_type].values():
            prms = [bati.prm for bati in bati_all_states if bati != None]
            for prm in prms:
                if prm not in fep_types[bati_type]:
                    fep_types[bati_type].append(prm)

    # DEBUG
    # for k,v in fep_types.iteritems():
        # print k
        # for v2 in v:
            # print v2

    # add reactive atoms from states that have bond==None to fep_reacting_atoms
    for atom_all_states in fep_changes["atoms"]:
        for atom in atom_all_states:
            if atom in fep_reacting_atoms:
                fep_reacting_atoms |= set(atom_all_states)




    ########################
    # Prepare the output
    ########################


    fep_l = {"atoms": [],
             "atom_types": [],
             "qcp_mass": [],
             "change_atoms": [],
             "change_charges": [],
             "soft_pairs": [],
             "off_diagonals": [],
             "bond_types": [],
             "change_bonds": [],
             "angle_types": [],
             "change_angles": [],
             "torsion_types": [],
             "change_torsions": [],
             "improper_types": [],
             "change_impropers": []}

    ####################
    # ATOMS
    # CHANGE_ATOMS
    # CHANGE_CHARGES
    ####################
    format_atoms = "{:<15} {:<10}     #  {:<15} {:<15} {:>3}"
    format_ch_atoms = "{:<10} " + " {:<12}"*num_evb_states + "    #  {:<}"
    format_ch_crgs = "{:<10} " + " {:12}"*num_evb_states + "    #  {:<10}"\
                   + " {:>12}"*(num_evb_states-1)
    format_qcp = "{:<15} {:<10}     #  {:<10}"

    fep_l["atoms"].append(format_atoms.format("#Q index", "PDB index",
                                              "St.1 PDB_ID", "St.1 LIB_ID", ""))

    tmp = ["#Q index"]
    tmp.extend(["Type st.{}".format(n+1) for n in range(num_evb_states)])
    tmp.append("St.1 PDB_ID")
    fep_l["change_atoms"].append(format_ch_atoms.format(*tmp))

    tmp = ["#Q index"]
    tmp.extend(["Charge st.{}".format(n+1) for n in range(num_evb_states)])
    tmp.append("St.1 PDB_ID")
    tmp.extend(["dq({}->{})".format(n+1, n+2) for n in range(num_evb_states-1)])
    fep_l["change_charges"].append(format_ch_crgs.format(*tmp))

    if fep_qcp_atoms:
        fep_l["qcp_mass"].append("[qcp_mass]")
        fep_l["qcp_mass"].append(format_qcp.format("#Q index", "Mass",
                                                   "St.1 PDB_ID"))

    for i, atom_all_states in enumerate(fep_changes["atoms"]):
        q_index = i + 1
        a = atom_all_states[0]
        pdb_id = "{}.{}".format(a.residue.index, a.name)
        lib_id = "{}.{}".format(a.residue.name, a.name)

        # atoms
        reacting_flag = " !" * bool([atom for atom in atom_all_states
                                     if atom in fep_reacting_atoms])
        fep_l["atoms"].append(format_atoms.format(q_index, "$"+pdb_id+"$",
                                                  pdb_id, lib_id,
                                                  reacting_flag))

        # change_atoms
        tmp = [q_index] + [a.prm.prm_id for a in atom_all_states] + [pdb_id]
        fep_l["change_atoms"].append(format_ch_atoms.format(*tmp))

        # charges
        crgs = [float(a.charge) for a in atom_all_states]
        tmp = [q_index] + crgs + [pdb_id] \
            + [crgs[n+1]-crgs[n] for n in range(num_evb_states-1)]
        fep_l["change_charges"].append(format_ch_crgs.format(*tmp))

        # qcp_atoms
        if atom_all_states in fep_qcp_atoms:
            fep_l["qcp_mass"].append(format_qcp.format(q_index,
                                                       "<FIX>", pdb_id))


    ###############
    # ATOM_TYPES
    ###############
    format_atypes = "{:<12} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}"
    if forcefield == "amber":
        fep_l["atom_types"].append(format_atypes.format("#Atom_type", "LJ_Rm",
                                                        "LJ_eps", "SP_Ci",
                                                        "SP_ai", "LJ_Rm",
                                                        "LJ_eps_14", "mass"))
    else:
        fep_l["atom_types"].append(format_atypes.format("#Atom_type", "LJ_A",
                                                        "LJ_B", "SP_Ci",
                                                        "SP_ai", "LJ_A_14",
                                                        "LJ_B_14", "mass"))

    fep_reacting_atoms_prms = [a.prm for a in fep_reacting_atoms]
    for prm in fep_types["atoms"]:
        sp_c = 1
        sp_a = 2.5
        if prm in fep_reacting_atoms_prms:
            sp_c = "<FIX>"
        if forcefield == "amber":
            lj1, lj2 = prm.lj_R, prm.lj_eps
            lj3, lj4 = lj1, round(lj2/1.2, 4)
        else:
            lj1, lj2 = prm.lj_A, prm.lj_B
            lj3, lj4 = round(lj1/(2**0.5), 4), round(lj2/(2**0.5), 4)
        fep_l["atom_types"].append(format_atypes.format(prm.prm_id, lj1, lj2,
                                                        sp_c, sp_a, lj3, lj4,
                                                        prm.mass))

    ###############
    # BOND_TYPES
    ###############
    format_hbonds = "{:<8}      {:>10}       {:>10}   # {}"
    format_mbonds = "{:<8} {:^10} {:^10} {:>10}   # {}"
    fep_l["bond_types"].append("## Harmonic format")
    fep_l["bond_types"].append(format_hbonds.format("#Index", "Fc",
                                                    "r0", "PRM_ID"))
    fep_l["bond_types"].append("## Morse format")
    fep_l["bond_types"].append(format_mbonds.format("#Index", "D", "alpha",
                                                    "r0", "PRM_ID"))

    for i, bond_type in enumerate(fep_types["bonds"]):
        b_index = i + 1
        if isinstance(bond_type, _FepPrmMorse):
            prm_id = "-".join(bond_type.harmonic_prm.prm_id.split())
            tmp = format_mbonds.format(b_index, "<FIX_D>", "<FIX_a>",
                                       "<FIX_r0>", prm_id)
            fep_l["bond_types"].append(tmp)
        else:
            prm_id = "-".join(bond_type.prm_id.split())
            tmp = format_hbonds.format(b_index, bond_type.fc, bond_type.r0,
                                       prm_id)
            fep_l["bond_types"].append(tmp)


    ###############
    # CHANGE_BONDS
    ###############
    format_bondch = "{:<10} {:<10} " + "{:^5} "*num_evb_states + "  # {}"
    tmp = ["#Atom1", "Atom2"]
    tmp.extend(["St.{}".format(n+1) for n in range(num_evb_states)])
    tmp.append("St.1 PDB_IDs")
    fep_l["change_bonds"].append(format_bondch.format(*tmp))

    for bond_key, bond_all_states in fep_changes["bonds"].iteritems():
        # bond_key == "PDB_ID1 PDB_ID2"
        prm_indexes = []
        for b in bond_all_states:
            if b == None:
                prm_indexes.append(0)
            else:
                btype_index = fep_types["bonds"].index(b.prm) + 1
                prm_indexes.append(btype_index)

        placeholders = ["${}$".format(a) for a in bond_key.split()]
        pdb_id = "-".join(bond_key.split())

        tmp = placeholders + prm_indexes + [pdb_id]
        fep_l["change_bonds"].append(format_bondch.format(*tmp))


    ###############
    # ANGLE_TYPES
    ###############
    format_angles = "{:<8} {:>10} {:>10}   # {}"
    fep_l["angle_types"].append(format_angles.format("#Index", "Fc",
                                                     "theta0", "PRM_ID"))
    for i, angle_type in enumerate(fep_types["angles"]):
        an_index = i + 1
        prm_id = "-".join(angle_type.prm_id.split())
        tmp = format_angles.format(an_index, angle_type.fc,
                                   angle_type.theta0, prm_id)
        fep_l["angle_types"].append(tmp)


    #################
    # CHANGE_ANGLES
    #################
    format_angch = "{:<10} {:<10} {:<10} " + "{:^5} "*num_evb_states + "  # {}"
    tmp = ["#Atom1", "Atom2", "Atom3"]
    tmp.extend(["St.{}".format(n+1) for n in range(num_evb_states)])
    tmp.append("St.1 PDB_IDs")
    fep_l["change_angles"].append(format_angch.format(*tmp))

    for angle_key, angle_all_states in fep_changes["angles"].iteritems():
        # angle_key == "PDB_ID1 PDB_ID2 PDB_ID3"
        prm_indexes = []
        for ang in angle_all_states:
            if ang == None:
                prm_indexes.append(0)
            else:
                atype_index = fep_types["angles"].index(ang.prm) + 1
                prm_indexes.append(atype_index)

        placeholders = ["${}$".format(a) for a in angle_key.split()]
        pdb_id = "-".join(angle_key.split())

        tmp = placeholders + prm_indexes + [pdb_id]
        fep_l["change_angles"].append(format_angch.format(*tmp))




    #################
    # TORSION_TYPES
    #################
    format_torsions = "{:<8} {:>10} {:>10} {:>10}   # {}"
    fep_l["torsion_types"].append(format_torsions.format("#Index", "Fc",
                                                         "mult", "psi0",
                                                         "PRM_ID"))
    tor_index = 1
    tor_indexes = []
    for i, torsion_type in enumerate(fep_types["torsions"]):
        prm_id = "-".join(torsion_type.prm_id.split())
        prm_indexes = []
        for fc, per, psi0, npath in torsion_type.get_prms():
            fc = fc/npath
            tmp = format_torsions.format(tor_index, fc, per, psi0, prm_id)
            fep_l["torsion_types"].append(tmp)
            prm_indexes.append(tor_index)
            tor_index += 1
        tor_indexes.append(prm_indexes)

    ###################
    # CHANGE_TORSIONS
    ###################
    format_torch = "{:<10} {:<10} {:<10} {:<10} " \
                 + "{:^5} "*num_evb_states + "  # {}"

    tmp = ["#Atom1", "Atom2", "Atom3", "Atom4"]
    tmp.extend(["St.{}".format(n+1) for n in range(num_evb_states)])
    tmp.append("St.1 PDB_IDs")
    fep_l["change_torsions"].append(format_torch.format(*tmp))

    for torsion_key, torsion_all_states in fep_changes["torsions"].iteritems():
        # torsion_key == "PDB_ID1 PDB_ID2 PDB_ID3 PDB_ID4"
        for state, tor in enumerate(torsion_all_states):
            if tor == None:
                continue

            for i in range(len(tor.prm.fcs)):
                tprm_index = fep_types["torsions"].index(tor.prm)
                ttype_index = tor_indexes[tprm_index][i]

                prm_indexes = [0 for _ in range(len(torsion_all_states))]
                prm_indexes[state] = ttype_index

                placeholders = ["${}$".format(t) for t in torsion_key.split()]
                pdb_id = "-".join(torsion_key.split())

                tmp = placeholders + prm_indexes + [pdb_id]
                fep_l["change_torsions"].append(format_torch.format(*tmp))



    #################
    # IMPROPER_TYPES
    #################
    format_impropers = "{:<8} {:>10} {:>10}   # {}"
    fep_l["improper_types"].append(format_impropers.format("#Index",
                                                           "Fc", "phi0",
                                                           "PRM_ID"))
    for i, improper_type in enumerate(fep_types["impropers"]):
        imp_index = i + 1
        prm_id = "-".join(improper_type.prm_id.split())
        tmp = format_impropers.format(imp_index, improper_type.fc,
                                      improper_type.phi0, prm_id)
        fep_l["improper_types"].append(tmp)


    ###################
    # CHANGE_IMPROPERS
    ###################
    format_impch = "{:<10} {:<10} {:<10} {:<10} " \
                 + "{:^5} "*num_evb_states + "  # {}"

    tmp = ["#Atom1", "Atom2", "Atom3", "Atom4"]
    tmp.extend(["St.{}".format(n+1) for n in range(num_evb_states)])
    tmp.append("St.1 PDB_IDs")
    fep_l["change_impropers"].append(format_impch.format(*tmp))

    for improper_key, improper_all_states in fep_changes["impropers"].iteritems():
        # improper_key == "PDB_ID1 PDB_ID2 PDB_ID3 PDB_ID4"
        prm_indexes = []
        for imp in improper_all_states:
            if imp == None:
                prm_indexes.append(0)
            else:
                itype_index = fep_types["impropers"].index(imp.prm) + 1
                prm_indexes.append(itype_index)

        placeholders = ["${}$".format(i) for i in improper_key.split()]
        pdb_id = "-".join(improper_key.split())

        tmp = placeholders + prm_indexes + [pdb_id]
        fep_l["change_impropers"].append(format_impch.format(*tmp))

    ##############
    # SOFT_PAIRS
    ##############
    for bond_key, bond_all_states in fep_changes["bonds"].iteritems():
        if None in bond_all_states:
            for state, bond in enumerate(bond_all_states):
                if bond == None:
                    continue
                atoms_in_state = [atom_all_states[state] for atom_all_states \
                                                        in fep_changes["atoms"]]
                a1_qindex = atoms_in_state.index(bond.atoms[0]) + 1
                a2_qindex = atoms_in_state.index(bond.atoms[1]) + 1
                fep_l["soft_pairs"].append("{:10} {:10}".format(a1_qindex,
                                                                a2_qindex))


    for k in fep_l.keys():
        fep_l[k] = "\n".join(fep_l[k])

    fepstr = """\
# Generated with Qtools, version {version}
# Date: {date}
# CWD: {cwd}
# CMDline: {cmd}
#

[FEP]
states {states}

[atoms]
{atoms}

[atom_types]
{atom_types}

[change_atoms]
{change_atoms}

[change_charges]
{change_charges}

[soft_pairs]
{soft_pairs}

[off_diagonals]
# State_i State_j  Atom1  Atom2  A_ij  mu_ij
#
## Example1, Hij=H12=0 (not known in advance)
## 1 2  13 14  0  0
## Example2, Hij=H12=C*exp(-mu * r_13_14)  (C=20.0, mu=0.45)
## 1 2  13 14  20.0  0.45
#
<FIX>

[bond_types]
{bond_types}

[change_bonds]
{change_bonds}

[angle_types]
{angle_types}

[change_angles]
{change_angles}

[torsion_types]
{torsion_types}

[change_torsions]
{change_torsions}

[improper_types]
{improper_types}

[change_impropers]
{change_impropers}

{qcp_mass}
""".format(states=num_evb_states, date=time.ctime(), cmd=" ".join(sys.argv),
           cwd=os.getcwd(), version=__version__, **fep_l)

    return fepstr

