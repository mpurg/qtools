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

import re
import logging
from collections import OrderedDict

from Qpyl.common import __version__, raise_or_log

logger = logging.getLogger(__name__)


class QLibError(Exception):
    pass

class QLib(object):
    """Class for reading and writing Q library files.

    Also supports parsing oplsaa_ffld and amber_lib.

    Args:
       ff_type (string):  Either 'oplsaa' or 'amber'
       ignore_errors (boolean):  Optional, default is False.
                                 If set to True, some non-vital
                                 exceptions are logged instead.
    """

    def __init__(self, ff_type, ignore_errors=False):
        self.ignore_errors = ignore_errors
        supported_ff = ['oplsaa', 'amber']
        ff_type = ff_type.lower()
        if ff_type not in supported_ff:
            raise QLibError("Force field type '{}' not supported. Use {}"
                            .format(ff_type, " or ".join(supported_ff)))

        self.ff_type = ff_type
        self.residue_dict = OrderedDict()


    def _add_residue(self, residue):
        """Internal method for adding residues"""
        # check for duplicates
        if residue.name in self.residue_dict:
            # check if they are the same...
            orig = self.residue_dict[residue.name].get_str()
            new = residue.get_str()
            if orig != new:
                raise_or_log("Duplicate library entries for residue '{}' "
                             "with different parameters"
                             .format(residue.name),
                             QLibError, logger, self.ignore_errors)
            else:
                logger.info("Duplicate library entry for residue '{}'"
                            .format(residue.name))

        self.residue_dict[residue.name] = residue


    def check_valid(self):
        """Call 'check_valid' on each _LibResidue object.

        Raises QLibError if something is not cool with the residues.
        (see _LibResidue for more info)
        """
        for residue in self.residue_dict.values():
            residue.check_valid()


    def read_lib(self, libfile):
        """Read and parse a Q library (.lib) file

        Add new residues to QLib.residue_dict as _LibResidue objects

        Args:
            libfile (string):  name/path of Q lib file
        """
        residues = []

        with open(libfile, 'r')  as lib:
            section = None
            for lnumber, line in enumerate(lib.readlines()):
                lnumber += 1
                # remove comments
                line = re.split("#|\*|\!", line, 1)[0].strip()
                if line == "":
                    continue
                if line.startswith("{"):
                    # get the 3 letter code and make a new object
                    resname = line.split("}")[0].strip("{} ").upper()
                    residues.append(_LibResidue(resname, self))
                    continue
                if line.startswith("["):
                    section = line.split("]")[0].strip("[] ").lower()
                    continue
                if not residues or not section:
                    raise QLibError("Line #{} in LIB file '{}' is not a "
                                    "comment and is not inside a {{XXX}} "
                                    "section and [XXX...X] subsection:\n{}"
                                    .format(lnumber, libfile, line))

                if residues:
                    residue = residues[-1]

                if section == "atoms":
                    try:
                        atom_name, atom_type, atom_charge = line.split()[1:4]
                        residue.atoms.append(_LibAtom(atom_name,
                                                      atom_type,
                                                      float(atom_charge),
                                                      residue))
                    except ValueError:
                        raise QLibError("Line #{} in LIB file '{}' couldn't"
                                        "be parsed (should look like this "
                                        "'aindex aname atype charge ...'):\n{}"
                                        .format(lnumber, libfile, line))

                elif section == "bonds":
                    try:
                        a1, a2 = line.split()
                    except ValueError:
                        raise QLibError("Line #{} in LIB file '{}' couldn't "
                                        "be parsed (should look like this "
                                        "'atom1 atom2'):\n{}"
                                        .format(lnumber, libfile, line))
                    else:
                        residue.bonds.append((a1, a2))


                elif section == "impropers":
                    try:
                        a1, a2, a3, a4 = line.split()
                    except ValueError:
                        raise QLibError("Line #{} in LIB file '{}' couldn't be"
                                        "parsed (should look like this "
                                        "'atom1 atom2 atom3 atom4'):\n{}"
                                        .format(lnumber, libfile, line))
                    else:
                        residue._add_improper(a2, [a1, a3, a4])


                elif section == "charge_groups":
                    if self.ff_type == "amber":
                        raise QLibError("'Charge_groups' section is not "
                                        "compatible with 'amber' forcefield")
                    cgrp = line.split()
                    residue.charge_groups.append(cgrp)


                elif section == "info":
                    try:
                        key, value = line.split()
                    except ValueError:
                        raise QLibError("Line #{} in LIB file '{}' couldn't be"
                                        "parsed (should look like this "
                                        "'keyword value'):\n{}"
                                        .format(lnumber, libfile, line))
                    else:
                        residue.info[key] = value

                elif section == "connections":
                    residue.connections.append(" ".join(line.split()))

                elif section == "build_rules":
                    residue.build_rules.append(" ".join(line.split()))
                else:
                    logger.warning("Unsupported section in '{}': {}"
                                   "".format(libfile, section))

        for residue in residues:
            self._add_residue(residue)

        if not residues:
            raise_or_log("No residues found", QLibError,
                         logger, self.ignore_errors)



    def read_amber_lib(self, libfile):
        """Read and parse an Amber library (.lib) file

        Add new residues to QLib.residue_dict as _LibResidue objects

        Args:
            libfile (string):  name/path of Amber lib file

        """

        if self.ff_type != "amber":
            raise QLibError("Function not supported with force field"
                            "'{}'".format(self.ff_type))

        residues = []
        section = None

        with open(libfile) as lib:
            for lnumber, line in enumerate(lib.readlines()):
                lnumber += 1
                line = line.strip()
                if not line:
                    continue
                if line.startswith("!!"):
                    continue # ignore comment

                if residues:
                    residue = residues[-1]

                if line[0] == "!":
                    lf = line.split()[0].split(".")
                    if len(lf) > 2:
                        rname = lf[1].upper()
                        if not residues or residue.name != rname:
                            residues.append(_LibResidue(rname, self))

                        section = lf[3]
                    else: raise QLibError("Line #{} in LIB file '{}' "
                                          "couldn't be parsed:\n{}"
                                          .format(lnumber, libfile, line))
                    continue

                if not section or not residues:
                    continue

                if section == "atoms":
                    lf = line.split()
                    try:
                        name, atype = lf[0].strip('"'), lf[1].strip('"')
                        charge = float(lf[7])
                        atype = atype.replace("*", "star")
                        residue.atoms.append(_LibAtom(name, atype,
                                                      charge,
                                                      residue))
                    except ValueError:
                        raise QLibError("Line #{} in LIB file '{}' "
                                        "couldn't be parsed:\n{}"
                                        .format(lnumber, libfile, line))

                elif section == "connectivity":
                    atomnames = [a.name for a in residue.atoms]
                    try:
                        ai1, ai2 = line.split()[:2]
                        a1 = atomnames[int(ai1)-1]
                        a2 = atomnames[int(ai2)-1]
                        if a1 not in atomnames or a2 not in atomnames:
                            raise QLibError("Undefined atom(s) {} and/or {} "
                                            "mentioned in the connectivity "
                                            "section of '{}', ln.{}."
                                            .format(a1, a2, libfile, lnumber))
                        residue.bonds.append((a1, a2))
                    except ValueError:
                        raise QLibError("Line #{} in LIB file '{}' "
                                        "couldn't be parsed:\n{}"
                                        .format(lnumber, libfile, line))

                elif section == "residueconnect":
                    atomnames = [a.name for a in residue.atoms]
                    try:
                        ai1, ai2 = line.split()[:2]
                        a1 = atomnames[int(ai1)-1]
                        a2 = atomnames[int(ai2)-1]
                        if a1 not in atomnames or a2 not in atomnames:
                            raise QLibError("Undefined atom(s) {} and/or {} "
                                            "mentioned in the residueconnect"
                                            " section of '{}', ln.{}."
                                            .format(a1, a2, libfile, lnumber))
                        residue.connections.append("head " + a1)
                        residue.connections.append("tail " + a2)
                    except ValueError:
                        raise QLibError("Line #{} in LIB file '{}' "
                                        "couldn't be parsed:\n{}"
                                        .format(lnumber, libfile, line))

        for residue in residues:
            self._add_residue(residue)
        if not residues:
            raise_or_log("No residues found", QLibError,
                         logger, self.ignore_errors)




    def read_mol2(self, mol2_file):
        """Read and parse a mol2 file.

        Add the residues to QLib.residue_dict as _LibResidue objects

        Args:
            mol2_file (string):  name/path of mol2 file

        """

        if self.ff_type != "amber":
            raise QLibError("Function not supported with "
                            "force field '{}'".format(self.ff_type))

        aindex, old_aindex = None, None
        rindex, old_rindex = None, None
        residues = []
        section = None
        for line in open(mol2_file, 'r').readlines():

            if line.startswith("@<TRIPOS>"):
                section = line.replace("@<TRIPOS>", "").strip()
                if section == "MOLECULE":
                    # lookup for bonds section
                    lookup_aindex = {}
                continue

            if section == "ATOM":

                if aindex != None:
                    old_aindex = aindex
                if rindex != None:
                    old_rindex = rindex

                lf = line.split()

                aindex, aname = int(lf[0]), lf[1]
                atype = lf[5]
                rindex = int(lf[6])
                rname = lf[7][0:4].upper()
                charge = float(lf[8])


                if old_aindex != None and aindex - old_aindex != 1:
                    logger.warning("Bad Mol2 format - atom "
                                   "index {} followed by {}"
                                   .format(old_aindex, aindex))

                if old_rindex == None or old_rindex != rindex:
                    residues.append(_LibResidue(rname, self))

                    if old_rindex and rindex - old_rindex != 1:
                        logger.warning("Bad Mol2 format - residue "
                                       "index {} followed by {}"
                                       .format(old_rindex, rindex))


                lib_residue = residues[-1]
                lib_residue.atoms.append(_LibAtom(aname, atype,
                                                  charge, lib_residue))
                lookup_aindex[lf[0]] = (aname, rindex, lib_residue)


            elif section == "BOND":

                lf = line.split()
                aname1, rindex1, residue1 = lookup_aindex[lf[1]]
                aname2, rindex2, residue2 = lookup_aindex[lf[2]]

                if rindex1 < rindex2:
                    residue1.connections.append("tail " + aname1)
                    residue2.connections.append("head " + aname2)
                elif rindex1 > rindex2:
                    residue1.connections.append("head " + aname1)
                    residue2.connections.append("tail " + aname2)
                else:
                    residue1.bonds.append((aname1, aname2))

        for residue in residues:
            self._add_residue(residue)
        if not residues:
            raise_or_log("No residues found", QLibError,
                         logger, self.ignore_errors)



    def read_prepin_impropers(self, prepin_file):
        """Read and parse an Amber prepin (.prepi) file

        NOTE: Extracts only improper definitions for residues already
              defined in QLib.residue_dict
              (usually obtained with read_amber_lib or read_mol2)
        """

        if self.ff_type != "amber":
            raise QLibError("Function not supported with "
                            "force field '{}'".format(self.ff_type))

        if len(self.residue_dict) == 0:
            raise_or_log("Cannot add impropers to empty library",
                         QLibError, logger, self.ignore_errors)

        # get all the impropers from the file
        section = None
        residue = None
        impropers = {}
        with open(prepin_file) as prepi:
            for line in prepi.readlines():
                lf = line.split()
                if not lf:
                    section = None
                    continue

                if len(lf) >= 2 and lf[1] == "INT":
                    residue = lf[0].upper()
                    if not residue in impropers:
                        impropers[residue] = []
                    continue

                if not residue:
                    continue

                if "DONE" in lf:
                    residue = None
                    continue

                if lf[0] == "IMPROPER":
                    section = "improper"
                    continue

                if section == "improper":
                    # example line (center atom is N)
                    # -M CA N H
                    imp = " ".join(line.split())
                    impropers[residue].append(imp)


        res_not_in_lib = set()
        # add them to the residues in the library
        for resname, imps in impropers.iteritems():
            for imp in imps:
                try:
                    # different naming convention for head and tail
                    imp = imp.replace("-M", "-C").replace("+M", "+N")
                    other_atoms = imp.split()
                    # third one is center
                    center_atom = other_atoms.pop(2)
                    self.residue_dict[resname]._add_improper(center_atom,
                                                             other_atoms)
                except KeyError:
                    res_not_in_lib.add(resname)
        for resname in res_not_in_lib:
            logger.warning("Prepi: Residue '{}' NOT found in library. "
                           "Impropers will be ignored. Check that residues "
                           "have the same names as in the library "
                           "(.lib, .mol2).".format(resname))


    def read_ffld(self, ffld_file, qstruct):
        """Read and parse a Macromodel's FFLD file for oplsaa parameters.

        Args:
            ffld_file (string):  path/name of ffld file
            qstruct (qstructure.QStruct):  object created with the same
                                           structure file as the ffld

        The second argument is a QStruct object created with
        the pdb or mol2 file that was used to create the ffld file.
        It's needed to map atoms and residue to their names in the structure.

        Note:
        The atoms in oplsaa have atom types defined such that the vdw types
        don't map 1 to 1 with bonding types ('symbols'). Also, atom type
        identifiers ('type') do not seem to uniquely map vdw-symbol
        combinations. Examples:

        Vdws and symbols map N:M
        type   vdw   symbol
        135    C1    CTH
        224    C1    CT1
        145    C4    CA
        351    C5    CA

        Type does not map them uniquely
        type   vdw   symbol
        135    C1    CT
        135    C1    CTH
        181    C1    CT

        One type can have different vdw
        type   vdw   symbol
        140    H1    HC
        140    H2    HC

        To prevent any clashes, this function creates Q atom types by
        combining 'symbol', 'vdw' and 'type':
         CA_C4_145
         CA_C5_351
        CTH_C1_135
         CT_C1_135
         CT_C1_181

        """

        if self.ff_type != "oplsaa":
            raise QLibError("Function not supported with "
                            "force field '{}'".format(self.ff_type))

        # keys are ffld atom names, values are tuples:
        # (StructAtom, LibResidue)
        lookup_aname = {}

        residues = []
        section = None

        for line in open(ffld_file, 'r').readlines():
            line = line.strip()
            if (line == "") or ("------" in line):
                continue
            elif line.startswith("atom   type  vdw  symbol"):
                section = "ATOMS"
                continue
            elif line.startswith("Stretch            k"):
                section = "BONDS"
                continue
            elif line.startswith("Bending                      k"):
                section = "ANGLES"
                continue
            elif line.startswith("proper Torsion"):
                section = "TORSIONS"
                continue
            elif line.startswith("improper Torsion"):
                section = "IMPROPERS"
                continue

            if section == "ATOMS":
#
#  C1      135  C1   CT      -0.0175   3.5000   0.0660 high   C: alkanes
#
                lf = line.split()
                name, type_, vdw, symbol = lf[0:4]
                charge, sigma, epsilon = map(float, lf[4:7])
                quality, comment = lf[7], lf[8:]


                aindex_struct = len(lookup_aname)
                atom_struct = qstruct.atoms[aindex_struct]
                residue_struct = atom_struct.residue

                # check if element from ffld matches the one in the structure
                # (just the first letters)
                if name[0].lower() != atom_struct.name[0].lower():
                    raise_or_log("Atom element mismatch, possible wrong "
                                 "order of atoms: '{}' (struct) '{}' (ffld)"
                                 .format(atom_struct.name, name),
                                 QLibError, logger, self.ignore_errors)

                # crate new library entry if first atom or different residue
                if not residues or \
                      residue_struct != qstruct.atoms[aindex_struct-1].residue:
                    residues.append(_LibResidue(residue_struct.name, self))

                residue = residues[-1]

                # append the atom to the residue
                atom_name = atom_struct.name
                atom_type = "{}_{}_{}".format(symbol, vdw, type_)
                residue.atoms.append(_LibAtom(atom_name, atom_type, charge,
                                              residue))

                lookup_aname[name] = (atom_struct, residue)

            elif section == "BONDS":
#
#  C1      H2      340.00000    1.09000   high      140  0   CT  -HC    ==> CT  -HC
#
                lf = line.split()
                atom1_struct, residue1 = lookup_aname[lf[0]]
                atom2_struct, residue2 = lookup_aname[lf[1]]

                atom_name1 = atom1_struct.name
                atom_name2 = atom2_struct.name
                rindex1 = atom1_struct.residue.index
                rindex2 = atom2_struct.residue.index

                if rindex1 < rindex2:
                    residue1.connections.append("tail " + atom_name1)
                    residue2.connections.append("head " + atom_name2)
                elif rindex1 > rindex2:
                    residue1.connections.append("head " + atom_name1)
                    residue2.connections.append("tail " + atom_name2)
                else:
                    residue1.bonds.append((atom_name1, atom_name2))

            elif section == "IMPROPERS":
#
#  C21     C22     C20     O19       2.200   high   aromatic atom
#
                lf = line.split()
                atoms = [lookup_aname[aname][0] for aname in lf[0:4]]
                center_atom = atoms.pop(2)  # third pos
                rindex = center_atom.residue.index

                atom_names = []
                for atom in atoms:
                    name = atom.name
                    if atom.residue.index > rindex:
                        name = "+" + name
                    elif atom.residue.index < rindex:
                        name = "-" + name
                    atom_names.append(name)

                atom_names.sort()
                residue = lookup_aname[lf[2]][1]
                residue._add_improper(center_atom.name, atom_names)



        for residue in residues:
            # add default charge group (all atoms)
            atom_names = [a.name for a in residue.atoms]
            residue.charge_groups.append(atom_names)
            # add residue to library
            self._add_residue(residue)

        if not residues:
            raise_or_log("No residues found", QLibError,
                         logger, self.ignore_errors)
        elif len(residues) > 1:
            logger.warning("Possible connections and improper "
                           "definitions of first and last residue "
                           "can be missing in the output!")



    def get_string(self):
        """
        Return the whole Q library in string format.

        To print specific residues, use
          self.residue_dict[name].get_str()
        """
        out = ""
        for residue in sorted(self.residue_dict.values(),
                              key=lambda x: x.name):
            out += residue.get_str()
            out += "*" + "-"*80 + "\n"
        return out

class _LibAtom(object):
    """Class containing library atom information.

    Args:
        name (string):   atom name (same as in structure)
        atype (string):  atom type (same as in parameter file)
        charge (float):  partial charge
        residue (_LibResidue):  parent residue object
        comment (string):  optional, comment
    """
    def __init__(self, name, atype, charge, residue, comment=None):
        self.name = name
        self.atom_type = atype
        self.charge = charge
        self.residue = residue
        self._comment = comment

    @property
    def comment(self):
        if self._comment:
            return " # {}".format(self._comment)
        else:
            return ""

    @comment.setter
    def comment(self, value):
        if self._comment:
            self._comment = "{} # {}".format(value, self._comment)
        else:
            self._comment = value

    def __repr__(self):
        return "_LibAtom(name={}, atype={}, charge={}, residue={})"\
               "".format(self.name, self.atom_type, self.charge,
                         self.residue.name)


class _LibResidue(object):
    """Class containing library residue entries.

    Provides direct access to properties of specific
      residues (atoms, charges, bonds, impropers...) and
      several useful functions for checking the
      integrity of the entry, rescaling charges and for
      printing it out.

    Args:
        resname (string):   three letter residue identifier (GLY, ARG...)
        library (QLib):   parent object
    """
    
    INFO_KEYS = ("solvent", "SYBILtype", "density")
    BUILD_RULES = ("torsion")
        
    def __init__(self, resname, library):
        self.name = resname.upper()
        self.atoms = []   # [ _LibAtom, _LibAtom... ]
        self.bonds = []   # [ ("CA","CB"), ("CA", "HA1"), ... ]
        self.impropers = []   # [ ("C1","C2","C3","C4"), ... ]  # C2 is center
        self.connections = []   # [ "head N", "tail C" ]
        self.charge_groups = []   # [ ["O1", "H1"], ["C1","H2","H3"], ... ]
        self.build_rules = []   # [ "torsion HE1 OE1 CD OE2 0", ... ]
        self.info = {}   # { "SYBILtype": "RESIDUE", "solvent": 1, ... }
        self.library = library

    def _add_improper(self, center_atom, other_atoms):
        """Add an explicit improper to the library entry. (private)

        Args:
            center_atom (string)
            other_atoms (list of strings)

        """
        #
        # Sort non-center atoms based on atom names (NOT atom types).
        # There is some ambiguity in improper definitions in Amber
        # due to sorting by atom types... Let's not propagate bad practice...
        # Also, sort head and tail atoms without + and - signs
        #
        imp = list(other_atoms)
        imp.sort(key=lambda x: x.strip("+-"))
        imp.insert(1, center_atom)  # center atom is in position 2 in Q

        if imp not in self.impropers:
            self.impropers.append(imp)
        else:
            logger.warning("Overwriting existing improper "
                           "'{}' in residue '{}'"
                           .format(imp, self.name))


    def check_valid(self):
        """Checks for duplicates, missing atoms, integer charges, ...

        Raises QLibError (or logs if QLib.ignore_errors is True).
        """

        names = [a.name for a in self.atoms]
        types = [a.atom_type for a in self.atoms]
        charges = [a.charge for a in self.atoms]
        # check for duplicates
        for i, name in enumerate(names):
            if name in names[i+1:]:
                raise QLibError("Duplicate atom name '{}' in residue '{}'"\
                                .format(name, self.name))
        # check charge groups for integer charges (and bad atom names)
        if self.charge_groups:
            for cg in self.charge_groups:

                net_charge = 0
                for aname in cg:
                    if aname not in names:
                        raise QLibError("Undefined atom {} in charge group"
                                        " '{}'".format(aname, " ".join(cg)))
                    net_charge += charges[names.index(aname)]

                if abs(net_charge - round(net_charge)) > 1e-7:
                    raise_or_log("Net charge of charge group '{}'"
                                 " in residue '{}' not integer: {}"\
                                 .format(" ".join(cg),
                                         self.name,
                                         net_charge),
                                 QLibError, logger, self.library.ignore_errors)
        # check the whole residue for integer charge
        else:
            net_charge = sum(charges)
            if abs(net_charge - round(net_charge)) > 1e-7:
                raise_or_log("Net charge of residue '{}' not integer: "
                             "{}".format(self.name, net_charge),
                             QLibError, logger, self.library.ignore_errors)

        # check bonds for bad atom names
        for bond_atoms in self.bonds:
            for aname in bond_atoms:
                if aname not in names:
                    raise QLibError("Undefined atom {} in bond '{}'"\
                                    .format(aname, " ".join(bond_atoms)))

        # check impropers for bad atom names
        for imp_atoms in self.impropers:
            for aname in imp_atoms:
                if aname not in names and aname not in ["-C", "+N"]:
                    raise QLibError("Undefined atom '{}' in improper '{}'"\
                                    .format(aname, " ".join(imp_atoms)))

        # check keywords in info section
        low_keys = (x.lower() for x in self.INFO_KEYS)
        for keyword in self.info:
            if keyword.lower() not in low_keys:
                raise QLibError("Keyword '{}' in [info] section not "
                                "not supported (residue '{}')."
                                "".format(keyword, self.name))

        # check build rules
        for build_rule in self.build_rules:
            br = build_rule.split()
            try:
                typ, anames, value = br[0], br[1:5], br[5]
            except IndexError:
                raise QLibError("Invalid build rule '{}'"
                                "".format(build_rule))
            if typ not in self.BUILD_RULES:
                raise QLibError("Build_rule '{}' not supported".format(typ))

            for aname in anames:
                if aname not in names:
                    raise QLibError("Undefined atom '{}' in build_rule '{}'"
                                    "".format(aname, build_rule))




    def rescale(self, atoms, threshold):
        """Rescale the charges of a group of atoms to nearest integer value.

        Args:
            atoms (list):  List of atom names
            threshold:  Maximum difference between sum and nearest integer.
                        If the difference is greater than this value,
                        QLibError, is raised.

        Used primarily in oplsaa for charge groups, or to fix rounding errors.

        According to this formula:
        (credits go to M.Repic)
        q_i = q_i_initial - weight(q_i) * diff

        where:
        weight(q_i) = abs(q_i) / sum([abs(q) for q in all_q])
        diff = sum(all_q) - round(sum(all_q))

        After rescaling, if the sum is not integer (rounding error), find
        the largest charge and adjust it to remove the excess charge.

        """


        # for some reason I feel obliged to check for stupidity
        if len(set(atoms)) < len(atoms):
            raise QLibError("Duplicate atom names in group '{}'"
                            "".format(" ".join(atoms)))

        # create atom_name: charge mapping for given atoms
        atom_dict = {a.name: round(a.charge, 6) for a in
                     self.atoms if a.name in atoms}

        # check if the given atoms actually exist
        for atom_name in atoms:
            if atom_name not in atom_dict:
                raise QLibError("Atom name '{}' not found in residue '{}'"
                                "".format(atom_name, self.name))

        # calculate absolute charges and sums
        sum_all_q = sum(atom_dict.values())
        sum_abs_all_q = sum([abs(q) for q in atom_dict.values()])
        target = round(sum_all_q)
        diff = sum_all_q - target
        if diff > threshold:
            raise QLibError("Difference between sum of charges and nearest "
                            "integer ({}) is greater than threshold"
                            "".format(diff))

        # rescale the charges
        for atom_name in atoms:
            charge_init = atom_dict[atom_name]
            weight = abs(charge_init) / sum_abs_all_q
            charge = charge_init - weight * diff
            atom_dict[atom_name] = round(charge, 6)

        sum_all_q_new = sum(atom_dict.values())

        # if the net charge is not integer (rounding errors)
        # find the atom with the largest absolute charge
        # that is not chemically equivalent to some other atom
        # and remove the excess charge

        excess = sum_all_q_new - target
        if abs(excess) > 1e-7:
            # only unique atoms, with abs charges
            atom_dict2 = {name: abs(charge) for name, charge in
                          atom_dict.iteritems() if \
                          atom_dict.values().count(charge) == 1}
            # maximum charge atom
            max_ch_atom = max(atom_dict2, key=lambda x: atom_dict2[x])

            atom_dict[max_ch_atom] -= excess

            logger.warning("Excess charge ({}) in group '{}' was "
                           "removed from atom {}."
                           .format(excess, " ".join(atoms), max_ch_atom))


        # change the charge in the actual atom definition
        for atom in self.atoms:
            if atom.name in atom_dict:
                ac_diff = atom_dict[atom.name] - atom.charge
                atom.comment = "{:10.6f} (dq={:+.6f}) {}".format(atom.charge,
                                                                 ac_diff,
                                                                 atom.comment)
                atom.charge = atom_dict[atom.name]




    def get_str(self):
        """Return the Q-lib formatted string for the residue."""

        infol, al, bl, il, cl, brl, col = [], [], [], [], [], [], []

        indent = "        "
        for k, v in self.info.iteritems():
            infol.append(indent + "{:30} {}".format(k, v))
        for i, atom in enumerate(self.atoms):
            al.append("    {:>5d}  {a.name:<5s}  {a.atom_type:<12s} "
                      "{a.charge:>10.6f}{a.comment}".format(i+1, a=atom))
        for bond in self.bonds:
            bl.append(indent + "{b[0]:<5s} {b[1]:s}".format(b=bond))
        for imp in self.impropers:
            tmp = " ".join("{:<5s}".format(a) for a in imp).rstrip()
            il.append(indent + tmp)
        for chgr in self.charge_groups:
            cl.append(indent + " ".join(chgr))
        for br in self.build_rules:
            brl.append(indent + br)
        for conn in self.connections:
            col.append(indent + conn)

        outl = OrderedDict((("info", infol),
                            ("atoms", al),
                            ("bonds", bl),
                            ("impropers", il),
                            ("build_rules", brl),
                            ("connections", col),
                            ("charge_groups", cl)))

        outstr = "{{{}}}\n".format(self.name)
        for section, lines in outl.iteritems():
            if lines:
                outstr += """\
    [{}]
{}
""".format(section, "\n".join(lines))
        return outstr

    def __repr__(self):
        return "_LibResidue({})".format(self.name)


