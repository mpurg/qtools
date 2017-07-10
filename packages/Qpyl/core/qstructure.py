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
from collections import namedtuple

from Qpyl.common import raise_or_log

logger = logging.getLogger(__name__)

PosVector = namedtuple("PosVector", ["x", "y", "z"])


_PLACEHOLDER_RE = re.compile("\$\S+\.\S+\$")
_COMMENTS_RE = re.compile(r"[#\!].*")

def find_placeholders(inputstring):
    """Find atom placeholders of the form $514.C3$

    It ignores comments (characters following # or !)
    
    See also QStruct.convert_placeholders
    """
    tmp = re.sub(_COMMENTS_RE, "", inputstring)
    return _PLACEHOLDER_RE.findall(tmp)


class QStructError(Exception):
    pass

class QStruct(object):
    """
    Class for processing the structure (coordinates)

    Args:
      filename (path to the structure file)
      filetype (type of structure file:  'pdb' or 'mol2')
      ignore_errors (boolean):  Optional, default is False.
                                If set to True, some non-vital
                                exceptions are logged instead.

    In contrast to QLib and QPrm, the 'read' methods in this
    class should not be called since the object should
    contain data from only one structure file.

    The structure data is stored in three lists:
      atoms, residues, molecules
    which contain _StructAtom, _StructResidue and _StructMolecule
    objects.

    """

    def __init__(self, filename, filetype, ignore_errors=False):
        self.ignore_errors = ignore_errors
        FILE_TYPES = {'pdb': self._read_pdb,
                      'mol2': self._read_mol2}
        self.filetype = filetype.lower()
        if self.filetype not in FILE_TYPES.keys():
            raise QStructError("Filetype {} not supported. Use {}"
                               .format(filetype,
                                       " or ".join(FILE_TYPES.keys())))

        self.atoms = []
        self.residues = []
        self.molecules = []
        self.filename = filename
        # TODO: some sort of lookup hashes if needed

        # run the parser function
        FILE_TYPES[self.filetype](filename)

        # check if we actually got something
        for t in ["atoms", "residues", "molecules"]:
            if len(self.__dict__[t]) == 0:
                raise QStructError("No {} found, check file '{}' and"
                                   " filetype '{}'".format(t,
                                                           self.filename,
                                                           self.filetype))




    def _read_mol2(self, mol2_file):
        """
        Read and parse a mol2 file for coordinates.

        Args:
            mol2_file (string):  name/path of file
        """

        molecule = None
        residue = None
        aindex, old_aindex = None, None
        section = None
        for line in open(mol2_file, 'r').readlines():

            if line.startswith("@<TRIPOS>"):
                section = line.replace("@<TRIPOS>", "").strip()
                if section == "MOLECULE":
                    if molecule != None:
                        self.molecules.append(molecule)

                    molecule = _StructMolecule(self)
                continue

            if section == "ATOM":
                if aindex != None:
                    old_aindex = aindex

                lf = line.split()

                aindex, aname = int(lf[0]), lf[1]
                x, y, z = map(float, lf[2:5])
                rindex = int(lf[6])
                rname = lf[7][0:4].upper()

                if old_aindex != None and aindex - old_aindex != 1:
                    raise_or_log("Bad Mol2 format - atom "
                                 "index {} followed by {}"
                                 .format(old_aindex, aindex),
                                 QStructError, logger, self.ignore_errors)

                if not residue or residue.index_struct != rindex:
                    if residue and rindex - residue.index_struct != 1:
                        raise_or_log("Bad Mol2 format - residue "
                                     "index {} followed by {}"
                                     .format(residue.index_struct, rindex),
                                     QStructError, logger, self.ignore_errors)
                    residue = _StructResidue(rindex, rname, molecule, self)
                    self.residues.append(residue)
                    molecule.add_residue(residue)

                atom = _StructAtom(aindex, aname, x, y, z, residue, self)
                self.atoms.append(atom)
                residue.add_atom(atom)

        # append last one after parsing
        if molecule != None and len(molecule.residues) > 0:
            self.molecules.append(molecule)




    def _read_pdb(self, pdb_file):
        """
        Read and parse a PDB file for coordinates.

        Args:
            pdb_file (string):  name/path of file

        """
        # make a new _StructMolecule object
        molecule = _StructMolecule(self)
        # parse the PDB file
        residue = None
        aindex, old_aindex = None, None
        for line in open(pdb_file, 'r').readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if aindex != None:
                    old_aindex = aindex

                aindex = int(line[6:12])
                if old_aindex != None and aindex - old_aindex != 1:
                    raise_or_log("Bad PDB format - atom "
                                 "index {} followed by {}"
                                 .format(old_aindex, aindex),
                                 QStructError, logger, self.ignore_errors)


                aname = line[12:17].strip()
                rname = line[17:20].strip().upper()
                rindex = int(line[22:26])
                x, y, z = map(float, (line[30:38], line[38:46], line[46:54]))

                if not residue or residue.index_struct != rindex:
                    if residue and rindex - residue.index_struct != 1:
                        raise_or_log("Bad PDB format - residue "
                                     "index {} followed by {}"
                                     .format(residue.index_struct, rindex),
                                     QStructError, logger, self.ignore_errors)

                    residue = _StructResidue(rindex, rname, molecule, self)
                    self.residues.append(residue)
                    molecule.add_residue(residue)

                if aname in [a.name for a in residue.atoms]:
                    raise_or_log("Bad PDB format - two atoms with same name "
                                 "({}) in residue {}.{}"
                                 "".format(aname, rname, rindex),
                                 QStructError, logger, self.ignore_errors)

                atom = _StructAtom(aindex, aname, x, y, z, residue, self)
                self.atoms.append(atom)
                residue.add_atom(atom)

            elif line.startswith("TER") or line.startswith("GAP"):
                self.molecules.append(molecule)
                residue = None
                molecule = _StructMolecule(self)

        # append last one if it didn't gave a TER/GAP
        if molecule != None and len(molecule.residues) > 0:
            self.molecules.append(molecule)


    def convert_placeholders(self, inputstring):
        """Convert atom placeholders ($514.C3$) to indexes.

        Placeholders are a combination of the residue id and
          atom name, encapsulated in $$ - $RESID.ATOM_NAME$
        In addition,there are some special values:
          $LAST.ID$ - id of last atom in the system

        Arguments:
          inputstring (string):  string with placeholders (input file contents)

        Returns:
          outputstring (string):  converted string

        """

        id_map = {"{}.{}".format(a.residue.index, a.name): str(a.index)
                                                         for a in self.atoms}
        last_id = "{}.{}".format(self.atoms[-1].residue.index,
                                 self.atoms[-1].name)

        outputstring = ""
        for line in inputstring.split("\n"):
            comment = ""
            if "#" in line:
                i = line.index("#")
                line, comment = line[:i], line[i:]

            c = find_placeholders(line)
            for pid in c:
                pid = pid.strip("$")
                pid2 = pid.replace("LAST.ID", last_id)
                try:
                    padding = (len(pid2)+2  - len(id_map[pid2])) * " "
                except KeyError:
                    raise QStructError("Atom '${}$' does not exist in the pdb "
                                       "structure.".format(pid2))
                line = re.sub("\$" + pid + "\$", id_map[pid2] + padding, line)

            outputstring += line + comment + "\n"
        return outputstring





class _StructAtom(object):
    """Contains structural information for an atom.

    Arguments:
      index_struct (int):  index as written in pdb or mol2
      name (string):  atom name
      x,y,z (float):  coordinates
      residue (_StructResidue):  parent residue object
      structure (_QStruct):  parent structure object

    Property 'index' (int) is the actual 1-based index of the atom
    in the atom list (as opposed to index_struct which was read from
    the file). It should correspond to the index in the generated topology.
    """
    def __init__(self, index_struct, name, x, y, z, residue, structure):
        self.index_struct = int(index_struct)
        self.name = name
        self.coordinates = PosVector(float(x), float(y), float(z))
        self.residue = residue
        self.structure = structure

    @property
    def index(self):
        return self.structure.atoms.index(self) + 1

    def __repr__(self):
        res = self.residue
        mol = res.molecule
        return "_StructAtom: {}.{}.{}".format(mol.index,
                                              res.index,
                                              self.index)


class _StructResidue(object):
    """Contains structural information for a residue.

    Arguments:
      index_struct (int):  index as written in pdb or mol2
      name (string):  residue name
      molecule (_StructMolecule):  parent molecule object
      structure (_QStruct):  parent structure object

    Property 'index' (int) is the actual 1-based index of the residue
    in the residue list (as opposed to index_struct which was read from
    the file). It should correspond to the index in the generated topology.
    """
    def __init__(self, index_struct, name, molecule, structure):
        self.atoms = []
        self.index_struct = int(index_struct)
        self.name = name
        self.molecule = molecule
        self.structure = structure

    @property
    def index(self):
        return self.structure.residues.index(self) + 1

    def add_atom(self, atom):
        self.atoms.append(atom)

    def __repr__(self):
        mol = self.molecule
        return "_StructResidue: {}.{}{}".format(mol.index,
                                                self.name,
                                                self.index)


class _StructMolecule(object):
    """Contains structural information for a molecule.

    Arguments:
      structure (_QStruct):  parent structure object

    Special property is 'index' (int). It is the actual
    1-based index of the molecule in the residue list (as it was appended).
    This should corresponds to the index in the generated topology.
    """
    def __init__(self, structure):
        self.residues = []
        self.structure = structure

    @property
    def index(self):
        return self.structure.molecules.index(self) + 1

    def add_residue(self, residue):
        self.residues.append(residue)

    def __repr__(self):
        return "_StructMolecule: {}".format(self.index)


