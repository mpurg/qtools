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
#
# This script reads in:
# MAP file, PDB structure file and .parm and .lib files produced by FFLD + q_ffld2q.py/mae2q (charges are usually calculated with some other approach like RESP)
# and outputs the Q FEP file (all atom, bond, angle, torsion and improper changes as well as soft_pairs and off diagonals)
#

# Example of a short MAP file strucure (the system shall contain 3 Q atoms on residues with indexes 395 and 514: their atom names are N2,C10 and H20 (reactive atom)
#   and 1 neighboring atom that is somehow connected to the reactive atom (by angle, torsion or improper) but not included in the Q region -  395.C4 )
#
# In the third and fourth (and fifth...) columns, the LIB IDs (resid.name) are specified for EVB State 1 and State 2 (and State 3...). 
# Comments start with # 
#
# Note: The libs and parms should be properly formatted  sections (mae2q .parm files do not contain these: [atom_types], [bonds], [angles]...)
# Note2: q_ffld2q.py creates these headers
#
"""
# PDB_ID  State 1 type   State 2 type  

q  395.N2     LFN.N2      LF2.N2
q  514.C10    DOP.C10     DO2.C10
q  514.H20    DOP.H20     LF2.H30

n  395.C4     LFN.C4      LF2.C4

"""

from qscripts_config import __version__, QScriptsConfig as QScfg

import sys
import os
import re
import time
import argparse

from Qpyl.common import backup_file


parser = argparse.ArgumentParser()
parser.add_argument('-s', nargs = 1, dest = 'pdb', help = 'pdb structure file', required=True)
parser.add_argument('-p', nargs = '+', dest = 'parms', help = '.parm files', required=True)
parser.add_argument('-l', nargs = '+', dest = 'libs', help = '.lib files', required=True)
parser.add_argument('-m', nargs = 1, dest = 'qmap', help = 'qmap file, see source for info', required=True)
parser.add_argument('-o', nargs = 1, dest = 'outfile', help = 'output file (default=generated.fep)', default=["generated.fep"])
parser.add_argument("--index", dest="pdb_index", action="store_true", help="Write out PDB indexes (8394) instead of placeholders ($593.N1$)", default=False)


if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

for k,v in vars(args).iteritems():
    if k in ["pdb", "parms", "libs", "qmap"]:
        for fn in v:
            if not os.path.lexists(fn):
                print "\n\nFATAL! File '%s' doesn't exist." % fn            
                sys.exit(1)

STATES = 0    # this constant changes :) no. of states is read from the map file
QATOM_MAP_FN = args.qmap[0]
PDB_FN       = args.pdb[0]
PARM_FNS     = args.parms
LIB_FNS      = args.libs
OUT_FN       = args.outfile[0]

pdb_dict = {}

fep_types = { "atoms" : [], "bonds" : [], "angles" : [], "torsions" : [], "impropers" : [] }
fep_changes = {"atoms" : [], "bonds" : [], "angles" : [], "torsions" : [], "impropers" : [] }


class AtomType(object):
    atypes_lookup = {}
    def __init__(self, atom_type, LJ_A, LJ_B, LJ_A_14, LJ_B_14, mass):
        # atom_type is a string like "CA" or "OW"
        # Lennard-Jones A and B
        self.atom_type = atom_type
        self.print_out = False
        self.LJ_A = LJ_A
        self.LJ_B = LJ_B
        self.LJ_A_14 = LJ_A_14
        self.LJ_B_14 = LJ_B_14
        # Soft_core Ci (change manually in the fep) and ai 
        self.SC_Ci = "1"
        self.SC_ai = "2.5"
        if float(mass) == 0.0:
            mass = "<FIX>"
        self.mass = mass

        AtomType.atypes_lookup[atom_type] = self

    def __repr__(self):
        return "AtomType (%s)" % (self.atom_type)

    @classmethod
    def find_by_atomtype(cls, atom_types):
        return AtomType.atypes_lookup[atom_type]


class BondType(object):
    atypes_lookup = {}
    def __init__(self, atom_types, fc, r):
        self.atom_types = atom_types   # list of atom_type strings
        self.fc = fc
        self.r = r
        self.index = None    # set when printing out types; get when printing changes
        self.print_out = False
        self.is_morse = False

        # bond_id = atom_types - sorted, "CA CB" instead of "CB CA", to prevent double entries
        bond_id = " ".join(sorted(atom_types))
        BondType.atypes_lookup[bond_id] = self

    def __repr__(self):
        return "BondType (%s)" % " ".join(self.atom_types)

    @classmethod
    def find_by_atomtypes(cls, atom_types):
        bond_id = " ".join(sorted(atom_types))
        return BondType.atypes_lookup[bond_id]


class AngleType(object):
    atypes_lookup = {}
    def __init__(self, atom_types, fc, theta):
        self.atom_types = atom_types
        self.fc = fc
        self.theta = theta
        self.index = None    # set when printing out types; get when printing changes
        self.print_out = False

        # angle_id = atom types - smaller of forward or reverse lists, "CA OS CT" not "CT OS CA"
        angle_id = " ".join( min(atom_types, list(reversed(atom_types))) )   # take the "smaller" of the normal and reverse lists to prevent double entries
        AngleType.atypes_lookup[angle_id] = self

    def __repr__(self):
        return "AngleType (%s)" % " ".join(self.atom_types)

    @classmethod
    def find_by_atomtypes(cls, atom_types):
        angle_id = " ".join( min(atom_types, list(reversed(atom_types))) )   
        return AngleType.atypes_lookup[angle_id]

class TorsionType(object):
    torsion_lookup = {}
    def __init__(self, atom_types, fc, rmult, psi):   # parameters == [ (fc,-1,psi), (fc,-2,psi), (fc, 3,psi) ]
        self.atom_types = atom_types
        self.fc = fc
        self.rmult = rmult    # -1,-2,3 or similar
        self.psi = psi
        self.index = None    # set when printing out types; get when printing changes
        self.print_out = False

# torsion_id = atom_types - smaller of forward or reverse lists, "CA CB CG OG1" instead of "OG1 CG CB CA"   
        torsion_id = " ".join( min(atom_types, list(reversed(atom_types))) )
# the lookup dictionary in torsions contains a dictionary for each atom_types combo which has kev:value pairs of self.rmult:self
        if not TorsionType.torsion_lookup.has_key(torsion_id):
            TorsionType.torsion_lookup[torsion_id] = {}
# remove decimal zeroes from rmult to prevent possible duplicates
# make it absolute for easier sorting and comparing later when printing the changes
        TorsionType.torsion_lookup[torsion_id][ abs(int(float(self.rmult))) ] = self       

    def __repr__(self):
        return "TorsionType (%s  %s)" % (" ".join(self.atom_types), self.rmult)

    @classmethod
    def find_by_atomtypes(cls, atom_types):
        torsion_id = " ".join( min(atom_types, list(reversed(atom_types))) )
        return TorsionType.torsion_lookup[torsion_id]    # this is a dictionary containing all torsion parameters for these atom_types (-1,-2,3)



class ImproperType(object):
    atypes_lookup = {}
    def __init__(self, atom_types, fc, psi):
        self.atom_types = atom_types  
        self.fc = fc
        self.psi = psi
        self.index = None    # set when printing out types; get when printing changes
        self.print_out = False

        # improper_id = atom types - sorted, with second one fixed, "CB CG OG1 OG2" instead of "OG1 CG CB OG2" or "CB CG OG2 OG1" or ...
        center_atom = atom_types.pop(1)   # sort all but the center atom
        atom_types.sort()
        atom_types.insert(1, center_atom)
        improper_id = " ".join(atom_types)     # take the "smaller" of the normal and reverse lists to prevent double entries
        ImproperType.atypes_lookup[improper_id] = self

    def __repr__(self):
        return "ImproperType (%s)" % " ".join(self.atom_types)

    @classmethod
    def find_by_atomtypes(cls, atom_types):
        center_atom = atom_types.pop(1)   # sort all but the center atom
        atom_types.sort()
        atom_types.insert(1, center_atom)
        improper_id = " ".join(atom_types)     # take the "smaller" of the normal and reverse lists to prevent double entries
        return ImproperType.atypes_lookup[improper_id]


class QAtom(object):
    pdbid_lookup = {}
    def __init__(self, pdb_index, pdb_id, lib_ids, isQ):
        self.pdb_index = pdb_index
        self.pdb_id = pdb_id
        QAtom.pdbid_lookup[pdb_id] = self
        self.isQ = isQ    # it could be a Q atom or one that only shares an angle or torsion with Q atoms
        self.lib_ids = lib_ids
        self.charges = []
        self.connections = [ [] for i in range(0,STATES) ]
        self.atom_types = []
        self.index = 0   
        self.is_reactive = False     # involved in at least one bond breaking or forming

    def __repr__(self):
        return "QAtom: " + self.pdb_id

    @classmethod
    def find_by_pdbid(cls, pdb_id):   # lookup an object by its pdb_id value
        return QAtom.pdbid_lookup[pdb_id]


class QBond(object):
    atoms_lookup = {}
    def __init__(self, atoms):
        self.bond_types = [ None for i in range(0,STATES) ]
        self.print_out = False

        # get pdb_indexes and sort the atoms accordingly 
        at_i = [int(a.pdb_index) for a in atoms]
        at_i = min( at_i, list(reversed(at_i)) )
        self.atoms = map(lambda x: [ a for a in atoms if int(a.pdb_index) == x ][0], at_i  )  # 2 QAtom instances 
        at_i = " ".join([str(ai) for ai in at_i])
        # add self to the lookup dict
        QBond.atoms_lookup[at_i] = self

    def __repr__(self):
        return "QBond (%s)" % ", ".join( [a.pdb_id for a in self.atoms] )

    @classmethod
    def find_by_atoms(cls, atoms):   # lookup an object that contains these atoms
        at_i = [int(a.pdb_index) for a in atoms]
        at_i = min( at_i, list(reversed(at_i)) )
        at_i = " ".join([str(ai) for ai in at_i])
        return QBond.atoms_lookup[at_i]

        
class QAngle(object):
    atoms_lookup = {}
    def __init__(self, atoms):
        self.angle_types = [ None for i in range(0,STATES) ]
        self.print_out = False

        # get pdb_indexes and sort the atoms accordingly 
        at_i = [int(a.pdb_index) for a in atoms]
        at_i = min( at_i, list(reversed(at_i)) )
        self.atoms = map(lambda x: [ a for a in atoms if int(a.pdb_index) == x ][0], at_i  )  # 2 QAtom instances 
        at_i = " ".join([str(ai) for ai in at_i])
        # add self to the lookup dict
        QAngle.atoms_lookup[at_i] = self

    def __repr__(self):
        return "QAngle (%s)" % ", ".join( [a.pdb_id for a in self.atoms] )

    @classmethod
    def find_by_atoms(cls, atoms):   # lookup an object that contains these atoms
        at_i = [int(a.pdb_index) for a in atoms]
        at_i = min( at_i, list(reversed(at_i)) )
        at_i = " ".join([str(ai) for ai in at_i])
        return QAngle.atoms_lookup[at_i]

# a Vtorsion can be defined with more than one function
# example:
# lfh.C11         lfh.C13         lfh.C14         lfh.H25         0.000     -1.000     0.000     1.000
# lfh.C11         lfh.C13         lfh.C14         lfh.H25         0.000     -2.000     180.000     1.000
# lfh.C11         lfh.C13         lfh.C14         lfh.H25         0.000    3.000     0.000     1.000
#
# when Q encounters a positive number in the 6th field (3.000) it knows this is the last parameter

class QTorsion(object):
    atoms_lookup = {}
    def __init__(self, atoms):
        self.atoms = atoms    # Four QAtom instances  
# Three torsiontypes for each Torsion can defined in Q but are not necessary
# Absolute of the "rmult" value (usually, -1,-2,3) is used here as the key   (update: 4th term is also used in newer versions of oplsaa)
        self.torsion_types = {1: [ None for i in range(0,STATES) ],
                              2: [ None for i in range(0,STATES) ],    
                              3: [ None for i in range(0,STATES) ],
                              4: [ None for i in range(0,STATES) ] }
        self.print_out = False

        # get pdb_indexes and sort the atoms accordingly 
        at_i = [int(a.pdb_index) for a in atoms]
        at_i = min( at_i, list(reversed(at_i)) )
        self.atoms = map(lambda x: [ a for a in atoms if int(a.pdb_index) == x ][0], at_i  )  # 4 QAtom instances 
        at_i = " ".join([str(ai) for ai in at_i])
        # add self to the lookup dict
        QTorsion.atoms_lookup[at_i] = self

    def __repr__(self):
        return "QTorsion (%s)" % ", ".join( [a.pdb_id for a in self.atoms] )

    @classmethod
    def find_by_atoms(cls, atoms):   # lookup an object that contains these atoms
        at_i = [int(a.pdb_index) for a in atoms]
        at_i = min( at_i, list(reversed(at_i)) )
        at_i = " ".join([str(ai) for ai in at_i])
        return QTorsion.atoms_lookup[at_i]


class QImproper(object):
    atoms_lookup = {}
    def __init__(self, atoms):
        self.improper_types = [ None for i in range(0,STATES) ]
        self.print_out = False

        self.atoms = atoms  # do not change the order, Q will fuck it up!   # 4 Qatom instances

        # get pdb_indexes and order the atoms accordingly 
        at_i = [int(a.pdb_index) for a in atoms]
        center_atom = at_i.pop(1)
        at_i.sort()
        at_i.insert(1, center_atom)
        at_i = " ".join([str(ai) for ai in at_i])
        # add self to the lookup dict
        QImproper.atoms_lookup[at_i] = self

    def __repr__(self):
        return "QImproper (%s)" % ", ".join( [a.pdb_id for a in self.atoms] )

    @classmethod
    def find_by_atoms(cls, atoms):   # lookup an object that contains these atoms
        at_i = [int(a.pdb_index) for a in atoms]
        center_atom = at_i.pop(1)
        at_i.sort()
        at_i.insert(1, center_atom)
        at_i = " ".join([str(ai) for ai in at_i])
        return QImproper.atoms_lookup[at_i]





# start parsing the libraries, parameters, pdb and qmap files (in that order)

# Save lib entries and parameters as dictionaries
q_library = { "impropers" : [] }

for lib_fn in LIB_FNS:
    with open(lib_fn, 'r') as lib:
        section=""
        resname=""
        for lnumber, line in enumerate(lib.readlines()):
            lnumber += 1
            line = re.split("#|\*|\!",line, 1)[0].strip()   # remove comments 
            if line == "":
                continue
            if line[0] == "{":
                resname = line.strip("{}")   # get the 3 letter code
                continue
            if line[0] == "[":
                section = line.strip("[]").upper()
                continue
            if not resname or not section:
                print "\n\nFATAL! Line #%s in LIB file '%s' is not a comment and is not inside a {XXX} section and [XXX...X] subsection:\n%s" % (lnumber, lib_fn,line)
                if section:
                    print "\nIs %s a lib or parm file? Double check the arguments." % (lib_fn)
                sys.exit(1)
            
            if section == "ATOMS":
                try:
                    atom_index,atom_name,atom_type,atom_charge = line.split()[0:4]
                except ValueError:
                    print "\n\nFATAL! Line %s in LIB file '%s' couldn't be parsed (should look like this 'aindex aname atype charge ...'):\n%s" % (lnumber, lib_fn,line)
                    sys.exit(1)
                if not q_library.has_key(resname):
                    q_library[resname] = {}
                q_library[resname][atom_name] = (atom_type, atom_charge, [])
            elif section == "BONDS":
                a1,a2 = line.split()
                try:
                    q_library[resname][a1][2].append(a2)    # add connectivity
                    q_library[resname][a2][2].append(a1)    # add connectivity 
                except KeyError:
                    print "\n\nFATAL! Undefined atom(s) %s and/or %s mentioned in the bonds section of '%s', ln.%s." % (a1,a2,lib_fn, lnumber)
                    sys.exit(1)
            elif section == "IMPROPERS":
                q_library["impropers"].append( [ resname + "." + name for name in line.split() ] )    # lib_ids


# create AtomType, BondType, AngleType, TorsionType and ImproperType objects and append them in fep_types
prms = ""
for parm_fn in PARM_FNS:
    with open(parm_fn, 'r') as parm:
        section=""
        for lnumber,line in enumerate(parm.readlines()):
            lnumber += 1
            line = re.split("#|\*|\!",line, 1)[0].strip()   # remove comments 
            if line == "":
                continue
            if line[0] == "[":
                section = line.split("]")[0].strip(" [").upper()         # it is apparently allowed to write text after the section identifier, so a simple strip isn't enough
                continue
            if not section:
                print "\n\nFATAL! Line %s in PARM file '%s' is not a comment and is not inside any section ([atom_types], [bonds], ...):\n%s" % (lnumber,parm_fn,line)
                sys.exit(1)
            if section == "ATOM_TYPES":
                parms = line.split()
                try:
                    atom_type = parms[0]
                    lj_a, lj_b, lj_a_14, lj_b_14, mass = parms[1], parms[3], parms[4], parms[5], parms[6]
                except IndexError:
                    print "\n\nFATAL! Could not parse line %s in [atom_types] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                    sys.exit(1)
                try:
                    AtomType.find_by_atomtype(atom_type)  # expecting an exception
                except KeyError:
                    fep_types["atoms"].append( AtomType(atom_type, lj_a, lj_b, lj_a_14, lj_b_14, mass) )
                else:
                    print "\n\nFATAL! Double VDW parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, atom_type)
                    sys.exit(1)
            elif section == "BONDS":
                parms = line.split()
                try:
                    atom_types = parms[0:2]
                    fc, r = parms[2], parms[3]
                except IndexError:
                    print "\n\nFATAL! Could not parse line %s in [bonds] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                    sys.exit(1)
                try:
                    BondType.find_by_atomtypes(atom_types)  # expecting an exception
                except KeyError:
                    fep_types["bonds"].append( BondType(atom_types, fc, r) )   # create new BondType
                else:
                    print "\n\nFATAL! Double BOND parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, " ".join(atom_types))
                    sys.exit(1)
            elif section == "ANGLES":
                parms = line.split()
                try:
                    atom_types = parms[0:3]
                    fc, theta = parms[3], parms[4]
                except IndexError:
                    print "\n\nFATAL! Could not parse line %s in [angles] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                    sys.exit(1)
                try:
                    AngleType.find_by_atomtypes(atom_types)  # expecting an exception
                except KeyError:
                    fep_types["angles"].append( AngleType(atom_types, fc, theta) )  # create new AngleType
                else:
                    print "\n\nFATAL! Double ANGLE parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, " ".join(atom_types))
                    sys.exit(1)
            elif section == "TORSIONS":
                parms = line.split()
                try:
                    atom_types = parms[0:4]
                    fc, rmult, psi = parms[4], parms[5], parms[6]
                except IndexError:
                    print "\n\nFATAL! Could not parse line %s in [torsions] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                    sys.exit(1)

                if not abs(int(float(rmult))) in [1,2,3,4]:
                    print "\n\nFATAL! Found torsion with unsupported mutliplier (not +- 1,2,3 or 4) in parm file '%s', ln.%s:\n%s" % (parm_fn, lnumber, " ".join(atom_types))
                    sys.exit(1)
                try:
                    tortypes = TorsionType.find_by_atomtypes(atom_types)  # this guy returns multiple tortypes in a dictionary or a KeyError if no torsion with this atom_types combo was defined yet
                except KeyError:
                    fep_types["torsions"].append(TorsionType(atom_types, fc, rmult, psi))    # create new TorsionType
                else:
                    if tortypes.has_key( abs(int(float(rmult))) ):
                        print "\n\nFATAL! Double TORSION parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, " ".join(atom_types))
                        sys.exit(1)
                    else:
                        fep_types["torsions"].append(TorsionType(atom_types, fc, rmult, psi))    # create new TorsionType
            elif section == "IMPROPERS":
                parms = line.split()
                try:
                    atom_types = parms[0:4]
                    fc, psi = parms[4], parms[5]
                except IndexError:
                    print "\n\nFATAL! Could not parse line %s in [torsions] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                    sys.exit(1)
                try:
                    ImproperType.find_by_atomtypes(atom_types)  # expecting an exception
                except KeyError:
                    fep_types["impropers"].append( ImproperType(atom_types, fc, psi) )  # create new AngleType
                else:
                    print "\n\nFATAL! Double IMPROPER parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, " ".join(atom_types))
                    sys.exit(1)
            elif section == "OPTIONS":
                pass   # ignore this section
            else:
                print "\n\nFATAL! Unknown section found in the parm file %s: %s" % (parm_fn, section)
                sys.exit(1)


# map atom indexes (values) to residue_id.atom_name tags (keys)
with open(PDB_FN, 'r') as pdb:
    for line in pdb.readlines():
        ri = line[20:26].strip() 
        an = line[11:17].strip()
        rian = ri + "." + an
        index = line[6:11].strip()
        pdb_dict[rian] = index


# read the QMAP file
# Lookup atom indexes from pdb_dict and check if library ids exist in the library
# Create QAtom objects and add them to the atom_list list
lib_ids_uniq = {}
with open(QATOM_MAP_FN, 'r') as qatom_map:
    for line in qatom_map.readlines():
        line = re.split("#|\*|\!",line, 1)[0].strip()   # remove comments
        if line == "":
            continue
        c = line.split()
        if c[0].lower() == "q":
            isQ = True
        elif c[0].lower() == "n":
            isQ = False
        else:
            print "\n\nFATAL! Lines in the MAP file should begin with either 'q' (qatom) or 'n' (not qatom)"
            sys.exit(1)
        pdb_id = c[1]
        try:
            QAtom.find_by_pdbid(pdb_id)
        except:
            pass
        else:
            print "\n\nFATAL! PDB IDs in the MAP file should be unique: %s" % pdb_id
            sys.exit(1)

        lib_ids = c[2:]
        STATES = len(lib_ids)
        for st, lib_id in enumerate(lib_ids):
            if lib_id in lib_ids_uniq.get(st, []):
                print "'\n\nFATAL! The library IDs in one EVB state should be unique (double '%s'), otherwise proper bonding can't be determined." % lib_id
                sys.exit(1)
            else:
                lib_ids_uniq[st] = lib_ids_uniq.get(st, [])
                lib_ids_uniq[st].append(lib_id)
            resname,name = lib_id.split(".")
            if not q_library.has_key(resname):
                print "\n\nFATAL! The library doesn't contain this residue: %s   (be careful of case)" % resname
                sys.exit(1)
            if not q_library[resname].has_key(name):
                print "\n\nFATAL! The library entry %s doesn't contain this atom: %s   (be careful of case)" % (resname, name)
                sys.exit(1)
        try:
            index = pdb_dict[pdb_id]
        except KeyError:
            print "\n\nFATAL! Atom '%s' in the MAP file is not in the PDB structure" % pdb_id
            sys.exit(1)
        a = QAtom(pdb_index=index, pdb_id=pdb_id, lib_ids=lib_ids, isQ=isQ)
        fep_changes["atoms"].append(a)

# Lookup the library - get charges, VDWs and connections for QAtom objects
for atom in fep_changes["atoms"]:
    for evb_state,lib_id in enumerate(atom.lib_ids):
        resname, name = lib_id.split(".")
        atom_type, charge, connections = q_library[resname][name] 
        atom.charges.append(charge)

        try:
            atype = AtomType.find_by_atomtype(atom_type)
        except KeyError:
            print "\n\nFATAL! Parameters for atom %s in resname %s (atom type: %s) were not found  (be careful of case)" % (name, resname, atom_type)
            sys.exit(1)

        atom.atom_types.append(atype)
        if atom.isQ:
            atype.print_out = True
                   
        for atom2_name in connections:        
            lib_id2 = resname + "." + atom2_name
            pdb_resindex = atom.pdb_id.split(".")[0]
            for atom2 in fep_changes["atoms"]:
                if atom2.lib_ids[evb_state] == lib_id2:         # locate the atom object that corresponds to this pdb_id
                    atom.connections[evb_state].append(atom2)   # append it to connections (in this state)
                    break

# Create QBond, Qangle and QTorsion and QImproper objects from QAtoms' connections
for evb_state in range(0,STATES):
    for atom in fep_changes["atoms"]:

        for at2_index,atom2 in enumerate(atom.connections[evb_state]):
            if atom2 == atom:
                continue
            atoms = [atom,atom2]
            try:
                bond = QBond.find_by_atoms(atoms)    # TODO: Make a function .get_or_make() that will get an existing or create a new QBond (same for all)
            except KeyError:
                bond = QBond(atoms)   # create new object
                fep_changes["bonds"].append(bond)

            atom_types = [a.atom_types[evb_state].atom_type for a in atoms]     #  atom_types_in_state_0 = [ a[0].atom_types[0].atom_type, a[1].atom_types[0].atom_type ]
            try:
                bond_type = BondType.find_by_atomtypes(atom_types)
            except KeyError:
                print "\n\nFATAL! Parameters for bond %s were not found  (be careful of case)" % bond
                sys.exit(1)
            else:
                bond.bond_types[evb_state] = bond_type


            for at3_index, atom3 in enumerate(atom.connections[evb_state]):   # iterate through all atom2 connections
                if atom3 in [atom,atom2]:
                    continue
                atoms = [atom2,atom,atom3]
                try:
                    angle = QAngle.find_by_atoms(atoms)
                except KeyError:
                    angle = QAngle(atoms)
                    fep_changes["angles"].append(angle)

                atom_types = [a.atom_types[evb_state].atom_type for a in atoms]     
                try:
                    angle_type = AngleType.find_by_atomtypes(atom_types)
                except KeyError:
                    print "\n\nFATAL! Parameters for angle '%s' were not found  (be careful of case)" % angle
                    sys.exit(1)
                else:
                    angle.angle_types[evb_state] = angle_type


                for at4_index, atom4 in enumerate(atom3.connections[evb_state]):
                    if atom4 in [atom,atom2,atom3]:
                        continue
                    atoms = [atom2,atom,atom3,atom4]
                    try:
                        torsion = QTorsion.find_by_atoms(atoms)
                    except KeyError:
                        torsion = QTorsion(atoms)
                        fep_changes["torsions"].append(torsion)

                    atom_types = [a.atom_types[evb_state].atom_type for a in atoms]     
                    try:
                        torsion_types = TorsionType.find_by_atomtypes(atom_types)    
                    except:
                        print "\n\nFATAL! Parameters for torsion '%s' were not found  (be careful of case)" % torsion
                        sys.exit(1)
                    else:
                        # torsion_types is a dictionary that looks like this: 
                        # { "1": TorsionType, "2":TorsionType, "3":TorsionType, "4":TorsionType}      (keys are abs(rmult))
                        for rmult, ttype in torsion_types.iteritems():
                            torsion.torsion_types[rmult][evb_state] = ttype

# replace lib_id (resname.atom_name) strings that define the atoms in the impropers in the library with QAtom objects (those that exist)
# if all four atoms have been replaced, create a new QImproper object
# 
# The order of atoms in an improper should not matter as long as the second atom is in second place, unfortunately, Q does not understand this.
#
#     For Q, C1 C2 C3 C4 is a different improper than C1 C2 C4 C3, therefore it is CRUCIAL to add it to the FEP file in the same order
#     as it is in the topology to prevent double entries. 
#     The following code segment is executed within a for loop that iterates over states from state=0 to state=N, therefore, 
#     the LIB_IDs of atoms from the first state will be processed first - creating new QImproper objects with the same atom order as in the topology (first state libs).
#     In other states, the improper will either be found if it exists in the first state (even if in diff. order - look at QImproper.find_by_atoms() code),
#     and filled with appropriate parameters, or it will be created (here the order doesn't matter, because it doesn't exist in the topology)
#     
        lib_id = atom.lib_ids[evb_state]
        for imp in q_library["impropers"]:
            if lib_id in imp:
                imp[imp.index(lib_id)] = atom

                atoms=[a for a in imp if isinstance(a,QAtom)]
                if len(atoms) == 4: 
                    try:
                        improper = QImproper.find_by_atoms(atoms)
                    except KeyError:
                        improper = QImproper(atoms)
                        fep_changes["impropers"].append(improper)

                    atom_types = [a.atom_types[evb_state].atom_type for a in atoms]
                    try:
                        improper_type = ImproperType.find_by_atomtypes(atom_types)
                    except KeyError:
                        print "\n\nFATAL! Parameters for improper '%s' were not found  (be careful of case)" % improper
                        sys.exit(1)
                    else:
                        improper.improper_types[evb_state] = improper_type
                    # revert back to lib_ids (the same improper could be reused in another state)
                    for i in range(4):
                        imp[i] = imp[i].lib_ids[evb_state]




# do some sorting
#

def custom_sort_criteria( bati ):
    return [ int(a.pdb_index) for a in bati.atoms ]

fep_changes["bonds"].sort(key=lambda x:custom_sort_criteria(x))
fep_changes["angles"].sort(key=lambda x:custom_sort_criteria(x))
fep_changes["torsions"].sort(key=lambda x:custom_sort_criteria(x))
fep_changes["impropers"].sort(key=lambda x:custom_sort_criteria(x))



# Set print_out=True to those parameters that should be printed and add indexes
#
# Also create a new BondType with is_morse=True in bonds that have None in bond_types (add it to fep_types too)
# ...making a new BondType to preserve the original harmonic in case some other bond uses it as well

def comp_str_floats(s1,s2):
    try:
        return abs( float(s1) - float(s2) ) < 0.00001
    except ValueError:
        return False

bindex = 1
for b in fep_changes["bonds"]:
    if None in b.bond_types:
        b.print_out = True
    else:
        bt1 = b.bond_types[0]
        for bt2 in b.bond_types[1:]:                       # check for differences between states
            if not comp_str_floats(bt1.fc,bt2.fc) or not comp_str_floats(bt1.r, bt2.r):
                b.print_out = True                         # if different, set the bond.print_out to True
    if b.print_out:                                        # if this, print out all its bond_types
        for bt in b.bond_types:
            if bt:
                if None in b.bond_types:     # set the is_morse and is_reactive flags
                    bt_new = BondType(list(bt.atom_types), bt.fc, bt.r)
                    bt_new.is_morse = True
                    bt_new.print_out = True
                    bt_new.index = bindex
                    bindex += 1
                    b.bond_types[ b.bond_types.index(bt) ] = bt_new
                    fep_types["bonds"].append(bt_new)
                    
                    b.atoms[0].is_reactive = True
                    b.atoms[1].is_reactive = True
                elif bt.index == None:                 # if a harmonic that hasn't been assigned an index... 
                    bt.print_out = True
                    bt.index = bindex
                    bindex += 1

aindex = 1
for ang in fep_changes["angles"]:
    if None in ang.angle_types:
        ang.print_out = True
    else:
        angt1 = ang.angle_types[0]
        for angt2 in ang.angle_types[1:]:
            if not comp_str_floats(angt1.fc, angt2.fc) or not comp_str_floats(angt1.theta, angt2.theta):
                ang.print_out = True
    if ang.print_out:
        for angt in ang.angle_types:
            if angt:
                angt.print_out = True
                if angt.index == None:
                    angt.index = aindex
                    aindex += 1

tindex=1
for tor in fep_changes["torsions"]:
# tor.torsion_types example for 2 evb states:   { 1:[TorsionType, TorsionType], 2:[TorsionType, None], 3:[TorsionType,TorsionType] } 
    for rmult, torsion_types in tor.torsion_types.iteritems():  
        if None in torsion_types:
            if not all( [ tort == None for tort in torsion_types ] ):   # if all of them are none, do nothing
                tor.print_out = True
        else:
            tort1 = torsion_types[0]
            for tort2 in torsion_types[1:]:
                if not comp_str_floats(tort1.fc, tort2.fc) or not comp_str_floats(tort1.psi, tort2.psi):
                    tor.print_out = True   
    if tor.print_out:
        torsion_types_all = zip(* zip(* sorted( tor.torsion_types.items() ) )[1] )
        for torsion_types_in_onestate in torsion_types_all:     # torsion_types_in_onestate is a list of TorsionTypes with different rmults, but same evb_state
            for tort in torsion_types_in_onestate:
                if tort:
                    tort.print_out = True
                    if tort.index == None:
                        tort.index = tindex
                        tindex += 1

iindex=1
for impr in fep_changes["impropers"]:
    if None in impr.improper_types:
        impr.print_out = True
    else:
        imprt1 = impr.improper_types[0]
        for imprt2 in impr.improper_types[1:]:
            if not comp_str_floats(imprt1.fc, imprt2.fc) or not comp_str_floats(imprt1.psi, imprt2.psi):
                impr.print_out = True
    if impr.print_out:
        for imprt in impr.improper_types:
            if imprt:
                imprt.print_out = True
                if imprt.index == None:
                    imprt.index = iindex
                    iindex += 1



# Create the fep file - nasty looking section...
#
fepstring = "[FEP]\n"
fepstring += "states %d\n\n" % STATES

fepstring += "[atoms]\n"
fepstring += "%-15s %-10s     #   %-10s  %-10s\n" % ("#Q index","PDB index", "PDB ID", "STATE 1 TYPE")
i = 0
for a in fep_changes["atoms"]:
    if a.isQ:
        i += 1
        a.index = i
        if args.pdb_index:
            fepstring += "%-15s %-10s     #   %-10s  %-10s %3s\n" % (a.index, a.pdb_index, a.pdb_id, a.lib_ids[0], " !" * a.is_reactive)
        else:
            fepstring += "%-15s %-10s     #   %-10s  %-10s %3s\n" % (a.index, ("$"+a.pdb_id+"$"), a.pdb_id, a.lib_ids[0], " !" * a.is_reactive)
    if a.is_reactive:
        for evb_state in range(0,STATES):
            a.atom_types[evb_state].SC_Ci = "<FIX>"

fepstring += "\n\n[atom_types]\n"
for at in fep_types["atoms"]:
    if at.print_out == False:
        continue
    fepstring += "%-10s %10s %10s %10s %10s %10s %10s %10s \n" % (at.atom_type, 
                                                                     at.LJ_A, 
                                                                     at.LJ_B, 
                                                                     at.SC_Ci,
                                                                     at.SC_ai, 
                                                                     at.LJ_A_14, 
                                                                     at.LJ_B_14, 
                                                                     at.mass)
                                                                     
fepstring += "\n\n[change_atoms]\n"
for a in fep_changes["atoms"]:
    if a.isQ == False:
        continue
    s = "%-10s" % a.index
    for st in range(0,STATES):
        s = s +  " %-12s" % a.atom_types[st].atom_type
    c = "     #    %-10s" % a.lib_ids[0]
    if a.is_reactive:
        c = c + " !"
    fepstring += s + c + "\n"


fepstring += "\n\n[change_charges]\n"
for a in fep_changes["atoms"]:
    if a.isQ == False:
        continue
    s = "%-5s" % a.index
    for st in range(0,STATES):
        s = s + " %10s" % a.charges[st]
    c = "     #    %-10s" % a.lib_ids[0]
    if STATES == 2:  # print difference in charges
        dq = float(a.charges[1]) - float(a.charges[0])
        c = c + " dq=%7.4f " % dq
    if a.is_reactive:
        c = c + " !"
    fepstring += s + c + "\n"

# Since the for loop is the same, get the soft_pairs and off_diagonals at the same time

soft_pairs = "\n\n[soft_pairs]\n"
off_diags = "\n\n[off_diagonals]\n"
for bond in fep_changes["bonds"]:
    if None in bond.bond_types:
        n = 0
        for state, btype in enumerate(bond.bond_types):
            if btype != None:
                for state2, btype in enumerate(bond.bond_types):
                    if btype == None:
                        st1,st2 = sorted( (state+1, state2+1) )
                        n += 1
                        if n > 1: 
                            off_diags += "#"
                        off_diags += "%-3s %-3s %5s %5s  0  0  # %-20s\n" % (st1, st2, bond.atoms[0].index, bond.atoms[1].index, "-".join( [ x.atom_types[state].atom_type for x in bond.atoms ] ) )
                soft_pairs += "%-5s %-5s   # %-20s\n" % (bond.atoms[0].index, bond.atoms[1].index, "-".join( [ x.atom_types[state].atom_type for x in bond.atoms ] ) ) 
fepstring += soft_pairs + "\n"
fepstring += off_diags + "\n"


fepstring += "\n\n[bond_types]\n"
for bt in sorted(fep_types["bonds"], key=lambda x: x.index):
    if not bt.print_out:
        continue
    s = "%-10s %10s %10s" % (bt.index, bt.fc, bt.r)
    com = "   # %-20s" % "-".join(bt.atom_types)
    if bt.is_morse:
        s = "%-10s %10s %10s %10s" % (bt.index, "<FIX_D>","<FIX_a>", "<FIX_r>" )
        com += " !"
    fepstring += s + com + "\n"


fepstring += "\n\n[change_bonds]\n"
for b in fep_changes["bonds"]:
    if not b.print_out:
        continue
    if args.pdb_index:
        s = "".join( [ "%-8s" % a.pdb_index for a in b.atoms ])
    else:
        s = "".join( [ "%-10s" % ("$"+a.pdb_id+"$") for a in b.atoms ])
    com = "   #    %-20s " % "-".join( [a.pdb_id for a in b.atoms] )
    for bt in b.bond_types:
        if bt:
            i = bt.index
            c = "-".join(bt.atom_types)
        else:
            i = 0
            c = "None"
        s = s + "%5d" % i
        com = com + " %-20s" % c
    fepstring += s + com + "\n"
    

fepstring += "\n\n[angle_types]\n"
for angt in sorted(fep_types["angles"], key=lambda x: x.index):
    if not angt.print_out:
        continue
    s = "%-5s %10s %10s" % (angt.index, angt.fc, angt.theta)
    com = "   # %-20s" % "-".join(angt.atom_types)
    fepstring += s + com + "\n"

fepstring += "\n\n[change_angles]\n"
for ang in fep_changes["angles"]:
    if not ang.print_out:
        continue
    if args.pdb_index:
        s = "".join( [ "%-8s" % a.pdb_index for a in ang.atoms ])
    else:
        s = "".join( [ "%-10s" % ("$"+a.pdb_id+"$") for a in ang.atoms ])
    com = "   #   %-30s" % "-".join( [a.pdb_id for a in ang.atoms] )
    for angt in ang.angle_types:
        if angt:
            i = angt.index  
            c = "-".join(angt.atom_types)
        else:
            i = 0
            c = "None"
        s = s + "%5d" % i
        com = com + " %-30s" % c
    fepstring += s + com + "\n"
    
fepstring += "\n\n[torsion_types]\n"
for tort in sorted(fep_types["torsions"], key=lambda x: x.index):
    if not tort.print_out:
        continue
    s = "%-5s %10s %10s %10s" % (tort.index, tort.fc, tort.rmult, tort.psi)
    com = "   # %-20s" % "-".join(tort.atom_types)
    fepstring += s + com + "\n"

fepstring += "\n\n[change_torsions]\n"
for tor in fep_changes["torsions"]:
    if not tor.print_out:
        continue
    for rmult,torsion_types in sorted(tor.torsion_types.items()):
        if all( [ tt == None for tt in torsion_types ] ):   # if it is None in all states, don't fepstring += it out
            continue
        if args.pdb_index:
            s = "".join( [ "%-8s" % a.pdb_index for a in tor.atoms ])
        else:
            s = "".join( [ "%-10s" % ("$"+a.pdb_id+"$") for a in tor.atoms ])
        com = "        # %-40s    " % "-".join( [a.pdb_id for a in tor.atoms] ) 
        for tort in torsion_types:   # in each state
            if tort:
                i = tort.index  
                c = "-".join(tort.atom_types)
            else:
                i = 0
                c = "None"
            s = s + "%5d" % i
            com = com + " %-35s" % c
        fepstring += s + com + "\n"
    
fepstring += "\n\n[improper_types]\n"
for imprt in sorted(fep_types["impropers"], key=lambda x: x.index):
    if not imprt.print_out:
        continue
    s = "%-5s %10s %10s" % (imprt.index, imprt.fc, imprt.psi)
    com = "   # %-20s" % "-".join(imprt.atom_types)
    fepstring += s + com + "\n"

fepstring += "\n\n[change_impropers]\n"
for impr in fep_changes["impropers"]:
    if not impr.print_out:
        continue
    if args.pdb_index:
        s = "".join( [ "%-8s" % a.pdb_index for a in impr.atoms ])
    else:
        s = "".join( [ "%-10s" % ("$"+a.pdb_id+"$") for a in impr.atoms ])
    com = "    #  %-40s"  % ("-".join( [a.pdb_id for a in impr.atoms]))
    for imprt in impr.improper_types:
        if imprt:
            i = imprt.index  
            c = "-".join(imprt.atom_types)
        else:
            i = 0
            c = "None"
        s = s + "%5d" % i
        com = com + " %-35s" % c
    fepstring += s + com + "\n"

fepstring += "\n\n\n#\n# Created with 'q_makefep.py' version '%s'\n# Date: %s\n#" % (__version__, time.ctime())

backup = backup_file(OUT_FN)
if backup:
    print "Backed up '%s' to '%s'" % (OUT_FN, backup)

open(OUT_FN, 'w').write(fepstring)
print "\nWrote fep file:  %s\n" % (OUT_FN)
print "Please inspect it for probable errors and add the soft pair and Morse parameters..."

# Debug - print out non-changing parameters
#
#for a in fep_changes["atoms"]:
#    print a, a.connections
#
#for k,v in fep_changes.iteritems():
#    try:
#        print k
#        for asd in v:
#            if not asd.print_out:
#                print asd
#    except:
#        pass
