#!/usr/bin/env python
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
This module implements an internal topology builder class QTopology.
QTopology creates a mapping between the system's structure (QStruct),
bonding patterns/charges (QLib), and the parameters (Qprm), allowing
evaluation of individual topological components of the system.
"""

from __future__ import absolute_import, unicode_literals, division
from Qpyl.core import qlibrary
from Qpyl.core import qparameter
from Qpyl.core import qstructure
from Qpyl.core import qpotential
from six.moves import range

class QTopologyError(Exception):
    pass

class QTopology(object):
    """
    Class for storing topology information.

    Args:
        qlib (qlibrary.QLib): library object
        qprm (qparameter.QPrm): parameter object
        qstruct (qstructure.QStruct): structure object

    This object is basically a mashup of library,
    parameter and structure data. It contains the lists
    of atoms, bonds, angles, torsions and impropers in
    the structure, along with their parameters.

    Typical usage:

    qlib = qlibrary.QLib("amber")
    qprm = qparameter.QPrm("amber")
    qstruct = qstructure.QStruct()
    qlib.read_lib(".../qamber14.lib")
    qprm.read_prm(".../qamber14.prm")
    qstruct.read_pdb(".../14u3.pdb")

    try:
        qtopo = QTopology(qlib, qprm, qstruct)
    except QTopologyError as e:
        print "Failed to make topology: " + str(e)

    for bond in qtopo.bonds:
        e, r = bond.calc()
        print(bond, bond.prm.fc, bond.prm.r0, e, r)

    """

    def __init__(self, qlib, qprm, qstruct):
        # do some type checking to prevent bad things from happening
        for arg, _type in ((qlib, qlibrary.QLib),
                           (qprm, qparameter.QPrm),
                           (qstruct, qstructure.QStruct)):
            if not isinstance(arg, _type):
                raise QTopologyError("{} not of type {}".format(arg, _type))

        if qlib.ff_type != qprm.ff_type:
            raise QTopologyError("QLib FF ({}) not "
                                 "compatible with QPrm FF ({})"
                                 "".format(qlib.ff_type, qprm.ff_type))
        self.qlib = qlib
        self.qprm = qprm
        self.qstruct = qstruct

        try:
            self.qlib.check_valid()  # check if lib entries are good
        except qlibrary.QLibError as e:
            raise QTopologyError(e)


        self.residues = []
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.torsions = []
        self.impropers = []

        for residue_struct in self.qstruct.residues:

            # see if it is defined in the library
            try:
                residue_lib = self.qlib.residue_dict[residue_struct.name]
            except KeyError:
                raise QTopologyError("Residue '{}' not found in library"
                                     .format(residue_struct.name))
            # create new object for the residue
            res_index = len(self.residues) + 1
            residue_top = _TopoResidue(res_index, residue_struct, residue_lib)
            self.residues.append(residue_top)

            # get atoms and bonds within the residue
            self._get_atoms(residue_top)
            self._get_bonds(residue_top)

	# get angles, torsions and impropers from the bonds
        self._get_angles_torsions()
        # get impropers (from the lib entries)
        self._get_impropers()


    def find_NB_pairs(self, atoms, atoms_surr=None, cutoff=None):
        """
        Returns non-bonded atom pairs between groups `atoms` and `atoms2`.

        Args:
            atoms (list): list of _TopoAtom objects or indexes (in top.atoms)
            atoms_surr (list): list of _TopoAtom objects or indexes
                               (in top.atoms), default=None
            cutoff (float): distance cut-off in Angstrom


        Returns:
            (nb_pairs): list of _TopoNonBondedPair objects


        The search for non-bonded atom pairs is performed against all atoms
        defined in `atoms`, as well as against all atom-atom combinations
        between the two groups `atoms` and `atoms_surr`. Pairs exclusively
        from `atoms2` are not included in the search.

        1-2 and 1-3 bonded pairs are excluded, while 1-4 bonded pairs are
        flagged accordingly.

        Note:
            The search for nonbonded atom pairs is done here and not during
            initialization because it's usually very expensive, and because 
            in most cases one is interested in a limited set of interacting
            atoms with a specific cutoff.

        Note2:
            The NB pairs are not stored in the topology object, nor are the
            references added to the atom objects (unlike Bonding objects).

        """


        raise NotImplementedError



    def _get_atoms(self, residue_top):
        """
        Creates _TopoAtom objects and adds them to
        _TopoResidue and QTopology.atoms
        """

        # iterate over the atoms in the library
        for atom in residue_top.lib.atoms:
            aname, atype, charge = atom.name, atom.atom_type, atom.charge
            # check if atoms exist in the structure
            try:
                atom_struct = [a for a in residue_top.struct.atoms
                               if a.name == aname][0]
            except IndexError:
                raise QTopologyError("Atom '{}' in residue '{}.{}'"
                                     " missing in the structure"
                                     .format(aname,
                                             residue_top.struct.index_struct,
                                             residue_top.struct.name))
            # check if atom parameters exist
            try:
                atom_prm = self.qprm.atom_types[atype]
            except KeyError:
                raise QTopologyError("Atom type '{}' not found!"
                                     .format(atype))

            # add new atom to list (and to the residue_top)
            atom_index = len(self.atoms) + 1
            a = _TopoAtom(atom_index, aname, charge, atom_prm,
                          atom_struct, residue_top)
            self.atoms.append(a)
            residue_top.add_atom(a)



    def _get_bonds(self, residue_top):
        """
        Creates _TopoBond objects and adds them to QTopology.bonds
        """

        # iterate over the bonds in the library
        for bond in residue_top.lib.bonds:
            # find the atom objects with those names
            atoms = [a for a in residue_top.atoms if a.name in bond]
            # find parameters
            atom_types = [a.prm.atom_type for a in atoms]
            prm_id = qparameter._PrmBond.get_id(atom_types)
            try:
                bond_prm = self.qprm.bonds[prm_id]
            except KeyError:
                raise QTopologyError("Bond type '{}' not found!"
                                     .format(prm_id))
            # create _TopoBond object
            self.bonds.append(_TopoBond(atoms, bond_prm))

        try:
            # -2 is assuming that the current one was just added (-1)
            prev_res = self.residues[-2]
        except IndexError:
            # first residue
            pass
            # don't check separate chains
        else:
            if residue_top.struct.molecule == prev_res.struct.molecule:
                for conn in residue_top.lib.connections:
                    for conn_prev in prev_res.lib.connections:
                        if "head" in conn and "tail" in conn_prev:

                            ahead = [a for a in residue_top.atoms if
                                     a.name == conn.split()[1]][0]

                            atail = [a for a in prev_res.atoms if
                                     a.name == conn_prev.split()[1]][0]

                            atoms = [atail, ahead]
                            atom_types = [a.prm.atom_type for a in atoms]
                            prm_id = qparameter._PrmBond.get_id(atom_types)
                            try:
                                bond_prm = self.qprm.bonds[prm_id]
                            except KeyError:
                                raise QTopologyError("Bond type '{}'"
                                                     "not found!"
                                                     .format(prm_id))
                            # create _TopoBond object
                            self.bonds.append(_TopoBond(atoms, bond_prm))




    def _get_angles_torsions(self):
        """
        Creates _TopoAngle and _TopoTorsion objects and
        adds them to QTopology.angles and QTopology.torsions
        """

        # to prevent backtracking
        processed_bonds = set()
        # iterate over all bonds and find the angles
        for bond1 in self.bonds:
            processed_angle_bonds = set()
            for bond2 in bond1.atoms[0].bonds + bond1.atoms[1].bonds:
                if bond2 == bond1 or bond2 in processed_bonds:
                    continue
                atoms1 = set(bond1.atoms)
                atoms2 = set(bond2.atoms)
                common_atom = atoms1 & atoms2
                side_atoms = atoms1 ^ atoms2
                angle_atoms = [side_atoms.pop(),
                               common_atom.pop(),
                               side_atoms.pop()]
                # find the angle parameter
                angle_atypes = [a.prm.atom_type for a in angle_atoms]
                prm_id = qparameter._PrmAngle.get_id(angle_atypes)
                try:
                    angle_prm = self.qprm.angles[prm_id]
                except KeyError:
                    raise QTopologyError("Angle type '{}' not found!"
                                         .format(prm_id))
                # create _TopoAngle object
                self.angles.append(_TopoAngle(angle_atoms, angle_prm))

                # find the torsions by looking at the bonds
                # of the angle's side atoms
                for side_atom_index in [0, 2]:
                    for bond3 in angle_atoms[side_atom_index].bonds:
                        if bond3 in processed_bonds or \
                           bond3 in processed_angle_bonds:
                            continue
                        try:
                            atom4 = [a for a in bond3.atoms
                                     if a not in angle_atoms][0]
                        except IndexError:
                            # both atoms are part of the angle
                            continue

                        if side_atom_index == 0:
                            torsion_atoms = [atom4] + angle_atoms
                        else:
                            torsion_atoms = angle_atoms + [atom4]

                        # TODO: QPrm.find_type() would be better
                        #
                        # find parameters
                        atom_types = [a.prm.atom_type for a in torsion_atoms]
                        prm_id = qparameter._PrmTorsion.get_id(atom_types)
                        try:
                            torsion_prm = self.qprm.torsions[prm_id]
                        except KeyError:
                            # see if generic parameters exist
                            gen_atypes = ["?"] + prm_id.split()[1:3] + ["?"]
                            prm_id_gen = qparameter._PrmTorsion.get_id(gen_atypes)
                            try:
                                torsion_prm = \
                                   self.qprm.generic_torsions[prm_id_gen]
                            except KeyError:
                                raise QTopologyError("Torsions type '{}' "
                                                     "for torsion '{}'"
                                                     "not found!"
                                                     .format(prm_id,
                                                     " ".join([a.name for a in
                                                      torsion_atoms])))
                        # create _TopoTorsion object
                        self.torsions.append(_TopoTorsion(torsion_atoms,
                                                          torsion_prm))
                # remove the 'angle' bond from the torsion search
                # (otherwise you get forward and reverse duplicates)
                processed_angle_bonds.add(bond2)

            # remove the bond from the search (to prevent backtracking)
            processed_bonds.add(bond1)



    def _get_impropers(self):
        # create impropers -
        # only those that are explicitly defined in the library
        for residue_index, residue in enumerate(self.residues):
            for improper in residue.lib.impropers:
                # find _TopoAtom-s involved
                atoms = []
                for aname in improper:
                    res = residue
                    # some impropers span to next or
                    # previous residues (-C, +N)
                    if "+" in aname:
                        if residue_index+1 == len(self.residues):
                            continue
                        res = self.residues[residue_index+1]
                    if "-" in aname:
                        if residue_index == 0:
                            continue
                        res = self.residues[residue_index-1]

                    # if separate chains, skip
                    if residue.struct.molecule != res.struct.molecule:
                        continue

                    aname = aname.strip("+-")

                    try:
                        atoms.append([a for a in res.atoms
                                      if a.name == aname][0])
                    except IndexError:
                        if not res.lib.connections:
                            # no connectivity between residues
                            # (end of protein - ligands or water)
                            continue
                        else:
                            raise QTopologyError("Bad improper '{}' between "
                                                 "residues '{}' and '{}'"
                                                 .format(" ".join(improper),
                                                         residue.index,
                                                         res.index))

                if len(atoms) != 4:
                    continue

                # find parameters
                other_atypes = [a.prm.atom_type for a in atoms]
                center_atype = other_atypes.pop(1)

                prm_id = qparameter._PrmImproper.get_id(center_atype, other_atypes)
                try:
                    improper_prm = self.qprm.impropers[prm_id]
                except KeyError:
                    improper_prm = None
                # see if general parameters exist - a lot of options here, example:
                # CA O2 CB CN
                #
                # Single wildcard:
                # ?  O2 CA CN
                # ?  O2 CA CB
                # ?  O2 CB CN
                if improper_prm is None:
                    for i in range(3):
                        ots = [other_atypes[i], other_atypes[(i+1)%3], "?"]
                        prm_id_gen = qparameter._PrmImproper.get_id(center_atype, ots)
                        try:
                            improper_prm = \
                                self.qprm.generic_impropers[prm_id_gen]
                            break
                        except KeyError:
                            improper_prm = None
                # Two wildcards:
                # ?  O2 CB ?
                # ?  O2 CA ?
                # ?  O2 CN ?
                if improper_prm is None:
                    for i in range(3):
                        otypes = [other_atypes[i], "?", "?"]
                        prm_id_gen = qparameter._PrmImproper.get_id(center_atype, otypes)
                        try:
                            improper_prm = \
                                self.qprm.generic_impropers[prm_id_gen]
                            break
                        except KeyError:
                            improper_prm = None

                if improper_prm is None:
                    raise QTopologyError("Improper type '{}' "
                                         "not found!"
                                         .format(prm_id))

                # create _TopoImproper object (same order as library)
                self.impropers.append(_TopoImproper(atoms, improper_prm))


################################################################################


class _TopoAtom(object):
    """
    Class containing topological information for an atom.

    Arguments:
      index (int):  topology index of atom (1-based)
      name (string):  atom name as defined in the library/structure
      charge (float):  charge as defined in the QLib library
      prm (_PrmAtom):  atom parameter as defined in QPrm
      struct (_StructAtom):  atom structure object (stuff from PDB)
      residue (_TopoResidue):  reference to its parent residue

    All these arguments are stored as object properties with the same name.
    Additionaly, bonds, angles, torsions and impropers are lists that contain
    references to _TopoBond, _TopoAngles etc. objects. These are filled in
    automatically when creating aforementioned objects.

    """


    def __init__(self, index, name, charge, prm, struct, residue):
        self.index = index
        self.name = name
        self.charge = charge
        self.prm = prm
        self.struct = struct
        self.residue = residue
        self.bonds = []
        self.angles = []
        self.torsions = []
        self.impropers = []
        self.bati_map = {_TopoBond: self.bonds,
                         _TopoAngle: self.angles,
                         _TopoTorsion: self.torsions,
                         _TopoImproper: self.impropers}

    def __repr__(self):
        return "_TopoAtom: {}_{}.{}".format(self.residue.name,
                                            self.residue.index,
                                            self.name)



    def add_ref(self, bond_angle_tor_imp):
        """
        Add bond, angle, torsion and improper references.

        See class docstring.
        """
        _type = type(bond_angle_tor_imp)
        try:
            self.bati_map[_type].append(bond_angle_tor_imp)
        except KeyError:
            raise TypeError("bond_angle_tor_imp of unsupported "
                            "type: {}".format(_type))


################################################################################


class _TopoResidue(object):
    """
    Class containing topological information for a residue.

    Arguments:
       index (int):  topology index of residue (1-based)
       struct (_StructResidue):  object with structure information
       lib (_LibResidue):  object with library information

    Besides the two properties above, it contains a list of its atoms
    as _TopoAtom objects in the 'atoms' property. This list is filled
    automatically as atoms are created with this residue passed in as an
    argument.

    """


    def __init__(self, index, struct, lib):
        self.index = index
        self.struct = struct
        self.lib = lib
        self.name = struct.name
        self.atoms = []

    def add_atom(self, atom):
        """Append a _TopoAtom object to the 'atoms' list"""
        self.atoms.append(atom)


################################################################################


class _TopoBond(object):
    """Contains topological information for a bond.

    Args:
        atoms (list of _TopoAtom): list of _TopoAtom objects
        prm (_PrmBond): bond parameter object

    """
    def __init__(self, atoms, prm):
        self.prm = prm
        self.atoms = sorted(atoms, key=lambda atom: atom.index)
        for atom in self.atoms:
            atom.add_ref(self)
    
    def calc(self, r=None):
        """Calculate bond distance and energy.
        
        Args:
            r (float, optional): define the bond distance instead of
                                 calculating it from the structure

        Returns tuple (E [kcal/mol], r [angstrom])
        """
        if not r:
            ac1, ac2 = [a.struct.coordinates for a in self.atoms]
            r = qpotential.bond_distance(ac1, ac2)

        e = qpotential.bond_energy(r, self.prm.fc, self.prm.r0)
        return (e,r)

    def __repr__(self):
        atoms_str = "-".join([a.name for a in self.atoms])
        return "{}: ({})".format(self.__class__.__name__, atoms_str)



################################################################################



class _TopoAngle(object):
    """Contains topological information for an angle.

    Args:
        atoms (list of _TopoAtom): list of _TopoAtom objects
        prm (_PrmAngle): angle parameter object

    """
    def __init__(self, atoms, prm):
        self.prm = prm
        atom_indexes = [(a.index, a) for a in atoms]
        atom_indexes = min(atom_indexes, list(reversed(atom_indexes)))
        self.atoms = [a for (i,a) in atom_indexes]
        for atom in self.atoms:
            atom.add_ref(self)

    def calc(self, theta=None):
        """Calculate angle and energy

        Args:
            theta (float, optional): define the angle instead of calculating it
                                     from the structure

        Returns tuple (E [kcal/mol], theta [degrees])
        """

        if theta is None:
            ac1, ac2, ac3 = [a.struct.coordinates for a in self.atoms]
            theta = qpotential.angle_angle(ac1, ac2, ac3)

        e = qpotential.angle_energy(theta,
                                    self.prm.fc,
                                    self.prm.theta0)
        return (e, theta)

    def __repr__(self):
        atoms_str = "-".join([a.name for a in self.atoms])
        return "{}: ({})".format(self.__class__.__name__, atoms_str)



################################################################################



class _TopoTorsion(object):
    """Contains topological information for a torsion.
    
    Args:
        atoms (list of _TopoAtom): list of _TopoAtom objects
        prm (_PrmTorsion): torsion parameter object

    """
    def __init__(self, atoms, prm):
        self.prm = prm
        atom_indexes = [(a.index, a) for a in atoms]
        atom_indexes = min(atom_indexes, list(reversed(atom_indexes)))
        self.atoms = [a for (i,a) in atom_indexes]
        for atom in self.atoms:
            atom.add_ref(self)
    
    def calc(self, phi=None):
        """Calculate torsion angle and energy

        Args:
            phi (float, optional): define the angle instead of calculating it
                                   from the structure

        Returns tuple (E [kcal/mol], phi [degrees])
        """

        if phi is None:
            ac1, ac2, ac3, ac4 = [a.struct.coordinates for a in self.atoms]
            phi = qpotential.torsion_angle(ac1, ac2, ac3, ac4)

        energy = 0
        for fc, multiplicity, phase, npaths in self.prm.get_prms():
            energy += qpotential.torsion_energy(phi,
                                                fc, multiplicity, 
                                                npaths, phase)
        return (energy, phi)

    @property
    def prm_full(self):
        """Return full parameter in case it is generic.

        Basically, make a copy of the generic parameter,
        but use actual atom-types instead of X's.
        """
        if self.prm.is_generic:
            atypes = [a.prm.prm_id for a in self.atoms]
            comment = "Generic: {}".format(self.prm.prm_id)
            full_prm = type(self.prm)(atypes, comment=comment)
            for p in self.prm.get_prms():
                full_prm.add_prm(*p)
            return full_prm
        else:
            return self.prm

    def __repr__(self):
        atoms_str = "-".join([a.name for a in self.atoms])
        return "{}: ({})".format(self.__class__.__name__, atoms_str)



################################################################################



class _TopoImproper(object):
    """Contains topological information for an improper.
    
    Args:
        atoms (list of _TopoAtom): list of _TopoAtom objects
        prm (_PrmImproper): improper parameter object

    """
    def __init__(self, atoms, prm):
        self.prm = prm
        # order is defined in the library
        self.atoms = atoms
        for atom in self.atoms:
            atom.add_ref(self)

    def calc(self, phi=None):
        """Calculate improper angle and energy

        Args:
            phi (float, optional): define the angle instead of calculating it
                                   from the structure

        Returns tuple (E [kcal/mol], phi [degrees])
        """

        if phi is None:
            ac1, ac2, ac3, ac4 = [a.struct.coordinates for a in self.atoms]
            phi = qpotential.improper_angle(ac1, ac2, ac3, ac4)
        e =  qpotential.improper_energy_periodic(phi,
                                                 self.prm.fc,
                                                 self.prm.multiplicity,
                                                 self.prm.phi0)
        return (e, phi)

    @property
    def prm_full(self):
        """Return full parameter in case it is generic.

        Basically, make a copy of the generic parameter,
        but use actual atom-types instead of X's.
        """
        if self.prm.is_generic:
            atypes = [a.prm.prm_id for a in self.atoms]
            center = atypes.pop(1)
            comment = "Generic: {}".format(self.prm.prm_id)
            full_prm = type(self.prm)(center, atypes, self.prm.fc,
                                      self.prm.phi0, self.prm.multiplicity,
                                      comment=comment)
            return full_prm
        else:
            return self.prm

    def __repr__(self):
        atoms_str = "-".join([a.name for a in self.atoms])
        return "{}: ({})".format(self.__class__.__name__, atoms_str)



################################################################################



class _TopoNonBondedPair(object):
    """Contains topological information for a non-bonded atom pair.
    
    Args:
        atoms (list of _TopoAtom): list of _TopoAtom objects
        is_14 (boolean): flag to indicate 1-4 bonded atoms

    """
    def __init__(self, atoms, is_14):
        self.is_14 = is_14
        self.atoms = sorted(atoms, key=lambda atom: atom.index)
    
    def calc_LJ(self, r=None):
        """Calculate distance and Lennard-Jones potential.
        
        Args:
            r (float, optional): define the distance instead of
                                 calculating it from the structure

        Returns tuple (E [kcal/mol], r [angstrom])
        """
        ac1, ac2 = [a.struct.coordinates for a in self.atoms]
        if not r:
            r = qpotential.bond_distance(ac1, ac2)

        if ac1.prm.lj_A is not None:
            lj_A = ac1.prm.lj_A * ac2.prm.lj_A
            lj_B = ac1.prm.lj_B * ac2.prm.lj_B
            e = qpotential.vdw_LJ_AB(r, lj_A, lj_B)
        else:
            lj_eps = (ac1.prm.lj_eps * ac2.prm.lj_eps)**0.5
            lj_R = ac1.prm.lj_R + ac2.prm.lj_R
            e = qpotential.vdw_LJ_epsR(r, lj_eps, lj_R)

        return (e,r)

    def calc_coulomb(self, r=None):
        """Calculate distance and Coulomb energy.
        
        Args:
            r (float, optional): define the distance instead of
                                 calculating it from the structure

        Returns tuple (E [kcal/mol], r [angstrom])
        """
        ac1, ac2 = [a.struct.coordinates for a in self.atoms]
        if not r:
            r = qpotential.bond_distance(ac1, ac2)
        e = qpotential.coulomb(r, ac1.charge, ac2.charge)
        return (e,r)

    def __repr__(self):
        atoms_str = "-".join([a.name for a in self.atoms])
        return "{}: ({})".format(self.__class__.__name__, atoms_str)

