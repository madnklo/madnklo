#!/usr/bin/env python
##########################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
##########################################################################################
"""Definition of all the classes and features relevant to the handling of 
higher order IR subtraction.
"""

import inspect
import itertools
import logging
import math
import sys
import os
import importlib
import shutil
import copy
import traceback

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(
        os.path.realpath(__file__)),
        os.path.pardir,
        os.path.pardir
    ))

from madgraph import MadGraph5Error, MG5DIR, InvalidCmd
import madgraph.core.base_objects as base_objects
import madgraph.fks.fks_common as fks_common
import madgraph.various.misc as misc
from madgraph.iolibs.files import cp, ln, mv
from bidict import bidict

logger = logging.getLogger('madgraph')
pjoin = os.path.join

#=========================================================================================
# Multinomial function
#=========================================================================================

def multinomial(i_s):
    """Compute the multinomial (i_1 + ... + i_n)!/(i_1! ... i_n!)."""

    itot = 0
    denom = 1
    for i in i_s:
        itot += i
        denom *= math.factorial(i)
    return math.factorial(itot) / denom

#=========================================================================================
# Ordered tuples
#=========================================================================================

def is_subset(A, B):
    """Check if the ordered tuple A is contained in the ordered tuple B,
    taking into account multiplicities.
    """

    A_iter = iter(A)
    B_iter = iter(B)
    try:
        A_last = next(A_iter)
    except StopIteration:
        return True
    try:
        B_last = next(B_iter)
    except StopIteration:
        return False
    while True:
        while (A_last > B_last):
            try:
                B_last = next(B_iter)
            except StopIteration:
                return False
        if A_last != B_last:
            return False
        try:
            A_last = next(A_iter)
        except StopIteration:
            return True
        try:
            B_last = next(B_iter)
        except StopIteration:
            return False

def is_disjoint(A, B):
    """Check if the ordered tuple A is disjoint from the ordered tuple B."""

    try:
        for Ai in A:
            B, B_iter = itertools.tee(B, 2)
            for Bi in B_iter:
                if Ai > Bi:
                    B.next()
                elif Bi == Ai:
                    return False
                else:
                    break
        return True
    except StopIteration: # B exhausted, no chance of intersection
        return True

class SortedUnion(object):
    """Iterator returning the sorted union of sorted iterators."""

    __end = object()

    def __init__(self, *A):

        self.A_iter = iter(())
        self.A_last = self.__end
        self.B_iter = iter(())
        self.B_last = self.__end
        if A:
            self.A_iter = iter(A[0])
            self.A_last = next(self.A_iter, self.__end)
            if len(A) > 1:
                if (len(A) == 2):
                    self.B_iter = iter(A[1])
                else:
                    self.B_iter = iter(SortedUnion(*A[1:]))
                self.B_last = next(self.B_iter, self.__end)

    def __iter__(self):

        return self

    def next(self):

        tmp = self.__end
        if self.A_last is not self.__end and self.B_last is not self.__end:
            # No set exhausted, return the first of the two
            if (self.A_last < self.B_last):
                tmp = self.A_last
                self.A_last = next(self.A_iter, self.__end)
            else:
                tmp = self.B_last
                self.B_last = next(self.B_iter, self.__end)
        elif self.A_last is self.__end and self.B_last is self.__end:
            # Both exhausted
            raise StopIteration
        elif self.A_last is self.__end:
            # A exhausted
            tmp = self.B_last
            self.B_last = next(self.B_iter, self.__end)
        else:
            # B exhausted
            tmp = self.A_last
            self.A_last = next(self.A_iter, self.__end)
        return tmp

def union(*A):
    """Merge ordered tuples."""

    return tuple(SortedUnion(*A))

class SortedDifference(object):
    """Iterator returning the difference of sorted iterators."""

    __end = object()

    def __init__(self, A, B):

        self.A_iter = iter(A)
        self.A_last = next(self.A_iter, self.__end)
        self.B_iter = iter(B)
        self.B_last = next(self.B_iter, self.__end)

    def __iter__(self):

        return self

    def next(self):

        while True:
            if self.A_last is self.__end:
                raise StopIteration
            if self.B_last is self.__end:
                tmp = self.A_last
                self.A_last = next(self.A_iter, self.__end)
                return tmp
            if self.A_last < self.B_last:
                tmp = self.A_last
                self.A_last = next(self.A_iter, self.__end)
                return tmp
            elif self.A_last == self.B_last:
                self.A_last = next(self.A_iter, self.__end)
                self.B_last = next(self.B_iter, self.__end)
            else:
                self.B_last = next(self.B_iter, self.__end)

def difference(A, B):
    """Difference of ordered tuples."""

    return tuple(SortedDifference(A, B))

#=========================================================================================
# SubtractionLeg
#=========================================================================================

class SubtractionLeg(tuple):
    """Leg object specialized for subtraction.
    A SubtractionLeg is (LEG_ID,PDG_ID,STATE) where
    - LEG_ID is the ID of this specific leg in the process
    - PDG_ID identifies the type of particle through the standard PDG codes
    - STATE is FINAL or INITIAL
    Initialization is performed using
    - SubtractionLeg(LEG_ID:int,PDG_ID:int,STATE:bool)
    - SubtractionLeg(Leg) using the Leg objects
    - SubtractionLeg(Leg_dict:dict) using a dict with keys  'number', 'id', 'state'
    """

    FINAL = base_objects.Leg.FINAL
    INITIAL = base_objects.Leg.INITIAL

    def __new__(cls, *args, **opts):
        """Initialize a SubtractionLeg object from various representations."""

        target = (0, SubtractionLeg.FINAL)
        # One argument passed
        if len(args) == 1:
            leg = args[0]
            # Initialization from Leg object
            if isinstance(leg, base_objects.Leg):
                target = (leg.get('number'), leg.get('id'), leg.get('state'))
            # Initialization from dictionary
            elif isinstance(leg, dict):
                target = (leg['number'], leg['id'], leg['state'])
            # Initialization from iterable
            elif (
                hasattr(leg, '__len__') and len(leg) == 3 and
                isinstance(leg[0], int) and isinstance(leg[1], int) and
                isinstance(leg[2], bool)
            ):
                target = leg
            # Invalid argument
            else:
                raise MadGraph5Error(
                    "SubtractionLeg cannot be initialized with "
                    "%s of type %s." % (leg, str(type(leg)))
                )
        # Initialize with 3 arguments
        elif (
            (len(args) == 3) and
            isinstance(args[0], int) and isinstance(args[1], int) and
            isinstance(args[2], bool)
        ):
            target = (args)
        # Invalid number of arguments
        else:
            raise MadGraph5Error(
                "SubtractionLeg cannot be initialized with argument %s." % str(args)
            )

        return super(SubtractionLeg, cls).__new__(cls, target)

    @property
    def n(self):

        return self[0]

    @property
    def pdg(self):

        return self[1]

    @property
    def state(self):

        return self[2]

    @staticmethod
    def state_str(state):

        if state == SubtractionLeg.FINAL:
            return 'f'
        return 'i'

    def __str__(self, print_n=True, print_pdg=False, print_state=False):
        """Return a string representation of this subtraction leg."""

        tmp = [str(self.n), str(self.pdg), SubtractionLeg.state_str(self.state)]
        if not print_state:
            del tmp[2]
        if not print_pdg:
            del tmp[1]
        if not print_n:
            del tmp[0]
        if tmp:
            if len(tmp) == 1:
                return tmp[0]
            else:
                return "(" + ", ".join(tmp) + ")"
        return ""

#=========================================================================================
# SubtractionLegSet
#=========================================================================================

class SubtractionLegSet(tuple):
    """Set of SubtractionLeg objects."""

    def __new__(cls, iterable=(), **opts):
        """Initialize set, trying to convert arguments into SubtractionLeg's."""    
        
        if isinstance(iterable, SubtractionLegSet):
            legs = tuple(SubtractionLeg(leg) for leg in iterable)
        else:
            legs = sorted(SubtractionLeg(leg) for leg in iterable)
        return super(SubtractionLegSet, cls).__new__(cls, legs)

    def __init__(self, iterable=(), **opts):
        """Initialize set, trying to convert arguments into SubtractionLeg's."""

        if isinstance(iterable, SubtractionLegSet):
            legs = tuple(SubtractionLeg(leg) for leg in iterable)
        else:
            legs = sorted(SubtractionLeg(leg) for leg in iterable)
        super(SubtractionLegSet, self).__init__(legs)

    def has_initial_state_leg(self):
        """Return True if this leg set has at least one initial state leg,
        False otherwise.
        """

        for leg in self:
            if leg.state == SubtractionLeg.INITIAL:
                return True
        return False

    def without_leg_numbers(self, discard_initial_leg_numbers=True):
        """Return a copy of this SubtractionLegSet without leg numbers."""
        
        if discard_initial_leg_numbers:
            return SubtractionLegSet((0, leg.pdg, leg.state) for leg in self)
        else:
            return SubtractionLegSet(
                (0, leg.pdg, leg.state) if leg.state==SubtractionLeg.FINAL else 
                (leg.n, leg.pdg, leg.state) for leg in self
            )

    def split_by_pdg_abs(self, pdgs, already_matched_PDGs):
        """Restructure the information about the SubtractionLegSet,
        grouping legs by absolute value of their pdg and their state.
        For each pdg abs, generate two sets of legs with opposite pdg sign,
        and count the positive/negative occurrences.
        This information is useful for the map_leg_numbers method.

        For example, calling
        split_by_pdg_abs(
            ((1, 2=u, I), (2, 1=d, F), (3, -2=u~, F), (4, 2=u, F), (5, 2=u, F), (6, 3=s, F), ),
            [2=u|u~, 3=s|s~],
            F )
        where I/F stands for initial/final yields
            legs_by_pdg = {
                2=u|u~: [((3, -2=u~, F), ), ((4, 2=u, F), (5, 2=u, F), )],
                3=s|s~: [(), ((6, 3=s, F), )],
            }
        as, for the specified pdgs 2=u|u~ and 3=s|s~, there are respectively
        1 negative / 2 positive and 0 negative / two positive legs.
        Note that the two sets of legs for a given pdg abs are sorted by leg count
        and not by positive/negative, and only the legs with state F are considered.
        The counts are expressed by
            pdg_counts = {(0, 1): [3=s|s~], (1, 2): [2=u|u~]}
        which means that the pdg abs values
        with 0 legs of one sign and 1 leg of the opposite are the list [3=s|s~],
        and the ones with 1 leg of one sign and 2 of the other are the list [2=u|u~].
        The function also returns the SubtractionLegSet without the legs
        that have a pdg/state which was specified, i.e.
            rest = ((1, 2=u, I), (2, 1=d, F), ).
        """

        legs_by_pdg = dict()
        pdg_counts = dict()
        rest = SubtractionLegSet(self)
        for pdg in set(map(abs, pdgs)):
            # Ignore PDGs that have already been matched
            if pdg in already_matched_PDGs:
                continue
            pos_initial = SubtractionLegSet(leg for leg in self
                                    if leg.pdg ==  pdg and leg.state == leg.INITIAL)
            pos_final = SubtractionLegSet(leg for leg in self
                                    if leg.pdg ==  pdg and leg.state == leg.FINAL)
            neg_initial = SubtractionLegSet(leg for leg in self
                                    if leg.pdg == -pdg and leg.state == leg.INITIAL)
            neg_final = SubtractionLegSet(leg for leg in self
                                    if leg.pdg == -pdg and leg.state == leg.FINAL)
            legs_by_pdg[pdg] = sorted((pos_initial,neg_initial), key=len)+\
                                sorted((pos_final,neg_final), key=len)
            count_key = tuple(map(len, legs_by_pdg[pdg]))
            pdg_counts.setdefault(count_key, [])
            pdg_counts[count_key].append(pdg)
            rest = SubtractionLegSet(difference(rest, union(pos_initial, pos_final,neg_initial,neg_final)))
        return legs_by_pdg, pdg_counts, rest


    def map_leg_numbers(self, target, already_matched_PDGs, equivalent_pdg_sets=()):
        """Find a correspondence between this SubtractionLegSet and another
        that only differs by leg numbers.
        The optional argument equivalent_pdg_sets allows to consider some particle species
        to be interchangeable in this matching. This is designed to identify processes
        where legs have the same quantum numbers up to some details which are irrelevant
        for subtraction. Typically this is the case of massless quarks, where one only
        keeps track of which quarks have the same flavor or are antiparticles.
        For instance
            (1, 1=d, F), (2, -1=d~, F)
        should match
            (7, 3=s, F), (5, -3=s~, F)    or    (4, -1=d~, F), (1, 1=d,    F),
        but not
            (4, 3=s, F), (5, -1=d~, F)    or    (7, 24=W+, F), (2, -24=W-, F).

        :param target: the target SubtractionLegSet
        :param equivalent_pdg_sets: a list of lists pdgs that should be considered
            equivalent in the matching (note that particles and antiparticles
            for these values will also be considered equivalent)
        :return: a dictionary of leg numbers that sends this set into the target one,
            if any; None otherwise.
        """

        if len(self) != len(target):
            return None

        s_rest = self
        t_rest = target
        leg_numbers_map = dict()

        new_already_matched_PDGs = dict(already_matched_PDGs)

        for pdgs in equivalent_pdg_sets:

            s_legs_by_pdgs, s_pdg_counts, s_rest = s_rest.split_by_pdg_abs(pdgs, already_matched_PDGs.keys())
            t_legs_by_pdgs, t_pdg_counts, t_rest = t_rest.split_by_pdg_abs(pdgs, already_matched_PDGs.values())
            for count in s_pdg_counts.keys():
                for i in range(len(s_pdg_counts[count])):
                    s_pdg = s_pdg_counts[count][i]
                    try:
                        t_pdg = t_pdg_counts[count][i]
                    except (KeyError, IndexError):
                        return None
                    if any(len(l1)!=len(l2) for l1,l2 in zip(s_legs_by_pdgs[s_pdg],t_legs_by_pdgs[t_pdg])):
                        raise MadGraph5Error(
                            "Inconsistent lists in SubtractionLegSet.map_leg_numbers, "
                            "SubtractionLegSet.split_by_pdg_abs is bugged.")
                    for l1, l2 in zip(s_legs_by_pdgs[s_pdg], t_legs_by_pdgs[t_pdg]):
                        leg_numbers_map.update(zip(
                            [leg.n for leg in l1],
                            [leg.n for leg in l2] ))
                        # Specify that now this matching leg pdgs (both positive and negative) are assigned
                        new_already_matched_PDGs.update(zip(
                            [leg.pdg for leg in l1],
                            [leg.pdg for leg in l2] ))
                        new_already_matched_PDGs.update(zip(
                            [-leg.pdg for leg in l1],
                            [-leg.pdg for leg in l2] ))

        for s_leg in s_rest:
            try:
                t_leg = next(leg for leg in t_rest if leg.state == s_leg.state and ( leg.pdg == s_leg.pdg if
                                s_leg.pdg not in already_matched_PDGs else already_matched_PDGs[s_leg.pdg]==leg.pdg) )
                leg_numbers_map[s_leg.n] = t_leg.n
                t_rest = SubtractionLegSet(difference(t_rest, (t_leg, )))
            except StopIteration:
                return None

        assert(not t_rest)

        # Now specify which template PDGs have now been assigned to particular PDG values
        already_matched_PDGs.update(new_already_matched_PDGs)

        return leg_numbers_map


#=========================================================================================
# SingularStructure
#=========================================================================================

class SingularStructure(object):
    """Object that represents a hierarchical structure of IR singularities."""


    def __init__(self, *args, **opts):
        """Initialize a hierarchical singular structure."""

        self.substructures = opts.get('substructures', [])
        if not isinstance(self.substructures, list):
            self.substructures = list(self.substructures)
        for substructure in self.substructures:
            assert isinstance(substructure, SingularStructure)
        self.legs = opts.get('legs', SubtractionLegSet())
        if not isinstance(self.legs, SubtractionLegSet):
            self.legs = SubtractionLegSet(self.legs)
        self.is_void = opts.get('is_void', False)
        assert isinstance(self.is_void, bool)
        for arg in args:
            if isinstance(arg, SingularStructure):
                self.substructures.append(arg)
            elif isinstance(arg, SubtractionLeg):
                self.legs = SubtractionLegSet(self.legs + (arg, ))
            elif isinstance(arg, base_objects.Leg):
                self.legs = SubtractionLegSet(self.legs + (SubtractionLeg(arg),))
            elif isinstance(arg, bool):
                self.is_void = self.is_void or arg
            else:
                raise MadGraph5Error(
                    "Invalid argument in SingularStructure.__init__: %s" % str(arg) )

    def get_copy(self):
        """Provide a modifiable copy of this singular structure."""

        return type(self)(
            legs=SubtractionLegSet(self.legs),
            substructures=[ss.get_copy() for ss in self.substructures],
            is_void=self.is_void )

    def __str__(self, print_n=True, print_pdg=False, print_state=False, sort=True):
        """Return a string representation of the singular structure."""

        if self.is_void:
            return 'Void structure'

        tmp_str = self.name() + "("
        sub_strs = ( sub.__str__(print_n, print_pdg, print_state, sort)
            for sub in self.substructures )
        if sort:
            sub_strs = sorted(sub_strs)
        tmp_str += ",".join(sub_strs)
        if self.substructures:
            tmp_str += ","
        tmp_str += ",".join(
            leg.__str__(print_n, print_pdg, print_state)
            for leg in self.legs
        )
        tmp_str += ")"
        return tmp_str

    def does_require_correlated_beam_convolution(self, subtraction_scheme_module):
        """ Checks whether this integrated counterterm requires a host contribution featuring
        correlated convolution of the beams (i.e. 'BS', 'VS', etc... contributions).
        Whether a correlated beam convolution is needed or not depends on the subtraction scheme
        module and this is why this decision is made through one of its functions."""

        return subtraction_scheme_module.exporter.does_require_correlated_beam_convolution(self)

#gl
    def does_require_torino_sub_BS(self, subtraction_scheme_module):
        """ Checks whether this integrated counterterm requires a host contribution featuring
        correlated convolution of the beams (i.e. 'BS', 'VS', etc... contributions).
        Whether a correlated beam convolution is needed or not depends on the subtraction scheme
        module and this is why this decision is made through one of its functions."""

        return subtraction_scheme_module.exporter.does_require_torino_sub_BS(self)

    def get_subtraction_prefactor(self):
        """Determine the prefactor due to the nested subtraction technique."""
        
        # If we are at the top level, we shall not include the factor
        # Indeed, we are integrating +R-C.
        if type(self) == SingularStructure:
            pref = 1
        else:
            pref = -1

        # Account for simplification of soft structures
        # TODO Check this works for collinears inside softs
        softs = []
        for sub in self.substructures:
            pref *= sub.get_subtraction_prefactor()
            if sub.name() == "S":
               softs += [len(sub.substructures) + len(sub.legs), ]
        pref *= multinomial(softs)
        
        return pref

    def get_beam_factorization_legs(self):
        """ Return the set of all the beam factorization legs involved in this singular
        structure. Note that beam factorization structures should always be at the top level
        and never nested."""
        res = set([])
        for struct in self.substructures:
            if isinstance(struct, BeamStructure):
                res |= set([l.n for l in struct.get_all_legs() ])
            else:
                res |=  struct.get_beam_factorization_legs()
        return res

    def annihilate(self):
        """When an elementary structure cannot be nested within this structure,
        remove this structure altogether by flagging it as void.
        """

        # Really clear structures and legs
        del self.substructures[:]
        self.legs = SubtractionLegSet()
        self.is_void = True

    def get_all_legs(self):
        """Return all legs involved in this singular structure."""

        all_legs = [ self.legs, ]
        for sub in self.substructures:
            all_legs.append(sub.get_all_legs())
        return SubtractionLegSet(union(*all_legs))

    def discard_leg_numbers(self, discard_initial_leg_numbers=True):
        """Set all leg numbers to zero, discarding this information forever."""

        self.legs = self.legs.without_leg_numbers(discard_initial_leg_numbers)
        for substructure in self.substructures:
            substructure.discard_leg_numbers(discard_initial_leg_numbers)

    def count_unresolved(self):
        """Count the number of unresolved legs."""

        total_unresolved = 0
        for substructure in self.substructures:
            total_unresolved += substructure.count_unresolved()
        return total_unresolved

    def unpack(self):
        """ Unpacks this nested singular structure (like     ((C(1,2),C(3,4)),)      ) into one that features only
        a single embedding singular structure (e.g. (C(1,2),C(3,4)) in that example) """

        ss = self
        while len(ss.substructures)>0 and ss.substructures[0].name()=='' and len(ss.substructures)==1:
            ss = ss.substructures[0]
        return ss

    def count_couplings(self):
        """Count the number of couplings covered by this structure."""

        total_couplings = 0
        for substructure in self.substructures:
            total_couplings += substructure.count_couplings()
        return total_couplings

    def get_canonical_representation(self, track_leg_numbers=True, sort=True):
        """Creates a canonical hashable representation of self."""
        
        canonical = dict()

        canonical['is_void'] = self.is_void
        if track_leg_numbers:
            canonical['legs'] = self.legs
        else:
            # Always track initial state leg numbers as they are always distinguishable
            canonical['legs'] = self.legs.without_leg_numbers(
                discard_initial_leg_numbers=False)
        canonical['name'] = self.name()
        canonical['substructures'] = tuple(
            structure.get_canonical_representation(track_leg_numbers)
            for structure in self.substructures)
        if sort:
            canonical['substructures'] = tuple(sorted(canonical['substructures']))
        return tuple(sorted(canonical.items()))

    def __hash__(self):

        return hash(self.get_canonical_representation())

    def __eq__(self, other, orderless=True):
        """Check if two singular structures are the same."""

        if not isinstance(other, SingularStructure):
            raise TypeError(
                "Comparing SingularStructure to %s (%s)" % (str(other), str(type(other))))
        self_can = self.get_canonical_representation(sort=orderless)
        other_can = other.get_canonical_representation(sort=orderless)
        return self_can == other_can

    def __ne__(self, other):

        return not (self == other)

    def name(self):

        return ""

    def decompose(self):
        """Decompose the singular structure as a flat list."""

        inner = list()
        for substructure in self.substructures:
            inner += substructure.decompose()
        if type(self) == SingularStructure:
            return inner
        foo = type(self)(legs=self.get_all_legs())
        inner.append(foo)
        return inner

    def act_here_needed(self, structure):
        """Return True if the counterterm obtained adding this limit
        on the given structure at the first level is needed for subtraction.
        """

        assert not self.substructures
        raise MadGraph5Error(
            "act_here_needed called in SingularStructure of unspecified type.")

    def non_overlapping_with(self, structure):
        """Determines if this limit does not affect any leg
        that appears within the given structure, at any level.
        """

        assert isinstance(structure, SingularStructure)
        assert not self.substructures
        # Check if the structure has any legs in common
        if not is_disjoint(self.legs, structure.legs):
            return False
        # Check if substructures have any leg in common
        for substructure in structure.substructures:
            if not self.non_overlapping_with(substructure):
                return False
        # All checks passed
        return True

    def act_on(self, structure):
        """Act with this elementary structure on another, non-elementary structure."""

        assert isinstance(structure, SingularStructure)
        assert not self.substructures

        # If the limit does not overlap with the existing structure at all,
        # just append it at the end
        if self.non_overlapping_with(structure):
            structure.substructures.append(self.get_copy())
            return

        # If the limit acts at least partly at the current level
        if not is_disjoint(self.legs, structure.legs):
            # If the limit acts at different levels,
            # it is not needed for subtraction
            if not is_subset(self.legs, structure.legs):
                structure.annihilate()
            # If the limit acts completely within this level,
            # it may be needed or not
            else:
                if self.act_here_needed(structure):
                    structure.substructures.append(self.get_copy())
                    structure.legs = SubtractionLegSet(difference(structure.legs, self.legs))
                else:
                    structure.annihilate()
            return

        # The limit acts deeper
        for substructure in structure.substructures:
            # Find the first substructure hit by the limit, even partially
            if not self.non_overlapping_with(substructure):
                self.act_on(substructure)
                # If the action on the substructure annihilated it,
                # then annihilate this structure as well
                if substructure.is_void:
                    structure.annihilate()
                return

    def nest(self, structure=None):

        assert isinstance(self, SingularStructure)
        if structure is None:
            structure = SingularStructure()
        substructures = self.substructures
        # Empty list of elementary structures to apply or invalid target structure: done
        if structure.is_void or not substructures:
            return structure
        # Apply last elementary structure and advance recursion
        most_structures = SingularStructure(substructures=substructures[:-1])
        substructures[-1].act_on(structure)
        # And now act using the rest of the elementary structures, if needed
        if structure.is_void:
            return structure
        return most_structures.nest(structure)

    @staticmethod
    def from_string(string, process, log=False):
        """Reconstruct a singular structure from a string representation,
        e.g. C(S(3),5,4)

        :return: The reconstructed singular structure if successful, else None.
        """

        singular_structures_name_dictionary = {
            '' : SingularStructure,
            'S': SoftStructure,
            'C': CollStructure,
            'F': BeamStructure
        }

        msg = "SingularStructure.from_string:"
        openpos = string.find('(')
        closepos = string.rfind(')')
        name = string[:openpos]
        if not singular_structures_name_dictionary.has_key(name):
            if log:
                logger.warning('%s %s: %s'%(msg, "no singular structure with name", name))
            return None
        legs = []
        substructures = []
        content = string[openpos+1:closepos]
        before_comma = ""
        while content:
            # print "Content is", content
            next_comma_pos = content.find(',')
            if next_comma_pos != -1:
                before_comma += content[:next_comma_pos]
            else:
                before_comma += content
                content = ""
            # print "Before comma:", before_comma
            content = content[next_comma_pos+1:]
            if '(' in before_comma or ')' in before_comma:
                if before_comma.count('(') == before_comma.count(')'):
                    substructure = SingularStructure.from_string(before_comma, process)
                    if substructure is None:
                        return None
                    substructures.append(substructure)
                    before_comma = ""
                else:
                    before_comma += ","
            else:
                try:
                    leg_n = int(before_comma)
                except:
                    if log:
                        logger.warning('%s %s %s'%("could not convert", before_comma, "to integer"))
                    return None
                for leg in process['legs']:
                    if leg['number'] == leg_n:
                        legs.append(leg)
                        before_comma = ""
                        break
                if before_comma:
                    if log:
                        logger.warning('%s %s %d %s'%(msg, "no leg number", leg_n, "found"))
                    return None
        if before_comma:
            if log:
                logger.warning("%s %s '%s'"%(msg, "bracket matching failed for:", before_comma[:-1]))
            return None
        return singular_structures_name_dictionary[name](
            substructures=substructures, legs=legs)

    def map_leg_numbers(self, target, equivalent_pdg_sets=(), already_matched_PDGs=None):
        """Find a correspondence between this SingularStructure and another
        that only differs by leg numbers.
        The optional argument equivalent_pdg_sets allows to consider some particle species
        to be interchangeable in this matching. This is designed to identify identical
        structures whose legs have the same quantum numbers up to details which are
        irrelevant for subtraction. Typically this is the case of massless quarks, where
        one only keeps track of which quarks have the same flavor or are antiparticles.
        For instance
            C(S((1, 21=g, F)), (2, -1=d~, F))
        should match
            C(S((3, 21=g, F)), (2, +3=s,  F))    or    C(S((7, 21=g, F)), (4, -1=d~, F))
        but not
            S(S((1, 21=g, F)), (2, -1=d~, F))    or    C(S((6, 21=g, F)), (4, 21=g,  F)).

        :param target: the target SingularStructure
        :param equivalent_pdg_sets: a list of lists pdgs that should be considered
            equivalent in the matching (note that particles and antiparticles
            for these values will also be considered equivalent)
        :return: a dictionary of leg numbers that sends this structure
            into the target one, if any; None otherwise.
        """

        # Set the map of already matched PDG if not already done:
        if not already_matched_PDGs:
            already_matched_PDGs = {}

        # 1) Match type
        if self.name() != target.name():
            return None

        # 2) Match legs
        leg_number_dict = self.legs.map_leg_numbers(target.legs, already_matched_PDGs, equivalent_pdg_sets)
        if leg_number_dict is None:
            return None

        # 3) Match substructures
        if len(self.substructures) != len(target.substructures):
            return None

        target_subs = [target_sub for target_sub in target.substructures]
        for sub in self.substructures:
            matching = None
            for i in range(len(target_subs)):
                sub_dict = sub.map_leg_numbers(target_subs[i], equivalent_pdg_sets, already_matched_PDGs=already_matched_PDGs)
                if sub_dict is not None:
                    matching = target_subs.pop(i)
                    leg_number_dict.update(sub_dict)
                    break
            if matching is None:
                return None

        return leg_number_dict

class SoftStructure(SingularStructure):

    def count_couplings(self):
        return self.count_unresolved()

    def count_unresolved(self):
        return len(self.get_all_legs())

    def name(self):
        return "S"

    def act_here_needed(self, structure):

        # Soft limits are not needed within soft structures
        if structure.name() == "S":
            return False

        # Soft limits are needed within collinear structures
        # if there is at least one hard particle left
        if structure.name() == "C":
            # If there is a hard particle at this level, the limit is needed
            if len(self.legs) < len(structure.legs):
                return True
            # If there is a hard collinear substructure, the limit is needed
            else:
                for substructure in structure.substructures:
                    if not substructure.name() == "S":
                        return True
            # If did not return yet, there is no hard particle left
            return False

        # By default, not skipping overlaps ensures correct subtraction
        return True

class CollStructure(SingularStructure):

    def count_couplings(self):
        return self.count_unresolved()

    def count_unresolved(self):
        return len(self.get_all_legs())-1

    def name(self):
        return "C"

    def act_here_needed(self, structure):
        return True

class BeamStructure(SingularStructure):

    # TODO this may only works at NLO...
    def count_couplings(self):
        return 1

    def count_unresolved(self):
        return 0

    def name(self):
        return "F"

    def act_here_needed(self, structure):
        return False

#=========================================================================================
# Currents
#=========================================================================================
class CurrentsBlock(tuple):
    """ This class implements the concept of a tuple of Currents which must be implemented
    at once in a single subtraction current class of the subtraction scheme because it
    corresponds to a non-factorisable combination of currents."""

    def __init__(self, *args, **opts):
        super(CurrentsBlock, self).__init__(*args, **opts)

    def to_single_current(self):

        if len(self)!=1:
            raise MadGraph5Error("A CurrentsBlock instance cannot be cast into a simple "+
                                 " current if it contains more than one current.")
        return self[0]

    def does_resolve_mother_spin_and_color(self):
        return any(c['resolve_mother_spin_and_color'] for c in self)

    def get_squared_orders(self):
        """ Return the aggregated squared orders for this currents block."""

        aggregated_squared_orders = {}
        for crt in self:
            for order, value in crt.get('squared_orders').items():
                if order in aggregated_squared_orders:
                    aggregated_squared_orders[order] += value
                else:
                    aggregated_squared_orders[order] = value
        return aggregated_squared_orders

    def __str__(self, *args, **opts):
        return '(%s)'%(', '.join('%s'%crt.__str__(*args, **opts) for crt in self))

    def get_key(self, *args, **opts):

        from madgraph.core.accessors import CompoundProcessKey
        return CompoundProcessKey(crt.get_key(*args, **opts) for crt in self)

#=========================================================================================
# Current
#=========================================================================================
class Current(base_objects.Process):

    def default_setup(self):
        """Default values for all properties specific to current,
        additional to Process.
        """

        super(Current, self).default_setup()
        # Assign a default value of n_loops to -1 which will be updated later
        # during current generation.
        self['n_loops'] = -1
        # And also by default assign no squared orders
        self['squared_orders'] = {}
        # If this current is directly connected to the underlying ME,
        # it will need to resolve spin and color of the parent leg
        # (i.e. specify the corresponding correlators).
        # If not, then it will be summed over.
        # WARNING defaulting to False might give wrong results instead of crashing
        # default to None?
        self['resolve_mother_spin_and_color'] = False
        self['singular_structure'] = SingularStructure()

        self['model'] = None

    def does_resolve_mother_spin_and_color(self):
        return self['resolve_mother_spin_and_color']

    def set(self,name,value, **opts):

        # There is no need to store a model for this class, so allow None
        if name=='model':
            base_objects.PhysicsObject.set(self, 'model', None, force=True)
        else:
            return super(Current, self).set(name, value, **opts)

    def count_unresolved(self):
        """Count the number of unresolved particles covered by this current."""

        return self['singular_structure'].count_unresolved()

    def count_couplings(self):
        """Count the number of couplings expected in this current."""

        return self['singular_structure'].count_couplings()

    def discard_leg_numbers(self, discard_initial_leg_numbers=True):
        """Discard all leg numbers in the singular_structure
        and in the parent_subtraction_leg,
        effectively getting the type of current
        without information about momenta for which it should be evaluated.
        """

        self['singular_structure'].discard_leg_numbers(discard_initial_leg_numbers)

    def __str__(
        self,
        print_n=True, print_pdg=True, print_state=True,
        print_loops=True, print_orders=True
    ):
        """Convert this current to a nice readable string."""

        readable_string = self['singular_structure'].__str__(
            print_n, print_pdg, print_state
        )
        n_loops = self.get('n_loops')
        if n_loops >= 0 and print_loops:
            readable_string += " @ " + str(n_loops) + " loop"
            if n_loops != 1:
                readable_string += "s"
        if self.get('orders') and print_orders:
            readable_string += " " + str(self.get('orders'))
        return readable_string

    def does_require_correlated_beam_convolution(self, *args, **opts):
        """ Checks whether this integrated counterterm requires a host contribution featuring
        correlated convolution of the beams (i.e. 'BS', 'VS', etc... contributions).
        This is the case only for integrated counterterms with disjoint structures having each
        at least one soft *BeamCurrents* that originated from a colorful subtraction currents 
        scheme with mappings such as ppToOneWalker that recoils the soft momentum equally 
        against both initial-state beams."""

        # For basic currents, this is never the case
        return False

#gl
    def does_require_torino_sub_BS(self, *args, **opts):
        """ Checks whether this integrated counterterm requires a host contribution featuring
        correlated convolution of the beams (i.e. 'BS', 'VS', etc... contributions).
        This is the case only for integrated counterterms with disjoint structures having each
        at least one soft *BeamCurrents* that originated from a colorful subtraction currents 
        scheme with mappings such as ppToOneWalker that recoils the soft momentum equally 
        against both initial-state beams."""

        # For basic currents, this is never the case
        return False

    def get_reduced_flavors(self, defining_flavors, IR_subtraction_module, routing_dict):
        """Given the defining flavors dictionary specifying the flavors of some leg numbers,
        add entries corresponding to the subtraction legs in this node and using the function
        get_parent_PDGs from the subtraction module."""

        get_particle = IR_subtraction_module.model.get_particle
        
        def get_parent(PDGs, is_initial):
            res = IR_subtraction_module.parent_PDGs_from_PDGs(PDGs)
            if len(res)!=1:
                str1 = "Multiple mother PDGs assignment "
                str2 = " is not yet supported by the function get_reduced_flavors"
                misc.sprint("get_parent called with PDGs = " + str(PDGs))
                raise MadGraph5Error(str1 + str(res) + str2)
            # Now flip back the identity of the parent PDG if in the initial state
            if is_initial:
                return get_particle(res[0]).get_anti_pdg_code()
            else:
                return res[0]

        # Go through the current's singular structure to find the disjoint unresolved sets
        ss = self.get('singular_structure')

        # In order to simplify the computation of the reduced flavours we unpack multiple nestings
        # like (C(1,2),C(3,4))  (used for local counterterms) ((C(1,2),C(3,4)),) (used for integrated CT)
        # into a structure that has a single embedding singular structure.
        ss = ss.unpack()

        for sub_ss in ss.substructures:
            all_legs = sub_ss.get_all_legs()

            # If all legs are unresolved, as it is the case in the pure soft configurations,
            # then nothing needs to be done for this substructure
            if sub_ss.count_unresolved() == len(all_legs):
                continue

            leg_numbers = [leg.n for leg in all_legs]

            try:
                # Swap the flavors of the initial states
                leg_flavors = [(
                    get_particle(defining_flavors[leg.n]).get_anti_pdg_code() if
                    leg.state == SubtractionLeg.INITIAL else get_particle(defining_flavors[leg.n]).get_pdg_code())
                    for leg in all_legs]

                if len(leg_numbers)>0:
                    mother_number = Counterterm.get_ancestor(frozenset(leg_numbers), routing_dict)
                    if mother_number not in defining_flavors:
                        defining_flavors[mother_number] = get_parent(leg_flavors,
                                       any(leg.state==SubtractionLeg.INITIAL for leg in all_legs) )
            except Exception as e:
                raise MadGraph5Error('\nError when computing reduced quantities for current:\n%s\n'%str(self)+
                                     'Defining flavors:\n%s\n'%str(defining_flavors)+
                                     'routing dict:\n%s\n'%str(routing_dict)+
                                     'Exception encountered:\n%s'%traceback.format_exc())

    def get_squared_orders(self):
        """ Necessary for common API with CurrentsBlock"""
        return self.get('squared_orders')

    def get_key(self, track_leg_numbers=False):
        """Return the ProcessKey associated to this current."""

        from madgraph.core.accessors import ProcessKey
        return ProcessKey(
                self,
                allowed_attributes = [
                    'singular_structure',
                    'n_loops',
                    'squared_orders',
                    'perturbation_couplings',
                    'resolve_mother_spin_and_color' ],
                track_leg_numbers = track_leg_numbers)

    def get_initial_state_leg_map(self, leg_numbers_map):
        """ Returns None if this current does not involve initial states, otherwise it returns 
        pair of SubtractionLeg instances in the form of the 2-tuple 
                             (mother_IS_leg, daughter_IS_leg)
        Notice that beam factorization across more than one (i.e. two) beams must be 
        factorized into two currents, so we should be guaranteed having to return at most
        one pair of IS legs here."""
        
    
        # First get the ancestors of all outermost leg numbers of this singular structures
        ancestors = Counterterm.get_ancestor( 
         [ leg.n for leg in self['singular_structure'].get_all_legs() ],  leg_numbers_map )
        
        # Filter initial state legs and make sure there is only one
        ancestors = [leg for leg in ancestors if leg.state==SubtractionLeg.INITIAL]
        
        # This current has no initial-state
        if len(ancestors) == 0:
            return None
        
        if len(ancestors)!=1:
            raise MadGraph5Error('Beam factorization current should involve'+
                ' exactly one ancestor initial state leg (the two beam factorization must be factorized).')
        initial_state_ancestor = ancestors[0]
        
        # Now find which external initial-state leg is connected to the initial
        # state ancestor identified.
        daughters = Counterterm.get_daughter_legs(
                                                   initial_state_ancestor, leg_numbers_map)
        # Filter initial state legs and make sure there is only one
        daughters = [leg for leg in daughters if leg.state==SubtractionLeg.INITIAL]
        if len(daughters)!=1:
            raise MadGraph5Error('Beam factorization current should involve'+
                ' exactly one daughter initial state leg (the two beam factorization must be factorized).')        
        initial_state_daughter = daughters[0]
        
        return (initial_state_ancestor, initial_state_daughter)

    def get_copy(self, copied_attributes=()):
        """Return a copy of this current with a deep-copy of its singular structure."""

        copied_current = super(Current, self).get_copy(
            tuple(attr for attr in copied_attributes if attr != 'singular_structure') )
        if 'singular_structure' in copied_attributes:
            copied_current['singular_structure'] = self['singular_structure'].get_copy()

        return copied_current

    @classmethod
    def get_type(cls):
        """ Return the type of current implemented by this instance."""

        return cls.__name__

    def get_number_of_couplings(self):
        """ Return a naive count of the number of couplings involved in this current."""
        return self.get('n_loops') + self['singular_structure'].count_couplings()

    def __eq__(self, other):
        """Compare two currents using their ProcessKey's."""

        return self.get_key().key_dict == other.get_key().key_dict

#=========================================================================================
# BeamCurrent
#=========================================================================================
class BeamCurrent(Current):
    """ A special class of current representing beam factorization terms 
    (aka PDF counterterms).

    If we have:
        PDF_CT(xi) = F_+(xi) + [F] \delta(xi-1)
    (which can also be explicitely written)
        PDF_CT(xi) = F(xi) + {F(xi)} \delta(xi-1) + [F] \delta(xi-1)

    Then this current represent the two pieces of the xi-dependent distribution F_+(xi),
    the first piece (F(xi)) with 'distribution_type' set to 'bulk' and the latter 
    ({F(xi)} \delta(xi-1)) with this attribute set to 'counterterm'.
    """
    
    def default_setup(self):
        """Default values for all properties specific to the BeamCurrent,
        additional to Process.
        """
        super(BeamCurrent, self).default_setup()
        self['n_loops'] = -1
        self['squared_orders'] = {}
        self['resolve_mother_spin_and_color'] = True
        self['singular_structure'] = SingularStructure()

        # Which piece of the plus distribution is intended to be returned
        # This can be 'bulk' or 'counterterm'.
        self['distribution_type'] = 'bulk'
        
        # List of the pdgs of the particles making up this beam
        self['beam_PDGs'] = []
        # Name of the beam factorization type that should be employed
        self['beam_type'] = None

    def does_require_correlated_beam_convolution(self, subtraction_scheme_module):
        """ Checks whether this integrated counterterm requires a host contribution featuring
        correlated convolution of the beams (i.e. 'BS', 'VS', etc... contributions).
        This is the case only for integrated counterterms with disjoint structures having each
        at least one soft *BeamCurrents* that originated from a colorful subtraction currents
        scheme with mappings such as ppToOneWalker that recoils the soft momentum equally
        against both initial-state beams. Given the dependency on the details of the implementation
        of the subbtraction scheme, this decision is made through a function of the subtraction scheme."""

        return self['singular_structure'].does_require_correlated_beam_convolution(subtraction_scheme_module)

#gl
    def does_require_torino_sub_BS(self, subtraction_scheme_module):
        """ Checks whether this integrated counterterm requires a host contribution featuring
        correlated convolution of the beams (i.e. 'BS', 'VS', etc... contributions).
        This is the case only for integrated counterterms with disjoint structures having each
        at least one soft *BeamCurrents* that originated from a colorful subtraction currents
        scheme with mappings such as ppToOneWalker that recoils the soft momentum equally
        against both initial-state beams. Given the dependency on the details of the implementation
        of the subbtraction scheme, this decision is made through a function of the subtraction scheme."""

        return self['singular_structure'].does_require_torino_sub_BS(subtraction_scheme_module)

    def get_key(self, track_leg_numbers=False):
        """Return the ProcessKey associated to this current."""

        from madgraph.core.accessors import ProcessKey
        return ProcessKey(
                self,
                allowed_attributes = [
                    'singular_structure',
                    'n_loops',
                    'squared_orders',
                    'perturbation_couplings',
                    'resolve_mother_spin_and_color',
                    'beam_PDGs',
                    'beam_type',
                    'distribution_type' ],
                track_leg_numbers = track_leg_numbers)
    
    def __str__(self, **opts):
        base_str = Current.__str__(self, **opts)
        # Denote "bulk" contributions embedded within colons ':%s:'
        if self['distribution_type']=='bulk':
            return ':%s:'%base_str
        elif self['distribution_type']=='counterterm':
            return '{%s}'%base_str
        return base_str

#=========================================================================================
# Integrated current
#=========================================================================================
class IntegratedCurrent(Current):
    """Class for integrated currents."""

    def __str__(
        self,
        print_n=True, print_pdg=True, print_state=True,
        print_loops=True, print_orders=True
    ):
        """Nice string representation of this integrated current."""

        res = super(IntegratedCurrent, self).__str__(
            print_n=print_n, print_pdg=print_pdg, print_state=print_state,
            print_loops=print_loops, print_orders=print_orders
        )
        # Denote integrated contributions embedded within '[%s]'
        return '[%s]' % res

class IntegratedBeamCurrent(IntegratedCurrent):
    """ A special class of current representing the integrated beam factorization terms 

    If we have:
        PDF_CT(xi) = F_+(xi) + [F] \delta(xi-1)
    (which can also be explicitely written)
        PDF_CT(xi) = F(xi) + :F(xi): \delta(xi-1) + [F] \delta(xi-1)

    Then this current returns the end-point of the distribution: [F] \delta(xi-1)
    """

    def default_setup(self):
        """Default values for all properties specific to the BeamCurrent,
        additional to Process.
        """
        super(IntegratedBeamCurrent, self).default_setup()
        self['n_loops'] = 0
        self['squared_orders'] = {}
        self['resolve_mother_spin_and_color'] = True
        self['singular_structure'] = SingularStructure()

        # List of the pdgs of the particles making up this beam
        self['beam_PDGs'] = []
        # Name of the beam factorization type that should be employed
        self['beam_type'] = None

    def get_key(self, track_leg_numbers=False):
        """Return the ProcessKey associated to this current."""

        from madgraph.core.accessors import ProcessKey
        return ProcessKey(
                self,
                allowed_attributes = [
                    'singular_structure',
                    'n_loops',
                    'squared_orders',
                    'perturbation_couplings',
                    'resolve_mother_spin_and_color',
                    'beam_PDGs',
                    'beam_type' ],
                track_leg_numbers = track_leg_numbers)

#=========================================================================================
# CountertermNode
#=========================================================================================

class CountertermNode(object):
    """Class representing a node in the tree of currents that make up a counterterm."""

    def __init__(self, *args, **opts):

        self.current = opts.get('current', Current())
        self.nodes = opts.get('nodes', [])
        return

    def __str__(self, lead="    ", tab="    "):

        tmp_str = lead + str(self.current) + "\n"
        for node in self.nodes:
            tmp_str += node.__str__(lead + tab, tab)
        return tmp_str

    def get_copy(self, copied_attributes=()):
        """Make sure that attributes args of the object returned
        and all its children can be modified without changing the original.
        """

        return CountertermNode(
            current=self.current.get_copy(copied_attributes),
            nodes=[node.get_copy(copied_attributes) for node in self.nodes]
        )

    def n_loops(self):
        """If number of loops are assigned,
        return the total number of loops in this node and its children.
        """

        result = self.current.get('n_loops')
        if result >= 0:
            for node in self.nodes:
                sub_n_loops = node.n_loops()
                if sub_n_loops >= 0:
                    result += sub_n_loops
                else:
                    raise MadGraph5Error(
                        "CountertermNode.n_loops: requested unassigned number of loops"
                    )
        return result

    def squared_orders(self):
        """Returns the total squared orders in the current of this node
        and its children.
        """

        result = self.current.get('squared_orders')
        for node in self.nodes:
            sub_squared_orders = node.squared_orders()
            for order, value in sub_squared_orders.items():
                try:
                    result[order] += value
                except KeyError:
                    result[order] = value
        return result

    def count_unresolved(self):
        """Count the number of unresolved particles covered by this counterterm node."""
        
        total_unresolved = self.current.count_unresolved()
        for node in self.nodes:
            total_unresolved += node.count_unresolved()
        return total_unresolved
    
    def get_number_of_couplings(self):
        """Count the number of couplings expected in the currents of this counterterm node."""
        
        total_couplings = self.current.get_number_of_couplings()
        for node in self.nodes:
            total_couplings += node.get_number_of_couplings()
        return total_couplings

    def find_leg(self, number):
        """Find the SubtractionLeg with number specified."""
        
        for subtraction_leg in self.current['singular_structure'].get_all_legs():
            if subtraction_leg.n == number:
                return subtraction_leg
        
        for node in self.nodes:
            subtraction_leg = node.find_leg(number)
            if not subtraction_leg is None:
                return subtraction_leg
        
        return None

    def get_all_currents(self):
        """Return a list of all currents involved in this node and its children."""

        currents = [self.current, ]
        for node in self.nodes:
            currents += node.get_all_currents()
        return currents

    def get_reduced_flavors(self, defining_flavors, IR_subtraction_module, routing_dict):
        """Given the defining flavors dictionary specifying the flavors of some leg numbers,
        add entries corresponding to the subtraction legs in this node and using the function
        get_parent_PDGs from the subtraction module."""

        for node in self.nodes:
            node.get_reduced_flavors(defining_flavors, IR_subtraction_module, routing_dict)  

        self.current.get_reduced_flavors(
                                    defining_flavors, IR_subtraction_module, routing_dict)

    def recursive_singular_structure(self, intermediate_leg_ns):
        """Reconstruct recursively the singular structure of this CountertermNode."""

        structure = self.current['singular_structure']
        legs = type(structure.legs)(
            leg for leg in structure.legs if leg.n not in intermediate_leg_ns
        )
        substructures = list(structure.substructures)
        substructures.extend([
            node.recursive_singular_structure(intermediate_leg_ns)
            for node in self.nodes
        ])
        return type(structure)(legs=legs, substructures=substructures)

    def distribute_loops(self, n_loops):
        """Split a counterterm in several ones according to the individual
        loop orders of its currents.
        Also, never reset the 'n_loops' attribute of a current that does not have it set to -1.
        """
        
        assert isinstance(n_loops, int)

        # Generate all combinations of nodes with total number of loops
        # less or equal to n_loops
        node_combinations = []
        # If this CountertermNode has no nodes, nothing to do
        if not self.nodes:
            node_combinations = [[], ]
        # Else construct all combinations recursively
        else:
            # Initialize node_combinations with the first node
            # at any loop number between 0 and n_loops
            first_without_loops = self.nodes[0]
            # Do not modify loop count of currents with n_loops != -1
            if self.nodes[0].current['n_loops'] != -1:
                node_combinations += [[self.nodes[0], ], ]
            else:
                for loop_n in range(n_loops + 1):
                    for first_with_loops in first_without_loops.distribute_loops(loop_n):
                        node_combinations += [[first_with_loops, ], ]

            # Add the other nodes one by one recursively,
            # taking care of never exceeding n_loops
            n_nodes = len(self.nodes)
            i_current = 1
            while i_current < n_nodes:
                new_node_combinations = []
                for combination in node_combinations:
                    combination_loops = sum(cur.n_loops() for cur in combination)
                    for new_loop_n in range(n_loops + 1 - combination_loops):
                        # Distribute the order of this node
                        # in all possible ways within the current itself
                        # (recursion step)
                        for ith_node in self.nodes[i_current].distribute_loops(new_loop_n):
                            new_node_combinations.append(combination + [ith_node, ] )
                node_combinations = new_node_combinations
                i_current += 1

        # This current should be evaluated with the number of loops
        # that is missing to reach exactly n_loops
        result = []
        for combination in node_combinations:
            combination_loops = 0
            for cur in combination:
                combination_loops += cur.current['n_loops']
            result.append(self.get_copy('n_loops'))
            result[-1].nodes = combination
            result[-1].current['n_loops'] = n_loops - combination_loops

        return result

    def distribute_orders(self, target_squared_orders):
        """Split a counterterm in several ones according to the individual
        coupling orders of its currents.
        """

        # For now this sets the orders recursively without creating copies
        # because there is only one choice
        # Eventually, it should create new instances, set orders there, and return the
        # new list of counterterm copies with orders
        # TODO actually create a new counterterm with orders set
        for node in self.nodes:
            node.distribute_orders(target_squared_orders)

        if isinstance(self.current, Current):
            self.current.set(
                'squared_orders', { 'QCD': self.current.get_number_of_couplings() * 2 } )
        else:
            # Here the squared order reduced process should be changed accordingly
            # and decreased for each squared order "sucked up" by the currents
            summed_orders = {}
            for current in self.get_all_currents():
                if isinstance(self.current, Current):
                    for key, value in current.get('squared_orders').items():
                        try:
                            summed_orders[key] += value
                        except KeyError:
                            summed_orders[key] = value

            sqo = self.current.get('squared_orders')
            for key, value in summed_orders.items():
                try:
                    sqo[key] -= value
                except KeyError:
                    raise MadGraph5Error(
                        "Subtraction currents have squared orders "
                        "absent from real-emission ME."
                    )

        return [self, ]

#=========================================================================================
# Counterterm
#=========================================================================================

class Counterterm(CountertermNode):
    """Class representing a tree of currents multiplying a matrix element."""

    def __init__(self, **opts):
        """Initialize a Counterterm using a dictionary. Valid entries are:

        process: The process associated to the reduced matrix element in the counterterm.
        current: An alternative name for process. If both are present, process is used.
        nodes: The nodes that represent the tree of currents of this counterterm.
        momenta_dict: A dictionary that specifies numbers for intermediate legs.
        prefactor: The global combinatorial prefactor carried by the counterterm.

        If no combinatorial prefactor is passed, Counterterm will attempt to compute it,
        forwarding options to the function get_counterterm() to speed up this calculation.
        See the documentation of that function for further information.
        """

        new_opts = {key: opts[key] for key in opts.keys() if key in ('current', 'nodes')}
        if opts.has_key('process'):
            new_opts['current'] = opts['process']
        elif not opts.has_key('current'):
            new_opts['current'] = base_objects.Process()
        super(Counterterm, self).__init__(**new_opts)
        self.momenta_dict = opts.get('momenta_dict', bidict())
        try:
            self.prefactor = opts['prefactor']
        except KeyError:
            self.prefactor = self.get_prefactor(**opts)

        # Both attributes below will be set via a call to `set_reduced_flavors_map`
        # during the generation of the counterterms
        # The dictionary below takes the form:
        #   { ( (reolved_initial_state_PDGs,), (resolved_final_state_PDGs, ) ) :
        #               ( (resolved_initial_state_PDGs,), (resolved_final_state_PDGs, ) ) }
        # In the future, the values can be upgraded to list of mapped reduced flavors.
        # But for now there can only be one at a time.
        # Also, an additional 'None' key is added which points to the reduced flavor
        # configuration of the defining process.
        self.reduced_flavors_map     = None

        # Local counterterms never require corelated beam convolutions
        self.requires_correlated_beam_convolution = False
        # Neither to they require non-factorisable ones
        self.n_non_factorisable_double_sided_convolution = 0
        #gl Local counterterms never require torino_sub_BS
        self.requires_torino_sub_BS = False

    @classmethod
    def get_ancestor(cls, particles, momentum_dict, stop_if_none_found=True):
        """Recursively explore the momentum dictionary to find an ancestor
        of the given particle set (i.e. these are leg numbers here).
        """

        try:
            return momentum_dict.inv[particles]
        except KeyError:
            # Find all keys that allow to to reduce the set of particles to find the ancestor for
            selected_keys = []
            for key in momentum_dict.inv.keys():
                if len(key) == 1 or key.isdisjoint(particles):
                    continue
                if key.issubset(particles):
                    selected_keys.append(key)

            for key in selected_keys:
                particles_for_next_step = particles.difference(key)
                particles_for_next_step = particles_for_next_step.union({momentum_dict.inv[key], })

                if len(particles_for_next_step) == 1:
                    if particles_for_next_step in momentum_dict.inv.keys():
                        return momentum_dict.inv[particles_for_next_step]
                    else:
                        return tuple(particles_for_next_step)[0]

                elif particles_for_next_step != particles:
                    outcome_for_this_branch = cls.get_ancestor(particles_for_next_step, momentum_dict, stop_if_none_found=False)
                    if outcome_for_this_branch is not None:
                        return outcome_for_this_branch

            if stop_if_none_found:
                raise KeyError(
                    "Could not find leg numbers " + str(particles) +
                    "in this momentum routing dictionary:\n" + str(momentum_dict) )
            else:
                return None

    def __str__(self, print_n=True, print_pdg=False, print_state=False, print_n_loops=True):

        suffix = ''
        n_loops_in_kernel = self.n_loops() - self.current.get('n_loops')
        n_loops_in_process = self.current.get('n_loops')
        if print_n_loops:
            if n_loops_in_kernel > 0 or n_loops_in_process > 0:
                suffix += ' @ %d(kernel)+%d(ME) loop%s' % (
                    n_loops_in_kernel, n_loops_in_process, 's' if n_loops_in_kernel + n_loops_in_process > 1 else '')

        return self.reconstruct_complete_singular_structure().__str__(
            print_n=print_n, print_pdg=print_pdg, print_state=print_state )+suffix

    def nice_string(self, lead="    ", tab="    "):
        
        CT_type = ''
        if type(self) == Counterterm:
            CT_type = '[local] '
        elif type(self) == IntegratedCounterterm:
            CT_type = '[integrated] '

        tmp_str  = lead + CT_type + self.process.nice_string(0, True, False)
        # From now on, only keep the : in the lead prefix
        lead = ':'.join(' '*len(s) for s in lead.split(':'))
        tmp_str += " ("
        tmp_str += " ".join(
            str(leg['number'])
            for leg in self.process['legs']
            if leg['state'] == SubtractionLeg.INITIAL
        )
        tmp_str += " > "
        tmp_str += " ".join(
            str(leg['number'])
            for leg in self.process['legs']
            if leg['state'] == SubtractionLeg.FINAL
        )
        tmp_str += ")"
        process_n_loops = self.current.get('n_loops')
        if process_n_loops >= 0:
            tmp_str += " @ " + str(process_n_loops) + " loop"
            if process_n_loops != 1:
                tmp_str += "s"
        tmp_str += "\n"
        if self.requires_correlated_beam_convolution is not None:
            if self.requires_correlated_beam_convolution:
                n_non_factorisable_double_sided_convolution = self.get_n_non_factorisable_double_sided_convolution()
                if n_non_factorisable_double_sided_convolution:
                    tmp_str += lead + "( supported by a non-factorisable double-sided beam convolution "+\
                                                    "(%d convolutions) )\n"%n_non_factorisable_double_sided_convolution
                else:
                    tmp_str += lead + "( supported by a correlated beam convolution )\n"
        for node in self.nodes:
            tmp_str += node.__str__(lead + tab, tab)
        tmp_str += lead + "Pseudoparticles: {"
        tmp_str += "; ".join(
            str(key) + ": (" +
            ",".join(str(n) for n in self.momenta_dict[key]) + ")"
            for key in self.momenta_dict
        )
        tmp_str += "}"

        return tmp_str

    def get_copy(self, copied_attributes=()):
        """Make sure that attributes args of the object returned
        and all its children can be modified without changing the original.
        """

        node = super(Counterterm, self).get_copy(copied_attributes)
        momenta_dict = self.momenta_dict
        if 'momenta_dict' in copied_attributes:
            momenta_dict = bidict(momenta_dict)
        CT_copy = type(self)(
            process=node.current, nodes=node.nodes,
            momenta_dict=momenta_dict, prefactor=self.prefactor
        )
        # Also forward the value of 'requires_correlated_beam_convolution' if it
        # was already set.
        CT_copy.requires_correlated_beam_convolution = self.requires_correlated_beam_convolution
        CT_copy.n_non_factorisable_double_sided_convolution = self.n_non_factorisable_double_sided_convolution
        #gl
        CT_copy.requires_torino_sub_BS = self.requires_torino_sub_BS

        return CT_copy

    @property
    def process(self):
        return self.current

    def find_leg(self, number):
        """Find the Leg or SubtractionLeg with number specified."""
        
        for leg in self.process.get('legs'):
            if leg.get('number') == number:
                return leg
        
        for node in self.nodes:
            subtraction_leg = node.find_leg(number)
            if not subtraction_leg is None:
                return subtraction_leg
        
        return None

    @classmethod
    def get_daughter_legs(cls, leg_number, momentum_dict):
        """Walk down the momentum dictionary tree to find the outermost daughter legs
        coming from the leg_number specified."""
        
        external_leg_numbers = []
        intermediate_leg_number = [leg_number, ]
        while len(intermediate_leg_number) > 0:
            next_leg = intermediate_leg_number.pop(0)
            # Check if this leg is external
            daughters = momentum_dict[next_leg]
            if daughters == frozenset([next_leg,]):
                external_leg_numbers.append(next_leg)
            else:
                intermediate_leg_number = list(daughters) + intermediate_leg_number
        return external_leg_numbers

    def get_daughter_pdgs(self, leg_number, state):
        """Walk down the tree of currents to find the pdgs 'attached'
        to a given leg number of the reduced process.
        """
        
        external_leg_numbers = self.get_daughter_legs(leg_number, self.momenta_dict)
        
        # Now that we have the external leg numbers, find the corresponding PDG
        pdgs = []
        for n in external_leg_numbers:
            leg = self.find_leg(n)
            if leg is None:
                raise MadGraph5Error(
                    "Could not find leg number %d in counterterm:\n%s"%(n,str(self))
                )
            if isinstance(leg, SubtractionLeg):
                if leg.state == state:
                    pdgs.append(leg.pdg)
            else:
                if leg['state'] == state:
                    pdgs.append(leg['id'])
        return pdgs

    def get_resolved_process_pdgs(self):
        """ Walks through the currents and obtain the pdgs list of resolved process."""

        reduced_initial_leg_numbers = [
            l.get('number') for l in self.process.get_initial_legs()
        ]
        all_initial_leg_pdgs = []
        for leg_number in reduced_initial_leg_numbers:
            all_initial_leg_pdgs.extend(
                self.get_daughter_pdgs(leg_number, SubtractionLeg.INITIAL)
            )

        reduced_final_leg_numbers = [ l.get('number') for l in self.process.get_final_legs() ]
        all_final_leg_pdgs= []
        for leg_number in reduced_final_leg_numbers:
            all_final_leg_pdgs.extend(
                self.get_daughter_pdgs(leg_number, SubtractionLeg.FINAL)
            )
        
        return ( tuple(all_initial_leg_pdgs), tuple(all_final_leg_pdgs) )

    def set_requires_correlated_beam_convolution(self, subtraction_scheme_module):
        """ Using the subtraction scheme module, this function sets whether this counterterm requires
        a correlated beam convolution or not and stores the result in the variable requires_correlated_beam_convolution.
        This information matters as it sets whether integrated counterterm requires a host contribution featuring
        correlated convolution of the beams (i.e. 'BS', 'VS', etc... contributions).
        This is typically the case for integrated counterterms with pure-soft *BeamCurrent* that originated
        from a colorful subtraction currents scheme with mappings such as ppToOneWalker that
        recoils the soft momentum equally against both initial-state beams. """

        # Fetch all integrated currents making up this integrated counterterm
        # In the current implementation, there can only be one for now.
        self.requires_correlated_beam_convolution = \
            self.get_integrated_current().does_require_correlated_beam_convolution(subtraction_scheme_module)

    def does_require_correlated_beam_convolution(self):
        """ Protected access to the variable requires_correlated_beam_convolution which must first be initialised."""

        if self.requires_correlated_beam_convolution is None:
            raise MadGraph5Error("The attribute 'requires_correlated_beam_convolution' of the counterterm:\n%s\n"%str(self)+
                "must be set using the IR subtraction scheme module with the function 'set_requires_correlated_beam_convolution'"+
                " before this information can be accessed.")

        return self.requires_correlated_beam_convolution

#gl
    def set_requires_torino_sub_BS(self, subtraction_scheme_module):

        # Fetch all integrated currents making up this integrated counterterm
        # In the current implementation, there can only be one for now.
        self.requires_torino_sub_BS = \
            self.get_integrated_current().does_require_torino_sub_BS(subtraction_scheme_module)

    def does_require_torino_sub_BS(self):
        """ Protected access to the variable requires_correlated_beam_convolution which must first be initialised."""

        if self.requires_torino_sub_BS is None:
            raise MadGraph5Error("The attribute 'requires_torino_sub_BS' of the counterterm:\n%s\n"%str(self)+
                "must be set using the IR subtraction scheme module with the function 'set_requires_torino_sub_BS'"+
                " before this information can be accessed.")

        return self.requires_torino_sub_BS

    def set_n_non_factorisable_double_sided_convolution(self, subtraction_scheme_module):
        """ Test if this counterterm involves a non-factorisable double-sided convolution which appears when one has
        a combination of a correlated convolution with single-sided ones (like for IF integrated collinear CT or PDF
        counterterms). For now it can only return 2 convolutions if non-factorisable; it may need to be extended for N^3LO."""

        # Fetch all integrated currents making up this integrated counterterm
        # In the current implementation, there can only be one for now.
        self.n_non_factorisable_double_sided_convolution = \
                            subtraction_scheme_module.exporter.get_n_non_factorisable_double_sided_convolution(self)

    def get_n_non_factorisable_double_sided_convolution(self):
        """ Counts how many non-factorisable beam convolutions this singular structure.
        At NNLO this can only be 0 or 2. For example ((C(2,4),S(5),),) will yield 2 but
        ((C(2,4),C(1,3)),) or ((S(2),C(4,5)),) would yield 0."""

        if self.n_non_factorisable_double_sided_convolution is None:
            raise MadGraph5Error("The attribute 'n_non_factorisable_double_sided_convolution' of the counterterm:\n%s\n"%str(self)+
                "must be set using the IR subtraction scheme module with the function 'set_n_non_factorisable_double_sided_convolution'"+
                " before this information can be accessed.")

        return self.n_non_factorisable_double_sided_convolution

    def reconstruct_complete_singular_structure(self):
        """Reconstruct the complete singular structure for this counterterm."""

        intermediate_leg_ns = frozenset(
            self.momenta_dict.inv[key] for key in self.momenta_dict.inv.keys()
            if len(key) > 1 )
        if len(self.nodes) > 0:
            return SingularStructure(
                substructures=[
                node.recursive_singular_structure(intermediate_leg_ns)
                for node in self.nodes ] )
        else:
            return SingularStructure()

    def get_prefactor(self, **opts):
        """Determine the overall prefactor of the counterterm
        associated with this singular structure.
        """
        pref = 1.

        # First get the subtraction prefactor.        
        complete_singular_structure = opts.get(
            'complete_singular_structure',
            self.reconstruct_complete_singular_structure()
        )
        
        pref *= complete_singular_structure.get_subtraction_prefactor()

        # Remove the final state symmetry factor from the reduced process,
        pref *= self.process.identical_particle_factor()

        # And enforce the final state symmetry factor of the resolve process.
        # For instance, for the limit C(5,6) of e+ e- > d d~ d d~
        # the net result is that we divide the counterterm by 4
        # because the reduced process e+ e- > d d~ g does
        # not have any identical particles in the final state.
        try:
            resolved_sym_factor = opts['resolved_process'].identical_particle_factor()
        except KeyError:
            # If the resolved process is not provided, we can reconstruct the information
            # by walking recursively through the nodes.
            final_pdgs = self.get_resolved_process_pdgs()[1]
            # If the resolved_process is provided, then construct the prefactor from it,
            # otherwise reconstruct it.
            resolved_sym_factor = 1.
            for final_pdg in set(final_pdgs):
                resolved_sym_factor *= final_pdgs.count(final_pdg)

        pref /= resolved_sym_factor

        return pref

    def is_singular(self):
        """Return whether this counterterm is just the implementation of the pure
        matrix element or if it has singular region (and the corresponding currents)
        attached to it.
        """

        return len(self.nodes) > 0
        
    def count_unresolved(self):
        """Count the number of unresolved particles covered by this counterterm."""
        
        total_unresolved = 0
        for node in self.nodes:
            total_unresolved += node.count_unresolved()
        return total_unresolved

    def get_number_of_couplings(self):
        """Count the number of couplings expected in the currents of this counterterm."""
        
        total_couplings = 0
        for node in self.nodes:
            total_couplings += node.get_number_of_couplings()
        return total_couplings

    def get_all_currents(self):
        """Return a list of all currents involved in this counterterm."""

        currents = []
        for node in self.nodes:
            currents += node.get_all_currents()
        return currents

    def get_currents_blocks(self):
        """ Return the list of (possible combination of) currents that are to be called at once when evaluating this
         counterterm. For example C(F,F) would yield just a list one element being a tuple of length one containing
         the current C(F,F) only. So -> [ ( C(F,F), ) ]
         Then (S(F), F1) in colorful would still return a list of a single element but this time being a tuple of length
         two containing the two currents S(F) and F1. So -> [ ( S(F), F1 ) ]
         Finally something the factorisable counterterm (F1, F2) would yield a list of two elements being the two
         length-one tuples [ ( F1, ), ( F2, ) ].
         """

        currents_blocks = []
        if self.get_n_non_factorisable_double_sided_convolution() >= 2:
            # In that case make sure that there is no nesting within this counterterm as this is not supported
            # for non-factorisable currents
            assert all(len(node.nodes)==0 for node in self.nodes)
            # This counterterm is non-factorisable, so we must aggregated the node into a single entity to export
            currents_blocks.append( CurrentsBlock(node.current for node in self.nodes) )
        else:
            currents_blocks.extend([ CurrentsBlock([current,]) for current in self.get_all_currents() ])

        return currents_blocks

    def get_beam_currents(self):
        """ Returns a list of dictionaries of the form
              {   'beam_one' : beam_current_1 # None for no convolution
                  'beam_two' : beam_current_2 # None for no convolution
              } """
        
        # Copy the beam currents from the reduced process to avoid border effects
        # when further post-processing it below.
        beam_currents = [dict(bf) for bf in self.current['beam_factorization']]

        # First fetch all currents making up this counterterm
        for current in self.get_all_currents():        
            if type(current) not in [BeamCurrent,]:
                continue

            # Now get all legs of its singular structure
            all_legs = current['singular_structure'].get_all_legs()
    
            if any(l.n==1 for l in all_legs):
                for bcs in beam_currents:
                    if bcs['beam_one'] is not None:
                        raise MadGraph5Error('The beam factorization currents from the reduced'+
                ' process and the ones in the counterterms must never apply to the same beam (#1).')
                    bcs['beam_one'] = current
            if any(l.n==2 for l in all_legs):
                for bcs in beam_currents:
                    if bcs['beam_two'] is not None:
                        raise MadGraph5Error('The beam factorization currents from the reduced'+
                ' process and the ones in the counterterms must never apply to the same beam (#2).')
                    bcs['beam_two'] = current
    
        return beam_currents

    def set_reduced_flavors_map(self, defining_process, mapped_processes, IR_subtraction):
        """Given the defining and mapped processes subject to this counterterm, 
        this function sets the instance attribute 'reduced_flavor_map' which can be 
        used later to easily access the corresponding reduced flavor."""

        if self.count_unresolved() > 1 and any(
            any(substruct.name()=='S' for substruct in current['singular_structure'].substructures)
                for current in self.get_all_currents() if isinstance(current,(
                                    IntegratedBeamCurrent,IntegratedCurrent,BeamCurrent))):
            raise MadGraph5Error("The function 'get_reduced_flavors' of the class Counterterm"+
                " does not yet support integrated NNLO soft structures such as: %s"%str(self))

        self.reduced_flavors_map = {}
        for i_proc, process in enumerate([defining_process,]+mapped_processes):
            number_to_flavors_map = { l['number'] : l['id'] for l in process.get('legs') }
            # Update the number_to_flavors_map dictionary by walking through the nodes
            for node in self.nodes:
                node.get_reduced_flavors(
                                 number_to_flavors_map, IR_subtraction, self.momenta_dict)
            reduced_process_leg_numbers = self.process.get_cached_initial_final_numbers()

            reduced_flavors = (
                tuple( number_to_flavors_map[n] for n in reduced_process_leg_numbers[0] ),
                tuple( number_to_flavors_map[n] for n in reduced_process_leg_numbers[1] )
            )

            self.reduced_flavors_map[process.get_cached_initial_final_pdgs()] = reduced_flavors
            # Also add a None key for the default reduced flavors corresponding to the 
            # resolved ones of the defining process.
            if i_proc==0:
                self.reduced_flavors_map[None] = reduced_flavors                

    def get_reduced_flavors(self, defining_flavors=None):
        """Given the defining flavors as a tuple (initial_state_PDGs, final_state_PDGs),
         returns in the same format the flavors corresponding to the flavor assignment of
         the reduced process  given the defining flavors.
         If the defining defining_flavors is set to None, then the default reduced_flavors
         will be returned, corresponding to the defining process that yielded this counterterm.
        """

        if self.reduced_flavors_map is None:
            raise MadGraph5Error("The function 'set_reduced_flavors_map' of the counterterm"+
                            " should be called before the function 'get_reduced_flavors'.")

        try:
            return self.reduced_flavors_map[defining_flavors]
        except KeyError:
            raise MadGraph5Error(
                "The reduced flavor configuration for the resolved configuration "+
                "'%s' was not found in the"%str(defining_flavors)+
                " reduced flavors map:\n%s\nof this counterterm:\n%s"%(
                                         str(self.reduced_flavors_map),self.nice_string()))

    def get_reduced_quantities(self, reduced_PS_dict, defining_flavors=None):
        """Given the PS *dictionary* providing the reduced kinematics returned
        by the mapping and the selected defining flavors of the resolved process,
        return a *list* of momenta corresponding to the reduced kinematics and
        a *list* of flavors corresponding to the flavor assignment to the reduced
        process."""
        
        from madgraph.integrator.phase_space_generators import LorentzVectorList
        # Here we want a reduced *List* of LorentzVectors, so as to be used directly
        # in the pass_cuts and observable functions. This is at variance with the 
        # integrated counterterm reduced kinematics which should be a *dictionary* with
        # keys aligned with the leg number of the reduced process and subtraction legs,
        # so that the integrated current evaluations can be performed.

        # First construct the reduced kinematic list
        reduced_PS = self.process.format_PS_point_for_ME_call(reduced_PS_dict)

        reduced_flavors = self.get_reduced_flavors(defining_flavors)
        
        return reduced_PS, reduced_flavors

    def distribute_loops(self, n_loops):
        """Split a counterterm in several ones according to the individual
        loop orders of its currents.
        Also, never reset the 'n_loops' attribute of a current that does not have it set to -1.
        """

        result = super(Counterterm,self).distribute_loops(n_loops)

        # Set the 'NLO_mode' and 'perturbation_couplings' to 'tree' and '[]' if n_loops=0
        for CT in result:
            if CT.process.get('n_loops') == 0:
                CT.process.set('NLO_mode', 'tree')
                CT.process.set('perturbation_couplings', [])

        return result

#=========================================================================================
# IntegratedCounterterm
#=========================================================================================

class IntegratedCounterterm(Counterterm):
    """A class for the integrated counterterm."""

    def __init__(self, *args, **opts):
        """ Similar instantiation as a normal counterterm except that one cannot know a-priori if
        this integrated counterterm will require correlated beam convolutions, as the specification
        of the subtraction scheme is necessary for this."""

        subtraction_scheme_module = opts.pop('subtraction_scheme_module', None)
        super(IntegratedCounterterm, self).__init__(*args, **opts)

        if subtraction_scheme_module is None:
            # Then this information cannot be constructed a priori and it will need to be set
            # later by calling the function set_
            self.requires_correlated_beam_convolution = None
            # Similarly for n_non_factorisable_double_sided_convolution
            self.n_non_factorisable_double_sided_convolution = None
            #gl
            self.requires_torino_sub_BS = None
        else:
            self.set_requires_correlated_beam_convolution(subtraction_scheme_module)
            self.set_n_non_factorisable_double_sided_convolution(subtraction_scheme_module)
            #gl
            self.set_requires_torino_sub_BS(subtraction_scheme_module)

    def __str__(self, print_n=True, print_pdg=False, print_state=False, print_n_loops = True):
        """ A nice short string for the structure of this integrated CT."""

        suffix = ''
        n_loops_in_kernel = self.n_loops() - self.current.get('n_loops')
        n_loops_in_process = self.current.get('n_loops')
        if print_n_loops:
            if n_loops_in_kernel > 0 or  n_loops_in_process > 0:
                suffix += ' @ %d(kernel)+%d(ME) loop%s' % (
                    n_loops_in_kernel, n_loops_in_process, 's' if n_loops_in_kernel + n_loops_in_process > 1 else '')

        reconstructed_ss = self.reconstruct_complete_singular_structure()
        
        # We deconstruct the reconstructed ss so as to be able to render the information about
        # wheter beam currents are local counterterms or bulk ones.
        res = []
        for i, node in enumerate(self.nodes):
            if isinstance(node.current, BeamCurrent) and node.current['distribution_type']=='bulk':
                template = ':%s:,'
            else:
                template = '%s,'
            res.append(template%(reconstructed_ss.substructures[i].__str__(
                                 print_n=print_n, print_pdg=print_pdg, print_state=print_state ))+suffix)

        return '[%s]'%(''.join(res))

    def get_integrated_current(self):
        """ This function only makes sense in the current context where we don't allow
        the factorization of the integrated counterterm into several integrated currents."""
        
        # The current way we use the IntegratedCounterterm structure is as follows:
        # a) It should never a nested death deeper than 1
        # b) The first element of the node corresponds to the genuine integrated counterterm, e.g. [C(1,2)]
        # c) If there are any, the subsequent elements in self.nodes correspond to PDF beam factorisation currents
        #    that have been merged into this counterterm because the convolution must be done analytically.
        if len(self.nodes) >= 1 and all(len(node.nodes) == 0 for node in self.nodes):
            assert isinstance(self.nodes[0].current, (IntegratedCurrent, BeamCurrent, IntegratedBeamCurrent))
            return self.nodes[0].current
        else:
            raise MadGraph5Error(
                """For now, MadNkLO only support simple integrated
                counterterms that consists of single current encompassing the full
                singular structure that must be analytically integrated over, not:\n%s"""%str(self))

    def merge_convolutions(self, subtraction_scheme_module):
        """ If this integrated counterterm implies multiple convolutions in the *same* convolution parameter
        then one of them will be performed analytically, thus inducing a new building-block counterterm implementing
        that composite object. Example:
            :[(C(1,4),)]:  with reduced process B x :F(1):
        Cannot be implemented as direct product of convolutions, instead it must correspond to novel composite object
        defined as follows:
            :[(C(1,4),F(1))]: with reduced process B x 1 (no convolution attached to this reduced process).

        """

        # First fetch what are the active beam currents of this counterterm
        integrated_current = self.get_integrated_current()

        if type(integrated_current) not in [BeamCurrent, IntegratedBeamCurrent]:
            return [self,]

        if all( (bf['beam_one'] is None and bf['beam_two'] is None) for bf in self.current['beam_factorization'] ):
            return [self,]

        # We can afford to hard-code that initial leg 1 corresponds to the first beam
        # while leg #2 would correspond to the second beam
        all_legs = integrated_current['singular_structure'].get_all_legs()
        beams_convolved_by_integrated_current = [
            beam_name for i_beam, beam_name in enumerate(['beam_one', 'beam_two']) if (
                integrated_current['singular_structure'].does_require_correlated_beam_convolution(subtraction_scheme_module) or
                self.get_n_non_factorisable_double_sided_convolution() or
                any(l.n == i_beam+1 for l in all_legs)
            )
        ]

        # Test if merging will be necessary
        if all( all(bf[beam_name] is None for bf in self.current['beam_factorization'])
                for beam_name in beams_convolved_by_integrated_current ):
            return [self,]

        merged_counterterms = []

        for bf in self.current['beam_factorization']:
            modified_bf = dict(bf)
            # The get_copy() below is very important and should make sure that the modification
            # of this copied counterterm will not have border effects.
            new_integrated_CT = self.get_copy(copied_attributes=tuple(["singular_structure",]))
            for beam_name in beams_convolved_by_integrated_current:
                if modified_bf[beam_name] is None:
                    continue
                # We must then remove this beam factorisation term from the reduced process and
                # Add it to the integrated current instead.
                new_integrated_CT.nodes.append(CountertermNode(current=bf[beam_name].get_copy()))
                # (new_integrated_CT.get_integrated_current()['singular_structure']).substructures.append(
                #                                                       modified_bf[beam_name].get('singular_structure'))
                modified_bf[beam_name] = None

            # Overwrite the beam factorization of the reduced process by the modified one where the
            # component merged with the integrated CT has been removed.
            new_integrated_CT.current.set('beam_factorization', [modified_bf,])

            merged_counterterms.append(new_integrated_CT)

        return merged_counterterms

    def get_beam_currents(self):
        """ Returns a list of dictionaries of the form
              {   'beam_one' : beam_current_1 # None for no convolution
                  'beam_two' : beam_current_2 # None for no convolution
      <optional>  'correlated_convolution' : beam_current # None for not correlated convolution
      <optional>  'non_factorisable_convolution' : <beam_current> # None for not correlated convolution
              }
           Note that whenever 'correlated_convolution' is not set to None, the other beam_currents
           will be None.
        """

        beam_currents = [dict(bf) for bf in self.current['beam_factorization']]

        # First fetch all integrated currents making up this integrated counterterm
        # In the current implementation, there can only be one for now.
        integrated_current = self.get_integrated_current()

        if type(integrated_current) not in [BeamCurrent, IntegratedBeamCurrent]:
            return beam_currents

        if self.does_require_correlated_beam_convolution():
            for bcs in beam_currents:
                if bcs['beam_one'] is not None or bcs['beam_two'] is not None:
                    raise MadGraph5Error('The beam factorization currents from the reduced'+
        ' process should all be None if the integrated counterterm necessitates a correlated'+
        ' beam convolution.')
                bcs['correlated_convolution'] = integrated_current
                bcs['torino_sub_BS'] = integrated_current   #gl
            return beam_currents

        if self.get_n_non_factorisable_double_sided_convolution()>=2:
            for bcs in beam_currents:
                if bcs['beam_one'] is not None or bcs['beam_two'] is not None:
                    raise MadGraph5Error('The beam factorization currents from the reduced' +
                                         ' process should all be None if the integrated counterterm necessitates a non-factorizable' +
                                         ' beam convolution.')
                bcs['non_factorisable_convolution'] = integrated_current
            return beam_currents

        # Now get all legs of its singular structure
        all_legs = integrated_current['singular_structure'].get_all_legs()

        if any(l.n==1 for l in all_legs):
            for bcs in beam_currents:
                if bcs['beam_one'] is not None:
                    raise MadGraph5Error('The beam factorization currents from the reduced'+
            ' process and the ones in the integrated CT must never apply to the same beam (#1).')
                bcs['beam_one'] = integrated_current
        elif any(l.n==2 for l in all_legs):
            for bcs in beam_currents:
                if bcs['beam_two'] is not None:
                    raise MadGraph5Error('The beam factorization currents from the reduced'+
            ' process and the ones in the integrated CT must never apply to the same beam (#2).')
                bcs['beam_two'] = integrated_current
        else:
            # If it does not contain any initial state then it must be the endpoint of an IntegratedBeamCounterterm
            # that has only final states but recoils against initial states in that particular subtraction scheme
            # The endpoint part of the distribution C(F,F) in colorful_pp would be an example of such integratedBeamCounterm
            # and we classify it here as a 'correlated_convolution' even though the rescaling xi variables supplied
            # would anyway be None. Note that we could also have decided to classify it as a 'non_factorisable_convolution',
            # it would have made no difference.
            for bcs in beam_currents:
                bcs['correlated_convolution'] = integrated_current
                bcs['torino_sub_BS'] = integrated_current   #gl

        return beam_currents

    def get_necessary_beam_convolutions(self):
        """ Returns a set of beam names ('beam_one' or 'beam_two') that must be active
        in the contribution that will host this counterterm"""
        
        necessary_beam_convolutions = set([])
        
        # First analyze the reduced process
        for bft in self.process['beam_factorization']:
            for beam_name, beam_current in bft.items():
                if not beam_current is None:
                    necessary_beam_convolutions.add(beam_name)

        # Then fetch all integrated currents making up this integrated counterterm
        # In the current implementation, there can only be one for now.
        integrated_current = self.get_integrated_current()
        if type(integrated_current) not in [ BeamCurrent, ]:
            return necessary_beam_convolutions
        
        # Check if this integrated counterterm necessitate a soft-recoil that demands
        # a correlated convolution of both beams
        if self.does_require_correlated_beam_convolution() or self.get_n_non_factorisable_double_sided_convolution()>=2:
            #print('subtraction.py -  does_require_correlated_beam_convolution : ' + str(self.does_require_correlated_beam_convolution()))
            #print('subtraction.py -  does_require_torino_sub_BS : ' + str(self.does_require_torino_sub_BS()))
            return set(['beam_one','beam_two'])
        
        # Now get all legs of its singular structure
        all_legs = integrated_current['singular_structure'].get_all_legs()
        # print('all_legs : ' + str(all_legs))

        if any(l.n==1 for l in all_legs):
            necessary_beam_convolutions.add('beam_one')
        if any(l.n==2 for l in all_legs):
            necessary_beam_convolutions.add('beam_two')
        # print('necessary_beam_convolutions : ' + str(necessary_beam_convolutions))
        return necessary_beam_convolutions

    def get_reduced_kinematics(self, input_reduced_PS):
        """Given the PS specifying the reduced kinematics corresponding to
        the current PS point thrown at the virtual contribution, return a *dictionary* of 
        momenta corresponding to the reduced kinematics."""

        from madgraph.integrator.phase_space_generators import LorentzVectorDict
        # Here we want a reduced *dictionary* of LorentzVectors with keys aligned with
        # the leg number of the reduced process and subtraction legs, so that the 
        # integrated current evaluations can be performed.
        # This is at variance with the local counterterms where a *list* of momenta
        # should be passed for the reduced kinematics so as to be used
        # in pass_cuts and the observable functions.
        reduced_PS = LorentzVectorDict()

        # First construct the reduced kinematic list, the leg numbers will not match the
        # keys of the reduced_PS_dict because in this case they will be consecutive since
        # they come from the probing of the virtual and not from any kind of mapping.

        leg_numbers = (n for sublist in self.process.get_cached_initial_final_numbers() 
                                                                          for n in sublist)                       
        if isinstance(input_reduced_PS, dict):
            for i, number in enumerate(leg_numbers):
                assert( i+1 in input_reduced_PS )
                reduced_PS[number] = input_reduced_PS[i+1]
        else:
            assert (len(input_reduced_PS) == len(leg_numbers))
            for i, number in enumerate(leg_numbers):
                reduced_PS[number] = input_reduced_PS[i]

        return reduced_PS


    def n_loops_in_host_contribution(self, host_will_have_correlated_beam_convolution, beams_with_convolutions):
        """ Returns the number of loops that the contribution hosting this integrated counterterm
        is supposed to have. Factorization contributions render this quantity more subtle than
        the naive n_loops() call to self.
        host_will_have_correlated_beam_convolution indicates if this counterterm will be hosted
        in a contribution of type BS, VS, etc...
        """

        # First account for all the loops of the building blocks
        host_contrib_n_loops  = self.n_loops()
        host_contrib_n_loops += self.current.get_n_loops_in_beam_factorization()

        # We must then add one loop for each unresolved parton of this integrated CT.
        host_contrib_n_loops += self.count_unresolved()

        # Similarly, the beamCurrent F terms have zero number of unresolved legs
        # but they must be placed together with the virtuals which have a loop count. We must
        # therefore increase the loop count by *at most* one for each of the two beams that has
        # at least one integratedBeamCurrent F structure attached
        beam_numbers_for_which_n_loops_is_already_incremented = set([])
        beam_currents = [ self.current['beam_factorization'][0][beam_name] for beam_name in ['beam_one', 'beam_two'] if
                          (self.current['beam_factorization'][0][beam_name] is not None) ]
        for current in (self.get_all_currents() + beam_currents):
            # We don't need to look recursively inside the singular structures since the beam
            # factorization ones are supposed to be at the top level since they factorize
            decomposed_singular_structure = current['singular_structure'].decompose()
            beam_numbers_in_currents = set(
                sum([[l.n for l in ss.legs] for ss in decomposed_singular_structure if ss.name() == 'F'], []))
            for beam_number in beam_numbers_in_currents:
                assert (beam_number in [1, 2]), "Inconsistent beam number found."
                if beam_number not in beam_numbers_for_which_n_loops_is_already_incremented:
                    beam_numbers_for_which_n_loops_is_already_incremented.add(beam_number)
                    host_contrib_n_loops += 1

        # Finally we should remove one loop count for *S contributions and 2 for double-sided non-correlated
        # beam convolutions.
        if host_will_have_correlated_beam_convolution:
            return host_contrib_n_loops-1
        else:
            return host_contrib_n_loops-len(beams_with_convolutions)

#=========================================================================================
# IRSubtraction
#=========================================================================================

class IRSubtraction(object):

    _allowed_model_names = ['sm', 'loop_sm', 'loopsm',
                            'simple_qcd','loop_qcd_qed_sm','hc_nlo_x0_ufo']

    def __init__(
        self, model,
        coupling_types=('QCD', ), beam_types=(None, None), subtraction_scheme=None):
        """Initialize a IR subtractions for a given model,ss
        correction order and type.
        """
        self.model          = model
        self.coupling_types = coupling_types
        self.beam_types     = beam_types

        if subtraction_scheme is None:
            self.subtraction_scheme_module = None
            self.soft_do_recoil_against_initial_states = False
        else:
            self.subtraction_scheme_module = SubtractionCurrentExporter.get_subtraction_scheme_module(subtraction_scheme)
            # Decide is soft recoil against initial states
            if self.subtraction_scheme_module.requires_soft_beam_factorization:
                self.soft_do_recoil_against_initial_states = True
            else:
                self.soft_do_recoil_against_initial_states = False
            
        # Map perturbed coupling orders to the corresponding interactions and particles.
        # The entries of the dictionary are
        # 'interactions', 'pert_particles' and 'soft_particles'
        self.IR_quantities_for_corrections_types = dict(
            (
                coupling,
                fks_common.find_pert_particles_interactions(
                    self.model,
                    pert_order=coupling
                )
            )
            for coupling in coupling_types
        )

    # =========================================================================================
    # Create and update momenta dictionaries
    # =========================================================================================
    @staticmethod
    def create_momenta_dict(process):
        """Create a new momenta dictionary for a given process."""

        momenta_dict = bidict()
        for leg in process['legs']:
            momenta_dict[leg['number']] = frozenset((leg['number'],))
        return momenta_dict

    @staticmethod
    def new_leg_number(momenta_dict):
        """Get a leg number that does not appear in momenta_dict."""

        all_leg_numbers = frozenset.union(
            frozenset.union(*momenta_dict.values()),
            frozenset(momenta_dict.keys()))
        return max(*all_leg_numbers) + 1

    @staticmethod
    def update_momenta_dict(momenta_dict, singular_structure, reduced_process=None):
        """Update momenta dictionary."""

        children = list(singular_structure.legs)
        parent_number = None
        for substructure in singular_structure.substructures:
            child = IRSubtraction.update_momenta_dict(momenta_dict, substructure, reduced_process)
            if child is not None:
                children.append(child)
        if singular_structure.name() in ["C", ]:
            parent_number = IRSubtraction.new_leg_number(momenta_dict)
            momenta_dict[parent_number] = frozenset(child.n for child in children)
        return parent_number

    def get_sectors(self, contrib_definition, process_map, counterterms=None, integrated_counterterms=None):
        """ Given a particular contribution definition, its process_map and the counterterms contained, this
        uses the subtraction module to build, for each defining process, the sectors to consider (or None if the subtraction
        does not require any) as well as the list of counterterms to consider for each of them.
        The sectors dictionary returned has the following format:

            sectors = {
                process_key: [
                    sector_1, sector_2, ....
                ]
            }

        where each sector_<i> is a dictionary with the following format:

            sector_<i> = {
                'sector'    :   sector_instance_or_identifier (exact format is up to the subtraction_scheme),
                'counterterms' : [list of counterterm indices (as per the ordering in self.counterterms of the integrand) to be considered],
                'integrated_counterterms' : [list of tuples (counterterm indices, input_mapping_index) for integrated CT to be considered]
            }

        'None' is returned for the sectors when not applicable for the particular subtraction scheme considered.

        """

        # If no subtraction module is loaded or if it doesn't require sectors, then return None
        if self.subtraction_scheme_module is None or self.subtraction_scheme_module.sector_generator is None:
            return None

        sectors = {}
        for process_key, (defining_process, mapped_processes) in process_map.items():
            sectors[process_key] = self.subtraction_scheme_module.sector_generator(
                contrib_definition, defining_process,
                None if counterterms is None else counterterms[process_key],
                None if integrated_counterterms is None else integrated_counterterms[process_key] )

        return sectors

    def can_be_IR_unresolved(self, PDG):
        """Check whether a particle given by its PDG can become unresolved
        and lead to singular behavior.
        """

        return any(
            (PDG in self.IR_quantities_for_corrections_types[coupling]['pert_particles'])
            for coupling in self.coupling_types )

    def parent_PDGs_from_PDGs(self, PDGs):
        """List all possible parent PDGs for a given set of children PDGs."""

        # WARNING: parent_PDGs_from_PDGs hardcoded to SM QCD (with sanity checks)
        # A generic way to build tree histories of PDGs
        # is using the model ref_dict_to1
        # misc.sprint(self.model.get('ref_dict_to1')[(-2,2)])
        # misc.sprint(self.model.get_particle(-2).get('spin'),
        #             self.model.get_particle(-2).get_color(),
        #             self.model.get_particle(-2).get('mass')=='zero'
        #             )

        if not any(
            self.model.get('name').lower().startswith(name)
            for name in self._allowed_model_names
        ):
            raise InvalidCmd(
                "parent_PDGs_from_PDGs is implemented for SM only, "
                "not in model %s." % self.model.get('name') )
        if any(order != 'QCD' for order in self.coupling_types):
            raise InvalidCmd(
                "The function parent_PDGs_from_PDGs is implemented for QCD only." )

        # Get parton flavors, eliminating gluons
        flavors = [pdg for pdg in PDGs if pdg != 21]

        # If all children were gluons, the only parent is a gluon
        if not flavors:
            return [21]

        # Consider last particle
        last = flavors.pop()
        # Look for a corresponding anti-particle
        for i in range(len(flavors)):
            anti = self.model.get_particle(last).get_anti_pdg_code()
            if flavors[i] == anti:
                # Eliminate it and start over
                flavors.pop(i)
                return self.parent_PDGs_from_PDGs(flavors)

        # If there was no anti-particle,
        # check if all other legs particles been emitted by the last one
        if self.parent_PDGs_from_PDGs(flavors) == [21]:
            return [last]

        # At this point, there is no valid parent: return empty list
        return []

    def parent_PDGs_from_legs(self, legs):
        """List all possible parent PDGs for a given set of legs."""

        # Determine if the legs contain an initial-state one
        initial_state = any(leg.state == leg.INITIAL for leg in legs)
        # Cross initial-state legs to final-state and get all PDGs
        pdgs = []
        for leg in legs:
            cur = leg.pdg
            if leg.state == leg.FINAL:
                pdgs.append(cur)
            else:
                pdgs.append(self.model.get_particle(cur).get_anti_pdg_code())

        # Get the parent PDG of this leg set
        parent_PDGs = self.parent_PDGs_from_PDGs(pdgs)

        # Cross back if the leg set referred to initial-state
        final_PDGs = []
        for parent_PDG in parent_PDGs:
            if initial_state:
                final_PDGs.append(
                    self.model.get_particle(parent_PDG).get_anti_pdg_code()
                )
            else:
                final_PDGs.append(parent_PDG)

        return final_PDGs
        
    def can_become_soft(self, legs):
        """Check whether a set of legs going simultaneously soft
        leads to singular behavior.
        """

        for pdg in self.parent_PDGs_from_legs(legs):
            particle = self.model.get_particle(pdg)
            if particle.get('spin') == 3 and particle.get('mass').lower() == 'zero':
                return True
        return False
    
    def can_become_collinear(self, legs):
        """Check whether a set of legs going collinear to each other
        leads to singular behavior.
        """
        
        for pdg in self.parent_PDGs_from_legs(legs):
            if self.can_be_IR_unresolved(pdg):
                return True
        return False        
    
    def get_all_elementary_structures(self, process, max_unresolved):
        """Generate all 'building blocks' structures relevant for the process
        passed in argument.
        """

        # Initialize variables
        elementary_structure_list = []
        # Eliminate particles that do not have a role in the subtraction
        legs = SubtractionLegSet(
            SubtractionLeg(leg) for leg in process.get('legs')
            if self.can_be_IR_unresolved(leg.get('id'))
        )
        fs_legs = SubtractionLegSet(
            leg for leg in legs if leg.state == SubtractionLeg.FINAL
        )
        is_legs = SubtractionLegSet(difference(legs, fs_legs))

        # Loop over number of unresolved particles
        for unresolved in range(1, max_unresolved + 1):
            # Get iterators at the start of the final-state list
            it = iter(fs_legs)
            soft_it, coll_final_it, coll_initial_it = itertools.tee(it, 3)
            # Final state particle sets going soft
            for soft_set in itertools.combinations(soft_it, unresolved):
                if self.can_become_soft(soft_set):
                    elementary_structure_list.append(SoftStructure(*soft_set))
            # Final state particle sets going collinear
            for coll_final_set in itertools.combinations(coll_final_it, unresolved + 1):
                if self.can_become_collinear(coll_final_set):
                    elementary_structure_list.append(CollStructure(*coll_final_set))
            # Initial-final collinear
            # For any final-state set with one less particle
            for coll_initial_set in itertools.combinations(coll_initial_it, unresolved):
                # Combine with all initial legs
                for coll_initial_leg in is_legs:
                    coll_set = (coll_initial_leg, ) + coll_initial_set
                    if self.can_become_collinear(coll_set):
                        elementary_structure_list.append(CollStructure(*coll_set))
                        
        # Finally add an elementary BeamStructure if the process acted on has factorization
        # Beam currents attached to it.
        beam_names_added = []
        for pbft in process.get('beam_factorization'):
            for beam_name in ['beam_one','beam_two']:
                if beam_name not in beam_names_added and (not pbft[beam_name] is None):
                    beam_names_added.append(beam_name)
                    elementary_structure_list.append(
                        BeamStructure(*pbft[beam_name]['singular_structure'].get_all_legs()))
        
        return elementary_structure_list

    def get_all_combinations(self, elementary_structures, max_unresolved, verbose=False):
        """Determine all combinations of elementary structures,
        applying simplification and discarding the ones that vanish.
        """

        # Duplicate elimination is ugly,
        # because more combinations than needed are generated.
        # If one were to construct the list as a product of binomials,
        # (1 - O1)(1 - O2) ...
        # adding one monomial at a time, duplicates would never appear from the beginning
        # and simplification might still occur at each stage.
        combos = [[[]]]
        strucs = [[SingularStructure()]]
        for n in range(len(elementary_structures)):
            if verbose:
                misc.sprint("Considering combinations of %d elementary structures" % (n+1))
            combos_n = []
            strucs_n = []
            n_filtered = 0
            n_void = 0
            n_duplicates = 0
            for el in elementary_structures:
                for combo in combos[-1]:
                    if el not in combo:
                        this_combo = combo + [el, ]
                        this_struc = SingularStructure(substructures=this_combo).nest()
                        if not this_struc.is_void:
                            if this_struc.count_unresolved() <= max_unresolved:
                                if this_struc not in strucs_n:
                                    combos_n.append(this_combo)
                                    strucs_n.append(this_struc)
                                else:
                                    n_duplicates += 1
                            else:
                                n_filtered += 1
                        else:
                            n_void += 1
            if verbose:
                misc.sprint(
                    "   valid: %d, void: %d, filtered: %d, duplicates: %d." %
                    (len(combos_n), n_void, n_filtered, n_duplicates) )
            combos.append(combos_n)
            strucs.append(strucs_n)
        return list(itertools.chain.from_iterable(strucs))

    def reduce_process(self, structure, process, momenta_dict):
        """Starting from a process, determine the reduced one in an unresolved limit.

        :param SingularStructure structure: singular structure which specifies the limit
        :param base_objects.Process process: the process to be reduced, will be modified!
        :param bidict momenta_dict: the momentum dictionary to be updated

        :return: the parent leg of the singular structure
        :rtype: SubtractionLeg
        """

        # Create the momenta dictionary and build the reduced process
        structure_legs = structure.get_all_legs()
        # Do not remove legs from the reduced process
        # for beam factorization contributions
        if structure.name() in ['F', ]:
            structure_leg_ns = []
        else:
            structure_leg_ns = [leg.n for leg in structure_legs]
        # Recurse into substructures
        for substructure in structure.substructures:
            subparent = self.reduce_process(substructure, process, momenta_dict)
            if subparent is not None:
                # Adding the parent
                structure_leg_ns.append(subparent.n)
                # Removing the children
                for n in [leg.n for leg in substructure.get_all_legs()]:
                    structure_leg_ns.remove(n)
        if structure.name() == "":
            return None

        # Determine parent leg and update momentum dictionary
        parent = None
        if structure.name() in ['C', ]:
            # Add entry to dictionary
            parent_number = IRSubtraction.new_leg_number(momenta_dict)
            momenta_dict[parent_number] = frozenset(structure_leg_ns)
            # Work out the complete SubtractionLeg for the parent
            parent_PDGs = self.parent_PDGs_from_legs(structure_legs)
            assert len(parent_PDGs) == 1
            parent_PDG = parent_PDGs[0]
            parent_state = SubtractionLeg.FINAL
            if structure_legs.has_initial_state_leg():
                parent_state = SubtractionLeg.INITIAL
            parent = SubtractionLeg(parent_number, parent_PDG, parent_state)
        elif structure.name() in ['S', ]:
            pass # No need to propagate parents

        elif structure.name() in ['F', ]:
            if len(structure.legs) != 1 or structure.legs[0].n not in [1, 2]:
                raise MadGraph5Error("Only beam factorizations attached to leg number 1 or" +
                                     " 2 are supported for now; not %s." % str(structure))
            beam_names = {1: 'beam_one', 2: 'beam_two'}
            # Change the distribution type of all beam_factorization terms of that beam
            # in the reduced process to be the identity (i.e. None).
            for bft in process.get('beam_factorization'):
                beam_name = beam_names[structure.legs[0].n]
                # It may be in the future that if there are several currents attached to this
                # beam we must only remove that particular one. However, it is likely that
                # any additional convolution between beam current will need to be carried out
                # analytically and treated as a singled overall counterterm.
                if not bft[beam_name] is None:
                    bft[beam_name] = None

        else:
            raise MadGraph5Error("Unrecognized structure of type : %s (name: %s)"%(str(type(structure)),structure.name()))

        # Adjust legs of the reduced process:
        # 1) remove children legs
        legs_to_remove = []
        for leg in process['legs']:
            if leg['number'] in structure_leg_ns:
                legs_to_remove.append(leg)
        for leg in legs_to_remove:
            process['legs'].remove(leg)
        # 2) add parent leg
        if parent:
            process['legs'].append(
                base_objects.Leg({
                    'number': parent.n,
                    'id': parent.pdg,
                    'state': parent.state}))
        # 3) re-sort preserving the initial state order
        rp_legs = process['legs']
        rp_legs.sort(key=lambda x: (x['state'], x['number']))
        if rp_legs[0]['number'] == 2:
            rp_legs[0], rp_legs[1] = rp_legs[1], rp_legs[0]

        # Return the parent leg if any (this is useful for the recursion)
        return parent

    def get_counterterm(self, structure, process):
        """Build the product of a set of currents and a matrix element
        that approximates the matrix element for process
        in the singular limit specified by structure.
        """

        assert type(structure) is SingularStructure
        assert isinstance(process, base_objects.Process)

        substructures = structure.substructures
        assert not structure.legs
        # Separate beam_factorization structures
        beam_structures = [s for s in substructures if s.name() == 'F']
        other_structures = [s for s in substructures if s.name() != 'F']
        other_structures = SingularStructure(substructures=other_structures)

        # Initialize a trivial momenta_dict
        momenta_dict = IRSubtraction.create_momenta_dict(process)
        # Get a modifiable copy of the reduced process
        reduced_process = process.get_copy(
            ['legs', 'n_loops', 'legs_with_decays', 'beam_factorization'])
        # We need a deep copy of the beam_factorization here
        reduced_process['beam_factorization'] = copy.deepcopy(
            reduced_process['beam_factorization'])
        # Empty legs_with_decays: it will be regenerated automatically when asked for
        reduced_process['legs_with_decays'][:] = []
        # The n_loops will be distributed later
        reduced_process.set('n_loops', -1)
        # TODO squared orders should also be reset and then adjusted later. Not done yet and optional for pure QCD corrections.
        all_bcf_combinations = {}
        # Split the reduced processes so that they have a fixed loop count in each of the two beam currents
        for bcf in process['beam_factorization']:
            key = ( None if bcf['beam_one'] is None else bcf['beam_one'].get('n_loops'),
                    None if bcf['beam_two'] is None else bcf['beam_two'].get('n_loops') )
            if key in all_bcf_combinations:
                all_bcf_combinations[key].append(bcf)
            else:
                all_bcf_combinations[key] = [bcf,]

        all_reduced_processes = {}
        for key, bcfs in all_bcf_combinations.items():
            # Get a modifiable copy of the reduced process
            reduced_process = process.get_copy(
                ['legs', 'n_loops', 'legs_with_decays', 'beam_factorization'])
            # We need a deep copy of the beam_factorization here
            reduced_process['beam_factorization'] = copy.deepcopy(bcfs)
            # Empty legs_with_decays: it will be regenerated automatically when asked for
            reduced_process['legs_with_decays'][:] = []
            # The n_loops will be distributed later
            reduced_process.set('n_loops', -1)
            all_reduced_processes[key] = reduced_process

        all_local_CT_generated = []

        for beams_loop_count, reduced_process in all_reduced_processes.items():

            # Build the nodes
            nodes = []

            # Add a node for each beam structure
            for beam_structure in beam_structures:
                # The beam factorization singular structure always have a single leg
                # with a number matching the beam number
                beam_number = beam_structure.legs[0].n
                assert (beam_number in [1, 2]), \
                    "BeamStructure is expected to have a single leg with number in [1, 2]."
                assert (beams_loop_count[beam_number-1] is not None), \
                    "There cannot be a PDF counterterm related to a process which has no beam factorisation current."
                current = BeamCurrent({
                    'beam_type': self.beam_types[beam_number - 1][0],
                    'beam_PDGs': self.beam_types[beam_number - 1][1],
                    'distribution_type': 'counterterm',
                    'singular_structure': SingularStructure(substructures=[beam_structure,]),
                    # loop counts have already been distributed when building the beam factorisation CTs
                    # so we should not distribute them again in the F(...) counterterms and therefore
                    # we do set it below to the correct value. Since this value won't be -1 when distributing loop, then
                    # the loop distribution process will keep the loop count of this current fixed.
                    'n_loops' : beams_loop_count[beam_number-1]
                })
                nodes.append(CountertermNode(current=current))
                # Modify the reduced process so as to remove the beam factorisation currents
                # like :F: which are now counterterm, i.e. [F]
                self.reduce_process(beam_structure, reduced_process, momenta_dict)

            # Add a single node for all other singular structures
            if other_structures != SingularStructure():
                current = Current({'singular_structure': other_structures})
                nodes.append(CountertermNode(current=current))

            # Determine reduced_process and momenta_dict
            self.reduce_process(other_structures, reduced_process, momenta_dict)

            # Set resolve_mother_spin_and_color to True for backward compatibility
            for subnode in nodes:
                subnode.current['resolve_mother_spin_and_color'] = True

            # Assemble the counterterm object and add it to the list of local CTs generated to return
            all_local_CT_generated.append( Counterterm(
                process=reduced_process,
                nodes=nodes,
                momenta_dict=momenta_dict,
                resolved_process=process,
                complete_singular_structure=structure
            ) )

        return all_local_CT_generated

    def get_integrated_counterterm(self, local_counterterm):
        """Given a local subtraction counterterm, return the corresponding integrated one.
        It has the same attributes, except that it contains only a single
        IntegratedCurrent, with the complete SingularStructure in it."""

        # First check that the local counterterm is singular, because if not then we
        # should of course not return any integrated counterterm.
        if not local_counterterm.is_singular():
            return []

        complete_singular_structure = local_counterterm.reconstruct_complete_singular_structure()

        # Useful to also have access to the flatten complete singular structure
        flatten_singular_structure = complete_singular_structure.decompose()

        reduced_process = local_counterterm.process.get_copy(
            ['legs', 'n_loops', 'legs_with_decays', 'squared_orders', 'beam_factorization'] )

        # The following sums all the loop numbers in all nodes and reduced process,
        # the latter of which must then be removed
        n_loops = local_counterterm.n_loops() - reduced_process.get('n_loops')

        # Retrieve the sum of squared orders in all nodes of the local counterterm
        # as well as this in the reduced process which should then be removed.
        squared_orders = local_counterterm.squared_orders()
        for order, value in reduced_process.get('squared_orders').items():
            try:
                squared_orders[order] -= value
            except KeyError:
                raise MadGraph5Error(
                    "Function squared_orders() of CountertermNode not working properly. "
                    "It should have at least the reduced process squared orders in it." )
        # filter zeros
        squared_orders = {order: value for order, value in squared_orders.items() if value!=0}

        ########
        #
        # TODO handle the mixed final-initial local counterterm case. This issue is more general, however and pertains
        # to the definition of the integrated countertem(S?) of any set of disjoint local counterterms.
        # For now only the case of a single initial state collinear local CT is supported.
        #
        ########

        # First tackle the special case of single initial-state collinear structure
        initial_state_legs = [leg for leg in complete_singular_structure.substructures[0].substructures[0].legs if
                                                                                leg.state == SubtractionLeg.INITIAL]

        # A list to hold all integrated CT generated here.
        integrated_CTs = []
        
        # Check if this local counterterm has a soft component that recoils symmetrically
        # against both beams.
        has_soft_symmetric_ISR_recoil = (self.soft_do_recoil_against_initial_states and 
                complete_singular_structure.does_require_correlated_beam_convolution(self.subtraction_scheme_module))

        # Handle the specific case of single initial-state pure collinear counterterm or 
        # soft recoil against IS.
        if (len(initial_state_legs) == 1 and not any(ss.name()=='S' 
                   for ss in flatten_singular_structure)) or has_soft_symmetric_ISR_recoil:
            
            if has_soft_symmetric_ISR_recoil:
                # In the case of soft recoiling against initial states, the beam_type should not matter.
                beam_type = None
                beam_PDGs = None
            else:  
                beam_number = initial_state_legs[0].n
                beam_names = {1:'beam_one', 2:'beam_two'}
                        
                assert (beam_number in [1,2]), "MadNkLO only supports initial state legs with number in [1,2]."
                beam_type = self.beam_types[beam_number-1][0] if self.beam_types[beam_number-1] else None
                beam_PDGs = self.beam_types[beam_number-1][1] if self.beam_types[beam_number-1] else None

            # if has_soft_symmetric_ISR_recoil:
            #     # In the case of soft recoiling against initial states, the beam_type should not matter.
            #     #gl
            #     beam_type = None
            #     beam_PDGs = None
            # else:  
            #     beam_number = initial_state_legs[0].n
            #     beam_names = {1:'beam_one', 2:'beam_two'}
                        
            #     assert (beam_number in [1,2]), "MadNkLO only supports initial state legs with number in [1,2]."
            #     beam_type = self.beam_types[beam_number-1][0] if self.beam_types[beam_number-1] else None
            #     beam_PDGs = self.beam_types[beam_number-1][1] if self.beam_types[beam_number-1] else None

            beam_current_options = {
                    'beam_type'     :   beam_type,
                    'beam_PDGs'     :   beam_PDGs,
                    'n_loops'       :   n_loops,
                    'squared_orders':   squared_orders,
                    'singular_structure': complete_singular_structure
                }
            # print('beam_current_options : ' + str(beam_current_options))

            # For now there will be no additional CT node, but for NNLO there may be because
            # the integrated ISR has to be combined with local FSR counterterms
            additional_CT_nodes = []

            # Handle the specific case of NLO single beam factorization counterterm
            ss_name = complete_singular_structure.substructures[0].substructures[0].name()
            # print('ss name : ' + str(ss_name))
            if ss_name == 'F':

                integrated_CTs.append(IntegratedCounterterm(
                    process     = reduced_process,
                    nodes       = [
                            CountertermNode(current=IntegratedBeamCurrent(beam_current_options)), 
                        ]+additional_CT_nodes,
                    momenta_dict= bidict(local_counterterm.momenta_dict),
                    prefactor   = -1. * local_counterterm.prefactor ))

            # Handle the specific case of ISR recoils, either because of ISR collinear or
            # soft counterterms in colorful and ppToOneWalker mapping.
            elif has_soft_symmetric_ISR_recoil or ss_name == 'C':
                # Now add the three pieces of the integrated ISR. Ultimately we may choose to instead
                # generate only the first piece ('bulk') and the other two via an iterative application
                # of the subtraction procedure. For now however, we will generate all three at once.
                
                # Starting with the integrated endpoint:
                integrated_CTs.append(IntegratedCounterterm(
                    process     = reduced_process,
                    nodes       = [
                            CountertermNode(current=IntegratedBeamCurrent(beam_current_options)), 
                        ]+additional_CT_nodes,
                    momenta_dict= bidict(local_counterterm.momenta_dict),
                    prefactor   = -1. * local_counterterm.prefactor ))
    
                # Then the bulk and counterterm integrated ISR pieces, which are both xi dependent.
                beam_current_options['distribution_type'] = 'bulk'
                integrated_CTs.append(IntegratedCounterterm(
                    process     = reduced_process,
                    nodes       = [
                            CountertermNode(current=BeamCurrent(beam_current_options)), 
                        ]+additional_CT_nodes,
                    momenta_dict= bidict(local_counterterm.momenta_dict),
                    prefactor   = -1. * local_counterterm.prefactor ))
                
                beam_current_options['distribution_type'] = 'counterterm'
                integrated_CTs.append(IntegratedCounterterm(
                    process     = reduced_process,
                    nodes       = [
                            CountertermNode(current=BeamCurrent(beam_current_options)), 
                        ]+additional_CT_nodes,
                    momenta_dict= bidict(local_counterterm.momenta_dict),
                    prefactor   = +1. * local_counterterm.prefactor ))
            
        else:
            # Here is the general solution chosen for arbitrary final state local CT
            integrated_current = IntegratedCurrent({
                'n_loops': n_loops,
                'squared_orders': squared_orders,
                'resolve_mother_spin_and_color': True,
                'singular_structure': complete_singular_structure })

            integrated_CTs.append(IntegratedCounterterm(
                process     = reduced_process,
                nodes       = [CountertermNode(current=integrated_current), ],
                momenta_dict= bidict(local_counterterm.momenta_dict),
                prefactor   = -1. * local_counterterm.prefactor ))

        return integrated_CTs

    @staticmethod
    def get_all_currents_blocks(counterterms,track_leg_numbers=False):
        """Deduce the list of currents needed to compute all counterterms given."""
        all_currents_blocks = []
        for counterterm in counterterms:
            for currents_block in counterterm.get_currents_blocks():
                # Remove duplicates already at this level
                if track_leg_numbers or currents_block not in all_currents_blocks:
                    all_currents_blocks.append(currents_block)
        return all_currents_blocks

    def get_all_counterterms(
        self, process, max_unresolved,
        extra_unresolved_in_combination=0,
        ignore_integrated_counterterms=False):
        """Generate all counterterms for the corrections specified in this module
        and the process given in argument."""
        
        elementary_structures = self.get_all_elementary_structures(
            process, max_unresolved)
        combinations = self.get_all_combinations(
            elementary_structures, max_unresolved + extra_unresolved_in_combination)

        target_n_loops = process['n_loops'] + process.get_n_loops_in_beam_factorization()

        all_counterterms = []
        all_integrated_counterterms = []
        for combination in combinations:
            for template_counterterm in self.get_counterterm(combination, process):

#               Uncomment the code below to remove all counterterms that have both soft and col
#                if any(
#                        (
#                            any(s.name()=='S' for s in c.get('singular_structure').decompose()) and
#                            any(s.name()=='C' for s in c.get('singular_structure').decompose())
#                        ) for c in template_counterterm.get_all_currents()
#                ):
#                    continue

#               Uncomment the code below to retain remove pure collinear counterterms
#                if any(
#                        (
#                            any(s.name()=='C' for s in c.get('singular_structure').decompose()) and True
#                            not (any(s.name()=='S' for s in c.get('singular_structure').decompose()))
#                        ) for c in template_counterterm.get_all_currents()
#                ):
#                    continue

                counterterms_with_loops = template_counterterm.distribute_loops(target_n_loops)


                for counterterm_with_loops in counterterms_with_loops:
                    # TODO
                    # For the time being, distribute_orders is given None instead of the squared orders
                    # because they should be retrieved from the process by looking at individual
                    # matrix elements
                    # That is, every process has a list of possible coupling orders assignations
                    # so we should loop over them
                    counterterm_with_loops_and_orders = counterterm_with_loops.distribute_orders(None)
                    for local_CT in counterterm_with_loops_and_orders:
                        all_counterterms.append(local_CT)
                        if not ignore_integrated_counterterms:
                            integrated_counterterms = self.get_integrated_counterterm(local_CT)
                            all_integrated_counterterms.extend(integrated_counterterms)

        # Now for each integrated counterterm generated, set the attribute 'requires_correlated_beam_convolution.'
        for integrated_CT in all_integrated_counterterms:
            integrated_CT.set_requires_correlated_beam_convolution(self.subtraction_scheme_module)
            integrated_CT.set_n_non_factorisable_double_sided_convolution(self.subtraction_scheme_module)
            #gl
            integrated_CT.set_requires_torino_sub_BS(self.subtraction_scheme_module)

        # Finally, if this integrated counterterm implies multiple convolutions in the *same* convolution parameter
        # then one of them will be performed analytically, thus inducing a new building-block counterterm implementing
        # that composite object. Example:
        #     :[(C(1,4),)]:  with reduced process B x :F(1):
        # Cannot be implemented as direct product of convolutions, instead it must correspond to novel composite object
        # defined as follows:
        #     :[(C(1,4),F(1))]: with reduced process B x 1 (no convolution attached to this reduced process).
        # Note however that we may have to create additional counterterms, because if the situation is:
        #     :[(C(1,4),)]:  with reduced process B x { :F(1)^(0): | :F(2)^(1): + :F(1)^(1): | :F(2)^(0): }
        # Then we have to split this into two counterterms:
        #     :[(C(1,4),F(1)^(0)]: with reduced process B x :F(2)^(1):
        #     :[(C(1,4),F(1)^(1)]: with reduced process B x :F(2)^(0):
        all_merged_integrated_counterterms = []
        for integrated_CT in all_integrated_counterterms:
            merged_integrated_counterterms = integrated_CT.merge_convolutions(self.subtraction_scheme_module)
            # We must refresh the attributes correlated beam and non factorisable convolutions.
            for merged_integrated_counterterm in merged_integrated_counterterms:
                merged_integrated_counterterm.set_requires_correlated_beam_convolution(self.subtraction_scheme_module)
                merged_integrated_counterterm.set_n_non_factorisable_double_sided_convolution(self.subtraction_scheme_module)
                #gl
                merged_integrated_counterterm.set_requires_torino_sub_BS(self.subtraction_scheme_module)
            all_merged_integrated_counterterms.extend(merged_integrated_counterterms)


        return all_counterterms, all_merged_integrated_counterterms

#=========================================================================================
# Subtraction current exporter
#=========================================================================================

class SubtractionCurrentExporter(object):
    """Class for mapping and exporting the subtraction currents to a given location
    and generate the corresponding accessors as well.
    """
    
    template_dir = pjoin(MG5DIR, 'madgraph', 'iolibs', 'template_files', 'subtraction')
    template_modules_path = 'madgraph.iolibs.template_files.subtraction'
    
    # The main module name in a process output
    main_module_name = 'SubtractionResources'

    @staticmethod
    def get_all_available_subtraction_schemes():
        """ Inspects the subtraction template directory to identify all available subtraction schemes defined."""

        subtraction_schemes_available = []
        for dir_name in os.listdir(pjoin(SubtractionCurrentExporter.template_dir,'subtraction_schemes')):
            if not os.path.exists(pjoin(SubtractionCurrentExporter.template_dir,
                                                        'subtraction_schemes', dir_name, '__init__.py')):
                continue
            subtraction_schemes_available.append(dir_name)
        return subtraction_schemes_available

    @staticmethod
    def get_subtraction_scheme_module(subtraction_scheme, root_path=None):
        """ Loads the subtractions scheme module specified.
        One can also decide to load this module from the specified root_path
        which can be a process output for instance."""

        # Attempt to import the subtraction scheme
        try:
            if root_path is None:
                subtraction_scheme_module = importlib.import_module('%s.%s.%s' % (
                    SubtractionCurrentExporter.template_modules_path, 'subtraction_schemes', subtraction_scheme))
                subtraction_utils = importlib.import_module('%s.%s.%s' % (
                    SubtractionCurrentExporter.template_modules_path, 'commons', 'utils'))
            else:
                sys.path.insert(0, root_path)
                subtraction_scheme_module = importlib.import_module('%s.%s.%s' % (
                         SubtractionCurrentExporter.main_module_name, 'subtraction_schemes', subtraction_scheme))
                subtraction_utils = importlib.import_module('%s.%s.%s' % (
                    SubtractionCurrentExporter.main_module_name, 'commons', 'utils'))
                sys.path.pop(0)
        except ImportError as e:
            raise MadGraph5Error("Specified subtraction module could not be found or loaded: %s. Error:\n%s"%(
                                                                                            subtraction_scheme, str(e)))

        # Now load this subtraction scheme (i.e. defining all static module variables)
        subtraction_scheme_module.load()
        for attr, value in subtraction_scheme_module.loaded_attributes.items():
            setattr(subtraction_scheme_module, attr, value)

        mandatory_attributes = ['all_subtraction_current_classes','__authors__','exporter',]
        missing_attributes = [ attr for attr in mandatory_attributes if not hasattr(subtraction_scheme_module,attr) ]
        if len(missing_attributes)>0:
            raise MadGraph5Error("The specified subtraction module '%s' is missing the following mandatory attributes: %s"%
                                                                        (subtraction_scheme, missing_attributes))

        # Sanity check that all these currents class implementations are valid
        for current_class in subtraction_scheme_module.all_subtraction_current_classes:
            if not (inspect.isclass(current_class) and hasattr(current_class, 'does_implement_these_currents') ):
                raise MadGraph5Error("The current '%s' does not implement the mandatory function 'does_implement_these_currents'."%
                                                                                                (current_class.__name__))

        # Add an attribute to the module which is a dictionary mapping a subtraction current unique identifier to the
        # class implementing it
        if not hasattr(subtraction_scheme_module, 'ordered_currents_identifiers'):
            subtraction_scheme_module.ordered_currents_identifiers = [
                current_class.__name__ for current_class in subtraction_scheme_module.all_subtraction_current_classes ]

        if not hasattr(subtraction_scheme_module, 'currents'):
            subtraction_scheme_module.currents = dict(
                (current_class.__name__, current_class)
                    for current_class in subtraction_scheme_module.all_subtraction_current_classes)

        # Add a default implementation for the currents used for debugging
        if 'DefaultCurrentImplementation' not in subtraction_scheme_module.ordered_currents_identifiers:
            # If there is a "default current" in the utils class,
            # presumably used for debugging only, then add this one at the very end of all_classes so that
            # it will be selected only if no other class matches.
            if hasattr(subtraction_utils, 'DefaultCurrentImplementation'):
                default_implementation_class = getattr( subtraction_utils, 'DefaultCurrentImplementation' )

                if (inspect.isclass(default_implementation_class) and
                    hasattr(default_implementation_class, 'does_implement_these_currents')):
                    subtraction_scheme_module.ordered_currents_identifiers.append('DefaultCurrentImplementation')
                    subtraction_scheme_module.currents['DefaultCurrentImplementation'] = default_implementation_class

        return subtraction_scheme_module

    def __init__(self, model, export_dir, subtraction_scheme):
        """Initialize the exporter.

        :param model: Model to be used when trying to match currents.
        :param export_dir: Target directory for currents to be exported to.
        :param subtraction_scheme: Lookup is performed within the subtraction scheme module specified
        """

        self.model       = model
        self.export_dir  = export_dir
        self.subtraction_scheme = subtraction_scheme

    def export(self, all_defining_currents_blocks):
        """Export the specified list of currents and return a list of accessors
        which contain the mapping information.
        """
        subtraction_scheme_module = SubtractionCurrentExporter.get_subtraction_scheme_module(self.subtraction_scheme)
        track_leg_numbers=subtraction_scheme_module.are_current_instances_for_specific_leg_numbers

        # Group all currents according to mapping rules.
        # The mapped currents is a dictionary of the form
        #         { (subtraction_scheme_name, current_identifier, instantiation_options_index),
        #                {'defining_currents': <...>,
        #                 'mapped_process_keys': [<...>],
        #                 'instantiation_options': instantiation_options
        #         }
        mapped_currents_blocks = {}
        
        # Look in all current implementation classes found
        # and find which one implements each current
        all_instantiation_options = []
        currents_with_default_implementation = []
        for currents_block in all_defining_currents_blocks:
            found_current_class = None
            for currents_identifier in subtraction_scheme_module.ordered_currents_identifiers:
                implementation_class = subtraction_scheme_module.currents[currents_identifier]

                #misc.sprint('Testing class %s for current: %s'%(currents_identifier, str(currents_block)))
                instantiation_options = implementation_class.does_implement_these_currents(currents_block, self.model )
                #misc.sprint('Result: %s'%str(instantiation_options))
                if instantiation_options is None:
                    continue
                if found_current_class is not None:
                    if currents_identifier == 'DefaultCurrentImplementation':
                        continue
                    logger.critical(
                        "%s found, %s already found" % (currents_identifier, found_current_class))
                    raise MadGraph5Error("Multiple implementations found for current %s." %currents_block)
                try:
                    instantiation_options_index = all_instantiation_options.index( instantiation_options )
                except ValueError:
                    all_instantiation_options.append(instantiation_options)
                    instantiation_options_index = len(all_instantiation_options)-1

                key = ( self.subtraction_scheme, currents_identifier, instantiation_options_index )

                if key in mapped_currents_blocks:
                    mapped_currents_blocks[key]['mapped_process_keys'].append(
                                                    currents_block.get_key(track_leg_numbers=track_leg_numbers) )
                else:
                    mapped_currents_blocks[key] = {
                        'defining_currents_block': currents_block,
                        'mapped_process_keys': [ currents_block.get_key(track_leg_numbers=track_leg_numbers), ],
                        'instantiation_options': instantiation_options }

                if currents_identifier == 'DefaultCurrentImplementation':
                    currents_with_default_implementation.append(currents_block)
                found_current_class = currents_identifier

            if found_current_class is None:
                raise MadGraph5Error(
                    "No implementation was found for current %s."%str(currents_block) )

        # Warn the user whenever DefaultCurrentImplementation is used
        # (it should never be used in production)
        if currents_with_default_implementation:
            currents_str = '\n'.join(
                ' > (%s) (type: (%s) )' % (', '.join(str(crt) for crt in crts), ', '.join(crt.get_type() for crt in crts) )
                for crts in currents_with_default_implementation )
            msg = """No implementation was found for the following subtraction currents:
%s
The class 'DefaultCurrentImplementation' will therefore be used for it
but results obtained in this way are very likely wrong 
and should be used for debugging only.""" % currents_str
            if __debug__:
                logger.critical(msg)
            else:
                raise MadGraph5Error(msg)

        # Export all the relevant resources to the specified export directory (process output)
        if self.export_dir is not None:
            subtraction_scheme_module.exporter.export(pjoin(self.export_dir, self.main_module_name), mapped_currents_blocks)

        return mapped_currents_blocks
        
#=========================================================================================
# Standalone main for debugging / standalone trials
#=========================================================================================

if __name__ == '__main__':
    misc.sprint("Put your standalone subtraction code here.")
