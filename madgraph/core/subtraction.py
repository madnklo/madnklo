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
    """Leg object specialized for subtraction."""

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

    def __str__(self, print_n = True, print_pdg = False, print_state = False):
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

    def without_leg_numbers(self):
        """Return a copy of this SubtractionLegSet without leg numbers."""

        return SubtractionLegSet((0, leg.pdg, leg.state) for leg in self)

#=========================================================================================
# SingularStructure
#=========================================================================================

class SingularStructure(object):
    """Object that represents a hierarchical structure of IR singularities."""

    def __init__(self, *args, **opts):
        """Initialize a hierarchical singular structure."""

        self.substructures = opts.get('substructures', [])
        self.legs = opts.get('legs', SubtractionLegSet())
        self.is_void = opts.get('is_void', False)
        for arg in args:
            if isinstance(arg, SingularStructure):
                self.substructures.append(arg)
            elif isinstance(arg, SubtractionLeg):
                self.legs = SubtractionLegSet(self.legs + (arg, ))
            elif isinstance(arg, bool):
                self.is_void = self.is_void or arg
            else:
                raise MadGraph5Error(
                    "Invalid argument in SingularStructure.__init__: %s" % str(arg)
                )

    def get_copy(self):
        """Provide a modifiable copy of this singular structure."""

        return type(self)(
            legs=SubtractionLegSet(self.legs),
            substructures=[ss.get_copy() for ss in self.substructures],
            is_void=self.is_void
        )

    def __str__(self, print_n=True, print_pdg=False, print_state=False):
        """Return a string representation of the singular structure."""

        if self.is_void:
            return 'Void structure'

        tmp_str = self.name() + "("
        tmp_str += ",".join(sorted(
            sub.__str__(print_n, print_pdg, print_state)
            for sub in self.substructures
        ))
        if self.substructures:
            tmp_str += ","
        tmp_str += ",".join(
            leg.__str__(print_n, print_pdg, print_state)
            for leg in self.legs
        )
        tmp_str += ")"
        return tmp_str

    def get_subtraction_prefactor(self):
        """Determine the prefactor due to the nested subtraction technique."""
        
        # If we are at the top level, we shall not include the factor
        # Indeed, we are integrating +R-C.
        if type(self) == SingularStructure:
            pref = 1
        else:
            pref = -1

        # Account for simplification of soft operators
        # TODO Check this works for collinears inside softs
        softs = []
        for sub in self.substructures:
            pref *= sub.get_subtraction_prefactor()
            if sub.name() == "S":
               softs += [len(sub.substructures) + len(sub.legs), ]
        pref *= multinomial(softs)
        
        return pref

    def annihilate(self):
        """When an operator cannot act on this structure,
        remove this structure altogether
        by flagging it as unapplicable to any other operator.
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

    def discard_leg_numbers(self):
        """Set all leg numbers to zero, discarding this information forever."""

        self.legs = self.legs.without_leg_numbers()
        for substructure in self.substructures:
            substructure.discard_leg_numbers()

    def count_unresolved(self):
        """Count the number of unresolved legs."""

        total_unresolved = 0
        for substructure in self.substructures:
            total_unresolved += substructure.count_unresolved()
        return total_unresolved

    def get_canonical_representation(self, track_leg_numbers=True):
        """Creates a canonical hashable representation of self."""
        
        canonical = dict()

        canonical['is_void'] = self.is_void
        if track_leg_numbers:
            canonical['legs'] = self.legs
        else:
            canonical['legs'] = self.legs.without_leg_numbers()
        canonical['name'] = self.name()
        canonical['substructures'] = tuple(
            structure.get_canonical_representation(track_leg_numbers)
            for structure in self.substructures
        )

        return tuple(sorted(canonical.items()))

    def __eq__(self, other):
        """Check if two singular structures are the same."""

        assert isinstance(other, SingularStructure)
        self_can = self.get_canonical_representation()
        other_can = other.get_canonical_representation()
        return self_can == other_can

    def name(self):

        return ""

class SoftStructure(SingularStructure):

    def count_unresolved(self):

        return len(self.get_all_legs())

    def name(self):

        return "S"

class CollStructure(SingularStructure):

    def count_unresolved(self):

        return len(self.get_all_legs())-1

    def name(self):

        return "C"

#=========================================================================================
# SingularOperator
#=========================================================================================

class SingularOperator(SubtractionLegSet):
    """Virtual base class for elementary singular operators."""

    def __new__(cls, *args, **opts):

        return super(SingularOperator, cls).__new__(cls, args, **opts)

    def __init__(self, *args, **opts):

        super(SingularOperator, self).__init__(args, **opts)

    def __str__(self):
        """Return a simple string representation of the singular operator."""

        return self.name() + super(SingularOperator, self).__str__()

    def name(self):
        """Symbol used to represent this operator within output."""

        raise MadGraph5Error(
            "name called in SingularOperator of unspecified type."
        )

    def get_structure(self):
        """Singular structure corresponding to this operator
        acting on hard particles.
        """

        raise MadGraph5Error(
            "get_structure called in SingularOperator of unspecified type."
        )

    def act_here_needed(self, structure):
        """Return True if the counterterm obtained acting with this operator
        on the given structure at the first level is needed for subtraction.
        """

        # WARNING Rules may not extend to non-QCD cases
        raise MadGraph5Error(
            "act_here_needed called in SingularOperator of unspecified type."
        )

    def non_overlapping_with(self, structure):
        """Determines if this operator does not affect any leg
        that appears within the given structure, at any level.
        """

        assert isinstance(structure, SingularStructure)

        if not is_disjoint(self, structure.legs):
            return False

        for substructure in structure.substructures:
            if not self.non_overlapping_with(substructure):
                return False

        return True

    def act_on(self, structure):
        """Act with an elementary operator on a non-elementary structure."""

        assert isinstance(structure, SingularStructure)

        # If the limit does not overlap with the existing structure at all,
        # just append it at the end
        if self.non_overlapping_with(structure):
            structure.substructures.append(self.get_structure())
            return

        # If the limit acts at least partly at the current level
        if not is_disjoint(self, structure.legs):
            # If the limit acts at different levels,
            # it is not needed for subtraction
            if not is_subset(self, structure.legs):
                structure.annihilate()
                return
            # If the limit acts completely within this level,
            # it may be needed or not
            else:
                if self.act_here_needed(structure):
                    structure.substructures.append(self.get_structure())
                    structure.legs = difference(structure.legs, self)
                    return
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

class SoftOperator(SingularOperator):
    """Object that represents a soft elementary singular operator."""

    def __init__(self, *args, **opts):

        super(SoftOperator, self).__init__(*args, **opts)

    def name(self):

        return "S"

    def get_structure(self):

        return SoftStructure(legs=self)

    def act_here_needed(self, structure):

        assert isinstance(structure, SingularStructure)

        # Soft limits are not needed within soft structures
        if structure.name() == "S":
            return False

        # Soft limits are needed within collinear structures
        # if there is at least one hard particle left
        if structure.name() == "C":
            # If there is a hard particle at this level, the limit is needed
            if len(self) < len(structure.legs):
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

class CollOperator(SingularOperator):
    """Object that represents a collinear elementary singular operator."""

    def name(self):

        return "C"

    def get_structure(self):

        return CollStructure(legs=self)

    def __init__(self, *args, **opts):

        super(CollOperator, self).__init__(*args, **opts)

    def act_here_needed(self, structure):

        assert isinstance(structure, SingularStructure)

        # Collinear limits are needed within soft _and_ collinear structures
        # WARNING: this statement has to be checked
        # By default, not skipping overlaps ensures correct subtraction
        return True

#=========================================================================================
# SingularOperatorList
#=========================================================================================

class SingularOperatorList(list):
    """List of subtraction operators."""

    def act_on(self, structure):
        """Act with a list of operators on a substructure."""

        # Empty list of operators to apply or invalid structure: done
        if structure.is_void or not self:
            return structure
        # Apply last operator and advance recursion
        all_but_last_operators = SingularOperatorList(self[:-1])
        # Act using the last operator
        self[-1].act_on(structure)
        # And now act using the rest of the operators, if needed
        if not structure.is_void:
            return all_but_last_operators.act_on(structure)

    def simplify(self):
        """Simplify a list of operators,
        returning the real structure of the corresponding singularities.
        """

        structure = SingularStructure()
        self.act_on(structure)
        return structure

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

    def count_unresolved(self):
        """Count the number of unresolved particles covered by this current."""

        return self['singular_structure'].count_unresolved()

    def discard_leg_numbers(self):
        """Discard all leg numbers in the singular_structure
        and in the parent_subtraction_leg,
        effectively getting the type of current
        without information about momenta for which it should be evaluated.
        """

        self['singular_structure'].discard_leg_numbers()

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

    def get_key(self):
        """Return the ProcessKey associated to this current."""

        from contributions import ProcessKey
        return ProcessKey(
                self,
                allowed_attributes = [
                    'singular_structure',
                    'n_loops',
                    'squared_orders',
                    'resolve_mother_spin_and_color' ] )

    def get_copy(self, copied_attributes=()):
        """Return a copy of this current with a deep-copy of its singular structure."""

        copied_current = super(Current, self).get_copy(
            tuple(attr for attr in copied_attributes if attr != 'singular_structure') )
        if 'singular_structure' in copied_attributes:
            copied_current['singular_structure'] = self['singular_structure'].get_copy()
        
        return copied_current

    def __eq__(self, other):
        """Compare two currents using their ProcessKey's."""

        return self.get_key().key_dict == other.get_key().key_dict

#=========================================================================================
# Integrated current
#=========================================================================================

class IntegratedCurrent(Current):
    """Class for integrated currents."""

    # TODO
    # For now, it behaves exactly as a local 4D current,
    # but it is conceptually different.

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
        return '[integrated] %s' % res

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

    def split_loops(self, n_loops):
        """Split a counterterm in several ones according to the individual
        loop orders of its currents.
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
            for loop_n in range(n_loops + 1):
                for first_with_loops in first_without_loops.split_loops(loop_n):
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
                        for ith_node in self.nodes[i_current].split_loops(new_loop_n):
                            new_node_combinations.append(
                                combination + [ith_node, ] )
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

    def split_orders(self, target_squared_orders):
        """Split a counterterm in several ones according to the individual
        coupling orders of its currents.
        """

        # For now this sets the orders recursively without creating copies
        # because there is only one choice
        # Eventually, it should create new instances, set orders there, and return the
        # new list of counterterm copies with orders
        # TODO actually create a new counterterm with orders set

        for node in self.nodes:
            node.split_orders(target_squared_orders)

        if isinstance(self.current, Current):
            self.current.set(
                'squared_orders',
                { 'QCD':
                        (self.current.get('n_loops') * 2 +
                         self.current[
                             'singular_structure'].count_unresolved() * 2) } )
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

    def __str__(self, print_n=True, print_pdg=False, print_state=False):

        return self.reconstruct_complete_singular_structure().__str__(
            print_n=print_n, print_pdg=print_pdg, print_state=print_state )

    def nice_string(self, lead="    ", tab="    "):

        tmp_str  = lead + self.process.nice_string(0, True, False)
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
        return type(self)(
            process=node.current, nodes=node.nodes,
            momenta_dict=momenta_dict, prefactor=self.prefactor
        )

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

    def get_daughter_pdgs(self, leg_number, state):
        """Walk down the tree of currents to find the pdgs 'attached'
        to a given leg number of the reduced process.
        """
        
        external_leg_numbers = []
        intermediate_leg_number = [leg_number, ]
        while len(intermediate_leg_number) > 0:
            next_leg = intermediate_leg_number.pop(0)
            # Check if this leg is final
            daughters = self.momenta_dict[next_leg]
            if  daughters == frozenset([next_leg,]):
                external_leg_numbers.append(next_leg)
            else:
                intermediate_leg_number = list(daughters) + intermediate_leg_number
        
        # Now that we have the external leg numbers, find the corresponding pdg
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

    def get_all_currents(self):
        """Return a list of all currents involved in this counterterm."""

        currents = []
        for node in self.nodes:
            currents += node.get_all_currents()
        return currents

    def get_reduced_flavors(self, defining_flavors=None):
        """Given the defining flavors corresponding to the resolved process, return a
         *list* of flavors corresponding to the flavor assignment of the reduced process 
         given the defining flavors"""

        # Now construct the reduced flavors list
        # /!\ In the case of integrated counterterm, this is not only relevant for the
        # flavour sensitive cuts in pass_cuts and for the observable but also for the
        # actual evalution of the *initial-final integrated* counterterm so as to 
        # decide which PDF to convolute against.
        if (defining_flavors is None) or True:
            # If no defining flavors are specified then simply use the flavors already
            # filled in the reduced process.
            reduced_flavors = self.process.get_initial_final_ids()
        else:
            ##############################################################################
            # --> TODO AND NECESSARY FOR NON-FLAVOR-BLIND OBSERVABLES
            #     For now always use the flavors of the defining process, which is
            #     fine for flavor blind observables
            ##############################################################################
            reduced_flavors = None

        return reduced_flavors

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
        reduced_PS = LorentzVectorList()
        
        # First construct the reduced kinematic list
        for leg in self.process.get_initial_legs():
            reduced_PS.append(reduced_PS_dict[leg.get('number')])
        for leg in self.process.get_final_legs():
            reduced_PS.append(reduced_PS_dict[leg.get('number')])

        reduced_flavors = self.get_reduced_flavors(defining_flavors)
        
        return reduced_PS, reduced_flavors

#=========================================================================================
# IntegratedCounterterm
#=========================================================================================

class IntegratedCounterterm(Counterterm):
    """A class for the integrated counterterm."""
    
    def get_reduced_kinematics(self, input_reduced_PS):
        """Given the PS *dictionary* providing the reduced kinematics corresponding to
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
        all_legs = self.process.get_initial_legs()+self.process.get_final_legs()
        if isinstance(input_reduced_PS, dict):
            for i, leg in enumerate(all_legs):
                assert( i+1 in input_reduced_PS )
                reduced_PS[leg.get('number')] = input_reduced_PS[i+1]
        else:
            n_legs = self.process.get_ninitial()
            n_legs += len(self.process.get_final_ids_after_decay())
            assert (len(input_reduced_PS) == n_legs)
            for i, leg in enumerate(all_legs):
                reduced_PS[leg.get('number')] = input_reduced_PS[i]

        return reduced_PS

#=========================================================================================
# IRSubtraction
#=========================================================================================

class IRSubtraction(object):

    def __init__(self, model, n_unresolved, coupling_types=('QCD', )):
        """Initialize a IR subtractions for a given model,
        correction order and type.
        """
        
        self.model = model
        self.n_unresolved = n_unresolved
        self.coupling_types = coupling_types
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
            for name in ['sm', 'loop_sm', 'loopsm', 'simple_qcd'] ):
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
        lead to singular behavior.
        """

        for pdg in self.parent_PDGs_from_legs(legs):
            particle = self.model.get_particle(pdg)
            if particle.get('spin') == 3 and particle.get('mass').lower() == 'zero':
                return True
        return False
    
    def can_become_collinear(self, legs):
        """Check whether a set of legs going collinear to each other
        lead to singular behavior.
        """
        
        for pdg in self.parent_PDGs_from_legs(legs):
            if self.can_be_IR_unresolved(pdg):
                return True
        return False        
    
    def get_all_elementary_operators(self, process):
        """Generate all 'building blocks' operators relevant for the process
        passed in argument.
        """
        
        elementary_operator_list = SingularOperatorList([])
        
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
        for unresolved in range(1, self.n_unresolved + 1):
            # Get iterators at the start of the final-state list
            it = iter(fs_legs)
            soft_it, coll_final_it, coll_initial_it = itertools.tee(it, 3)
            # Final state particle sets going soft
            for soft_set in itertools.combinations(soft_it, unresolved):
                if self.can_become_soft(soft_set):
                    elementary_operator_list.append(SoftOperator(*soft_set))
            # Final state particle sets going collinear
            for coll_final_set in itertools.combinations(coll_final_it, unresolved + 1):
                if self.can_become_collinear(coll_final_set):
                    elementary_operator_list.append(CollOperator(*coll_final_set))
            # Initial-final collinear
            # For any final-state set with one less particle
            for coll_initial_set in itertools.combinations(coll_initial_it, unresolved):
                # Combine with all initial legs
                for coll_initial_leg in is_legs:
                    coll_set = (coll_initial_leg, ) + coll_initial_set
                    if self.can_become_collinear(coll_set):
                        elementary_operator_list.append(CollOperator(*coll_set))
        return SingularOperatorList(elementary_operator_list)


    def get_all_combinations(
        self, elementary_operators,
        max_unresolved=None, verbose=False
    ):
        """Determine all combinations of elementary operators,
        applying simplification and discarding the ones that vanish.
        """

        if max_unresolved is None:
            unresolved = self.n_unresolved
        else:
            unresolved = max_unresolved

        combos = [[SingularOperatorList()]]
        strucs = [[SingularStructure()]]
        for n in range(len(elementary_operators)):
            if verbose:
                misc.sprint("Considering combinations of %d operators" % (n+1))
            combos_n = []
            strucs_n = []
            n_filtered = 0
            n_void = 0
            for op in elementary_operators:
                for combo in combos[-1]:
                    if op not in combo:
                        this_combo = SingularOperatorList(combo + [op, ])
                        this_struc = this_combo.simplify()
                        if not this_struc.is_void:
                            if this_struc.count_unresolved() <= unresolved:
                                combos_n.append(this_combo)
                                strucs_n.append(this_struc)
                            else:
                                n_filtered += 1
                        else:
                            n_void += 1
            if verbose:
                misc.sprint(
                    "   valid: %d, void: %d, filtered: %d." %
                    (len(combos_n), n_void, n_filtered) )
            combos.append(combos_n)
            strucs.append(strucs_n)
        return list(itertools.chain.from_iterable(strucs))

    def get_counterterm(
        self,
        structure,
        process,
        momenta_dict_so_far=None
    ):
        """Build the product of a set of currents and a matrix element
        that approximates the matrix element for process
        in the singular limit specified by structure.
        Also build the integrated counterterm that cancels the contribution of the local
        subtraction counterterm inclusively over the splitting phase-space in d-dimension.
        """

        # 1. Initialize variables

        assert isinstance(structure, SingularStructure)
        assert isinstance(process, base_objects.Process)

        reduced_process = process

        # If no momenta dictionary was passed
        if not momenta_dict_so_far:
            # Initialize it with process legs
            momenta_dict_so_far = bidict()
            for leg in process['legs']:
                # Check that legs are numbered progressively
                # from 1 to len(process['legs']),
                # else a more elaborate treatment of indices is needed
                assert leg['number'] == len(momenta_dict_so_far) + 1
                momenta_dict_so_far[leg['number']] = frozenset((leg['number'], ))

            # The squared orders of the reduced process will be set correctly later
            reduced_process = reduced_process.get_copy(
                ['legs', 'n_loops', 'legs_with_decays']
            )
            # Empty legs_with_decays: it will be regenerated automatically when asked for
            reduced_process['legs_with_decays'][:] = []
            # TODO If resetting n_loops, what about orders?
            # The n_loops will be distributed later
            reduced_process.set('n_loops', -1)

        nodes = []

        # 2. Recursively look into substructures

        current_args = set(structure.legs)
        for substructure in structure.substructures:
            node = self.get_counterterm(
                substructure,
                reduced_process,
                momenta_dict_so_far
            )
            current_structure = node.current['singular_structure']
            current_legs = current_structure.get_all_legs()
            current_leg_ns = frozenset(
                leg.n for leg in current_legs
            )
            # Replace collinear substructures with their parent
            if current_structure.name() == "C":
                # The parent has already been generated by recursion,
                # retrieve its number
                parent_index = momenta_dict_so_far.inv[current_leg_ns]
                # Build the SubtractionLeg that appears in the current arguments
                parent_PDGs = self.parent_PDGs_from_legs(current_legs)
                assert len(parent_PDGs) == 1
                parent_PDG = parent_PDGs[0]
                parent_state = SubtractionLeg.FINAL
                if current_legs.has_initial_state_leg():
                    parent_state = SubtractionLeg.INITIAL
                current_args.add(
                    SubtractionLeg(parent_index, parent_PDG, parent_state)
                )
                # Eliminate soft sub-nodes without losing their children
                for subnode in node.nodes:
                    if subnode.current['singular_structure'].name() == "S":
                        node.nodes += subnode.nodes
                        node.nodes.remove(subnode)
            # Replace soft structures with their flattened versions
            elif current_structure.name() == "S":
                current_args.add(current_structure)
            # Other structures need to be implemented
            else:
                raise MadGraph5Error(
                    "Unrecognized current of type %s" %
                    str(type(current_structure))
                )
            # Add this node
            nodes.append(node)

        # If this is the outermost level,
        # the recursion was all that needed to be done
        if type(structure) is SingularStructure:
            for subnode in nodes:
                subnode.current['resolve_mother_spin_and_color'] = True
            return Counterterm(
                process=reduced_process,
                nodes=nodes,
                momenta_dict=momenta_dict_so_far,
                resolved_process=process,
                complete_singular_structure=structure
            )

        # 3. Else build the current and update
        #    the reduced process as well as the dictionary
        current_type = type(structure)(*current_args)
        current = Current({
            'singular_structure': current_type })
        structure_legs = current_type.get_all_legs()
        structure_leg_ns = frozenset(leg.n for leg in structure_legs)
        parent = None
        if structure.name() == "C":
            # Add entry to dictionary
            parent_index = len(momenta_dict_so_far) + 1
            momenta_dict_so_far[parent_index] = structure_leg_ns
            # Work out the complete SubtractionLeg for the parent
            parent_PDGs = self.parent_PDGs_from_legs(structure_legs)
            assert len(parent_PDGs) == 1
            parent_PDG = parent_PDGs[0]
            parent_state = SubtractionLeg.FINAL
            if structure_legs.has_initial_state_leg():
                parent_state = SubtractionLeg.INITIAL
            parent = SubtractionLeg(parent_index, parent_PDG, parent_state)
        elif structure.name() == "S":
            # No parent
            pass
        else:
            raise MadGraph5Error(
                "Building unrecognized current of type %s" %
                str(type(structure)) )
        # Remove legs of this structure
        legs_to_remove = []
        for leg in reduced_process['legs']:
            if leg['number'] in structure_leg_ns:
                legs_to_remove.append(leg)
        for leg in legs_to_remove:
            reduced_process['legs'].remove(leg)
        # Add parent of this structure
        if parent:
            reduced_process['legs'].append(
                base_objects.Leg({
                    'number': parent.n,
                    'id': parent.pdg,
                    'state': parent.state }) )
        # Sort preserving the initial state order
        rp_legs = reduced_process['legs']
        rp_legs.sort(key=lambda x: (x['state'], x['number']))
        if rp_legs[0]['number'] == 2:
            rp_legs[0], rp_legs[1] = rp_legs[1], rp_legs[0]

        # Finally return the counterterm node
        return CountertermNode(current=current, nodes=nodes)

    def get_integrated_counterterm(self, local_counterterm):
        """Given a local subtraction counterterm, return the corresponding integrated one.
        It has the same attributes, except that it contains only a single
        IntegratedCurrent, with the complete SingularStructure in it."""

        # First check that the local counterterm is singular, because if not then we
        # should of course not return any integrated counterterm.
        if not local_counterterm.is_singular():
            return None

        complete_singular_structure = local_counterterm. \
            reconstruct_complete_singular_structure()

        reduced_process = local_counterterm.process.get_copy(
            ['legs', 'n_loops', 'legs_with_decays'] )

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

        ########
        # TODO
        #
        # Modify the reduced process leg numbers so as to follow the order of appearance in the list of legs
        # (because this is how the PS point will be provided in the virtual contribution).
        #
        # Modify the bi-dictionary so as to match the modification of the leg numbers above in the reduced process.
        # What leg numbers are chosen for the unresolved legs is irrelevant, as long as I can deduce the parent leg
        # number from the singular structure and the bidictionary momentum map
        #
        ########
        
        
        integrated_current = IntegratedCurrent({
            'n_loops': n_loops,
            'squared_orders': squared_orders,
            'resolve_mother_spin_and_color': True,
            'singular_structure': complete_singular_structure })

        return IntegratedCounterterm(
            process=reduced_process,
            nodes=[CountertermNode(current=integrated_current), ],
            momenta_dict=bidict(local_counterterm.momenta_dict),
            prefactor=-1. * local_counterterm.prefactor )

    @staticmethod
    def get_all_currents(counterterms):
        """Deduce the list of currents needed to compute all counterterms given."""

        all_currents = []
        for counterterm in counterterms:
            for current in counterterm.get_all_currents():
                # Remove duplicates already at this level
                if current not in all_currents:
                    all_currents.append(current)
        return all_currents

    def get_all_counterterms(self, process):
        """Generate all counterterms for the corrections specified in this module
        and the process given in argument."""

        elementary_operators = self.get_all_elementary_operators(process)

        combinations = self.get_all_combinations(elementary_operators)
        all_counterterms = []
        all_integrated_counterterms = []
        for combination in combinations:
            template_counterterm = self.get_counterterm(combination, process)
            template_integrated_counterterm = \
                                   self.get_integrated_counterterm(template_counterterm)
            counterterms_with_loops = template_counterterm.split_loops(process['n_loops'])
            # TODO
            # For the time being, split_loops is given None instead of the squared orders
            # because they should be retrieved from the process by looking at individual
            # matrix elements
            # That is, every process has a list of possible coupling orders assignations
            # so we should loop over them
            for counterterm_with_loops in counterterms_with_loops:
                all_counterterms.extend( counterterm_with_loops.split_orders(None) )
            # Now also distribute the template integrated counterterm if it is not None
            if not template_integrated_counterterm is None:
                integrated_counterterms_with_loops = template_integrated_counterterm.split_loops(
                    process['n_loops'] )
                for integrated_counterterm_with_loops in integrated_counterterms_with_loops:
                    all_integrated_counterterms.extend(
                        integrated_counterterm_with_loops.split_orders(None) )

        return all_counterterms, all_integrated_counterterms

#=========================================================================================
# Subtraction current exporter
#=========================================================================================

class SubtractionCurrentExporter(object):
    """Class for mapping and exporting the subtraction currents to a given location
    and generate the corresponding accessors as well.
    """
    
    template_dir = pjoin(MG5DIR,'madgraph','iolibs','template_files','subtraction')
    template_modules_path = 'madgraph.iolibs.template_files.subtraction'
    
    # The main module name is not meant to be changed. If changed, then beware that the
    # import statements of the implementation of the subtraction currents must be 
    # updated accordingly
    main_module_name = 'SubtractionCurrents'

    def __init__(self, model, export_dir=None):
        """Initialize the exporter with a model and target export directory."""

        self.model      = model
        self.export_dir = export_dir

    def collect_modules(self, modules_path=[]):
        """Return a list of modules to load, from a given starting location."""
        
        base_path = pjoin(self.template_dir, pjoin(*modules_path))
        collected_modules = []

        for python_file in misc.glob(pjoin(base_path, '*.py')):
            python_file_name = os.path.basename(python_file)
            if python_file_name == '__init__.py':
                continue
            relative_module_path = '%s.%s'%('.'.join(modules_path), python_file_name[:-3])
            absolute_module_path = '%s.%s'%(self.template_modules_path, relative_module_path)
            collected_modules.append(
                (relative_module_path, importlib.import_module(absolute_module_path)) )

        for dir_name in os.listdir(base_path):
            if os.path.isdir(dir_name) and os.path.isfile(pjoin(dir_name, '__init__.py')):
                collected_modules.extend(self.collect_modules(modules_path+[dir_name, ]))
        
        return collected_modules

    def export(self, currents):
        """Export the specified list of currents and return a list of accessors
        which contain the mapping information.
        """
        
        subtraction_utils_module_path = '%s.%s'%(
            self.template_modules_path,'subtraction_current_implementations_utils' )
        subtraction_utils = importlib.import_module(subtraction_utils_module_path)
        
        if not self.export_dir is None:
            # First copy the base files to export_dir
            if not os.path.isdir(pjoin(self.export_dir, self.main_module_name)):
                os.mkdir(pjoin(self.export_dir, self.main_module_name))

            for file in ['__init__.py','subtraction_current_implementations_utils.py']:
                if not os.path.isfile(
                    pjoin(self.export_dir, self.main_module_name, file)):
                    cp(
                        pjoin(self.template_dir,file),
                        pjoin(self.export_dir,self.main_module_name, file))

        # Now load all modules specified in the templates
        # and identify the current implementation classes
        all_classes = []
        for dir_name in os.listdir(self.template_dir):
            if not os.path.isfile(pjoin(self.template_dir, dir_name, '__init__.py')):
                continue
            all_modules = self.collect_modules([dir_name])
            for (module_path, module) in all_modules:
                for class_name in dir(module):
                    implementation_class = getattr(module, class_name)
                    if (inspect.isclass(implementation_class) and
                        hasattr(implementation_class, 'does_implement_this_current') ):
                        all_classes.append(
                            (dir_name, module_path, class_name, implementation_class) )
        
        # If there is a "default current" in the subtraction_current_implementations_utils class,
        # presumably used for debugging only, then add this one at the very end of all_classes so that
        # it will be selected only if no other class matches.
        if hasattr(subtraction_utils,'DefaultCurrentImplementation'):
            default_implementation_class = getattr(
                subtraction_utils, 'DefaultCurrentImplementation' )
            if (inspect.isclass(default_implementation_class) and
                hasattr(default_implementation_class, 'does_implement_this_current') ):
                all_classes.append((
                    '', 'subtraction_current_implementations_utils',
                    'DefaultCurrentImplementation', default_implementation_class ) )
        
        # Save directories that will need to be exported
        directories_to_export = set([])
        
        # Group all currents according to mapping rules.
        # The mapped currents is a dictionary of the form
        #         { (module_path, class_name, instantiation_options_index), 
        #                {'defining_current': <...>,
        #                 'mapped_process_keys': [<...>],
        #                 'instantiation_options': instantiation_options
        #         }
        mapped_currents = {}
        
        # Look in all current implementation classes found
        # and find which one implements each current
        all_instantiation_options = []
        currents_with_default_implementation = []
        for current in currents:
            found_current_class = False
            for (dir_name, module_path, class_name, implementation_class) in all_classes:
                instantiation_options = implementation_class.does_implement_this_current(
                    current, self.model )
                if instantiation_options is None:
                    continue
                try:
                    instantiation_options_index = all_instantiation_options.index(
                        instantiation_options )
                except ValueError:
                    all_instantiation_options.append(instantiation_options)
                    instantiation_options_index = len(all_instantiation_options)-1
                # If no export is asked for,
                # load directly from the subtraction template directory
                if not self.export_dir is None:
                    main_module = self.main_module_name
                else:
                    main_module = 'subtraction'                    
                key = (
                    '%s.%s'%(main_module, module_path),
                    class_name, instantiation_options_index )
                if key in mapped_currents:
                    mapped_currents[key]['mapped_process_keys'].append(current.get_key())
                else:
                    if dir_name != '':
                        directories_to_export.add(dir_name)
                    mapped_currents[key]={'defining_current': current,
                                          'mapped_process_keys': [current.get_key()],
                                          'instantiation_options': instantiation_options}
                if class_name == 'DefaultCurrentImplementation':
                    currents_with_default_implementation.append(current)
                found_current_class = True
                break
            if not found_current_class:
                raise MadGraph5Error("No implementation was found for current %s."%str(current))
        
        # Warn the user whenever DefaultCurrentImplementation is used
        # (it should never be used in production)
        if currents_with_default_implementation:
            currents_str = '\n'.join(
                ' > %s' % str(crt)
                for crt in currents_with_default_implementation )
            msg = """No implementation was found for the following subtraction currents:
                %s
                The class 'DefaultCurrentImplementation' will therefore be used for it
                but results obtained in this way are very likely wrong 
                and should be used for debugging only.""" % currents_str
            if __debug__:
                logger.critical(msg)
            else:
                raise MadGraph5Error(msg)
        
        # Now copy all the relevant directories
        if not self.export_dir is None:
            for directory_to_export in directories_to_export:
                dir_path = pjoin(self.export_dir, self.main_module_name, directory_to_export)
                def ignore_function(d, files):
                    return [
                        f for f in files
                        if (
                            os.path.isfile(pjoin(d, f)) and
                            f.split('.')[-1] in ['pyc','pyo','swp'] ) ]
                if not os.path.isdir(dir_path):
                    shutil.copytree(
                        pjoin(self.template_dir, directory_to_export), dir_path,
                        ignore=ignore_function )
        
        return mapped_currents
        
#=========================================================================================
# Standalone main for debugging / standalone trials
#=========================================================================================

if __name__ == '__main__':
    misc.sprint("Put your standalone subtraction code here.")
