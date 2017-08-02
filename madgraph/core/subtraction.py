#!/usr/bin/env python
################################################################################
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
################################################################################
"""Definition of all the classes and features relevant to the handling of 
higher order IR subtraction.
"""

import copy
import itertools
import collections
import logging
import math
import sys
import os
if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(
        os.path.realpath(__file__)),
        os.path.pardir,
        os.path.pardir
    ))
from madgraph import MadGraph5Error, MG5DIR, InvalidCmd
import madgraph.various.misc as misc 
import madgraph.core.base_objects as base_objects
import madgraph.fks.fks_common as fks_common

logger = logging.getLogger('madgraph')
pjoin = os.path.join

#===============================================================================
# SubtractionLeg
#===============================================================================
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

#===============================================================================
# SubtractionLegSet
#===============================================================================
class SubtractionLegSet(frozenset):
    """Set of SubtractionLeg objects."""

    def __new__(cls, *args, **opts):
        """Initialize set, trying to convert arguments into SubtractionLeg's."""    
        
        if not args:
            return super(SubtractionLegSet, cls).__new__(cls)

        if isinstance(args[0], collections.Iterable) and not isinstance(args[0],(dict, SubtractionLeg)):
            legs = args[0]
        else:
            legs = args

        return super(SubtractionLegSet, cls).__new__(cls, [SubtractionLeg(leg) for leg in legs])


    def __init__(self, *args, **opts):
        """Initialize set, trying to convert arguments into SubtractionLeg's."""

        if not args:
            super(SubtractionLegSet, self).__init__()
            return

        if isinstance(args[0], collections.Iterable) and not isinstance(args[0],(dict, SubtractionLeg)):
            legs = args[0]
        else:
            legs = args

        super(SubtractionLegSet, self).__init__([SubtractionLeg(leg) for leg in legs])
        return

    def has_initial_state_leg(self):
        """Returns True if this leg set has at least one initial state leg,
        False otherwise.
        """

        for leg in self:
            if leg.state == SubtractionLeg.INITIAL:
                return True
        return False

#===============================================================================
# SingularStructure
#===============================================================================
class SingularStructure(object):
    """Object that represents a hierarchical structure of IR singularities."""

    def __init__(self, *args, **opts):
        """Initialize a hierarchical singular structure."""

        if (
            args and
            isinstance(args[0], collections.Iterable) and
            not isinstance(args[0], (dict, SubtractionLeg) )
        ):
            components = args[0]
        else:
            components = args

        # If this Structure is annihilated, any operator acting on it would
        # not apply and this structure evaluated to False.
        self.is_annihilated = False

        # Type check
        for a in components:
            if not isinstance(a, (SingularStructure, SubtractionLeg)):
                raise MadGraph5Error(
                        "SubtractionOperator initialized with invalid argument "
                        "%s of type %s." % (a, str(type(a)))
                )

        # List of substructures that this SubtractionOperators acts on
        self.substructures = [
            a for a in components if isinstance(a, SingularStructure)
        ]
        # Set of simple legs this SubtractionOperators acts on
        self.legs = SubtractionLegSet(
            # HACK without tuple and expansion fails for empty components
            *(a for a in components if isinstance(a, SubtractionLeg))
        )

    def __str__(self):
        """Return a string representation of the singular structure."""

        if self.is_annihilated:
            return 'NULL'

        tmp_str = self.name() + "("
        tmp_str += ",".join(sorted(str(sub) for sub in self.substructures))
        if self.substructures:
            tmp_str += ","
        tmp_str += ",".join(sorted(str(leg.n) for leg in self.legs))
        tmp_str += ")"
        return tmp_str

    def is_void(self):

        return self.is_annihilated

    def annihilate(self):
        """When an operator cannot act on this structure,
        remove this structure altogether
        by flagging it as unapplicable to any other operator.
        """

        # Really clear structures and legs
        del self.substructures[:]
        self.legs = SubtractionLegSet()
        self.is_annihilated = True

    def get_all_legs(self):
        """Return all legs involved in this singular structure."""

        all_legs = set(sub.get_all_legs() for sub in self.substructures)
        all_legs.add(self.legs)
        return SubtractionLegSet().union(*all_legs)

    def discard_leg_numbers(self):
        """Set all leg numbers to zero, discarding this information forever."""

        self.legs = SubtractionLegSet(
                SubtractionLeg(0, leg.pdg, leg.state)
                for leg in self.legs
        )
        for substructure in self.substructures:
            substructure.discard_leg_numbers()
        return

    def name(self):

        return ""

class SoftStructure(SingularStructure):

    def name(self):

        return "S"

class CollStructure(SingularStructure):

    def name(self):

        return "C"

#===============================================================================
# SingularOperator
#===============================================================================
class SingularOperator(SubtractionLegSet):
    """Virtual base class for elementary singular operators."""

    def __init__(self, *args, **opts):

        super(SingularOperator, self).__init__(*args, **opts)

    def __str__(self):
        """Return a simple string representation of the singular operator."""

        return self.name() + str(sorted((leg.n for leg in self)))

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
            "structure called in SingularOperator of unspecified type."
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

        if not self.isdisjoint(structure.legs):
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
        if not self.isdisjoint(structure.legs):
            # If the limit acts at different levels,
            # it is not needed for subtraction
            if not self.issubset(structure.legs):
                structure.annihilate()
                return
            # If the limit acts completely within this level,
            # it may be needed or not
            else:
                if self.act_here_needed(structure):
                    structure.substructures.append(self.get_structure())
                    structure.legs = structure.legs.difference(self)
                    return
                else:
                    structure.annihilate()
                    return

        # The limit acts deeper
        for substructure in structure.substructures:
            # Find the first substructure hit by the limit, even partially
            if not self.non_overlapping_with(substructure):

                self.act_on(substructure)
                # If the action on the substructure annihilated it, then annihilate this 
                # structure as well
                if substructure.is_void():
                    structure.annihilate()
                return

class SoftOperator(SingularOperator):
    """Object that represents a soft elementary singular operator."""

    def __init__(self, *args, **opts):

        super(SoftOperator, self).__init__(*args, **opts)
        return

    def name(self):

        return "S"

    def get_structure(self):

        return SoftStructure(self)

    def act_here_needed(self, structure):

        assert isinstance(structure, SingularStructure)

        # Soft limits are not needed within soft structures
        if isinstance(structure, SoftStructure):
            return False

        # Soft limits are needed within collinear structures
        # if there is at least one hard particle left
        if isinstance(structure, CollStructure):
            # If there is a hard particle at this level, the limit is needed
            if len(self) < len(structure.legs):
                return True
            # If there is a hard collinear substructure, the limit is needed
            else:
                for substructure in structure.substructures:
                    if not isinstance(substructure, SoftStructure):
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

        return CollStructure(self)

    def __init__(self, *args, **opts):

        super(CollOperator, self).__init__(*args, **opts)
        return

    def act_here_needed(self, structure):

        assert isinstance(structure, SingularStructure)

        # Collinear limits are needed within soft structures
        # WARNING: this statement has to be checked
        # at NNLO only S(C(i,j)), at N3LO we first get S(i,C(j,k))
        if isinstance(structure, SoftStructure):
            return True

        # Collinear limits are always needed within collinear structures
        if isinstance(structure, CollStructure):
            return True

        # By default, not skipping overlaps ensures correct subtraction
        return True

#===============================================================================
# SingularOperatorList
#===============================================================================
class SingularOperatorList(list):
    """List of subtraction operators."""

    def __str__(self):
        """Print this singular operator list."""

        return "[" + ", ".join([str(leg) for leg in self]) + "]"
    
    def act_on(self, structure):
        """Act with a list of operators on a substructure."""

        # Empty list of operators to apply or invalid structure: done
        if structure.is_void() or not self:
            return structure
        # Apply last operator and advance recursion
        all_but_last_operators = SingularOperatorList(self[:-1])
        # Act using the last operator
        self[-1].act_on(structure)
        # And now act using the rest of the operators, if needed
        if not structure.is_void():
            return all_but_last_operators.act_on(structure)

    def simplify(self):
        """Simplify a list of operators,
        returning the real structure of the corresponding singularities.
        """

        structure = SingularStructure()
        self.act_on(structure)
        return structure

#===============================================================================
# Current
#===============================================================================
class Current(base_objects.Process):

    def default_setup(self):
        """Default values for all properties specific to current,
        additional to Process.
        """

        super(Current, self).default_setup()
        self['singular_structure'] = SingularStructure()
        self['parent_subtraction_leg'] = None

        return

    def discard_leg_numbers(self):
        """Discard all leg numbers in the singular_structure
        and in the parent_subtraction_leg,
        effectively getting the type of current
        without information about momenta for which it should be evaluated.
        """

        self['singular_structure'].discard_leg_numbers()
        self['parent_subtraction_leg'] = SubtractionLeg(
                0,
                self['parent_subtraction_leg'].pdg,
                self['parent_subtraction_leg'].state
        )
        return

#===============================================================================
# IRSubtraction
#===============================================================================
class IRSubtraction(object):
    
    def __init__(self, model, correction_order='NLO', correction_types=['QCD']):
        """Initialize a IR subtractions for a given model,
        correction order and type.
        """
        
        self.model = model
        self.correction_order = correction_order
        self.correction_types = correction_types
        # Map perturbed coupling orders to the corresponding relevant interactions and particles.
        # The values of the dictionary are 'interactions', 'pert_particles' and 'soft_particles'
        self.IR_quantities_for_corrections_types = dict(
            (order, fks_common.find_pert_particles_interactions(self.model, pert_order = order))
            for order in correction_types)

    def can_be_IR_unresolved(self, PDG):
        """Check whether a particle given by its PDG can become unresolved
        and lead to singular behavior.
        """

        return any(
            (PDG in self.IR_quantities_for_corrections_types[order]['pert_particles'])
            for order in self.correction_types
        )

    def parent_PDGs(self, legs):
        """List all possible parent PDGs to a given set of legs."""

        # BALDY: parent_PDGs hardcoded to SM QCD (with some sanity checks)
        # A generic way to build histories of possible combination of external legs is using the 
        # model ref_dict_to1
        # misc.sprint(self.model.get('ref_dict_to1')[(-2,2)])
        # misc.sprint(self.model.get_particle(-2).get('spin'),
        #             self.model.get_particle(-2).get_color(),
        #             self.model.get_particle(-2).get('mass')=='zero'
        #             )

        if not self.model.get('name').startswith('sm'):
            raise InvalidCmd(
                "The function parent_PDGs is implemented for the SM only."
            )
        if any(order!='QCD' for order in self.correction_types):
            raise InvalidCmd(
                "The function parent_PDGs is implemented for QCD only."
            )

        # Get parton flavors, eliminating gluons
        flavored_legs = [leg for leg in legs if leg.pdg != 21]

        # If all daughters were gluons, the only parent is a gluon
        if not flavored_legs:
            return [21]

        # Consider last particle
        last_leg = flavored_legs.pop()
        ll_state = last_leg.state
        ll_id    = last_leg.pdg
        # Look for a corresponding anti-particle
        for leg in range(len(flavored_legs)):
            cur_state = flavored_legs[leg].state
            cur_id    = flavored_legs[leg].pdg
            if (
                        (cur_state == ll_state and cur_id == -ll_id)
                        or
                        (cur_state != ll_state and cur_id == ll_id)
            ):
                # Eliminate it and start over
                flavored_legs.pop(leg)
                return self.parent_PDGs(flavored_legs)

        # If there was no anti-particle,
        # check if all other legs have been emitted by the last one
        if self.parent_PDGs(flavored_legs) == [21]:
            # Return the 'initial-state PDG' of the particle
            if ll_state == SubtractionLeg.INITIAL:
                return [ll_id]
            else:
                return [-ll_id]

        # At this point, there is no valid parent: return empty list
        return []
        
    def can_become_soft(self, legs):
        """Check whether a bunch of legs going simultaneously soft 
        lead to singular behavior.
        """

        for pdg in self.parent_PDGs(legs):
            particle = self.model.get_particle(pdg)
            if (particle.get('spin')==3 and particle.get('mass')=='zero'):
                return True
        return False
    
    def can_become_collinear(self, legs):
        """Check whether a bunch of legs going collinear to each other 
        lead to singular behavior.
        """
        
        for pdg in self.parent_PDGs(legs):
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
        is_legs = SubtractionLegSet(legs.difference(fs_legs))

        # Loop over number of unresolved particles
        for unresolved in range(1, self.correction_order.count('N')+1):
            # Get iterators at the start of the final-state list
            it = iter(fs_legs)
            soft_it, coll_final_it, coll_initial_it = itertools.tee(it, 3)
            # Final state particle sets going soft
            for soft_set in itertools.combinations(soft_it, unresolved):
                if self.can_become_soft(soft_set):
                    elementary_operator_list.append(SoftOperator(soft_set))
            # Final state particle sets going collinear
            for coll_final_set in itertools.combinations(coll_final_it, unresolved + 1):
                if self.can_become_collinear(coll_final_set):
                    elementary_operator_list.append(CollOperator(coll_final_set))
            # Initial-final collinear
            # For any final-state set with one less particle
            for coll_initial_set in itertools.combinations(coll_initial_it, unresolved):
                # Combine with all initial legs
                for coll_initial_leg in is_legs:
                    coll_set = coll_initial_set + (coll_initial_leg, )
                    if self.can_become_collinear(coll_set):
                        elementary_operator_list.append(CollOperator(coll_set))

        return SingularOperatorList(elementary_operator_list)

    def get_all_raw_combinations(self, elementary_operators):
        """Determine all combinations of elementary operators,
        without applying any simplification.
        """

        combos = []
        for nop in range(len(elementary_operators) + 1):
            combos += [
                SingularOperatorList(combo)
                for combo in itertools.combinations(
                        iter(elementary_operators),
                        nop
                )
            ]
        return combos

    def get_all_combinations(self, elementary_operators):
        """Determine all combinations of elementary operators,
        applying simplification and discarding the ones that vanish.
        """

        return [
            simple_combo
            for simple_combo in [
                combo.simplify()
                for combo in self.get_all_raw_combinations(elementary_operators)
                ]
            if not simple_combo.is_void()
        ]

    def count_unresolved(self, structure):
        """Count the number of unresolved particles in some structure."""

        assert(type(structure) is SingularStructure)
        number_of_unresolved_particles = 0
        for sub in structure.substructures:
            if isinstance(sub, SoftStructure):
                number_of_unresolved_particles += len(sub.get_all_legs())
            elif isinstance(sub, CollStructure):
                number_of_unresolved_particles += len(sub.get_all_legs()) - 1
            else:
                raise MadGraph5Error(
                        "Unrecognised structure of type %s"
                        "in IRSubtraction.count_unresolved" %
                        str(type(structure))
                )
        return number_of_unresolved_particles

    def filter_combinations(self, combinations):
        """Filter out combinations with too many unresolved particles
        for the order that has been set.
        """

        max_unresolved = self.correction_order.count('N')
        return [
            combo
            for combo in combinations
            if self.count_unresolved(combo) <= max_unresolved
        ]

    def get_elementary_currents(self, structure):
        """Convert this structure to a list of elementary currents
        at the present level, without the corresponding momenta.
        """

        print "get_elementary_currents called with argument", structure

        currents = []

        # Handle lists of structures as a special case
        if not isinstance(structure, SingularStructure):
            for real_structure in structure:
                for current in self.get_elementary_currents(real_structure):
                    if current not in currents:
                        currents.append(current)
            return currents
        # Treat non-overlapping structures separately, avoiding recursion
        if type(structure) is SingularStructure:
            for real_structure in structure.substructures:
                for current in self.get_elementary_currents(real_structure):
                    if current not in currents:
                        currents.append(current)
            return currents

        # 1. Build a list of possible entries for any position in the current

        # All legs have to be entries, add them discarding particle numbers
        possible_entries_with_parent = [
            [(SubtractionLeg(0, leg.pdg, leg.state), )*2, ]
            for leg in structure.legs
        ]
        # Loop over substructures
        for substructure in structure.substructures:
            sub_currents = self.get_elementary_currents(substructure)
            print "sub_currents are", sub_currents
            # Recursively add sub_currents as needed currents
            for sub_current in sub_currents:
                if sub_current not in currents:
                    currents.append(sub_current)
            if isinstance(substructure, SoftStructure):
                # If it is a soft within something else (should be collinear),
                # preserve the structure because C(i,...,S(j,...))
                # is an elementary current with an expression of its own
                # Note: take the singular_structure from the current,
                # so that internal simplifications like
                #   S(C(i,j,...)) -> S(parent of i,j,...)
                # have already happened (and leg numbers were discarded)
                possible_entries_with_parent.append([
                    (
                        current['parent_subtraction_leg'],
                        current['singular_structure']
                    ) for current in sub_currents
                ])
            else:
                # else just add the parent of the substructure
                possible_entries_with_parent += [[
                    (
                        current['parent_subtraction_leg'],
                        current['parent_subtraction_leg']
                    ) for current in sub_currents
                ],]

        print "Possible entries with parent:", possible_entries_with_parent

        # 2. Expand it, creating all tuples where the i-th element
        #    comes from the possible_entries for position i

        entry_sets_with_parent = [[],]
        while possible_entries_with_parent:
            print "entry_sets_with_parent is", entry_sets_with_parent
            print "possible_entries_with_parent is", possible_entries_with_parent
            last_entries_with_parent = possible_entries_with_parent.pop()
            new_entry_sets_with_parent = []
            for entry_set_with_parent in entry_sets_with_parent:
                for last_entry_with_parent in last_entries_with_parent:
                    new_entry_set_with_parent = entry_set_with_parent + [last_entry_with_parent, ]
                    print "adding", new_entry_set_with_parent, "to entry sets"
                    new_entry_sets_with_parent.append(
                            new_entry_set_with_parent
                    )
            entry_sets_with_parent = new_entry_sets_with_parent

        print "Entry sets with parent:", entry_sets_with_parent

        # 3. Now create the current by adding the parent_subtraction_leg
        #    and converting the list of entries in the appropriate object

        for entry_set_with_parent in entry_sets_with_parent:
            daughters = SubtractionLegSet(
                entry_with_parent[0]
                for entry_with_parent in entry_set_with_parent
            )
            singular_structure = type(structure)(
                entry_with_parent[1]
                for entry_with_parent in entry_set_with_parent
            )
            parent_pdgs  = self.parent_PDGs(daughters)
            parent_state = SubtractionLeg.FINAL
            if daughters.has_initial_state_leg():
                parent_state = SubtractionLeg.INITIAL
            for parent_pdg in parent_pdgs:
                parent_subtraction_leg = SubtractionLeg(0, parent_pdg, parent_state)
                current = Current({
                    'parent_subtraction_leg': parent_subtraction_leg,
                    'singular_structure': singular_structure
                })
                if current not in currents:
                    currents.append(current)
                else:
                    print "current", current, "was already found"

        print "returning from get_elementary_currents with argument", structure
        return currents

    # # TODO Fix this
    # def reduced_me(self, structure, process):
    #     """
    #     Get the reduced matrix element after the action of operators
    #     :param unresolved_structure: Simplified structure of unresolved particles
    #     :param process: Particle PDGs for the considered matrix element
    #     :return: List of momenta and corresponding particle PDGs
    #     """
    #     # TODO add spin and color correlations
    #     all_unres = [unresolved_particles(group) for group in
    #                  structure]
    #     new_enum = set(enumerate(process))
    #     # For every group of unresolved particles
    #     for iur in range(structure.substructures):
    #         # Erase all unresolved particles from the reduced process
    #         for erase in all_unres[iur]:
    #             new_enum = set(
    #                     [prtcl for prtcl in new_enum if prtcl[0] != erase])
    #         # Add back the parent if the group was a collinear one
    #         if structure.substructures[iur] == coll:
    #             tmpparent = unresolved_flavor(all_unres[iur], process)
    #             new_enum.add((frozenset(all_unres[iur]), tmpparent))
    #     return new_enum

#===============================================================================
# Standalone main for debugging / standalone trials
#===============================================================================
if __name__ == '__main__':
    misc.sprint("Put your standalone subtraction code here.")
