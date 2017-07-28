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

    def __new__(self, *args, **opts):
        """Initialize a SubtractionLeg object from various representations."""

        target = (0, 0, SubtractionLeg.FINAL)
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
                hasattr(leg, '__len__') and
                len(leg) == 3 and isinstance(leg[0], int) and
                isinstance(leg[1], int) and isinstance(leg[2], bool)
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
                (len(args) == 3) and isinstance(args[0], int) and
                isinstance(args[1], int) and isinstance(args[2], bool)
        ):
            target = (args)
        # Invalid number of arguments
        else:
            raise MadGraph5Error(
                    "SubtractionLeg cannot be initialized with argument %s."
                    % args
            )

        return super(SubtractionLeg, self).__new__(self, target)

    def n(self):
        return self[0]

    def pdg(self):
        return self[1]

    def state(self):
        return self[2]

#===============================================================================
# SubtractionLegSet
#===============================================================================
class SubtractionLegSet(frozenset):
    """Set of SubtractionLeg objects."""

    def __new__(self, *legs):
        """Initialize set, trying to convert arguments into SubtractionLeg's."""

        return super(SubtractionLegSet, self).__new__(
                self,
                frozenset(SubtractionLeg(leg) for leg in legs)
        )

#===============================================================================
# SingularStructure
#===============================================================================
class SingularStructure(object):
    """Object that represents a hierarchical structure of IR singularities."""

    def __init__(self, *args, **opts):
        """Initialize a hierarchical singular structure."""

        # Type check
        for arg in args:
            if not isinstance(arg, (SingularStructure, SubtractionLeg)):
                raise MadGraph5Error(
                        "SubtractionOperator initialized with invalid argument "
                        "%s of type %s." % (arg, str(type(arg)))
                )

        # List of substructures that this SubtractionOperators acts on
        self.substructures = [
            arg for arg in args if isinstance(arg, SingularStructure)
        ]
        # Set of simple legs this SubtractionOperators acts on
        self.legs = SubtractionLegSet(
            # HACK without tuple and expansion fails for empty args
            *(arg for arg in args if isinstance(arg, SubtractionLeg))
        )

    def __str__(self):
        """Return a string representation of the singular structure."""
        tmp_str = self.name() + "("
        tmp_str += ",".join(sorted(str(sub) for sub in self.substructures))
        if self.substructures:
            tmp_str += ","
        tmp_str += ",".join(sorted(str(leg.n()) for leg in self.legs))
        tmp_str += ")"
        return tmp_str

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

    def __str__(self):
        """Return a simple string representation of the singular operator."""
        return self.name() + str(sorted((leg.n() for leg in self)))

    def name(self):
        """Symbol used to represent this operator within output."""
        raise MadGraph5Error(
            "name called in SingularOperator of unspecified type."
        )

    def structure(self):
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
            structure.substructures.append(self.structure())
            return structure

        # If the limit acts at least partly at the current level
        if not self.isdisjoint(structure.legs):
            # If the limit acts at different levels,
            # it is not needed for subtraction
            if not self.issubset(structure.legs):
                structure = None
                return structure
            # If the limit acts completely within this level,
            # it may be needed or not
            else:
                if self.act_here_needed(structure):
                    structure.substructures.append(self.structure())
                    structure.legs = structure.legs.difference(self)
                    return structure
                else:
                    structure = None
                    return structure

        # The limit acts deeper
        for substructure in structure.substructures:
            # Find the first substructure hit by the limit, even partially
            if not self.non_overlapping_with(substructure):
                # Act on that substructure and propagate None if needed
                if self.act_on(substructure) is None:
                    structure = None
                return structure

class SoftOperator(SingularOperator):
    """Object that represents a soft elementary singular operator."""

    def name(self):
        return "S"

    def structure(self):
        return SoftStructure(*self)

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

    def structure(self):
        return CollStructure(*self)

    def act_here_needed(self, structure):

        assert isinstance(structure, SingularStructure)

        # Collinear limits are not needed within soft structures
        if isinstance(structure, SoftStructure):
            return False

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
        if (structure is None) or (not self):
            return structure
        # Apply last operator and advance recursion
        most = SingularOperatorList(self[:-1])
        last = self[-1]
        return most.act_on(last.act_on(structure))

    def simplify(self):
        """Simplify a list of operators,
        returning the real structure of the corresponding singularities.
        """

        return self.act_on(SingularStructure())

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
        flavored_legs = [leg for leg in legs if leg.get('id') != 21]

        # If all daughters were gluons, the only parent is a gluon
        if not flavored_legs:
            return [21]

        # Consider last particle
        last_leg = flavored_legs.pop()
        ll_state = last_leg.get('state')
        ll_id    = last_leg.get('id')
        # Look for a corresponding anti-particle
        for leg in range(len(flavored_legs)):
            cur_state = flavored_legs[leg].get('state')
            cur_id    = flavored_legs[leg].get('id')
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
        
        # Eliminate particles that do not take place in the subtraction
        legs = [
            l for l in process.get('legs')
            if self.can_be_IR_unresolved(l.get('id'))
        ]
        fs_legs = [l for l in legs if l.get('state') == SubtractionLeg.FINAL]
        is_legs = [l for l in legs if l.get('state') == SubtractionLeg.INITIAL]

        # Loop over number of unresolved particles
        for unresolved in range(1, self.correction_order.count('N')+1):
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
                    coll_set = coll_initial_set + (coll_initial_leg, )
                    if self.can_become_collinear(coll_set):
                        elementary_operator_list.append(CollOperator(*coll_set))

        return SingularOperatorList(elementary_operator_list)

    def get_all_combos(self, elementary_operators):
        """Determine all combinations of elementary operators."""

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

    # # TODO Fix this
    # def reduced_me(unresolved_structure, process):
    #     """
    #     Get the reduced matrix element after the action of operators
    #     :param unresolved_structure: Simplified structure of unresolved particles
    #     :param process: Particle PDGs for the considered matrix element
    #     :return: List of momenta and corresponding particle PDGs
    #     """
    #     # TODO add spin and color correlations
    #     all_unres = [unresolved_particles(group) for group in
    #                  unresolved_structure]
    #     new_enum = set(enumerate(process))
    #     # For every group of unresolved particles
    #     for iur in range(len(all_unres)):
    #         # Erase all unresolved particles from the reduced process
    #         for erase in all_unres[iur]:
    #             new_enum = set(
    #                     [prtcl for prtcl in new_enum if prtcl[0] != erase])
    #         # Add back the parent if the group was a collinear one
    #         if unresolved_structure[iur][0] == coll:
    #             tmpparent = unresolved_flavor(all_unres[iur], process)
    #             new_enum.add((frozenset(all_unres[iur]), tmpparent))
    #     return new_enum

    # # TODO Fix this
    # def current(unresolved_group, process):
    #     """
    #     Get the current for unresolved particles
    #     :param unresolved_group: Simplified structure of unresolved particles
    #     :param process: Flavors for the considered matrix element
    #     :return: To be decided
    #     """
    #     # TODO implement currents
    #     buffer = []
    #     subgroups = set()
    #     for subgroup in unresolved_group[1]:
    #         # Particle is not part of a nested limit
    #         if isinstance(subgroup, int):
    #             subgroups.add((subgroup, process[subgroup]))
    #         # Soft subgroup in a collinear limit
    #         elif unresolved_group[0] == coll and subgroup[0] == soft:
    #             # Soft groups can contain no further subgroups
    #             subgroups.add(current(subgroup, process)[0])
    #         # Iterated limit decomposition
    #         else:
    #             buffer += current(subgroup, process)
    #             ups = unresolved_particles(subgroup)
    #             if len(ups) == 1:
    #                 subgroups.add((
    #                     tuple(ups)[0],
    #                     unresolved_flavor(ups, process)
    #                 ))
    #             else:
    #                 subgroups.add((
    #                     frozenset(ups),
    #                     unresolved_flavor(ups, process)
    #                 ))
    #     buffer = [(unresolved_group[0], frozenset(subgroups)), ] + buffer
    #     return buffer
    #
    # # TODO Fix this
    # def currents(unresolved_structures, process):
    #     dbllist = [current(cur, process) for cur in unresolved_structures]
    #     return list(itertools.chain.from_iterable(dbllist))

# ===============================================================================
# Standalone main for debugging / standalone trials
# ===============================================================================
if __name__ == '__main__':
    misc.sprint("Put your standalone subtraction code here.")
