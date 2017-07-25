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
""" Definition of all the classes and features relevant to the handling of 
higher order IR subtraction."""

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
class SubtractionLeg(base_objects.Leg):
    """ Leg object specialized for subtraction. """

    def __init__(self, *args, **opts):
        super(SubtractionLeg, self).__init__(*args, **opts)

#===============================================================================
# SubtractionLegSet
#===============================================================================
class SubtractionLegSet(set):
    """ Set of SubtractionLeg objects. """

    def __init__(self, *args, **opts):
        super(SubtractionLegSet, self).__init__(*args, **opts)

#===============================================================================
# SubtractionOperator 
#===============================================================================
class SubtractionOperator(object):
    """ Object that represents a generic subtraction operator
    or a hierarchical structure of IR singularities. """

    def __init__(self, *args, **opts):
        """ Initialize a subtraction operator with legs and operator sets. """
        
        # Set of simple legs this SubtractionOperators acts on
        self.legs = SubtractionLegSet(
            [a for a in args if isinstance(a, SubtractionLeg)]
        )
        # List of substructures that this SubtractionOperators acts on
        self.substructures = SubtractionOperatorList(
            [a for a in args if isinstance(a, SubtractionOperator)]
        )
    
    def is_elementary(self):
        """ Return whether this SubtractionOperator is elementary,
         in the sense that it is not an empty structure
         and it hasn't been applied on yet. """
        return (self.legs and not self.substructures)

    def act_on(self, structure):
        """ Act with an elementary operator on a non-elementary structure. """

        # Sanity check
        if (not self.is_elementary()) or structure.is_elementary():
            raise Exception("You can only apply elementary operators \
            to non-elementary structures.")

        # Count the common legs between the current operator
        # and the structure at the highest level
        overlap = len(self.legs.intersection(structure.legs))

        # No overlap, we must look into substructure to find a match
        if overlap == 0:
            for substructure in structure.substructures:
                # Check if the operator acts entirely on a deeper level
                # If so, act and return
                result = self.act_on(substructure)
                if not result is None:
                    return structure
            # If no match was found within substructures,
            # the operator has no legs in common with the whole structure
            # and we should simply add this operator to substructures.
            # In other terms, non-overlapping limits are additively combined.
            structure.substructures.append(self)

        # Partial overlaps are always set to zero
        elif overlap != len(self.legs):
            structure = None
        
        # Full overlap, combine the affected legs into a substructure
        else:
            structure.substructure.append(self)
            structure.legs = structure.legs.difference(self.legs)

        return structure

class SoftOperator(SubtractionOperator):
    def __init__(self, *args, **opts):
        """ Initialize a soft operator with all relevant arguments. """
        super(SoftOperator, self).__init__(self, *args, **opts)
        
class CollOperator(SubtractionOperator):
    def __init__(self, *args, **opts):
        """ Initialize a collinear operator with all relevant arguments. """
        super(CollOperator, self).__init__(self, *args, **opts) 

#===============================================================================
# SubtractionOperatorList
#===============================================================================
class SubtractionOperatorList(list):
    """ List of subtraction operators. """
    
    def __init__(self, *args, **opts):
        super(SubtractionOperatorList, self).__init__(*args, **opts)

#===============================================================================
# Standalone main for debugging / standalone trials 
#===============================================================================
if __name__ == '__main__':
    misc.sprint("Put your standalone subtraction code here.")

#===============================================================================
# IRSubtraction
#===============================================================================
class IRSubtraction(object):
    
    def __init__(self, model, correction_order='NLO', correction_types=['QCD']):
        """ Initialize a IR subtractions for a given model,
        correction order and type. """
        
        self.model = model
        self.correction_order = correction_order
        self.correction_types = correction_types
        # Map perturbed coupling orders to the corresponding relevant interactions and particles.
        # The values of the dictionary are 'interactions', 'pert_particles' and 'soft_particles'
        self.IR_quantities_for_corrections_types = dict(
            (order, fks_common.find_pert_particles_interactions(self.model, pert_order = order))
            for order in correction_types)
        
        pass

    def can_be_IR_unresolved(self, PDG):
        """ Checks whether a particle given by its PDG can become unresolved 
        and lead to singular behavior. """
        
        return any(
            (PDG in self.IR_quantities_for_corrections_types[order])
            for order in self.correction_types
        )

    def parent_PDGs(self, legs):
        """ List all possible parent PDGs to a given set of legs."""

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
        """ Check whether a bunch of legs going simultaneously soft 
        lead to singular behavior. """

        for pdg in self.parent_PDGs(legs):
            particle = self.model.get_particle(pdg)
            if (particle.get('spin')==3 and particle.get('mass')=='zero'):
                return True
        return False
    
    def can_become_collinear(self, legs):
        """ Check whether a bunch of legs going collinear to each other 
        lead to singular behavior. """
        
        for pdg in self.parent_PDGs(legs):
            if self.can_be_IR_unresolved(pdg):
                return True
        return False        
    
    def get_all_elementary_operators(self, process):
        """ Generate all the 'building blocks' operator relevant for the process
        passed in argument. """
        
        elementary_operator_list = SubtractionOperatorList([])
        
        # Eliminate particles that do not take place in the subtraction
        all_legs = process.get('legs')
        active_legs = [l for l in all_legs if self.can_be_IR_unresolved(l.get('id'))]
        fs_legs = [l for l in active_legs if l.get('state') == SubtractionLeg.FINAL]
        is_legs = [l for l in active_legs if l.get('state') == SubtractionLeg.INITIAL]

        # Loop over number of unresolved particles
        for unresolved in range(1, self.correction_order.count('N')):
             # Get iterators at the start of the final-state list
             it = iter(fs_legs)
             soft_it, coll_final_it, coll_initial_it = itertools.tee(it, 3)
             # Final state particle sets going soft
             for soft_set in itertools.combinations(soft_it, unresolved):
                 if any(self.can_become_soft(p) for p in self.parent_PDGs(soft_set)):
                     elementary_operator_list.append(SoftOperator(soft_set))
             # Final state particle sets going collinear
             for coll_final_set in itertools.combinations(coll_final_it, unresolved + 1):
                 if any(self.can_become_collinear(p) for p in self.parent_PDGs(coll_final_set)):
                     elementary_operator_list.append(CollOperator(coll_final_set))
             # Initial-final collinear
             # For any final-state set with one less particle
             for coll_initial_set in itertools.combinations(coll_initial_it, unresolved):
                 # Combine with all initial legs
                 for coll_initial_leg in is_legs:
                     coll_set = coll_initial_set
                     coll_set.append(coll_initial_leg)
                     if any(self.can_become_collinear(p) for p in self.parent_PDGs(coll_set)):
                        elementary_operator_list.append(CollOperator(coll_set))

        return elementary_operator_list
        