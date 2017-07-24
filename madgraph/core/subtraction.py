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
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir))
from madgraph import MadGraph5Error, MG5DIR, InvalidCmd
import madgraph.various.misc as misc 
import madgraph.core.base_objects as base_objects
import madgraph.fks.fks_common as fks_common

logger = logging.getLogger('madgraph')
pjoin = os.path.join

#===============================================================================
# SubtractionOperator 
#===============================================================================
class SubtractionOperator(object):
    """ Subtraction operator definition"""

    def __init__(self, *args, **opts):
        """ Initialize a subtraction operator with legs and operator sets."""
        
        # Ordered list of operators that this SubtractionOperators acts on.
        self.legs       = SLegSet([a for a in args if isinstance(a, SLeg)])
        self.operators  = SLegSet([a for a in args if isinstance(a, SLeg)])
    
    def is_elementary(self):
        """ Return whether this Operator is elementary in the sense that it is not
        an empty structure and it hasn't been applied on yet."""
        return (self.legs and not self.operators)

    def act_on(self, operand):
        """ Core of the algorithm defining the action of an elementary operator on 
        a non-elementary structure."""
        
        if (not self.is_elementary()) or operand.is_elementary():
            raise "You can only apply elementary operators to act on non-elementary structures."

        # How many legs are common between the current operators and the operand at the highest level
        overlap = len(self.legs.intersection(operands.legs))

        # No overlap, we must look into the operators of the operand to find a match.
        if overlap == 0:
            for operator in self.operand.operators:
                # Investigate the action of operators within the operand (recursion starts)
                result = self.act_on(operator)
                if not result is None:
                    return result
            # If no match was found also in the operator list of the operand, this means that the
            # legs are no where to be found in the operand and we should simply add on this operator
            # to the list of operators. In other terms, non overlapping limits are additively combined.
            operand.operators.append(self)
            return True

        # Partial overlap always set to zero. 
        elif overlap != len(self.legs):
            return None
        
        # Full overlap, combine the legs into this current operator.
        else:
            operand.operators.append(self)
            # Remove the legs from the 
            operand.legs = operand.legs.difference(self.legs)
            return True

class SoftOperator(SubtractionOperator):
    def __init__(self, *args, **opts):
        """ Initialize a soft Subtraction operator with all relevant arguments."""
        super(SoftOperator, self).__init__(self, *args, **opts)
        
class CollOperator(SubtractionOperator):
    def __init__(self, *args, **opts):
        """ Initialize a Collinear operator with all relevant arguments."""
        super(CollOperator, self).__init__(self, *args, **opts) 

#===============================================================================
# List of subtraction operators with dedicated handling
#===============================================================================
class SubtractionOperatorList(list):
    """ Set of subtraction operators."""
    
    def __init__(self, *args, **opts):
        res = super(SubtractionOperatorList, self).__init__(*args, **opts)
        #self.sort()
        return res

    def append(self, op):
        """ Add op to current list while respecting the ordr of operators."""
        self.append(op)
        #self.sort()

    def sort(self):
        """ Chose a useful ordering"""
        # Check if really necessary in the future:
#        self.sort(lambda key: len(k.legs))

#   The following two functions could eventually be used to implement 
#   simplifications between different operators.
#    def act_on(self, other):
#    def simplify():

#===============================================================================
# SLeg 
#===============================================================================

class SLeg(base_objects.Leg):
    """ A subtraction Leg object."""

    def __init__(self, *args, **ops):
        super(SubtractionLeg,self).__init__(*args, **opts)


#===============================================================================
# Standalone main for debugging / standalone trials 
#===============================================================================
if __name__ == '__main__':
    misc.sprint("Put your standalone subtraction code here.")


class IRSubtraction(object):
    
    def __init__(self, model, correction_order='NLO', correction_types=['QCD']):
        """ Initialize a subtracter for a given model, correction order and type."""
        
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
        """ Checks whether a particle given by its PDG can become unresolved and lead 
        to singular behavior."""
        
        return any( (PDG in self.IR_quantities_for_corrections_types[order])
                                                    for order in self.correction_types)

    def parent_PDGs(self, legs):
        """ List all possible parent PDGs to a given set of PDGs specified."""
        
        # A generic way to build histories of possible combination of external legs is using the 
        # model ref_dict_to1
        #misc.sprint(self.model.get('ref_dict_to1')[(-2,2)])
        #misc.sprint(self.model.get_particle(-2).get('spin'),
        #            self.model.get_particle(-2).get_color(),
        #            self.model.get_particle(-2).get('mass')=='zero'
        #            )

        if not self.model.get('name').startswith('sm'):
            raise InvalidCmd("The function parent PDG is implemented for the SM only.")             
        
        if any(order!='QCD' in self.correction_types):
            raise InvalidCmd("The function parent PDG is implemented for QCD corrections only.")             
    
        flavs = [leg.get('id') for leg in legs if abs(leg.get('id')) != 21]

        if not flavs:
            return [21]
        
        last = flavs.pop()
        try:
            flavs.remove(-last)
            return parent(flavs)
        except ValueError:
            if parent(flavs) == 21:
                return [last]
            else:
                return None
        
        return None

    def can_become_soft(self, legs):
        """ Check whether all the PDGs specified can become soft together."""

        for pdg in self.parent_PDGs(legs):
            particle = self.model.get_particle(pdg)
            if (particle.get('spin')==3 and particle.get('mass')=='zero'):
                return True
        return False
    
    def can_become_collinear(self, legs):
        """ Test whether all the PDGs specified can be IR divergent when all collinear to one-another."""
        
        for pdg in self.parent_PDGs(legs):
            if self.can_be_IR_unresolved(pdg):
                return True
        return False        
    
    def get_all_elementary_operators(self, process):
        """ Generate all the 'building blocks' operator relevant for the process
        passed in argument."""
        
        elementary_operator_list = SubtractionOperatorList([])
        
        # Eliminate particles that do not take place in the subtraction
        active_legs = [l for l in process.get('legs') if self.can_be_IR_unresolved(l.get('id'))]
    
        # Loop over number of unresolved particles
        for unresolved in range(1, self.correction_order.count('N')):
             # Get iterators at the start of the final-state list
             it = iter(enum)
             it.next()
             it.next()
             sit, cit, ifit = itertools.tee(it, 3)
             # Final state particle sets (with gluon parent) going soft
#             for sset in itertools.combinations(sit, unresolved):
#                 if parent([x[-1] for x in sset]) == gluon:
#                     foo = frozenset(x[0] for x in sset)
#                     bolist += [(soft, foo), ]
#             # Final state particle sets (with one parent parton) going collinear
#             for cset in itertools.combinations(cit, unresolved + 1):
#                 if parent([x[-1] for x in cset]) is not None:
#                     foo = frozenset(x[0] for x in cset)
#                     bolist += [(coll, foo), ]
#             # Initial-final collinear
#             # The particle going into the hard process has to be a single parton
#             for cset in itertools.combinations(ifit, unresolved):
#                 if parent([-enum[0][-1]] + [x[-1] for x in cset]) is not None:
#                     foo = frozenset([0, ] + [x[0] for x in cset])
#                     bolist += [(coll, foo), ]
#                 if parent([-enum[1][-1]] + [x[-1] for x in cset]) is not None:
#                     foo = frozenset([1, ] + [x[0] for x in cset])
#                     bolist += [(coll, foo), ]
#         return bolist
        
        