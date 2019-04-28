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
"""Implementation of NLO type of integrated currents for colorful_pp."""

import os
import sys
import math

import madgraph.core.subtraction as subtraction
import madgraph.integrator.phase_space_generators as PS_utils
from madgraph.core.base_objects import EpsilonExpansion
import madgraph.various.misc as misc

import commons.utils as utils
import commons.QCD_local_currents as currents

from integrated_current_expressions import HE

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError

# Mother function grouping functionalities common to all integrated FF NLO QCD currents
class integrated_NLO_FF_QCD_current(utils.IntegratedCurrent, currents.QCDCurrent):
    """ Just a template class for all Final-Final NLO QCD local subtraction current."""
    
    @classmethod
    def common_does_implement_this_current(cls, current, QCD_squared_order=None, n_loops=None):
        """ General class of checks common to all currents inheriting from this class."""
    
        # Make sure it is an integrated subtraction counterterm and not a local one.
        if not isinstance(current, subtraction.IntegratedCurrent):
            return None

        # Also make sure this is not a beam factorization current as these are implemented
        # by the currents in beam_factorization.py
        if isinstance(current, (subtraction.BeamCurrent, subtraction.IntegratedBeamCurrent)):
            return None   

        squared_orders = current.get('squared_orders')

        # First check that we indeed have a pure NLO QCD current
        if squared_orders['QCD'] != 2 or \
           any(squared_orders[order] != 0 for order in squared_orders if order!='QCD'):
            return None

        # Now check that it is tree-level
        if current.get('n_loops')!=0:
            return None

        # Make sure we don't need to sum over the quantum number of the mother leg
        if not current.get('resolve_mother_spin_and_color'):
            return None

        # Integrated counterterm have an inert singular structure that contains
        # exactly one actual singular structure
        if current.get('singular_structure').name()!='':
            return None
        if len(current.get('singular_structure').substructures)!=1:
            return None
        singular_structure = current.get('singular_structure').substructures[0].substructures[0]
        
        # Finally check that the singular structure and PDG matches 
        
        # Check that there is at most one substructure (it's NLO here)
        for substructure in singular_structure.substructures:
            if len(substructure.substructures)>0:
                return None

        return {}

class integrated_NLO_FF_QCD_softcollinear_gq(integrated_NLO_FF_QCD_current):
    """ Implements the NLO_FF_QCD_collinear_gq current."""

    # The soft-collinear integrated counterterm has been accounted for completely in the 
    # soft integrated counterterm
    is_zero = True
    
    def __init__(self, *args, **opts):
        
        # Make sure it is initialized with the proper set of options and remove them
        # before calling the mother constructor
        if 'color_charge' not in opts:
            raise CurrentImplementationError(
                "The current '%s' must be instantiated with "%self.__class__.__name__+
                                                " a 'color_charge' option specified.")
        color_charge = opts.pop('color_charge')
        
        super(integrated_NLO_FF_QCD_softcollinear_gq, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False
        # At this state color_charge is the string of the argument to retrieve ('CA' or 'CF')
        # And now that the mother constructor is called, the group factors have been initialized
        # and we can retrieve them.
        self.color_charge = getattr(self, color_charge)

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Returns None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning 
        a dictionary, it specifies potential options that must passed upon instantiating
        this implementation for the current given in argument. """
        
        # Check the general properties common to FF NLO QCD
        if super(integrated_NLO_FF_QCD_softcollinear_gq, cls).common_does_implement_this_current(
                                                                current, model) is None:
            return None

        # Finally check that the singular structure and PDG matches
        singular_structure = current.get('singular_structure').substructures[0].substructures[0]
        
        # It main structure should be of collinear type
        if singular_structure.name()!='C':
            return None

        # It should have only one leg left, the other one being in the nested soft structure
        if len(singular_structure.legs)!=1:
            return None

        # Make sure legs are final and massless
        for leg in singular_structure.legs:
            if model.get_particle(leg.pdg).get('mass').upper()!='ZERO':
                return None

        # The leg not soft must be quark or a gluon      
        if not abs(singular_structure.legs[0].pdg) in [21,]+range(1,7):
            return None

        # It should have exactly one nested structures
        if len(singular_structure.substructures)!=1:
            return None
        
        sub_singular_structure = singular_structure.substructures[0]
        
        # Make sure this substructure is soft
        if sub_singular_structure.name()!='S':
            return None
        
        # Make sure it contains a single soft leg
        if len(sub_singular_structure.legs)!=1:
            return None
        
        soft_leg = sub_singular_structure.legs[0]
        
        # Make sure the soft leg is massless final and a gluon
        if model.get_particle(soft_leg.pdg).get('mass').upper()!='ZERO':
            return None
        if soft_leg.pdg != 21:
            return None
        
        # We now know that this current is implemented here. We return
        # the specific color charge to instantiate this kernel with,
        # in the form of a the name of the group factor to retrieve upon
        # initialization.
        if singular_structure.legs[0].pdg == 21:
            # This is a 'g > g g' soft-collinear splitting
            color_charge = 'CA'
        else:
            # This is a 'q > g g' soft-collinear splitting
            color_charge = 'CF'

        return {'color_charge': color_charge}

    def evaluate_integrated_current(self,  current, 
                                            PS_point, 
                                            reduced_process = None,
                                            leg_numbers_map = None,
                                            hel_config = None,
                                            compute_poles = False,
                                            **opts
                                     ):

        """ Now evalaute the current and return the corresponding instance of
        SubtractionCurrentResult. See documentation of the mother function for more details."""

        if not hel_config is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                            "%s does not support helicity assignment."%self.__class__.__name__) 
        if leg_numbers_map is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                      "%s requires the leg_number_map."%self.__class__.__name__)
        if reduced_process is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                "%s requires the reduced_process."%self.__class__.__name__)

        result = utils.SubtractionCurrentResult()

        ss = current.get('singular_structure').substructures[0].substructures[0]
        
        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        children_numbers = (ss.legs[0].n, ss.substructures[0].legs[0].n) 
        parent_number    = leg_numbers_map.inv[frozenset(children_numbers)]

        p12      = PS_point[parent_number]
        Q        = sum([PS_point[l.get('number')] for l in reduced_process.get_initial_legs()])
        Q_square = Q.square()
        y12      = 2.*Q.dot(p12) / Q_square

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [ None ],
            'values'              : { (0,0): { } }
          }
        )

        # The soft-collinear integrated counterterm has been accounted for completely in the 
        # soft integrated counterterm
        value = EpsilonExpansion({ 0 : 0., -1 : 0., -2 : 0.})
        
        logMuQ = math.log(mu_r**2/Q_square)
        
        prefactor = EpsilonExpansion({ 0 : 1., 1 : logMuQ, 2 : 0.5*logMuQ**2 })
        prefactor *= self.SEpsilon
        
        # Now add the normalization factors
        value *= prefactor*(alpha_s / (2.*math.pi))*self.CF
        value.truncate(min_power=-2, max_power=0)
        
        # Now register the value in the evaluation
        evaluation['values'][(0,0)] = value.to_human_readable_dict()       

        # And add it to the results
        result.add_result(
            evaluation, 
            hel_config=hel_config, 
            squared_orders=tuple(sorted(current.get('squared_orders').items()))
        )
 
        return result


class integrated_NLO_QCD_soft_gluon(integrated_NLO_FF_QCD_current):
    """ Implements the soft gluon Eikonel current.
    See Eq.4.12-4.13 of ref. https://arxiv.org/pdf/0903.1218.pdf"""

    def __init__(self, *args, **opts):
        misc.sprint("")
        misc.sprint("=======================================================================================")
        misc.sprint("Warning: the integrated Final-final soft no longer matches the definition of the locals")
        misc.sprint("=======================================================================================")
        misc.sprint("")
        super(integrated_NLO_QCD_soft_gluon, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False
        
    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Returns None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning 
        a dictionary, it specifies potential options that must passed upon instantiating
        this implementation for the current given in argument. """

        # Check the general properties common to FF NLO QCD
        if super(integrated_NLO_QCD_soft_gluon, cls).common_does_implement_this_current(
                                                                current, model) is None:
            return None
        
        # Finally check that the singular structure and PDG matches
        singular_structure = current.get('singular_structure').substructures[0].substructures[0]

        # It should be a soft type of structure
        if singular_structure.name()!='S':
            return None

        # It should not have nested structures
        if singular_structure.substructures:
            return None

        # It should consist in exactly one legs going soft
        if len(singular_structure.legs)!=1:
            return None

        for leg in singular_structure.legs:
            if abs(leg.pdg) not in [21]:
                return None
            if model.get_particle(leg.pdg).get('mass').upper()!='ZERO':
                return None
        
        # We now know that this current is implemented here, so we return
        # an empty dictionary which could potentially have contained specific
        # options to specify upon instantiating this class.
        return {}
    
    def evaluate_integrated_current(self,  current, 
                                            PS_point, 
                                            reduced_process = None,
                                            leg_numbers_map = None,
                                            hel_config = None,
                                            compute_poles = False,
                                            **opts
                                     ):
        """ Evaluates this current and return the corresponding instance of
        SubtractionCurrentResult. See documentation of the mother function for more details."""
        
        if not hel_config is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                            "%s does not support helicity assignment."%self.__class__.__name__) 

        if leg_numbers_map is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                      "%s requires the leg_number_map."%self.__class__.__name__)
        
        if reduced_process is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                      "%s requires a reduced_process."%self.__class__.__name__)
        
        result = utils.SubtractionCurrentResult()

        ss = current.get('singular_structure').substructures[0].substructures[0]
        
        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        soft_leg_number = ss.legs[0].n
        # Use the momenta map, in case it has been remapped.
        # Although for the soft current it's typically not the case
        soft_leg_number   = leg_numbers_map.inv[frozenset([soft_leg_number,])]
        
        Q        = sum([PS_point[l.get('number')] for l in reduced_process.get_initial_legs()])
        Q_square = Q.square()
        
        # Now find all colored leg numbers in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color')==1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [ ],
            'values'              : { }
          })
        
        logMuQ = math.log(mu_r**2/Q_square)
        
        prefactor = EpsilonExpansion({ 0 : 1., 1 : logMuQ, 2 : 0.5*logMuQ**2 })
        prefactor *= self.SEpsilon
        
        # Now add the normalization factors
        prefactor *= (alpha_s / (2.*math.pi))
        prefactor.truncate(min_power=-2, max_power=2)

        #Virtuality cut in the integration
        y_0=0.5

        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        for i, a in enumerate(all_colored_parton_numbers):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in all_colored_parton_numbers[i+1:]:
                evaluation['color_correlations'].append( ((a, b), ) )
                # We multiply by a factor 2. because we symmetrized the sum below
                value = prefactor*2.
                pa = PS_point[a]
                pb = PS_point[b]
                Y = (pa.dot(pb) * Q_square) / (2. * Q.dot(pa) * Q.dot(pb))
                finite_part = HE.SoftFF_Finite_Gabor_DIVJAC_NOD0(y_0,Y)
                value *= EpsilonExpansion({
                    0   : finite_part,
                    -1  : math.log(Y),
                    -2  : 0.
                })
                # Truncate expansion so as to keep only relevant terms
                value.truncate(min_power=-2, max_power=0)
                evaluation['values'][(0,color_correlation_index)] = value.to_human_readable_dict()
                color_correlation_index += 1
        
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items()))
        )
 
        return result
