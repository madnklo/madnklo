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
"""Implementation of NLO type of currents."""

import os
import sys
import math

import madgraph.core.subtraction as subtraction
import madgraph.integrator.phase_space_generators as PS_utils
import madgraph.various.misc as misc

pjoin = os.path.join

# The subtraction_current_implementations_utils module will be located
# at an arbitrary location but identical to the one of this file, so we
# add it here
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(root_path, os.path.pardir))
import subtraction_current_implementations_utils as utils

CurrentImplementationError = utils.CurrentImplementationError
CS_utils = utils.CS_utils

# Shorthands
TR = utils.TR
CF = utils.CF
NC = utils.NC

class NLO_FF_QCD_collinear_qqx(utils.VirtualCurrentImplementation):
    """ Implemen ts the GluonToQQbar NLO collinear current."""

    def __init__(self, *args, **opts):
        super(NLO_FF_QCD_collinear_qqx, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Returns None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning 
        a dictionary, it specifies potential options that must passed upon instanciating
        the current implementation for the current given in argument. """

        # For debugging accept all currents
        # return {} 

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
        
        # Finally check that the singular structure and PDG matches
        singular_structure = current.get('singular_structure')

        # It should be a collinear type of structure
        if singular_structure.name()!='C':
            return None

        # It should not have nested structures
        if singular_structure.substructures:
            return None

        # It should consist in exactly two legs going collinear
        if len(singular_structure.legs)!=2:
            return None

        for leg in singular_structure.legs:
            if leg.state != subtraction.SubtractionLeg.FINAL:
                return None
            if abs(leg.pdg) not in range(1,7):
                return None
            if model.get_particle(leg.pdg).get('mass').upper()!='ZERO':
                return None
        
        # We now know that this current is implemented here, so we return
        # an empty dictionary which could potentially have contained specific
        # options to specify upon instantiating this class.
        return {}

    def get_cache_and_result_key(self,  current, 
                                        PS_point,
                                        reduced_process=None,
                                        leg_numbers_map=None,
                                        hel_config = None,
                                        mapping_variables = {},
                                        **opts
                      ):
        """ If this subtraction current depends on more than just the PS point and
        alpha_s, then complement the cache key here with the additional necessary information. """
        

        cache_key, result_key = super(NLO_FF_QCD_collinear_qqx, self).get_cache_and_result_key(
            current, PS_point,
            reduced_process=reduced_process, leg_numbers_map=leg_numbers_map, 
            hel_config=hel_config, mapping_variables=mapping_variables)
        
        # cache_key['another_call_specifier'] = 'set_it_here'

        return cache_key, result_key

    def evaluate_subtraction_current(self,  current, 
                                            PS_point, 
                                            reduced_process = None,
                                            leg_numbers_map = None,
                                            hel_config = None,
                                            mapping_variables = {}
                                     ):
        """ Now evalaute the current and return the corresponding instance of
        SubtractionCurrentResult. See documentation of the mother function for more details.
        Implementation according to Eq.12 of https://arxiv.org/pdf/hep-ph/9810389.pdf."""

        if not hel_config is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                            "%s does not support helicity assignment."%self.__class__.__name__) 

        if leg_numbers_map is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                      "%s requires the leg_number_map."%self.__class__.__name__)

        result = utils.SubtractionCurrentResult()

        ss = current.get('singular_structure')
        
        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        children_numbers = tuple(leg.n for leg in ss.legs)
        parent_number    = leg_numbers_map.inv[children_numbers]
        
        kin_variables = CS_utils.get_massless_collinear_CS_variables(
                PS_point, parent_number, children_numbers, mapping_variables=mapping_variables)
        z  = kin_variables['z%d'%ss.legs[0].n]
        kT_vec = kin_variables['kt%d'%ss.legs[0].n]
#        misc.sprint(z,kT_vec,kT_vec.square())
        kT_vec = kT_vec/math.sqrt(abs(kT_vec.square()))

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None, (parent_number,(tuple(kT_vec),))],
            'color_correlations'  : [ None ],
            'values'              : { (0,0): { 'finite' : None },
                                      (1,0): { 'finite' : None },
                                  }
          }
        )

        evaluation['values'][(0,0)]['finite'] = -TR
        evaluation['values'][(1,0)]['finite'] = 4.0*z*(1-z)
        
        # Inefficient, we should add this to the kinematic variables also generated in the mapping!
        q_sum = PS_utils.LorentzVector(4)
        for leg in ss.legs:
            q_sum += PS_utils.LorentzVector(PS_point[leg.n])
        s12 = q_sum.square()
        
        # Now add the normalization factors
        norm = 4.0*math.pi*alpha_s*(2.0/s12)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm
        
        result.add_result(evaluation, 
                          hel_config=hel_config, 
                          squared_orders=tuple(sorted(current.get('squared_orders').items()))
                         )
 
        return result
