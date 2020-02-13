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
"""Implementation of NNLO colorful_pp local currents."""

import math

import madgraph.integrator.vectors as vectors
import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings
from bidict import bidict

import commons.utils as utils
from commons.universal_kernels import AltarelliParisiKernels, SoftKernels

import commons.general_current as general_current
import colorful_pp_config
import variables as kernel_variables

from madgraph.core.base_objects import EpsilonExpansion

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

# =========================================================================================
# NNNLO soft currents
# =========================================================================================

class QCD_S_FgFgFg(general_current.GeneralCurrent):
    """ Soft g(final) g(final) g(final) current"""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    soft_structure = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, 21, sub.SubtractionLeg.FINAL),            
        )
    )
    structure = sub.SingularStructure(substructures=(soft_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 6},
                'singular_structure'                : structure
            }),
        ),
    ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.SoftStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(12, 21, sub.SubtractionLeg.FINAL),                    
                )
            ),)),
            'mapping'               : colorful_pp_config.soft_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({}),
            'variables'             : None,
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_parton_numbers = sorted(colored_partons.keys())

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']
        higher_PS_point = all_steps_info[0]['higher_PS_point']

        colored_partons_momenta = {leg_number: lower_PS_point[leg_number] for leg_number in colored_parton_numbers}
        soft_momenta = [ higher_PS_point[leg_number] for leg_number in
                                        all_steps_info[0]['bundles_info'][0]['final_state_children'] ]

        # All contributions from the soft counterterms are color correlated so we should remove
        # the default value None from the list of color correlations
        evaluation['color_correlations'] = []

        color_correlation_index = 0
        for color_correlator, weight in SoftKernels.gg(
                self, colored_partons_momenta, soft_momenta, colored_parton_numbers):
            evaluation['color_correlations'].append(color_correlator)
            evaluation['values'][(0, color_correlation_index, 0, 0)] = {'finite': weight[0]}
            color_correlation_index += 1

        return evaluation

# =========================================================================================
# NNNLO collinear currents
# =========================================================================================

class QCD_C_FqFqxFqpFqpx(general_current.GeneralCurrent):
    """ Quadruple final-collinear current g(final) > q(final) q~(final) q'(final) q'~(final)"""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),            
            sub.SubtractionLeg(12, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(13, -2, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 6},
                'singular_structure'                : structure
            }),
        ),
    ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),            
                    sub.SubtractionLeg(12, +2, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(13, -2, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12,13))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this FFFF counterterm given the supplied variables. """

        # Retrieve the collinear variable z's
        zs = all_steps_info[0]['variables'][0]['zs'][:4]
        # WARNING: the colorful_pp_FFn_variables do not yet specify how kTs must be computed for a quadruple collinear.
        # This will need to be specified in that function 'colorful_pp_FFn_variables'
        kTs = all_steps_info[0]['variables'][0]['kTs'][:4]
        s_invariants = all_steps_info[0]['variables'][0]['ss']
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        higher_PS_point = all_steps_info[0]['higher_PS_point']
        pC = sum(higher_PS_point[child] for child in all_steps_info[0]['bundles_info'][0]['children'])
        # Of course, this "propagator" structure may be more complex for the quadrilinear
        propagator = 1./(pC.square())**2

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqqpqp(self, zs, s_invariants, kTs):
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': weight[0]*propagator}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0, 0)] = {'finite': weight[0]*propagator}

        return evaluation

class QCD_C_FgFgFgFg(general_current.GeneralCurrent):
    """ Quadruple final-collinear current g(final) > g(final) g(final) g(final) g(final)"""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),            
            sub.SubtractionLeg(12, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(13, 21, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 6},
                'singular_structure'                : structure
            }),
        ),
    ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),            
                    sub.SubtractionLeg(12, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(13, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12,13))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this FFFF counterterm given the supplied variables. """

        # Retrieve the collinear variable z's
        zs = all_steps_info[0]['variables'][0]['zs'][:4]
        # WARNING: the colorful_pp_FFn_variables do not yet specify how kTs must be computed for a quadruple collinear.
        # This will need to be specified in that function 'colorful_pp_FFn_variables'
        kTs = all_steps_info[0]['variables'][0]['kTs'][:4]
        s_invariants = all_steps_info[0]['variables'][0]['ss']
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        higher_PS_point = all_steps_info[0]['higher_PS_point']
        pC = sum(higher_PS_point[child] for child in all_steps_info[0]['bundles_info'][0]['children'])
        # Of course, this "propagator" structure may be more complex for the quadrilinear
        propagator = 1./(pC.square())**2

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_gggg(self, zs, s_invariants, kTs):
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': weight[0]*propagator}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0, 0)] = {'finite': weight[0]*propagator}

        return evaluation
