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
# NNLO initial-collinear currents
# =========================================================================================

class QCD_C_IqFqpFqpx(general_current.GeneralCurrent):
    """ Triple initial-collinear current q(initial) > q(initial) q'(final) q'(final)"""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, -2, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(11, +2, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(12, -2, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """


        # Retrieve the collinear variable x
        x_a, x_r, x_s = all_steps_info[0]['variables'][0]['xs'][:3]
        kT_a, kT_r, kT_s = all_steps_info[0]['variables'][0]['kTs'][:3]
        s_invariants = all_steps_info[0]['variables'][0]['ss']
        s_ar, s_as, s_rs = s_invariants['ss'][(0,1)], s_invariants['ss'][(0,2)], s_invariants['ss'][(1,2)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorization formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= 1.

        higher_PS_point = all_steps_info[0]['higher_PS_point']
        pC = higher_PS_point[all_steps_info[0]['bundles_info'][0]['children'][0]]-\
                sum(higher_PS_point[child] for child in all_steps_info[0]['bundles_info'][0]['children'][1:])
        propagator = 1./(pC.square())**2

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqpqp(self,
                1./x_a, -x_r/x_a, -x_s/x_a, -s_ar, -s_as, s_rs, kT_a, kT_r, kT_s
            ):
            complete_weight = weight*initial_state_crossing_factor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]*propagator}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0, 0)] = {'finite': complete_weight[0]*propagator}

        return evaluation

class QCD_C_IqFqFqx(general_current.GeneralCurrent):
    """ Triple initial-collinear current q(initial) > q(initial) q(final) q(final)"""


    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, -1, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(12, -1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """


        # Retrieve the collinear variable x
        x_a, x_r, x_s = all_steps_info[0]['variables'][0]['xs'][:3]
        kT_a, kT_r, kT_s = all_steps_info[0]['variables'][0]['kTs'][:3]
        s_invariants = all_steps_info[0]['variables'][0]['ss']
        s_ar, s_as, s_rs = s_invariants['ss'][(0,1)], s_invariants['ss'][(0,2)], s_invariants['ss'][(1,2)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorization formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= 1.

        higher_PS_point = all_steps_info[0]['higher_PS_point']
        pC = higher_PS_point[all_steps_info[0]['bundles_info'][0]['children'][0]]-\
                sum(higher_PS_point[child] for child in all_steps_info[0]['bundles_info'][0]['children'][1:])
        propagator = 1./(pC.square())**2

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqq(self,
                1./x_a, -x_r/x_a, -x_s/x_a, -s_ar, -s_as, s_rs, kT_a, kT_r, kT_s
            ):
            complete_weight = weight * initial_state_crossing_factor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]*propagator}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {
                    'finite': complete_weight[0]*propagator}

        return evaluation

class QCD_C_IqFgFg(general_current.GeneralCurrent):
    """ Triple initial-collinear current q(initial) > q(initial) g(final) g(final)"""


    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, 21, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(12, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """


        # Retrieve the collinear variable x
        x_a, x_r, x_s = all_steps_info[0]['variables'][0]['xs'][:3]
        kT_a, kT_r, kT_s = all_steps_info[0]['variables'][0]['kTs'][:3]
        s_invariants = all_steps_info[0]['variables'][0]['ss']
        s_ar, s_as, s_rs = s_invariants['ss'][(0,1)], s_invariants['ss'][(0,2)], s_invariants['ss'][(1,2)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorization formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= 1.

        higher_PS_point = all_steps_info[0]['higher_PS_point']
        pC = higher_PS_point[all_steps_info[0]['bundles_info'][0]['children'][0]]-\
                sum(higher_PS_point[child] for child in all_steps_info[0]['bundles_info'][0]['children'][1:])
        propagator = 1./(pC.square())**2

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qgg(self,
                1./x_a, -x_r/x_a, -x_s/x_a, -s_ar, -s_as, s_rs, kT_a, kT_r, kT_s
            ):
            complete_weight = weight * initial_state_crossing_factor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]*propagator}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {
                    'finite': complete_weight[0]*propagator}

        return evaluation

# =========================================================================================
# NNLO final-collinear currents
# =========================================================================================

class QCD_C_FqFqpFqpx(general_current.GeneralCurrent):
    """ collinear FSR tree-level current. q(final) > q(final) q'(final) q' (final) """

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, -2, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                    sub.SubtractionLeg(11, +2, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(12, -2, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """

        # Retrieve the collinear variable z's
        z_i, z_r, z_s = all_steps_info[0]['variables'][0]['zs'][:3]
        kTs = all_steps_info[0]['variables'][0]['kTs'][:3]
        # This current does not involve KTs.
        kT_i, kT_r, kT_s = None, None, None
        s_invariants = all_steps_info[0]['variables'][0]['ss']
        s_ir, s_is, s_rs = s_invariants['ss'][(0, 1)], s_invariants['ss'][(0, 2)], s_invariants['ss'][(1, 2)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        higher_PS_point = all_steps_info[0]['higher_PS_point']
        pC = sum(higher_PS_point[child] for child in all_steps_info[0]['bundles_info'][0]['children'])
        propagator = 1./(pC.square())**2

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqpqp(self,
                z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s
            ):
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': weight[0]*propagator}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0, 0)] = {'finite': weight[0]*propagator}

        return evaluation

class QCD_C_FqFqFqx(general_current.GeneralCurrent):
    """ collinear FSR tree-level current. q(final) > q(final) q(final) q(final) """

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, -1, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                    sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(12, -1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : colorful_pp_config.generalised_cuts,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """

        # Retrieve the collinear variable z's
        z_i, z_r, z_s = all_steps_info[0]['variables'][0]['zs'][:3]
        kTs = all_steps_info[0]['variables'][0]['kTs'][:3]
        # This current does not involve KTs.
        kT_i, kT_r, kT_s = None, None, None
        s_invariants = all_steps_info[0]['variables'][0]['ss']
        s_ir, s_is, s_rs = s_invariants['ss'][(0, 1)], s_invariants['ss'][(0, 2)], s_invariants['ss'][(1, 2)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        higher_PS_point = all_steps_info[0]['higher_PS_point']
        pC = sum(higher_PS_point[child] for child in all_steps_info[0]['bundles_info'][0]['children'])
        propagator = 1./(pC.square())**2

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqq(self,
                z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s
            ):
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': weight[0]*propagator}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0, 0)] = {'finite': weight[0]*propagator}

        return evaluation

class QCD_C_FqFgFg(general_current.GeneralCurrent):
    """ collinear FSR tree-level current. q(final) > q(final) g(final) g(final)"""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, 21, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(12, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : colorful_pp_config.generalised_cuts,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """

        # Retrieve the collinear variable z's
        z_i, z_r, z_s = all_steps_info[0]['variables'][0]['zs'][:3]
        kTs = all_steps_info[0]['variables'][0]['kTs'][:3]
        # This current does not involve KTs.
        kT_i, kT_r, kT_s = None, None, None
        s_invariants = all_steps_info[0]['variables'][0]['ss']
        s_ir, s_is, s_rs = s_invariants['ss'][(0, 1)], s_invariants['ss'][(0, 2)], s_invariants['ss'][(1, 2)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        higher_PS_point = all_steps_info[0]['higher_PS_point']
        pC = sum(higher_PS_point[child] for child in all_steps_info[0]['bundles_info'][0]['children'])
        propagator = 1./(pC.square())**2

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qgg(self,
                z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s
            ):
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': weight[0]*propagator}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0, 0)] = {'finite': weight[0]*propagator}

        return evaluation

# =========================================================================================
# NNLO soft currents
# =========================================================================================

class QCD_S_FqFqx(general_current.GeneralCurrent):
    """ Soft q(final) q(final) current"""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    soft_structure = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(soft_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
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
        for color_correlator, weight in SoftKernels.qqx(
                self, colored_partons_momenta, soft_momenta, colored_parton_numbers):
            evaluation['color_correlations'].append(color_correlator)
            evaluation['values'][(0, color_correlation_index, 0, 0)] = {'finite': weight[0]}
            color_correlation_index += 1

        return evaluation

class QCD_S_FgFg(general_current.GeneralCurrent):
    """ Soft g(final) g(final) current"""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    soft_structure = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(soft_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
# NNLO soft-collinear currents
# =========================================================================================

class QCD_S_FqFqx_C_IgFqFqx(general_current.GeneralCurrent):
    """ Soft-collinear NNLO current g(initial) > g(initial) q(final) q~(final)"""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    soft_structure = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, -1, sub.SubtractionLeg.FINAL),
        )
    )
    coll_structure_g = sub.CollStructure(
        substructures=tuple([soft_structure,]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure_g,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
                'singular_structure'                : structure
            }),
        ),
    ]

    # The q > q g g soft-collinear configuration factorizes the color factor CA
    color_charge = 'CA'

    # An now the mapping rules (which will also be the one used in the QCD_S_FgFg_C_IqFqFqx as it inherits from this one)
    mapping_rules = [
        {
            'singular_structure'    : structure.get_copy(),
            'mapping'               : colorful_pp_config.initial_soft_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_IFF_softFF_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """


        # Retrieve the collinear variable x
        x_a, x_r, x_s = all_steps_info[0]['variables'][0]['xs'][:3]
        kT_a, kT_r, kT_s = all_steps_info[0]['variables'][0]['kTs'][:3]
        s_invariants = all_steps_info[0]['variables'][0]['ss']
        s_ar, s_as, s_rs = s_invariants['ss'][(0,1)], s_invariants['ss'][(0,2)], s_invariants['ss'][(1,2)]
        s_a_rs = s_ar + s_as

        color_charge = getattr(self, self.color_charge)

        lower_PS_point = all_steps_info[0]['lower_PS_point']
        higher_PS_point = all_steps_info[0]['higher_PS_point']

        p_a_tilde = lower_PS_point[global_variables['overall_parents'][0]]
        # The soft momenta are the last two children given that the ordering follows the ascending order given
        # by the leg number assigned in the matching structures
        soft_momenta_summed = higher_PS_point[global_variables['overall_children'][0][1]]+ \
                              higher_PS_point[global_variables['overall_children'][0][2]]

        propagator = 1.0 / ((2. * p_a_tilde.dot(soft_momenta_summed)) * soft_momenta_summed.square())

        kernel = 2. * ( 1. / (x_r + x_s) - (((s_ar * x_s - s_as * x_r) ** 2) / (s_a_rs * s_rs * ((x_r + x_s) ** 2))) )

        evaluation['values'][(0, 0, 0, 0)] = {'finite': color_charge * propagator * self.TR * kernel }

        return evaluation

class QCD_S_FqFqx_C_IqFqFqx(QCD_S_FqFqx_C_IgFqFqx):
    """ Soft-collinear NNLO current q(initial) > q(initial) q(final) q~(final) and q(initial) > q(initial) q'(final) q~'(final)"""

    # Now define the matching singular structures
    soft_structure = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, -1, sub.SubtractionLeg.FINAL),
        )
    )
    coll_structure_q = sub.CollStructure(
        substructures=tuple([soft_structure,]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
        )
    )
    coll_structure_qp = sub.CollStructure(
        substructures=tuple([soft_structure,]),
        legs=(
            sub.SubtractionLeg(10, +2, sub.SubtractionLeg.INITIAL),
        )
    )

    defining_currents = [
        # The tuple below indicates that any of the two currents can be matched to.
        tuple([
            sub.Current({
                'resolve_mother_spin_and_color' : True,
                'n_loops'                       : 0,
                'squared_orders'                : {'QCD': 4},
                'singular_structure'            : sub.SingularStructure(substructures=(coll_structure_q,))
            }),
            sub.Current({
                'resolve_mother_spin_and_color' : True,
                'n_loops'                       : 0,
                'squared_orders'                : {'QCD': 4},
                'singular_structure'            : sub.SingularStructure(substructures=(coll_structure_qp,))
            }),
        ]),
    ]

    # The q > q q q soft-collinear configuration factorizes the color factor CF
    color_charge = 'CF'

    # The actual kernel evaluation is then the same as the mother class QCD_S_FgFg_C_IgFqFqx and therefore need
    # not be repeated here.

#=========================================================================================
#
# Nested currents
#
# These are more advanced currents that require dedicated implementations of the cuts,
# variables and class attributes
#
#=========================================================================================

class QCD_C_FqFqx_C_IqpFqFqx(general_current.GeneralCurrent):
    """ Nested FF (q_qx) collinear within IFF (q_q'q'x)."""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    sub_coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -2, sub.SubtractionLeg.FINAL),
        )
    )
    # Match both the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(sub_coll_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                    sub.SubtractionLeg(10, +2, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, -2, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({1001:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL)]),
        },
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
            'momenta_dict'          : bidict({-1: frozenset((1001, 1))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        kT_IF = all_steps_info[1]['variables'][0]['kTs'][0]
        x_IF  = all_steps_info[1]['variables'][0]['xs'][0]
        s_a_rs = all_steps_info[1]['variables'][0]['ss'][(0,1)]

        p_a_hat = all_steps_info[1]['higher_PS_point'][
            all_steps_info[1]['bundles_info'][0]['initial_state_children'][0]
        ]
        p_rs_hat = all_steps_info[1]['higher_PS_point'][
            all_steps_info[1]['bundles_info'][0]['final_state_children'][0]
        ]

        parent = all_steps_info[1]['bundles_info'][0]['parent']

        #misc.sprint(s_rs, s_a_rs)
        #misc.sprint(z_FF, x_IF)
        #misc.sprint(kT_FF, kT_IF)
        #misc.sprint(p_a_hat, parent)

        # We must include here propagator factors, but no correction for symmetry
        # or averaging factor is necessary in this particular pure-quark kernel
        initial_state_crossing = -1.0
        prefactor = initial_state_crossing*(1./(s_rs*s_a_rs))
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_q_qpqp(self, z_FF, 1./x_IF, kT_FF, -p_a_hat, p_rs_hat):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

class QCD_C_FqFqx_C_IqFqFqx(QCD_C_FqFqx_C_IqpFqFqx):
    """ Nested FF (q_qx) collinear within IFF (q_qqx).
    Identical to IFF (q_q'q'x) except for flavor differences in the mapping rules which are anyway irrelevant. """

    # Now define the matching singular structures
    sub_coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )
    # Match both the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(sub_coll_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({1001:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL)]),
        },
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
            'momenta_dict'          : bidict({-1: frozenset((1001, 1))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

def QCD_C_IqFqx_C_IqFqFqx_global_IFF_collinearIF_variables(reduced_process, all_steps, global_info):
    """ Calls the IF variables by forcing the parent momentum for the IF variable computation to be the
    overall parent momentum of the C(C(IF),F) structure.
    Needed for QCD_C_IqFqx_C_IqFqFqx where the hatted xs use Q_hat.
    """

    Q_original = global_info['Q']
    p_a = all_steps[0]['higher_PS_point'][
            all_steps[0]['bundles_info'][0]['initial_state_children'][0]
          ]
    p_b = Q_original - p_a
    p_a_hat = all_steps[1]['higher_PS_point'][
            all_steps[1]['bundles_info'][0]['initial_state_children'][0]
          ]
    p_a_tilde = all_steps[1]['lower_PS_point'][ global_info['overall_parents'][0] ]
    IF_variables = kernel_variables.colorful_pp_IF1_variables(
        all_steps[1]['higher_PS_point'],
        p_a_tilde,
        tuple( list(all_steps[1]['bundles_info'][0]['initial_state_children']) +
               list(all_steps[1]['bundles_info'][0]['final_state_children']) ),
        Q = Q_original,
        Q_hat = p_a_hat + p_b
    )[0]

    return IF_variables

class QCD_C_IqFqx_C_IqFqFqx(general_current.GeneralCurrent):
    """ Nested IF (q_qx) collinear within IFF (q_qqx)."""

    # TO ADJUST

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = staticmethod(QCD_C_IqFqx_C_IqFqFqx_global_IFF_collinearIF_variables)

    # Now define the matching singular structures
    sub_coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
        )
    )
    # This counterterm will be used if any of the structures of the list below matches
    structure = [
        # Match both the case of the initial state being a quark and/or antiquark
        sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(sub_coll_structure,),
            legs=(sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),)
        ),)),
        sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(sub_coll_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL),)
        ),)),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({1001:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
            'is_cut'                : colorful_pp_config.generalised_cuts,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL)]),
        },
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
            'momenta_dict'          : bidict({-1: frozenset((1001, 1))}),
            # needs special set of variables where the xs are produced using Q_hat (which is p_a_hat + p_b in this case)
            'variables'             : None,
            'is_cut'                : colorful_pp_config.generalised_cuts,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """

        kT_s  = all_steps_info[0]['variables'][0]['kTs'][1]
        x_a    = all_steps_info[0]['variables'][0]['xs'][0]
        s_as   = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        kT_rH  = global_variables['kTs'][0]
        x_aH   = global_variables['xs'][0]
        x_rH   = global_variables['xs'][1]
        s_r_as = global_variables['ss'][(0,1)]

        p_as_hat = all_steps_info[1]['higher_PS_point'][
            all_steps_info[1]['bundles_info'][0]['initial_state_children'][0]
        ]
        p_r_hat = all_steps_info[1]['higher_PS_point'][
            all_steps_info[1]['bundles_info'][0]['final_state_children'][0]
        ]

        #kT_rH  = all_steps_info[1]['variables'][0]['kTs'][0]
        #x_aH   = all_steps_info[1]['variables'][0]['xs'][0]
        #x_rH   = 1.-x_aH
        #s_r_as = all_steps_info[1]['variables'][0]['ss'][(0,1)]
        #p_b = global_variables['Q'] + all_steps_info[0]['bundles_info'][0]['cut_inputs']['pA']
        #QH = p_as_hat + p_b
        #x_aH = 1. - p_r_hat.dot(QH)/p_as_hat.dot(QH)

        parent = all_steps_info[1]['bundles_info'][0]['parent']

        # We must include here propagator factors, but no correction for symmetry
        # or averaging factor is necessary in this particular pure-quark kernel
        initial_state_crossing = 1.0
        prefactor = initial_state_crossing*(1./(s_as*s_r_as))
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_q_qpqp(self, 1./x_a, -x_rH/x_aH, kT_s, p_r_hat, -p_as_hat):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

class QCD_S_FqFqx_C_FqFqx(general_current.GeneralCurrent):
    """ Nested soft FF (q_qx) limit within collinear FF (q_qx) limit."""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    sub_coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -2, sub.SubtractionLeg.FINAL),
        )
    )
    # This counterterm will be used if any of the structures of the list below matches
    # Match both the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.SoftStructure(
            substructures=(sub_coll_structure,),
            legs=tuple([])
        ),))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
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
                    sub.SubtractionLeg(10, +2, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, -2, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({1001:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.SoftStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.soft_mapping,
            # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
            'momenta_dict'          : bidict({}),
            'variables'             : None,
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./s_rs

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qq(self, z_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ( (parent, (spin_correlation_vector,) ), ) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a collinear type of splitting kernel, which *does* need to know about the reduced process
        Should be specialised by the daughter class if not dummy
        """

        new_evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations' : [ None, ],
            'color_correlations': [ ],
            'reduced_kinematics': evaluation['reduced_kinematics'],
            'values': { }
        })

        overall_lower_PS_point = all_steps_info[-1]['lower_PS_point']
        soft_leg_number = all_steps_info[-1]['bundles_info'][0]['final_state_children'][0]
        pr = all_steps_info[-1]['higher_PS_point'][soft_leg_number]
        colored_parton_numbers = sorted(colored_partons.keys())

        for i, a in enumerate(colored_parton_numbers):
            for b in colored_parton_numbers[i:]:
                # Write the eikonal for that pair
                if a!=b:
                    mult_factor = 1.
                else:
                    mult_factor = 1./2.

                pi = overall_lower_PS_point[a]
                pk = overall_lower_PS_point[b]
                composite_weight = EpsilonExpansion({'finite': 0.})
                for (sc, cc, rk), coll_weight in evaluation['values'].items():
                    if evaluation['spin_correlations'][sc] is None:
                        # We *subtract* here the contribution because the non-spin-correlated contributin is -g^{\mu\nu}
                        composite_weight -= SoftKernels.eikonal_g(self, pi, pk, pr, spin_corr_vector=None)*EpsilonExpansion(coll_weight)
                    else:
                        # Normally the collinear current should have built spin-correlations with the leg number corresponding
                        # to the soft one in the context of this soft current
                        assert len(evaluation['spin_correlations'][sc])==1
                        parent_number, spin_corr_vecs = evaluation['spin_correlations'][sc][0]
                        assert soft_leg_number==parent_number
                        assert len(spin_corr_vecs)==1
                        spin_corr_vec = spin_corr_vecs[0]
                        composite_weight += SoftKernels.eikonal_g(self, pi, pk, pr, spin_corr_vector=spin_corr_vec)*EpsilonExpansion(coll_weight)
                new_evaluation['color_correlations'].append( ((a, b), ) )
                new_evaluation['values'][(0,len(new_evaluation['color_correlations'])-1,0)] = composite_weight*mult_factor

        return new_evaluation

def QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx_global_IFF_softFF_variables(reduced_process, all_steps, global_info):
    """ Calls the IF variables by forcing the parent momentum for the IF variable computation to be the
    overall parent momentum of the C(S(C(FF)),I) structure.
    """

    p_a_tilde = all_steps[1]['lower_PS_point'][ global_info['overall_parents'][0] ]
#    misc.sprint(p_a_tilde)
#    misc.sprint(all_steps[1]['higher_PS_point'])
#    misc.sprint(tuple( [ global_info['leg_numbers_map'][1], ] +
#               list(all_steps[1]['bundles_info'][0]['final_state_children']) ))
    IF_variables = kernel_variables.colorful_pp_IFF_softFF_variables(
        all_steps[1]['higher_PS_point'],
        p_a_tilde,
        tuple( [ global_info['leg_numbers_map'][1], ] +
               list(all_steps[1]['bundles_info'][0]['final_state_children']) ),
        Q=global_info['Q']
    )[0]
    IF_variables['p_a_tilde'] = p_a_tilde
    return IF_variables

class QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx(general_current.GeneralCurrent):
    """ Nested soft FF (q_qx) limit within collinear FF (q_qx) limit with collinear limit IFF (q' q_qx)."""

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # In order to build the IF variables using the initial parent momentum which is absent from any mapping
    # structure, we use the global variables
    variables = staticmethod(QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx_global_IFF_softFF_variables)

    # Now define the matching singular structures
    sub_coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -2, sub.SubtractionLeg.FINAL),
        )
    )
    soft_structure = sub.SoftStructure(
            substructures=(sub_coll_structure,),
            legs=tuple([])
    )
    # This counterterm will be used if any of the structures of the list below matches
    # Match the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(soft_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
                'singular_structure'                : structure
            }),
        ),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, +2, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, -2, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({1001:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),]),
        },
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.SoftStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.soft_mapping,
            # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
            # The momenta dictionary below is irrelevant for the soft_mapping used above will make sure that
            # it applies the necessary relabelling of the final-state leg 33 to the parent -1 which will be
            # used by the reduced ME called with it.
            'momenta_dict'          : bidict({-1: frozenset((1,))}),
            'variables'             : None,
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        kT_IF = global_variables['kTs'][0]
        x_IF  = global_variables['xs'][0]
        s_a_rs = global_variables['ss'][(0,1)]

        p_a_tilde = global_variables['p_a_tilde']

        p_rs_hat = all_steps_info[0]['lower_PS_point'][
            all_steps_info[0]['bundles_info'][0]['parent']
        ]

#        misc.sprint(s_rs,s_a_rs)
#        misc.sprint(p_a_tilde,p_rs_hat,p_a_tilde.dot(p_rs_hat))
#        misc.sprint(p_a_tilde,kT_FF,p_a_tilde.dot(kT_FF))
#        misc.sprint(kT_FF, kT_FF.square())
#        misc.sprint(x_IF)
#        misc.sprint(z_FF,(1.-z_FF))
        evaluation['values'][(0,0,0)] = EpsilonExpansion({'finite':
            (2./(s_rs*s_a_rs))*self.TR*self.CF*(
                1./(1.-x_IF) + z_FF * (1. - z_FF) * ((2.*p_a_tilde.dot(kT_FF))**2)/(kT_FF.square()*(2.*p_a_tilde.dot(p_rs_hat)))
            )
        })
        return evaluation

class QCD_S_FqFqx_C_FqFqx_C_IqFqFqx(QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx):
    """ Nested soft FF (q_qx) limit within collinear FF (q_qx) limit with collinear limit IFF (q q_qx).
    Same as QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx except for different flavors in the matching structures.
    """

    # Now define the matching singular structures
    sub_coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )
    soft_structure = sub.SoftStructure(
            substructures=(sub_coll_structure,),
            legs=tuple([])
    )
    # This counterterm will be used if any of the structures of the list below matches
    # Match the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(soft_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color'     : True,
                'n_loops'                           : 0,
                'squared_orders'                    : {'QCD': 4},
                'singular_structure'                : structure
            }),
        ),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({1001:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),]),
        },
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.SoftStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.soft_mapping,
            # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
            # The momenta dictionary below is irrelevant for the soft_mapping used above will make sure that
            # it applies the necessary relabelling of the final-state leg 33 to the parent -1 which will be
            # used by the reduced ME called with it.
            'momenta_dict'          : bidict({-1: frozenset((1,))}),
            'variables'             : None,
            'is_cut'                : None,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]
