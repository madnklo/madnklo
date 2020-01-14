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

# =========================================================================================
# NNLO final-collinear currents
# =========================================================================================

# =========================================================================================
# NNLO soft currents
# =========================================================================================

# =========================================================================================
# NNLO soft-collinear currents
# =========================================================================================


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
    # This counterterm will be used if any of the structures of the list below matches
    structure = [
        # Match both the case of the initial state being a quark and/or antiquark
        sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(sub_coll_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),)),
    ]

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
            'is_cut'                : colorful_pp_config.generalised_cuts,
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
            'is_cut'                : colorful_pp_config.generalised_cuts,
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
    structure = [
        # Match both the case of the initial state being a quark and/or antiquark
        sub.SingularStructure(substructures=(sub.SoftStructure(
            substructures=(sub_coll_structure,),
            legs=tuple([])
        ),)),
    ]

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
            'is_cut'                : colorful_pp_config.generalised_cuts,
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
            'is_cut'                : colorful_pp_config.generalised_cuts,
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
    overall parent momentum of the C(S(C(FF)),I) structure."""

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
    structure = [
        # Match the case of the initial state being a quark and/or antiquark
        sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(soft_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),)),
    ]

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
            'is_cut'                : colorful_pp_config.generalised_cuts,
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
            'is_cut'                : colorful_pp_config.generalised_cuts,
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