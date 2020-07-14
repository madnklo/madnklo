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
"""Implementation of NLO distributed soft currents."""

import math
from bidict import bidict

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings

import commons.utils as utils
import commons.general_current as general_current

import distributed_soft_config as config
import variables as kernel_variables

from commons.universal_kernels import AltarelliParisiKernels, SoftKernels
from distributed_softs_kernels import AltarelliParisiKernels_soft, SoftKernels_soft

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# NLO final-collinear currents
#=========================================================================================

# class QCD_C_FqFqx(general_current.GeneralCurrent):
#     """ FF C(q_qx)"""
#
#     is_zero = True  # is this neaded is called from NLO local currents in test IR limits
#
#     # Enable the flag below to debug this current
#     DEBUG = False
#
#     divide_by_jacobian = config.divide_by_jacobian_local # = False
#
#     # We should not need global variables for this current
#     variables = None
#
#     # Now define the matching singular structures
#     coll_structure = sub.CollStructure(
#         substructures=tuple([]),
#         legs=(
#             sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
#             sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
#         )
#     )
#
#     # This counterterm will be used if any of the current of the list below matches
#     currents = [
#         sub.Current({
#             'resolve_mother_spin_and_color'     : True,
#             'n_loops'                           : 0,
#             'squared_orders'                    : {'QCD': 4}, # should this be 2 ?
#             'singular_structure'                : sub.SingularStructure(substructures=(coll_structure,)),
#         }),
#     ]
#
#     # The defining currents correspond to a currents block composed of several lists of template currents
#     # for matching each currents of the target currents block (in an unordered fashion
#     defining_currents = [ currents, ]
#
#     # An now the mapping rules
#     mapping_rules = [
#         # {
#         #     'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
#         #         substructures=tuple([]),
#         #         legs=(
#         #             sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
#         #             sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
#         #         )
#         #     ),)),
#         #     'mapping'               : config.final_coll_mapping,
#         #     # Intermediate legs should be strictly superior to a 1000
#         #     'momenta_dict'          : bidict({-1:frozenset((10,11))}),
#         #     'variables'             : general_current.CompoundVariables(kernel_variables.distributed_soft_FFn_variables),
#         #     'is_cut'                : config.generalised_cuts,
#         #     'reduced_recoilers'     : config.get_final_state_recoilers,
#         #     'additional_recoilers'  : sub.SubtractionLegSet([]),
#         # },
#     ]
#
#     def kernel(self, evaluation, all_steps_info, global_variables): # fine
#         """ Evaluate this counterterm given the variables provided. """
#
#         # print("local kernel qqx")
#
#         kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
#         z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
#         s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]
#
#         parent = all_steps_info[0]['bundles_info'][0]['parent']
#         # p_rs_hat = all_steps_info[0]['lower_PS_point'][parent]
#
#         # We must include here propagator factors
#         prefactor = 1. /(s_rs * config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True))
#
#         for spin_correlation_vector, weight in AltarelliParisiKernels_soft.P_qq(self, z_FF, kT_FF):
#             complete_weight = weight * prefactor
#             if spin_correlation_vector is None:
#                 evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
#             else:
#                 evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
#                 evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}
#
#         return evaluation

class QCD_C_FqFqpFqpx(general_current.GeneralCurrent):
    """ collinear FSR tree-level current. q(final) > q(final) q'(final) q' (final) """

    # is_zero = True

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = config.divide_by_jacobian_local  # = False

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
            'mapping'               : config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.distributed_soft_FFn_variables),
            'is_cut'                : config.generalised_cuts,
            'reduced_recoilers'     : config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):  # fine
        """ Evaluate this counterterm given the variables provided. """

        # import ipdb
        # ipdb.set_trace()

        print("local kernel Cqqpqpx")
        print(all_steps_info[0]['bundles_info'][0]['final_state_children'])


        # Retrieve the collinear variable z's
        z_1, z_2, z_3 = all_steps_info[0]['variables'][0]['zs'][:3] # 2,3 are same flavour 1 different flavour
        # This current does not involve KTs. :: So are set to None
        # kTs = list()
        #
        # for i in all_steps_info[0]['variables'][0]['kTs']:
        #     kTs.append(all_steps_info[0]['variables'][0]['kTs'][i])

        kT_1, kT_2, kT_3 = None, None, None # kTs

        s_12, s_13, s_23 = all_steps_info[0]['variables'][0]['ss'][(0, 1)], all_steps_info[0]['variables'][0]['ss'][(0, 2)], all_steps_info[0]['variables'][0]['ss'][(1, 2)]
        # parent = all_steps_info[0]['bundles_info'][0]['parent']

        higher_PS_point = all_steps_info[0]['higher_PS_point']

        parent = all_steps_info[0]['bundles_info'][0]['parent']
        # p_rs_hat = all_steps_info[0]['lower_PS_point'][parent]

        # print (all_steps_info)
        # print (global_variables)
        #
        # print (all_steps_info[0]['bundles_info'][0]['final_state_children'])
        # print (higher_PS_point)
        # print (z_i, z_r, z_s)
        # print (s_ir, s_is, s_rs)

        # We must include here propagator factors
        pC = sum(higher_PS_point[child] for child in all_steps_info[0]['bundles_info'][0]['final_state_children'])
        propagator = 1./((pC.square())**2 * config.Jacobian_determinant(all_steps_info, global_variables,all_mapped_masses_are_zero=True))
        # propagator = 2. / ((pC.square()) ** 2 * config.Jacobian_determinant(all_steps_info, global_variables,all_mapped_masses_are_zero=True)) # wrong
        # prefactor = 1. / ((s_ir + s_is + s_rs)**2) # * config.Jacobian_determinant(all_steps_info, global_variables,all_mapped_masses_are_zero=True))
        # use propagator instead ??

        # print ((pC.square())**2, (s_12 + s_13 + s_23)**2)

        for spin_correlation_vector, weight in AltarelliParisiKernels_soft.P_qqxq_sub (
                self, z_1=z_2, z_2=z_3, z_3=z_1, s_12=s_23, s_13=s_12, s_23=s_13, kT_1=kT_2, kT_2=kT_3, kT_3=kT_1):
            complete_weight = weight * propagator
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append(((parent, spin_correlation_vector),))
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a soft type of splitting kernel, which *does* need to know about the reduced process.
        Should be specialised by the daughter class if not dummy
        """

        print("local soft kernel Cqqpqpx")
        print(all_steps_info[0]['bundles_info'][0]['final_state_children'])

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']
        higher_PS_point = all_steps_info[0]['higher_PS_point']

        pS_0 = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][0]]  # first final state children
        pS_1 = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][1]]  # second final state children
        pS_2 = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][2]]  # third final state children

        # Normalization factors
        norm = 2. / config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True) # factor of two to make limits match where from ??

        colinear_parent_leg = all_steps_info[0]['bundles_info'][0]['parent']

        color_correlation_index = 1

        for i, a in enumerate(colored_partons):
            if a == colinear_parent_leg:
                continue
            pa = lower_PS_point[a]  # mapped momentum
            # s_tilde = 2.0 * pa.dot(lower_PS_point[colinear_parent_leg])
            eikonal = SoftKernels_soft.eikonal_qqx(self, pS_1, pS_2, pS_0, pa)  # changed to the eikonal for distributed softs
            evaluation['color_correlations'].append(((colinear_parent_leg, a),))
            evaluation['values'][(0, color_correlation_index, 0, 0)] = {'finite': norm * (eikonal)}
            color_correlation_index += 1

        return evaluation

class QCD_C_FqFqFqx(general_current.GeneralCurrent):
    """ collinear FSR tree-level current. q(final) > q(final) q(final) q(final) """

    is_zero = True  # remove when implemented

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = config.divide_by_jacobian_local  # = False

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
            'mapping'               : config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.distributed_soft_FFn_variables),
            'is_cut'                : config.generalised_cuts,
            'reduced_recoilers'     : config.get_final_state_recoilers,
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

    divide_by_jacobian = config.divide_by_jacobian_local  # = False

    is_zero = True # remove when implemented

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
            'mapping'               : config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11,12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.distributed_soft_FFn_variables),
            'is_cut'                : config.generalised_cuts,
            'reduced_recoilers'     : config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """


        return evaluation

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a soft type of splitting kernel, which *does* need to know about the reduced process.
        Should be specialised by the daughter class if not dummy
        """

        """Important note about the IF CS:
        - in this scheme we want to use "tilded" mapped momenta for the dipole legs in eikonals.
           This is explicitly implemented in the soft local current
        - this implies that the correct form for the local C(ir)S(r) taken as the collinear limit of the eikonals is
        1/ (p_r + p_i_tilde)^2 (1-z_r)/z_r where z_r = p_r.Q/(p_r+p_i_tilde).Q
        - Specializing to the case where the collinear partner of the soft particle is an initial state particle (i = a ), we have
        p_a_tilde = xi p_a and 2p_a.Q = Q^2 so that the soft-collinear takes the form
        1/(p_r+xi p_a)^2 * xi/y_rQ where y_rQ is the usual Hungarian variable
        this simplifies to
        1/(p_r+p_a)^2 * 1/y_rQ which is exactly the soft collinear as computed *without* tilded variables
        (i.e. exactly eq.5.29 of 0903.1218)

        As a result we use exactly the same way of evaluating the counterterms as in honest-to-god colorful.
        """

        # print("local soft kernel gg")

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']
        higher_PS_point = all_steps_info[0]['higher_PS_point']

        # pS_0 = higher_PS_point[
        #     all_steps_info[0]['bundles_info'][0]['final_state_children'][0]]  # first final state children
        # pS_1 = higher_PS_point[
        #     all_steps_info[0]['bundles_info'][0]['final_state_children'][1]]  # second final state children
        # # here both colinears are gluons
        #
        # # Normalization factors
        # norm = -1. / config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True)
        #
        # colinear_parent_leg = all_steps_info[0]['bundles_info'][0]['parent']
        #
        # color_correlation_index = 1
        #
        # for i, a in enumerate(colored_partons):
        #     if a == colinear_parent_leg:
        #         continue
        #     pa = lower_PS_point[a]  # mapped momentum :: should be higher PS point??
        #     s_tilde = 2.0 * pa.dot(lower_PS_point[colinear_parent_leg])
        #     # pa = higher_PS_point[a]
        #     eikonal_0 = SoftKernels_soft.eikonal_dipole_soft_mod(pS_1, pa, pS_0,
        #                                                          s_tilde)  # changed to the eikonal for distributed softs
        #     # in Lionetti's thesis (pS_1 = p_j, pa = p_k, pS_0 = p_i)
        #     # p_i soft momentum
        #     eikonal_1 = SoftKernels_soft.eikonal_dipole_soft_mod(pS_0, pa, pS_1,
        #                                                          s_tilde)  # changed to the eikonal for distributed softs
        #     # in Lionetti's thesis (pS_0 = p_j, pa = p_k, pS_1 = p_i)
        #     # p_i soft momentum
        #
        #     evaluation['color_correlations'].append(((colinear_parent_leg, a),))
        #     evaluation['values'][(0, color_correlation_index, 0, 0)] = {'finite': norm * (eikonal_0 + eikonal_1)}
        #     color_correlation_index += 1

        return evaluation


class QCD_C_FqFqx_C_FqFqx(general_current.GeneralCurrent):
    """ Counterterm for two quark antiquark pairs going colinear"""

    # need this counterterm ??

    # Enable the flag below to debug this current
    DEBUG = False

    is_zero = True # remove when implemented

    divide_by_jacobian = config.divide_by_jacobian_local # = False

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures :: find correct form of the singular structure
    coll_structure_1 = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, -1, sub.SubtractionLeg.FINAL),
        )
    )

    coll_structure_2 = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(3, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(4, -2, sub.SubtractionLeg.FINAL),
        )
    )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 4},
            'singular_structure'                : sub.SingularStructure(substructures=(coll_structure_1,coll_structure_2)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion )
    defining_currents = [ currents, ]

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
            'mapping'               : config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.distributed_soft_FFn_variables),
            'is_cut'                : config.generalised_cuts,
            'reduced_recoilers'     : config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables): # fine
        """ Evaluate this counterterm given the variables provided. """


        return evaluation

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a soft type of splitting kernel, which *does* need to know about the reduced process.
        Should be specialised by the daughter class if not dummy
        """

        """Important note about the IF CS:
        - in this scheme we want to use "tilded" mapped momenta for the dipole legs in eikonals.
           This is explicitly implemented in the soft local current
        - this implies that the correct form for the local C(ir)S(r) taken as the collinear limit of the eikonals is
        1/ (p_r + p_i_tilde)^2 (1-z_r)/z_r where z_r = p_r.Q/(p_r+p_i_tilde).Q
        - Specializing to the case where the collinear partner of the soft particle is an initial state particle (i = a ), we have
        p_a_tilde = xi p_a and 2p_a.Q = Q^2 so that the soft-collinear takes the form
        1/(p_r+xi p_a)^2 * xi/y_rQ where y_rQ is the usual Hungarian variable
        this simplifies to
        1/(p_r+p_a)^2 * 1/y_rQ which is exactly the soft collinear as computed *without* tilded variables
        (i.e. exactly eq.5.29 of 0903.1218)

        As a result we use exactly the same way of evaluating the counterterms as in honest-to-god colorful.
        """

        # print("local soft kernel gg")

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']
        higher_PS_point = all_steps_info[0]['higher_PS_point']

        # pS_0 = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][0]]  # first final state children
        # pS_1 = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][1]]  # second final state children
        # # here both colinears are gluons
        #
        # # Normalization factors
        # norm = -1./config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True)
        #
        # colinear_parent_leg = all_steps_info[0]['bundles_info'][0]['parent']
        #
        # color_correlation_index = 1
        #
        # for i, a in enumerate(colored_partons):
        #     if a == colinear_parent_leg:
        #         continue
        #     pa = lower_PS_point[a] # mapped momentum :: should be higher PS point??
        #     s_tilde = 2.0 * pa.dot(lower_PS_point[colinear_parent_leg])
        #     # pa = higher_PS_point[a]
        #     eikonal_0 = SoftKernels_soft.eikonal_dipole_soft_mod(pS_1, pa, pS_0,s_tilde) # changed to the eikonal for distributed softs
        #     # in Lionetti's thesis (pS_1 = p_j, pa = p_k, pS_0 = p_i)
        #     # p_i soft momentum
        #     eikonal_1 = SoftKernels_soft.eikonal_dipole_soft_mod(pS_0, pa,pS_1,s_tilde)  # changed to the eikonal for distributed softs
        #     # in Lionetti's thesis (pS_0 = p_j, pa = p_k, pS_1 = p_i)
        #     # p_i soft momentum
        #
        #     evaluation['color_correlations'].append(((colinear_parent_leg, a),))
        #     evaluation['values'][(0, color_correlation_index, 0, 0)] = {'finite': norm * (eikonal_0 + eikonal_1)}
        #     color_correlation_index += 1

        return evaluation

#=========================================================================================
# NNLO final-colinear-colinear currents
#=========================================================================================

class QCD_C_FqFqx_C_FpFqpFqpx(general_current.GeneralCurrent):
    """ Nested FF (q_qx) collinear within IFF (q_q'q'x)."""

    # is_zero = True

    divide_by_jacobian = config.divide_by_jacobian_local  # = False

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
            legs=(sub.SubtractionLeg(12, +1, sub.SubtractionLeg.FINAL),)
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
            'mapping'               : config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({1001:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.distributed_soft_FFn_variables),
            'is_cut'                : config.generalised_cuts,
            'reduced_recoilers'     : config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([sub.SubtractionLeg(12, +1, sub.SubtractionLeg.FINAL)]), # sub.SubtractionLegSet([]) #
        },
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(12, +1, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : config.final_coll_mapping,
            # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
            'momenta_dict'          : bidict({-1: frozenset((1001, 12))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.distributed_soft_FFn_variables),
            'is_cut'                : config.generalised_cuts,
            'reduced_recoilers'     : config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """

        print("local kernel Cqqx Cqqpqpx")
        print(all_steps_info[0]['bundles_info'][0]['final_state_children'])
        print(all_steps_info[1]['bundles_info'][0]['final_state_children'])

        # import ipdb
        # ipdb.set_trace()

        kT_1 = all_steps_info[0]['variables'][0]['kTs'][(0, (1,))]
        # kT_2 = all_steps_info[0]['variables'][0]['kTs'][(1, (0,))]
        z_1  = all_steps_info[0]['variables'][0]['zs'][0]
        s_12  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        kT_hat_12 = all_steps_info[1]['variables'][0]['kTs'][(1, (0,))]
        z_hat_12  = all_steps_info[1]['variables'][0]['zs'][1] # parent of double colinear
        z_hat_3   = all_steps_info[1]['variables'][0]['zs'][0]
        s_hat_123 = all_steps_info[1]['variables'][0]['ss'][(0,1)]

        p_rs_hat = all_steps_info[1]['higher_PS_point'][
            all_steps_info[1]['bundles_info'][0]['final_state_children'][0]
        ]

        parent = all_steps_info[1]['bundles_info'][0]['parent']


        # higher_PS_point = all_steps_info[1]['higher_PS_point']
        #
        #
        # children = all_steps_info[1]['bundles_info'][0]['final_state_children']
        # all_p_fs = [higher_PS_point[child] for child in children]
        # zs = []
        # n = config.ref_n  # now with massless reference vector :: before massive Q
        # for p_fs in all_p_fs:
        #     zs.append(p_fs.dot(n))
        # normalisation = sum(zs)
        # for i in range(len(zs)):
        #     zs[i] /= normalisation
        # print(zs)
        # print(z_hat_12, z_hat_3)
        # We must include here propagator factors, but no correction for symmetry
        # or averaging factor is necessary in this particular pure-quark kernel

        prefactor = 1.0/(s_12 * s_hat_123 * config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True))
        # prefactor *= -1.0 # factor tomach limit ?? :: wrong
        for spin_correlation_vector, weight in AltarelliParisiKernels_soft.P_qpqp_q_sub(self, z_1, 1.0 - z_1, z_hat_12, 1.0 - z_hat_12, kT_1, -kT_1, kT_hat_12, -kT_hat_12):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a soft type of splitting kernel, which *does* need to know about the reduced process.
        Should be specialised by the daughter class if not dummy
        """

        print("local soft kernel Cqqx Cqqpqpx")
        print(all_steps_info[0]['bundles_info'][0]['final_state_children'])
        print(all_steps_info[1]['bundles_info'][0]['final_state_children'])

        # import ipdb
        # ipdb.set_trace()

        kT_1 = all_steps_info[0]['variables'][0]['kTs'][(0, (1,))]
        kT_2 = all_steps_info[0]['variables'][0]['kTs'][(1, (0,))]
        z_1  = all_steps_info[0]['variables'][0]['zs'][0]
        z_2 = all_steps_info[0]['variables'][0]['zs'][1]
        s_12  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        Q = global_variables['Q']

        # n = config.ref_n

        # kT_hat_12 = all_steps_info[1]['variables'][0]['kTs'][(1, (0,))]
        # z_hat_12  = all_steps_info[1]['variables'][0]['zs'][1]
        # s_hat_123 = all_steps_info[1]['variables'][0]['ss'][(0,1)]

        # p_rs_hat = all_steps_info[1]['higher_PS_point'][
        #     all_steps_info[1]['bundles_info'][0]['final_state_children'][0]
        # ]

        # parent = all_steps_info[1]['bundles_info'][0]['parent']

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[1]['lower_PS_point']
        higher_PS_point = all_steps_info[1]['higher_PS_point']

        pS_0 = higher_PS_point[all_steps_info[1]['bundles_info'][0]['final_state_children'][0]]  # first final state children
        pS_1 = higher_PS_point[all_steps_info[1]['bundles_info'][0]['final_state_children'][1]]  # second final state children

        # Normalization factors
        norm = 2. / config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True)
        # for counter counterterm ?? :: -1.
        colinear_parent_leg = all_steps_info[1]['bundles_info'][0]['parent']

        # gluon = global_variables['leg_numbers_map'][1001]  # leg number of the gluon

        color_correlation_index = 1

        for i, a in enumerate(colored_partons):
            if a == colinear_parent_leg:
                continue
            pa = lower_PS_point[a]  # mapped momentum
            # pa = higher_PS_point[a] # wrong
            # s_tilde = 2.0 * pa.dot(lower_PS_point[colinear_parent_leg])
            if all_steps_info[1]['bundles_info'][0]['final_state_children'][0] == 1001:  # first is the gluon
                # ipdb.set_trace()
                p_12 = all_steps_info[0]['higher_PS_point'][all_steps_info[0]['bundles_info'][0]['final_state_children'][0]] + \
                       all_steps_info[0]['higher_PS_point'][all_steps_info[0]['bundles_info'][0]['final_state_children'][1]]
                # norm *=  Q.dot(p_12) / Q.dot(pS_0) # pS_0 is p_hat_12
                # print(Q.dot(p_12) / Q.dot(pS_0))
                # norm *= n.dot(p_12) / n.dot(pS_0)  # pS_0 is p_hat_12
                # print(n.dot(p_12) / n.dot(pS_0))
                eikonal = SoftKernels_soft.eikonal_qqx_soft_colinear(
                    self, z_1, z_2, s_12, kT_1, kT_2, p_hat_12=pS_0, p_hat_i=pS_1, p_hat_j=pa)  # changed to the eikonal for distributed softs

            if all_steps_info[1]['bundles_info'][0]['final_state_children'][1] == 1001:  # second is the gluon
                # ipdb.set_trace()
                p_12 = all_steps_info[0]['higher_PS_point'][all_steps_info[0]['bundles_info'][0]['final_state_children'][0]] + \
                       all_steps_info[0]['higher_PS_point'][all_steps_info[0]['bundles_info'][0]['final_state_children'][1]]
                # norm *= Q.dot(p_12) / Q.dot(pS_1) # pS_1 is p_hat_12
                # print(Q.dot(p_12) / Q.dot(pS_1))
                # norm *= n.dot(p_12) / n.dot(pS_1)  # pS_1 is p_hat_12
                # print(n.dot(p_12) / n.dot(pS_1))
                eikonal = SoftKernels_soft.eikonal_qqx_soft_colinear(
                    self, z_1, z_2, s_12, kT_1, kT_2, p_hat_12=pS_1,p_hat_i=pS_0, p_hat_j=pa)  # changed to the eikonal for distributed softs

            evaluation['color_correlations'].append(((colinear_parent_leg, a),))
            evaluation['values'][(0, color_correlation_index, 0, 0)] = {'finite': norm * (eikonal)}
            color_correlation_index += 1

        return evaluation

#=========================================================================================
# NNLO final-soft currents
#=========================================================================================

class QCD_S_FqFqx(general_current.GeneralCurrent):
    """ Soft q(final) q(final) current"""

    # In the distributed soft subtraction scheme, the soft current is split among the collinear sectors.
    is_zero = True

    divide_by_jacobian = config.divide_by_jacobian_local # = False

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
    mapping_rules = [ ]


class QCD_S_FgFg(general_current.GeneralCurrent):
    """ Soft g(final) g(final) current"""

    # In the distributed soft subtraction scheme, the soft current is split among the collinear sectors.
    is_zero = True

    divide_by_jacobian = config.divide_by_jacobian_local # = False

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
    mapping_rules = [ ]

#=========================================================================================
# NLO final-soft-collinear currents
#=========================================================================================

class QCD_S_FqFqx_C_FqFqx(general_current.GeneralCurrent):
    """ Nested soft FF (q_qx) limit within collinear FF (q_qx) limit."""

    is_zero = True # not implemented

    divide_by_jacobian = config.divide_by_jacobian_local  # = False

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
            'mapping'               : config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({1001:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.distributed_soft_FFn_variables),
            'is_cut'                : config.generalised_cuts,
            'reduced_recoilers'     : config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.SoftStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : config.final_coll_mapping,
            # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
            'momenta_dict'          : bidict({}),
            'variables'             : None,
            'is_cut'                : None,
            'reduced_recoilers'     : config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this I(FF) counterterm given the supplied variables. """

        return evaluation

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a collinear type of splitting kernel, which *does* need to know about the reduced process
        Should be specialised by the daughter class if not dummy
        """
        return evaluation

class QCD_S_FqFqx_C_FqFqpFqpx(general_current.GeneralCurrent):
    """ Soft-collinear NNLO current q(initial) > q(initial) q(final) q~(final) and q(initial) > q(initial) q'(final) q~'(final)"""

    # In the distributed soft subtraction scheme, the soft current is split among the collinear sectors.
    is_zero = True

    divide_by_jacobian = config.divide_by_jacobian_local  # = False

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
    coll_structure_q = sub.CollStructure(
        substructures=tuple([soft_structure,]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
        )
    )
    coll_structure_qp = sub.CollStructure(
        substructures=tuple([soft_structure,]),
        legs=(
            sub.SubtractionLeg(10, +2, sub.SubtractionLeg.FINAL),
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

class QCD_S_FqFqx_C_FgFqFqx(general_current.GeneralCurrent):

    # In the distributed soft subtraction scheme, the soft current is split among the collinear sectors.
    is_zero = True

    divide_by_jacobian = config.divide_by_jacobian_local  # = False

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
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
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


class QCD_S_FqFqx_C_FqFqx_C_FqFqpFqpx(general_current.GeneralCurrent):
    """ Nested soft FF (q_qx) limit within collinear FF (q_qx) limit with collinear limit IFF (q' q_qx)."""

    is_zero = True

    divide_by_jacobian = config.divide_by_jacobian_local  # = False

    # In order to build the IF variables using the initial parent momentum which is absent from any mapping
    # structure, we use the global variables
    variables = None # staticmethod(QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx_global_IFF_softFF_variables)

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
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL),)
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
        # {
        #     'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
        #         substructures=tuple([]),
        #         legs=(
        #             sub.SubtractionLeg(10, +2, sub.SubtractionLeg.FINAL),
        #             sub.SubtractionLeg(11, -2, sub.SubtractionLeg.FINAL),
        #         )
        #     ),)),
        #     'mapping'               : config.final_coll_mapping,
        #     # Intermediate legs should be strictly superior to a 1000
        #     'momenta_dict'          : bidict({1001:frozenset((10,11))}),
        #     'variables'             : general_current.CompoundVariables(kernel_variables.distributed_soft_FFn_variables),
        #     'is_cut'                : None,
        #     'reduced_recoilers'     : config.get_final_state_recoilers,
        #     'additional_recoilers'  : sub.SubtractionLegSet([sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL),]),
        # },
        # {
        #     'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
        #         substructures=tuple([]),
        #         legs=(
        #             sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
        #         )
        #     ),)),
        #     'mapping'               : config.final_coll_mapping,
        #     # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
        #     # The momenta dictionary below is irrelevant for the soft_mapping used above will make sure that
        #     # it applies the necessary relabelling of the final-state leg 33 to the parent -1 which will be
        #     # used by the reduced ME called with it.
        #     'momenta_dict'          : bidict({-1: frozenset((1,))}),
        #     'variables'             : None,
        #     'is_cut'                : None,
        #     'reduced_recoilers'     : config.get_final_state_recoilers,
        #     'additional_recoilers'  : sub.SubtractionLegSet([]),
        # },
    ]

    # def kernel(self, evaluation, all_steps_info, global_variables):
    #     """ Evaluate this I(FF) counterterm given the supplied variables. """
    #
    #     # import ipdb
    #     # ipdb.set_trace()
    #
    #     kT_1 = all_steps_info[0]['variables'][0]['kTs'][(0, (1,))]
    #     # kT_2 = all_steps_info[0]['variables'][0]['kTs'][(1, (0,))]
    #     z_1  = all_steps_info[0]['variables'][0]['zs'][0]
    #     s_12  = all_steps_info[0]['variables'][0]['ss'][(0,1)]
    #
    #     kT_hat_12 = all_steps_info[1]['variables'][0]['kTs'][(1, (0,))]
    #     z_hat_12  = all_steps_info[1]['variables'][0]['zs'][1]
    #     s_hat_123 = all_steps_info[1]['variables'][0]['ss'][(0,1)]
    #
    #     p_rs_hat = all_steps_info[1]['higher_PS_point'][
    #         all_steps_info[1]['bundles_info'][0]['final_state_children'][0]
    #     ]
    #
    #     parent = all_steps_info[1]['bundles_info'][0]['parent']
    #
    #     # We must include here propagator factors, but no correction for symmetry
    #     # or averaging factor is necessary in this particular pure-quark kernel
    #
    #     prefactor = 1.0/(s_12 * s_hat_123) # add divide by Jacobian
    #     for spin_correlation_vector, weight in AltarelliParisiKernels_soft.P_S_qpqp_q(self, z_1, 1.0 - z_1, z_hat_12, 1.0 - z_hat_12, kT_1, -kT_1, kT_hat_12, -kT_hat_12):
    #         complete_weight = weight * prefactor
    #         if spin_correlation_vector is None:
    #             evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
    #         else:
    #             evaluation['spin_correlations'].append((parent, spin_correlation_vector))
    #             evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0)] = {'finite': complete_weight[0]}
    #
    #     return evaluation























