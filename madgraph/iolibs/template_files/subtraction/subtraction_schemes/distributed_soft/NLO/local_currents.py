"""Implementation of NLO distributed softs local currents."""

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

#Different Splitting kernels

subtract_colinear = False # Set to True to subtract the pure colinear counterterms


#=========================================================================================
# NLO final collinears
#=========================================================================================

class QCD_C_FqFg(general_current.GeneralCurrent):
    """ FF C(q_g)"""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = config.divide_by_jacobian_local # = False

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
    # Two structures for quark and antiquark (+1 or -1) should deal with both cases
    coll_structure_q = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL), # gluon
        )
    )
    coll_structure_qx = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, -1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL), # gluon
        )
    )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(coll_structure_q,)),
        }),
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(coll_structure_qx,)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
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

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        # print("local kernel qg")

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))] # reads out the information from all_steps_info
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']
        # p_rs_hat = all_steps_info[0]['lower_PS_point'][parent]

        gluon = global_variables['leg_numbers_map'][11]  # leg number of the gluon

        # We must include here propagator factors
        prefactor = 1./(s_rs * config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True))

        if all_steps_info[0]['bundles_info'][0]['final_state_children'][0] == gluon:  # first is the gluon
            # The P_gq kernel is the one that matches the z definitions above.
            for spin_correlation_vector, weight in AltarelliParisiKernels_soft.P_gq(self, z_FF, kT_FF): # Spin correlation vectoror = None, kt
                complete_weight = weight * prefactor
                if spin_correlation_vector is None:
                    evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
                else:
                    evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                    evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

            if subtract_colinear:
                for spin_correlation_vector, weight in AltarelliParisiKernels.P_gq(self, z_FF,kT_FF):  # Spin correlation vectoror = None, kt
                    complete_weight = weight * prefactor
                    if spin_correlation_vector is None:
                        evaluation['values'][(0, 0, 0, 0)]['finite'] -= complete_weight[0]
                    else:
                        evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)]['finite'] -= complete_weight[0]

        if all_steps_info[0]['bundles_info'][0]['final_state_children'][1] == gluon:  # second is the gluon
            # The P_gq kernel is the one that matches the z definitions above.
            for spin_correlation_vector, weight in AltarelliParisiKernels_soft.P_qg(self, z_FF, kT_FF): # Spin correlation vectoror = None, kt
                complete_weight = weight * prefactor
                if spin_correlation_vector is None:
                    evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
                else:
                    evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                    evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

            if subtract_colinear:
                for spin_correlation_vector, weight in AltarelliParisiKernels.P_qg(self, z_FF,kT_FF):  # Spin correlation vectoror = None, kt
                    complete_weight = weight * prefactor
                    if spin_correlation_vector is None:
                        evaluation['values'][(0, 0, 0, 0)]['finite'] -= complete_weight[0]
                    else:
                        evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)]['finite'] -= complete_weight[0]

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

        # print("local soft kernel qg")

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']
        higher_PS_point = all_steps_info[0]['higher_PS_point']

        gluon = global_variables['leg_numbers_map'][11] # leg number of the gluon

        pS_0 = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][0]] #first final state children
        pS_1 = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][1]] #second final state children

        # Normalization factors
        norm = -1./config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True) # change this ??

        # All contributions from the soft counterterms are color correlated so we should remove
        # the default value None from the list of color correlations
        #evaluation['color_correlations'] = []

        colinear_parent_leg = all_steps_info[0]['bundles_info'][0]['parent'] # leg number of the parent of the two colinears

        color_correlation_index = 1

        #Loop over all final state hard partons exept the parent of the colinears

        for i, a in enumerate(colored_partons): # i is unused
            if a == colinear_parent_leg:
                continue
            pa = lower_PS_point[a] # mapped momentum :: should be higher PS point??
            s_tilde = 2.0 * pa.dot(lower_PS_point[colinear_parent_leg])
            # pa = higher_PS_point[a]
            if all_steps_info[0]['bundles_info'][0]['final_state_children'][0] == gluon: #first is the gluon
                eikonal = SoftKernels_soft.eikonal_dipole_soft_mod(pS_1, pa, pS_0, s_tilde) # changed to the eikonal for distributed softs
                # in Lionetti's thesis (pS_1 = p_j, pa = p_k, pS_0 = p_i)
                # p_i soft momentum

            if all_steps_info[0]['bundles_info'][0]['final_state_children'][1] == gluon:  # second is the gluon
                eikonal = SoftKernels_soft.eikonal_dipole_soft_mod(pS_0, pa,pS_1,s_tilde)  # changed to the eikonal for distributed softs
                # in Lionetti's thesis (pS_0 = p_j, pa = p_k, pS_1 = p_i)
                # p_i soft momentum

            evaluation['color_correlations'].append(((colinear_parent_leg, a),))
            evaluation['values'][(0, color_correlation_index, 0, 0)] = {'finite': norm * eikonal}
            color_correlation_index += 1

        return evaluation

class QCD_C_FqFqx(general_current.GeneralCurrent):
    """ FF C(q_qx)"""

    # is_zero = True # remove later

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = config.divide_by_jacobian_local # = False

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(coll_structure,)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion
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

        print("local kernel qqx")
        print(all_steps_info[0]['bundles_info'][0]['final_state_children'])

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']
        # p_rs_hat = all_steps_info[0]['lower_PS_point'][parent]

        # We must include here propagator factors
        prefactor = 1./(s_rs * config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True))

        for spin_correlation_vector, weight in AltarelliParisiKernels_soft.P_qq(self, z_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        if subtract_colinear:
            for spin_correlation_vector, weight in AltarelliParisiKernels.P_qq(self, z_FF, kT_FF):
                complete_weight = weight * prefactor
                if spin_correlation_vector is None:
                    evaluation['values'][(0, 0, 0, 0)]['finite'] -= complete_weight[0]
                else:
                    evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)]['finite'] -= complete_weight[0]

        return evaluation


class QCD_C_FgFg(general_current.GeneralCurrent):
    """ FF C(g_g)"""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = config.divide_by_jacobian_local # = False

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
    )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(coll_structure,)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
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

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        # print("local kernel gg")

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']
        # p_rs_hat = all_steps_info[0]['lower_PS_point'][parent]

        # We must include here propagator factors
        prefactor = 1./(s_rs * config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True)) # divergent behavior of the counterterm

        for spin_correlation_vector, weight in AltarelliParisiKernels_soft.P_gg(self, z_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        if subtract_colinear:
            for spin_correlation_vector, weight in AltarelliParisiKernels.P_gg(self, z_FF, kT_FF):
                complete_weight = weight * prefactor
                if spin_correlation_vector is None:
                    evaluation['values'][(0, 0, 0, 0)]['finite'] -= complete_weight[0]
                else:
                    evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)]['finite'] -= complete_weight[0]

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

        pS_0 = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][0]]  # first final state children
        pS_1 = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][1]]  # second final state children
        # here both colinears are gluons

        # Normalization factors
        norm = -1./config.Jacobian_determinant(all_steps_info, global_variables, all_mapped_masses_are_zero=True)

        colinear_parent_leg = all_steps_info[0]['bundles_info'][0]['parent']

        color_correlation_index = 1

        for i, a in enumerate(colored_partons):
            if a == colinear_parent_leg:
                continue
            pa = lower_PS_point[a] # mapped momentum :: should be higher PS point??
            s_tilde = 2.0 * pa.dot(lower_PS_point[colinear_parent_leg])
            # pa = higher_PS_point[a]
            eikonal_0 = SoftKernels_soft.eikonal_dipole_soft_mod(pS_1, pa, pS_0,s_tilde) # changed to the eikonal for distributed softs
            # in Lionetti's thesis (pS_1 = p_j, pa = p_k, pS_0 = p_i)
            # p_i soft momentum
            eikonal_1 = SoftKernels_soft.eikonal_dipole_soft_mod(pS_0, pa,pS_1,s_tilde)  # changed to the eikonal for distributed softs
            # in Lionetti's thesis (pS_0 = p_j, pa = p_k, pS_1 = p_i)
            # p_i soft momentum

            evaluation['color_correlations'].append(((colinear_parent_leg, a),))
            evaluation['values'][(0, color_correlation_index, 0, 0)] = {'finite': norm * (eikonal_0 + eikonal_1)}
            color_correlation_index += 1

        return evaluation


#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_S_g(general_current.GeneralCurrent):
    """ F S(g). This current is zero in the distributed soft subtraction scheme."""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = config.divide_by_jacobian_local

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    soft_structure = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
    )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_structure,)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [ ]

    # In the distributed soft subtraction scheme, the soft current is split among the collinear sectors.
    is_zero = True


#=========================================================================================
# NLO final soft-collinear currents
#=========================================================================================

class QCD_CS_FgFq(general_current.GeneralCurrent): # changed to general current
    """ FF C(S(g),q)
    NLO tree-level (final) soft-collinear currents. The momenta used in this current are the
    mapped momenta from the soft mapping."""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = config.divide_by_jacobian_local  # = False; Do not need this since False is default

    # We should not need global variables for this current
    variables = None # Default value

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) ) # soft gluon leg 11
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(11,  +1, sub.SubtractionLeg.FINAL), ) ) # gluon collinear tu quark leg 10

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_q,)), # colinear or soft_colinear?
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = []

    # In the distributed soft subtraction scheme, the soft current is split among the collinear sectors.
    is_zero = True

class QCD_CS_FgFg(general_current.GeneralCurrent):
    """ FF C(S(g),g)
    NLO tree-level (final) soft-collinear currents. The momenta used in this current are the
    mapped momenta from the soft mapping."""

    divide_by_jacobian = config.divide_by_jacobian_local

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )

    soft_coll_structure_g = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL), ) )


    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_g,)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = []

    # In the distributed soft subtraction scheme, the soft current is split among the collinear sectors.
    is_zero = True
















