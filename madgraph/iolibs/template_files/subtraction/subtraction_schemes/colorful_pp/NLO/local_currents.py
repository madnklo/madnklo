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
"""Implementation of NLO colorful_pp local currents."""

import math
from bidict import bidict

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings

import commons.utils as utils
import commons.general_current as general_current
import commons.factors_and_cuts as factors_and_cuts

import colorful_pp_config
import variables as kernel_variables

from commons.universal_kernels import AltarelliParisiKernels, SoftKernels

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# NLO final collinears
#=========================================================================================

class QCD_C_FqFqx(general_current.GeneralCurrent):
    """ FF C(q_qx)"""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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

    # Te defining currents correspond to a currents block composed of several lists of template currents
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
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : colorful_pp_config.generalised_cuts,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']
        p_rs_hat = all_steps_info[0]['lower_PS_point'][parent]

        #misc.sprint(s_rs,)
        #misc.sprint(z_FF)
        #misc.sprint(kT_FF)
        #misc.sprint(p_rs_hat, parent)

        # We must include here propagator factors
        prefactor = 1./s_rs
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qq(self, z_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

class QCD_C_FgFg(general_current.GeneralCurrent):
    """ FF C(g_g)"""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion
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
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : colorful_pp_config.generalised_cuts,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']
        p_rs_hat = all_steps_info[0]['lower_PS_point'][parent]

        # We must include here propagator factors
        prefactor = 1./s_rs
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_gg(self, z_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

class QCD_C_FqFg(general_current.GeneralCurrent):
    """ FF C(q_g)"""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
    coll_structure_q = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
    )
    coll_structure_qx = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, -1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
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

    # Te defining currents correspond to a currents block composed of several lists of template currents
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
            'mapping'               : colorful_pp_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : colorful_pp_config.generalised_cuts,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']
        p_rs_hat = all_steps_info[0]['lower_PS_point'][parent]

        #misc.sprint(s_rs,)
        #misc.sprint(z_FF)
        #misc.sprint(kT_FF)
        #misc.sprint(p_rs_hat, parent)

        # We must include here propagator factors
        prefactor = 1./s_rs
        # The P_gq kernel is the one that matches the z definitions above.
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_gq(self, z_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_S_g(general_current.GeneralCurrent):
    """ F S(q_g)"""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.SoftStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.soft_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({}),
            'variables'             : None,
            'is_cut'                : colorful_pp_config.generalised_cuts,
            # Note that even though the soft mapping recoils against the final state, it still needs to
            # receive the list of *final state* as the recoilers, because those are the ones which will absorb
            # the transverse momentum recoil.
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

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

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']
        higher_PS_point = all_steps_info[0]['higher_PS_point']

        pS = higher_PS_point[all_steps_info[0]['bundles_info'][0]['final_state_children'][0]]

        # Normalization factors
        norm = -1./2.

        # All contributions from the soft counterterms are color correlated so we should remove
        # the default value None from the list of color correlations
        evaluation['color_correlations'] = []

        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        for i, a in enumerate(colored_partons):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in colored_partons[i:]:
                # Write the eikonal for that pair
                if a != b:
                    mult_factor = 2.
                else:
                    mult_factor = 1.
                pa = lower_PS_point[a]
                pb = lower_PS_point[b]
                eikonal = SoftKernels.eikonal_dipole(pa, pb, pS)
                evaluation['color_correlations'].append(((a, b),))
                evaluation['values'][(0, color_correlation_index, 0, 0)] = {
                    'finite': norm * mult_factor * eikonal}
                color_correlation_index += 1

        return evaluation

#=========================================================================================
# NLO final soft-collinear currents
#=========================================================================================

class QCD_CS_FgFg(general_current.GeneralCurrent):
    """ FF C(S(g),g)
    NLO tree-level (final) soft-collinear currents. The momenta used in this current are the
    mapped momenta from the soft mapping."""

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

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # For the gg soft-collinear the color prefactor is CA
    color_factor = "CA"

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.SoftStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.soft_mapping,
            # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
            # The momenta dictionary below is irrelevant for the soft_mapping used above, but it will make sure that
            # it applies the necessary relabelling of the final-state leg 11 to the parent -1 which will be
            # used by the reduced ME called with it.
            'momenta_dict'          : bidict({-1: frozenset((11,))}),
            'variables'             : None,
            'is_cut'                : colorful_pp_config.generalised_cuts,
            # Note that even though the soft mapping recoils against the final state, it still needs to
            # receive the list of *final state* as the recoilers, because those are the ones which will absorb
            # the transverse momentum recoil.
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        soft_leg_number = all_steps_info[0]['bundles_info'][0]['final_state_children'][0]
        p_soft = all_steps_info[0]['higher_PS_point'][soft_leg_number]

        # We could use a custom "global_variables" functional instead of computing the variable directly here
        # but this would be overkill in this case, we choose instead to directly reconstruct them here.
        collinear_parent = global_variables['overall_parents'][0]
        Q = global_variables['Q']
        p_rs_hat = all_steps_info[0]['lower_PS_point'][collinear_parent]

        z = p_soft.dot(Q) / p_rs_hat.dot(Q)

        # We must include here propagator factors
        prefactor = 1./(p_rs_hat+p_soft).square()

        # For now, we hard-code the soft-collinear current here and do not fetch it from universal kernels
        # as this is a sub limit
        color_factor = getattr(self, self.color_factor)
        evaluation['values'][(0, 0, 0, 0)] = {'finite' : prefactor * color_factor * 2.*(1.-z) / z }

        return evaluation

class QCD_CS_FgFq(QCD_CS_FgFg):
    """ FF C(S(g),q)
    NLO tree-level (final) soft-collinear currents. The momenta used in this current are the
    mapped momenta from the soft mapping."""

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(11,  +1, sub.SubtractionLeg.FINAL), ) )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_q,)),
        }),
    ]

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # For the qg soft-collinear the color prefactor is CF
    color_factor = "CF"

    # The rest of the class is identical to that of the g g soft-collinear so we do not need to specify anything further.

#=========================================================================================
# NLO initial-collinear currents
#=========================================================================================


#=========================================================================================
# NLO soft initial-collinear currents
#=========================================================================================

