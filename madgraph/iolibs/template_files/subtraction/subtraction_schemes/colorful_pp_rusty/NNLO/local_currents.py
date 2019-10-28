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

import commons.QCD_local_currents as currents
import colorful_pp_config
import variables as kernel_variables

from madgraph.core.base_objects import EpsilonExpansion

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# NNLO initial-collinear currents
#=========================================================================================

class QCD_initial_collinear_0_qqpqp(currents.GeneralQCDLocalCurrent):
    """qg collinear FSR tree-level current. q(final) > q(final) q'(final) qbar'(final) """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, -2, sub.SubtractionLeg.FINAL),
        )),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                #substructures=tuple([]), ???
                legs=(
	            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL),
        	    sub.SubtractionLeg(1, +2, sub.SubtractionLeg.FINAL),
     	 	    sub.SubtractionLeg(2, -2, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            'momenta_dict'          : bidict({-1:frozenset((0,1,2))}),
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
            'is_cut'                : colorful_pp_config.initial_coll_cut,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]


    def kernel(self, evaluation, all_steps_info, global_variables):

        kT_a = all_steps_info[0]['variables'][0]['kTs'][0]
        kT_r = all_steps_info[0]['variables'][0]['kTs'][1]
        kT_s = all_steps_info[0]['variables'][0]['kTs'][2]

        x_a = all_steps_info[0]['variables'][0]['xs'][0]
        x_r = all_steps_info[0]['variables'][0]['xs'][1]
        x_s = all_steps_info[0]['variables'][0]['xs'][2]

        s_ar = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        s_as = all_steps_info[0]['variables'][0]['ss'][(0,2)]
        s_rs = all_steps_info[0]['variables'][0]['ss'][(1,2)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./(s_rs - s_ar - s_as)**2
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqpqp(self,
                1./x_a, -x_r/x_a, -x_s/x_a, -s_ar, -s_as, s_rs, kT_a, kT_r, kT_s
            ):
            complete_weight = weight*prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


class QCD_initial_collinear_0_qqq(currents.GeneralQCDLocalCurrent):
    """qg collinear FSR tree-level current. q(final) > q(final) q(final) qbar(final) """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, -1, sub.SubtractionLeg.FINAL),
        )),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                #substructures=tuple([]), ???
                legs=(
	            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL),
        	    sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL),
     	 	    sub.SubtractionLeg(2, -1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            'momenta_dict'          : bidict({-1:frozenset((0,1,2))}),
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
            'is_cut'                : colorful_pp_config.initial_coll_cut,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]


    def kernel(self, evaluation, all_steps_info, global_variables):

        kT_a = all_steps_info[0]['variables'][0]['kTs'][0]
        kT_r = all_steps_info[0]['variables'][0]['kTs'][1]
        kT_s = all_steps_info[0]['variables'][0]['kTs'][2]

        x_a = all_steps_info[0]['variables'][0]['xs'][0]
        x_r = all_steps_info[0]['variables'][0]['xs'][1]
        x_s = all_steps_info[0]['variables'][0]['xs'][2]

        s_ar = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        s_as = all_steps_info[0]['variables'][0]['ss'][(0,2)]
        s_rs = all_steps_info[0]['variables'][0]['ss'][(1,2)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./(s_rs - s_ar - s_as)**2
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqq(self,
                1./x_a, -x_r/x_a, -x_s/x_a, -s_ar, -s_as, s_rs, kT_a, kT_r, kT_s
            ):
            complete_weight = weight*prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


class QCD_initial_collinear_0_qgg(currents.GeneralQCDLocalCurrent):
    """qg collinear FSR tree-level current. q(final) > q(final) g(final) g(final) """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, 21, sub.SubtractionLeg.FINAL),
        )),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                #substructures=tuple([]), ???
                legs=(
	            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL),
        	    sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL),
     	 	    sub.SubtractionLeg(2, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.initial_coll_mapping,
            'momenta_dict'          : bidict({-1:frozenset((0,1,2))}),
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
            'is_cut'                : colorful_pp_config.initial_coll_cut,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]


    def kernel(self, evaluation, all_steps_info, global_variables):

        kT_a = all_steps_info[0]['variables'][0]['kTs'][0]
        kT_r = all_steps_info[0]['variables'][0]['kTs'][1]
        kT_s = all_steps_info[0]['variables'][0]['kTs'][2]

        x_a = all_steps_info[0]['variables'][0]['xs'][0]
        x_r = all_steps_info[0]['variables'][0]['xs'][1]
        x_s = all_steps_info[0]['variables'][0]['xs'][2]

        s_ar = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        s_as = all_steps_info[0]['variables'][0]['ss'][(0,2)]
        s_rs = all_steps_info[0]['variables'][0]['ss'][(1,2)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./(s_rs - s_ar - s_as)**2
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qgg(self,
                1./x_a, -x_r/x_a, -x_s/x_a, -s_ar, -s_as, s_rs, kT_a, kT_r, kT_s
            ):
            complete_weight = weight*prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


#=========================================================================================
# NNLO final-collinear currents
#=========================================================================================
class QCD_final_collinear_0_qqpqp(currents.GeneralQCDLocalCurrent):
    """qg collinear FSR tree-level current. q(final) > q(final) q'(final) qbar'(final) """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(1, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, -2, sub.SubtractionLeg.FINAL),
        )),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                #substructures=tuple([]), ???
                legs=(
	            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        	    sub.SubtractionLeg(1, +2, sub.SubtractionLeg.FINAL),
     	 	    sub.SubtractionLeg(2, -2, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            'momenta_dict'          : bidict({-1:frozenset((0,1,2))}),
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : colorful_pp_config.final_coll_cut,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]


    def kernel(self, evaluation, all_steps_info, global_variables):

        kT_i = all_steps_info[0]['variables'][0]['kTs'][(0,(1,2))]
        kT_r = all_steps_info[0]['variables'][0]['kTs'][(1,(0,2))]
        kT_s = all_steps_info[0]['variables'][0]['kTs'][(2,(0,1))]

        z_i = all_steps_info[0]['variables'][0]['zs'][0]
        z_r = all_steps_info[0]['variables'][0]['zs'][1]
        z_s = all_steps_info[0]['variables'][0]['zs'][2]

        s_ir = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        s_is = all_steps_info[0]['variables'][0]['ss'][(0,2)]
        s_rs = all_steps_info[0]['variables'][0]['ss'][(1,2)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./(s_ir + s_is + s_rs)**2
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqpqp(self,
                z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s
            ):
            complete_weight = weight*prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


class QCD_final_collinear_0_qqq(currents.GeneralQCDLocalCurrent):
    """qg collinear FSR tree-level current. q(final) > q(final) q(final) qbar(final) """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, -1, sub.SubtractionLeg.FINAL),
        )),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                #substructures=tuple([]), ???
                legs=(
	            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        	    sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL),
     	 	    sub.SubtractionLeg(2, -1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            'momenta_dict'          : bidict({-1:frozenset((0,1,2))}),
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : colorful_pp_config.final_coll_cut,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]


    def kernel(self, evaluation, all_steps_info, global_variables):

        kT_i = all_steps_info[0]['variables'][0]['kTs'][(0,(1,2))]
        kT_r = all_steps_info[0]['variables'][0]['kTs'][(1,(0,2))]
        kT_s = all_steps_info[0]['variables'][0]['kTs'][(2,(0,1))]

        z_i = all_steps_info[0]['variables'][0]['zs'][0]
        z_r = all_steps_info[0]['variables'][0]['zs'][1]
        z_s = all_steps_info[0]['variables'][0]['zs'][2]

        s_ir = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        s_is = all_steps_info[0]['variables'][0]['ss'][(0,2)]
        s_rs = all_steps_info[0]['variables'][0]['ss'][(1,2)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./(s_ir + s_is + s_rs)**2
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqq(self,
                z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s
            ):
            complete_weight = weight*prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


class QCD_final_collinear_0_qgg(currents.GeneralQCDLocalCurrent):
    """qg collinear FSR tree-level current. q(final) > q(final) g(final) g(final) """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, 21, sub.SubtractionLeg.FINAL),
        )),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
                #substructures=tuple([]), ???
                legs=(
	            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        	    sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL),
     	 	    sub.SubtractionLeg(2, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.final_coll_mapping,
            'momenta_dict'          : bidict({-1:frozenset((0,1,2))}),
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
            'is_cut'                : colorful_pp_config.final_coll_cut,
            'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]


    def kernel(self, evaluation, all_steps_info, global_variables):

        kT_i = all_steps_info[0]['variables'][0]['kTs'][(0,(1,2))]
        kT_r = all_steps_info[0]['variables'][0]['kTs'][(1,(0,2))]
        kT_s = all_steps_info[0]['variables'][0]['kTs'][(2,(0,1))]

        z_i = all_steps_info[0]['variables'][0]['zs'][0]
        z_r = all_steps_info[0]['variables'][0]['zs'][1]
        z_s = all_steps_info[0]['variables'][0]['zs'][2]

        s_ir = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        s_is = all_steps_info[0]['variables'][0]['ss'][(0,2)]
        s_rs = all_steps_info[0]['variables'][0]['ss'][(1,2)]

        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./(s_ir + s_is + s_rs)**2
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qgg(self,
                z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s
            ):
            complete_weight = weight*prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


# =========================================================================================
# NNLO soft currents
# =========================================================================================
class QCD_S_FqFqx(currents.GeneralQCDLocalCurrent):
    """
    quark-antiquark double soft
    """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    # This counterterm will be used if any of the structures of the list below matches
    structure = [
        sub.SingularStructure(substructures=(sub.SoftStructure(
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
            )
        ),)),
    ]

    # And now the mapping rules
    mapping_rules = [
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.SoftStructure(
                #substructures=tuple([]), ???
                legs=(
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.soft_mapping,
            'momenta_dict'          : bidict({}),
            'variables'             : None,
            #'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables), ??
            'is_cut'                : colorful_pp_config.generalised_cuts,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):

        new_evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations' : [ None, ],
            'color_correlations': [ ],
            'reduced_kinematics': evaluation['reduced_kinematics'],
            'values': { }
        })

        overall_lower_PS_point = all_steps_info[-1]['lower_PS_point']
        soft_leg_number_1 = all_steps_info[-1]['bundles_info'][0]['final_state_children'][0]
        soft_leg_number_2 = all_steps_info[-1]['bundles_info'][0]['final_state_children'][1]
        pr = all_steps_info[-1]['higher_PS_point'][soft_leg_number_1]
        ps = all_steps_info[-1]['higher_PS_point'][soft_leg_number_2]
        colored_parton_numbers = sorted(colored_partons.keys())

	color_correlators_added = {}
	color_correlation_max_index = 0

	overall_factor = 1./(2.*pr.dot(ps))**2

        for i, a in enumerate(colored_parton_numbers):
            for b in colored_parton_numbers[i:]:
                if a!=b:
                    mult_factor = 2.
                else:
                    mult_factor = 1.

                pi = overall_lower_PS_point[a]
                pk = overall_lower_PS_point[b]
		color_correlator = (
                         ( ( (a,-1,a), ), ( (b,-1,b), ) ),
                        )
		the_kernel = SoftKernels.eikonal_qqx(self,pi,pk,pr,ps)*overall_factor*mult_factor
                if color_correlator in color_correlators_added:
                    color_correlation_index = color_correlators_added[color_correlator]
                    new_evaluation['values'][(0, color_correlation_index,0)]['finite'] += the_kernel[0]
                else:
                    new_evaluation['color_correlations'].append(color_correlator)
                    color_correlation_index = color_correlation_max_index
                    color_correlators_added[color_correlator] = color_correlation_max_index
                    color_correlation_max_index += 1
                    new_evaluation['values'][(0, color_correlation_index,0)] = {
								'finite': the_kernel[0]
								}

        return new_evaluation


class QCD_S_FgFg(currents.GeneralQCDLocalCurrent):
    """ Nested soft FF (g_g) limit."""

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    # This counterterm will be used if any of the structures of the list below matches
    structure = [
        sub.SingularStructure(substructures=(sub.SoftStructure(
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
            )
        ),)),
    ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure': sub.SingularStructure(substructures=(sub.SoftStructure(
                #substructures=tuple([]), ???
                legs=(
                    sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : colorful_pp_config.soft_mapping,
            'momenta_dict'          : bidict({}),
            'variables'             : None,
            #'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables), ??
            'is_cut'                : colorful_pp_config.generalised_cuts,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    #copied from Valentine's S_FgFg template
    def create_CataniGrazzini_correlator(self, (i,j),(k,l)):
        """ Returns the correlator of Catani-Grazzini (Eq.113 of hep-ph/9908523v1)
                <M| ( T^-1_i \dot T^-1_j ) * ( T^-1_k \dot T^-1_l ) | M > 
                
            converted into MadGraph's conventions:
            
              ( (a,-1,a),(b,-2,b) ) , ( (c,-1,c),(d,-2,d) ) --> T^-1_a T^-2_b T^-1_c T^-2_d
            """

        # It is important to never commute two color operators acting on the same index, so we 
	# must choose carefully which index we pick to carry the gluon index '-2' of the first 
	# connection. This can be either 'k' or 'l'.
        if j!=k and j!=l:
            # If all indices are different, we can pick either k or l, it is irrelevant
            index1, index2, index3, index4 = i, k, j, l
        elif j==k and j!=l:
            # If j is equal to k, we must pick l
            index1, index2, index3, index4 = i, l, j, k
        elif j==l and j!=k:
            # If j is equal to l, we must pick k
            index1, index2, index3, index4 = i, k, j, l
        elif j==l and j==k:
            # If j is equal to both l and k, then agin it doesn't matter and we can pick k
            index1, index2, index3, index4 = i, k, j, l

        # The sorting according to the first index of each tuple of each of the two convention 
	# is to match Madgraph's convention for sorting color connection in the color correlators 
	# definition
        return (
            tuple(sorted([ (index1,-1,index1), (index2,-2,index2) ], key = lambda el: el[0], reverse=True)),
            tuple(sorted([ (index3,-1,index3), (index4,-2,index4) ], key = lambda el: el[0], reverse=True))
        )

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """
        Should be specialised by the daughter class if not dummy
        """

        new_evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations' : [ None, ],
            'color_correlations': [ ],
            'reduced_kinematics': evaluation['reduced_kinematics'],
            'values': { }
        })

        overall_lower_PS_point = all_steps_info[-1]['lower_PS_point']
        soft_leg_number_1 = all_steps_info[-1]['bundles_info'][0]['final_state_children'][0]
        soft_leg_number_2 = all_steps_info[-1]['bundles_info'][0]['final_state_children'][1]
        pr = all_steps_info[-1]['higher_PS_point'][soft_leg_number_1]
        ps = all_steps_info[-1]['higher_PS_point'][soft_leg_number_2]
        colored_parton_numbers = sorted(colored_partons.keys())

	color_correlators_added = {}
	color_correlation_max_index = 0

        for i, a in enumerate(colored_parton_numbers):
            for b in colored_parton_numbers[i:]:
                if a!=b:
                    mult_factor_1 = 2.
                else:
                    mult_factor_1 = 1.

                pi = overall_lower_PS_point[a]
                pk = overall_lower_PS_point[b]
                lvl1_weight = SoftKernels.eikonal_g(self, pi, pk, pr, spin_corr_vector=None)
		#Implement the non-abelian piece
		non_abelian_correlator = (
                         ( ( (a,-1,a), ), ( (b,-1,b), ) ),
                        )
                non_abelian_kernel = SoftKernels.eikonal_2g(self, pi, pk, pr, ps)*(-self.CA/4.)*mult_factor_1
                if non_abelian_correlator in color_correlators_added:
                    color_correlation_index = color_correlators_added[non_abelian_correlator]
                    new_evaluation['values'][(0, color_correlation_index,0)]['finite'] += non_abelian_kernel[0]
                else:
                    new_evaluation['color_correlations'].append(non_abelian_correlator)
                    color_correlation_index = color_correlation_max_index
                    color_correlators_added[non_abelian_correlator] = color_correlation_max_index
                    color_correlation_max_index += 1
                    new_evaluation['values'][(0, color_correlation_index,0)] = {
								'finite': non_abelian_kernel[0]
								}

                for j, c in enumerate(colored_parton_numbers):
                    for d in colored_parton_numbers[j:]:
                        if c!=d:
                            mult_factor_2 = 2.
                        else:
                            mult_factor_2 = 1.
                        pj = overall_lower_PS_point[c]
                        pl = overall_lower_PS_point[d]
                        lvl2_weight = SoftKernels.eikonal_g(self, pj, pl, ps, spin_corr_vector=None)
			# Implement the abelian piece
                        abelian_correlator_A = ( self.create_CataniGrazzini_correlator((a,b),(c,d)), )
                        abelian_correlator_B = ( self.create_CataniGrazzini_correlator((c,d),(a,b)), )
			abelian_kernel = lvl1_weight*lvl2_weight*mult_factor_1*mult_factor_2/8.

                        for correlator in [abelian_correlator_A, abelian_correlator_B]:
                            if correlator in color_correlators_added:
                                color_correlation_index = color_correlators_added[correlator]
                                new_evaluation['values'][(0, color_correlation_index,0)]['finite'] += abelian_kernel[0]
                            else:
                                new_evaluation['color_correlations'].append(correlator)
                                color_correlation_index = color_correlation_max_index
                                color_correlators_added[correlator] = color_correlation_max_index
                                color_correlation_max_index += 1
                                new_evaluation['values'][(0, color_correlation_index,0)] = {
								'finite': abelian_kernel[0]
								}

        return new_evaluation


# =========================================================================================
# NNLO soft-collinear currents
# =========================================================================================
#Three implementations differing only in singular structure and color charge. These were unified in the colorful_pp implementation and they should be merged once again.
class QCD_initial_soft_collinear_0_qqpqp(currents.GeneralQCDLocalCurrent):
    """q' qbar' soft collinear ISR tree-level current. q(initial) > q(initial_after_emission) Soft(q'(final) qbar' (final)) """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian


    soft_structure = sub.SoftStructure(
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )
    soft_coll_structure_q2 = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(0, +2, sub.SubtractionLeg.INITIAL), ) )
    structure_q2 = sub.SingularStructure(substructures=(soft_coll_structure_q2, ))

    structure = [structure_q2]


    mapping_rules = [
        {
            'singular_structure'    : structure_q2,
            'mapping'               : colorful_pp_config.initial_soft_coll_mapping,
            'momenta_dict'          : bidict({-1:frozenset((0,10,11))}),
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFF_softFF_variables), 
            'is_cut'                : colorful_pp_config.soft_cut,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        kT_a = all_steps_info[0]['variables'][0]['kTs'][0]
        kT_r = all_steps_info[0]['variables'][0]['kTs'][1]
        kT_s = all_steps_info[0]['variables'][0]['kTs'][2]

        x_a = all_steps_info[0]['variables'][0]['xs'][0]
        x_r = all_steps_info[0]['variables'][0]['xs'][1]
        x_s = all_steps_info[0]['variables'][0]['xs'][2]

        s_ar = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        s_as = all_steps_info[0]['variables'][0]['ss'][(0,2)]
        s_rs = all_steps_info[0]['variables'][0]['ss'][(1,2)]
        s_a_rs = s_ar + s_as

        kernel = 2./s_a_rs/s_rs * ( 1. / (x_r + x_s) - (((s_ar * x_s - s_as * x_r) ** 2) / (s_a_rs * s_rs * ((x_r + x_s) ** 2))) )

        evaluation['values'][(0, 0, 0)] = {'finite': self.CF * self.TR * kernel }
        return evaluation


class QCD_initial_soft_collinear_0_qqq(currents.GeneralQCDLocalCurrent):
    """q' qbar' soft collinear ISR tree-level current. q(initial) > q(initial_after_emission) Soft(q'(final) qbar' (final)) """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian


    soft_structure = sub.SoftStructure(
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )
    soft_coll_structure_q1 = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL), ) )
    structure_q1 = sub.SingularStructure(substructures=(soft_coll_structure_q1, ))

    structure = [structure_q1]


    mapping_rules = [
        {
            'singular_structure'    : structure_q1,
            'mapping'               : colorful_pp_config.initial_soft_coll_mapping,
            'momenta_dict'          : bidict({-1:frozenset((0,10,11))}),
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFF_softFF_variables), 
            'is_cut'                : colorful_pp_config.soft_cut,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        kT_a = all_steps_info[0]['variables'][0]['kTs'][0]
        kT_r = all_steps_info[0]['variables'][0]['kTs'][1]
        kT_s = all_steps_info[0]['variables'][0]['kTs'][2]

        x_a = all_steps_info[0]['variables'][0]['xs'][0]
        x_r = all_steps_info[0]['variables'][0]['xs'][1]
        x_s = all_steps_info[0]['variables'][0]['xs'][2]

        s_ar = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        s_as = all_steps_info[0]['variables'][0]['ss'][(0,2)]
        s_rs = all_steps_info[0]['variables'][0]['ss'][(1,2)]
        s_a_rs = s_ar + s_as

        kernel = 2./s_a_rs/s_rs * ( 1. / (x_r + x_s) - (((s_ar * x_s - s_as * x_r) ** 2) / (s_a_rs * s_rs * ((x_r + x_s) ** 2))) )

        evaluation['values'][(0, 0, 0)] = {'finite': self.CF * self.TR * kernel }
        return evaluation


class QCD_initial_soft_collinear_0_gqq(currents.GeneralQCDLocalCurrent):
    """q' qbar' soft collinear ISR tree-level current. g(initial) > g(initial_after_emission) Soft(q'(final) qbar' (final)) """

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian


    soft_structure = sub.SoftStructure(
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )
    soft_coll_structure_g = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.INITIAL), ) )
    structure_g = sub.SingularStructure(substructures=(soft_coll_structure_g, ))

    structure = [structure_g]


    mapping_rules = [
        {
            'singular_structure'    : structure_g,
            'mapping'               : colorful_pp_config.initial_soft_coll_mapping,
            'momenta_dict'          : bidict({-1:frozenset((0,10,11))}),
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFF_softFF_variables), 
            'is_cut'                : colorful_pp_config.soft_cut,
            'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        kT_a = all_steps_info[0]['variables'][0]['kTs'][0]
        kT_r = all_steps_info[0]['variables'][0]['kTs'][1]
        kT_s = all_steps_info[0]['variables'][0]['kTs'][2]

        x_a = all_steps_info[0]['variables'][0]['xs'][0]
        x_r = all_steps_info[0]['variables'][0]['xs'][1]
        x_s = all_steps_info[0]['variables'][0]['xs'][2]

        s_ar = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        s_as = all_steps_info[0]['variables'][0]['ss'][(0,2)]
        s_rs = all_steps_info[0]['variables'][0]['ss'][(1,2)]
        s_a_rs = s_ar + s_as

        kernel = 2./s_a_rs/s_rs * ( 1. / (x_r + x_s) - (((s_ar * x_s - s_as * x_r) ** 2) / (s_a_rs * s_rs * ((x_r + x_s) ** 2))) )

        evaluation['values'][(0, 0, 0)] = {'finite': self.CA * self.TR * kernel }
        return evaluation


#=========================================================================================
#
# Nested currents
#
# These are more advanced currents that require dedicated implementations of the cuts,
# variables and class attributes
#
#=========================================================================================

class QCD_C_FqFqx_C_IqpFqFqx(currents.GeneralQCDLocalCurrent):
    """ Nested FF (q_qx) collinear within IFF (q_q'q'x)."""

    squared_orders = {'QCD': 4}
    n_loops = 0
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
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
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
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
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

#Should be joined with the class above since the subtraction is same for qpqqx and qqqx. Unable to do so because of additional_recoilers in mapping_rules.
class QCD_C_FqFqx_C_IqFqFqx(currents.GeneralQCDLocalCurrent):
    """ Nested FF (q_qx) collinear within IFF (q_qqx)."""

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    sub_coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
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
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
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
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
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

#Needed for QCD_C_IqFqx_C_IqFqFqx where the hatted xs use Q_hat
def QCD_C_IqFqx_C_IqFqFqx_global_IFF_collinearIF_variables(all_steps, global_info):
    """ Calls the IF variables by forcing the parent momentum for the IF variable computation to be the
    overall parent momentum of the C(C(IF),F) structure."""

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


class QCD_C_IqFqx_C_IqFqFqx(currents.GeneralQCDLocalCurrent):
    """ Nested IF (q_qx) collinear within IFF (q_qqx)."""

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
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
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



class QCD_S_FqFqx_C_FqFqx(currents.GeneralQCDLocalCurrent):
    """ Nested soft FF (q_qx) limit within collinear FF (q_qx) limit."""

    squared_orders = {'QCD': 4}
    n_loops = 0
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
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
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


def QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx_global_IFF_softFF_variables(all_steps, global_info):
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


class QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx(currents.GeneralQCDLocalCurrent):
    """ Nested soft FF (q_qx) limit within collinear FF (q_qx) limit with collinear limit IFF (q' q_qx)."""

    squared_orders = {'QCD': 4}
    n_loops = 0
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
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
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

#Inherits everything from the different-flavor class. Structure and mapping_rules redefined to match same-flavor case.
class QCD_S_FqFqx_C_FqFqx_C_IqFqFqx(currents.GeneralQCDLocalCurrent):
    """ Nested soft FF (q_qx) limit within collinear FF (q_qx) limit with collinear limit IFF (q' q_qx)."""

    squared_orders = {'QCD': 4}
    n_loops = 0
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # In order to build the IF variables using the initial parent momentum which is absent from any mapping
    # structure, we use the global variables
    variables = staticmethod(QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx_global_IFF_softFF_variables)

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
    structure = [
        # Match the case of the initial state being a quark and/or antiquark
        sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(soft_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),)),
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
            'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
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

        evaluation['values'][(0,0,0)] = EpsilonExpansion({'finite':
            (2./(s_rs*s_a_rs))*self.TR*self.CF*(
                1./(1.-x_IF) + z_FF * (1. - z_FF) * ((2.*p_a_tilde.dot(kT_FF))**2)/(kT_FF.square()*(2.*p_a_tilde.dot(p_rs_hat)))
            )
        })
        return evaluation

