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
"""Implementation of NLO local currents for the Torino subtraction scheme
https://arxiv.org/abs/1806.09570
"""
# Ez
import logging
logger1 = logging.getLogger('madgraph')

import math
from bidict import bidict

import madgraph.core.subtraction as sub
from madgraph.core.base_objects import EpsilonExpansion
import madgraph.integrator.mappings as mappings

import commons.utils as utils
import commons.general_current as general_current
import dipole_current as dipole_current
import madgraph.integrator.vectors as vectors


import torino_config
import variables as kernel_variables

from commons.universal_kernels import AltarelliParisiKernels, SoftKernels

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError


def compensate_sector_wgt(PS_point, global_variables, CT_type):
    """returns W_ab^X / W_ab, where X=S,C,SC for CT_type = 1,2,3.
     This is done in order to have the correct sector function together with the CT
     """
    return global_variables['sector_info'][0](PS_point, [], sector_type=CT_type) / \
            global_variables['sector_info'][0](PS_point, [])


#=========================================================================================
# NLO final collinears
#=========================================================================================


class QCD_TRN_C_FgFg(general_current.GeneralCurrent):
    """ FF C(g_g)"""

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
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
            'mapping'               : torino_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_FFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][1]
        x_FF  = all_steps_info[0]['variables'][0]['xs'][0]
        x_oth = all_steps_info[0]['variables'][0]['xs'][1]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./s_rs #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 2)
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_gg(self, x_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


class QCD_TRN_C_FqFqx(general_current.GeneralCurrent):
    """ FF C(q_qx)"""

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
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
    # for matching each currents of the target currents block (in an unordered fashion)
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
            'mapping'               : torino_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_FFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][1]
        x_FF  = all_steps_info[0]['variables'][0]['xs'][0]
        x_oth = all_steps_info[0]['variables'][0]['xs'][1]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./s_rs #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 2)
        # include the soft_collinear counterterm here, as in the torino paper
        # (see the definition of 'hard-collinear' splitting function there)
        soft_col = 0.
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qq(self, x_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                complete_weight += soft_col * (- prefactor)
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation



class QCD_TRN_C_FgFq(general_current.GeneralCurrent):
    """ FF C(q_g)"""

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
    coll_structure_q = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
        )
    )
    coll_structure_qx = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
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
                    sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : torino_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_FFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][0]
        x_FF  = all_steps_info[0]['variables'][0]['xs'][1]
        x_oth = all_steps_info[0]['variables'][0]['xs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./s_rs #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 2)
        # include the soft_collinear counterterm here, as in the torino paper
        # (see the definition of 'hard-collinear' splitting function there)
#        soft_col = EpsilonExpansion({0: self.CF * 2 * x_FF / x_oth, 1: 0.})
        soft_col = 0
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qg(self, x_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                complete_weight += soft_col * (- prefactor)
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_TRN_S_g(dipole_current.DipoleCurrent):
    """ F S(q_g)"""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = torino_config.divide_by_jacobian

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
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.SoftStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : torino_config.soft_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({}),
            'variables'             : None,
            'is_cut'                : torino_config.generalised_cuts,
            # Note that even though the soft mapping recoils against the final state, it still needs to
            # receive the list of *final state* as the recoilers, because those are the ones which will absorb
            # the transverse momentum recoil.
            'reduced_recoilers'     : torino_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a soft type of splitting kernel, which *does* need to know about the reduced process.
        Should be specialised by the daughter class if not dummy
        """

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

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
                    continue # MZ as long as we just care for massless partons, this is irrelevant.
                    # note that this has to be taken care of for the case of massive one, in
                    # particular how to identify the extra recoiler for the momentum mapping
                    mult_factor = 1.

                # Retrieve the reduced kinematics as well as the soft momentum
                lower_PS_point = all_steps_info[(a, b)][0]['lower_PS_point']
                higher_PS_point = all_steps_info[(a, b)][0]['higher_PS_point']

                sector_prefactor = compensate_sector_wgt(higher_PS_point, global_variables, 1)

                ###pS = higher_PS_point[all_steps_info[(a, b)][0]['bundles_info'][0]['final_state_children'][0]]
                # the above expression cannot be used any more since we use a mapping structure of collinear
                #  type also for this current
                pS = higher_PS_point[ self.map_leg_number(
                                    global_variables['leg_numbers_map'],
                                    self.soft_structure.legs[0].n,
                                    global_variables['overall_parents'])]

                pa = higher_PS_point[a]
                pb = higher_PS_point[b]
                # some special behaviour may be needed in the case
                # either a or b are one of the particles defining
                # the sector, in particular for the SC limit.

                # At the moment, we do not implement anything special
                eikonal = SoftKernels.eikonal_dipole(pa, pb, pS)
                evaluation['color_correlations'].append(((a, b),))
                evaluation['values'][(0, color_correlation_index, 0, color_correlation_index)] = {
                    'finite': norm * mult_factor * eikonal * sector_prefactor}
                evaluation['reduced_kinematics'].append(('Dip %d-%d' % (a, b), lower_PS_point))
                color_correlation_index += 1

        return evaluation

#=========================================================================================
# NLO final soft-collinear currents
#=========================================================================================


class QCD_TRN_CS_FgFq(general_current.GeneralCurrent):
    """ FF C(S(g),q)
    NLO tree-level (final) soft-collinear currents. Returns zero, since we have already subtracted the
    soft-col singularity inside the splitting functions"""

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(11,  +1, sub.SubtractionLeg.FINAL), ) )
#    soft_coll_structure_qx = sub.CollStructure(
#        substructures=(soft_structure, ),
#        legs=(sub.SubtractionLeg(11,  -1, sub.SubtractionLeg.FINAL), ) )


#    is_zero = True

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_q,)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=[sub.SoftStructure(
                               legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )],
                legs=(sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),)
            ),)),
            'mapping'               : torino_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_FFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        x_soft  = all_steps_info[0]['variables'][0]['xs'][0]
        x_oth = all_steps_info[0]['variables'][0]['xs'][1]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        prefactor = 1./s_rs #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 3)
        # include the soft_collinear counterterm here, as in the torino paper
        # (see the definition of 'hard-collinear' splitting function there)
        soft_col = EpsilonExpansion({0: self.CF * 2 * x_oth / x_soft, 1: 0.})
        complete_weight = soft_col * prefactor
        evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


class QCD_TRN_CS_FgFg(general_current.GeneralCurrent):
    """ FF C(S(g),g)
    NLO tree-level (final) soft-collinear currents. """

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
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=[sub.SoftStructure(
                    legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )],
                legs=(sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),)
            ),)),
            'mapping'               : torino_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_FFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]


    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        x_soft  = all_steps_info[0]['variables'][0]['xs'][0]
        x_oth = all_steps_info[0]['variables'][0]['xs'][1]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        prefactor = 1./s_rs #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 3)
        # include the soft_collinear counterterm here, as in the torino paper
        # (see the definition of 'hard-collinear' splitting function there)
        soft_col = EpsilonExpansion({0: self.CA * 2 * x_oth / x_soft, 1: 0.})
        complete_weight = soft_col * prefactor
        evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


#=========================================================================================
# NLO initial-collinear currents
#=========================================================================================

class QCD_TRN_C_IgFg(general_current.GeneralCurrent):
    """ IF C(g_g) """

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
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
                    sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : torino_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_IFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][1]
        x_FF  = all_steps_info[0]['variables'][0]['xs'][0]
        x_oth = all_steps_info[0]['variables'][0]['xs'][1]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./s_rs/x_FF #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 2)
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_gg(self, x_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


class QCD_TRN_C_IqFq(general_current.GeneralCurrent):
    """ IF C(q_q) """

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
    coll_structure_q = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
        )
    )
    coll_structure_qx = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, -1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color': True,
            'n_loops': 0,
            'squared_orders': {'QCD': 2},
            'singular_structure': sub.SingularStructure(substructures=(coll_structure_q,)),
        }),
        sub.Current({
            'resolve_mother_spin_and_color': True,
            'n_loops': 0,
            'squared_orders': {'QCD': 2},
            'singular_structure': sub.SingularStructure(substructures=(coll_structure_qx,)),
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
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : torino_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_IFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][1]
        x_FF  = all_steps_info[0]['variables'][0]['xs'][0]
        x_oth = all_steps_info[0]['variables'][0]['xs'][1]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./s_rs/x_FF #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 2)
        # include the soft_collinear counterterm here, as in the torino paper
        # (see the definition of 'hard-collinear' splitting function there)
        soft_col = 0.
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_gq(self, x_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                complete_weight += soft_col * (- prefactor)
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

class QCD_TRN_C_IgFq(general_current.GeneralCurrent):
    """ IF C(g_q) """

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
    coll_structure_q = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
        )
    )
    coll_structure_qx = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
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
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : torino_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_IFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][0]
        x_FF  = all_steps_info[0]['variables'][0]['xs'][0]
        x_oth = all_steps_info[0]['variables'][0]['xs'][1]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

        prefactor = 1./s_rs/x_FF #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 2)
        # include the soft_collinear counterterm here, as in the torino paper
        # (see the definition of 'hard-collinear' splitting function there)
#        soft_col = EpsilonExpansion({0: self.CF * 2 * x_FF / x_oth, 1: 0.})
        soft_col = 0.
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qq(self, x_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                complete_weight += soft_col * (- prefactor)
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

class QCD_TRN_C_IqFg(general_current.GeneralCurrent):
    """ IF C(q_g) """

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
    coll_structure_q = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
    )
    coll_structure_qx = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, -1, sub.SubtractionLeg.INITIAL),
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

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : torino_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_IFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][0]
# Ez
#        x_FF  = all_steps_info[0]['variables'][0]['xs'][1]
#        x_oth = all_steps_info[0]['variables'][0]['xs'][0]
        x_FF  = all_steps_info[0]['variables'][0]['xs'][0]
        x_oth = all_steps_info[0]['variables'][0]['xs'][1]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]
        parent = all_steps_info[0]['bundles_info'][0]['parent']

# Ez
#        logger1.info("")
#        logger1.info("local_currents.QCD_TRN_C_IqFg.kernel")
#        logger1.info("x_FF: "+str(x_FF))
#        logger1.info("x_oth: "+str(x_oth)+"   1-x_FF: "+str(1-x_FF))
#        logger1.info("all_steps_info:\n"+str(all_steps_info))

        prefactor = 1./s_rs/x_FF #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 2)
        # include the soft_collinear counterterm here, as in the torino paper
        # (see the definition of 'hard-collinear' splitting function there)
#        soft_col = EpsilonExpansion({0: self.CF * 2 * x_FF / x_oth, 1: 0.})
        soft_col = 0.
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qg(self, x_FF, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                complete_weight += soft_col * (- prefactor)
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        # Ez
        #logger1.info("QCD_TRN_C_IqFg.kernel")
        #logger1.info("prefactor: " + str(prefactor) + "  x_FF: " + str(x_FF) + "  x_oth: " + str(x_oth))
        #logger1.info("complete_weight[0]: " + str(prefactor * self.CF * 2 * x_FF / x_oth))

        return evaluation

#=========================================================================================
# NLO soft initial-collinear currents
#=========================================================================================

class QCD_TRN_CS_IqFg(general_current.GeneralCurrent):
    """ FI C(S(g),q)
    NLO tree-level (final) soft-collinear currents. Returns zero, since we have already subtracted the
    soft-col singularity inside the splitting functions"""

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(11,  +1, sub.SubtractionLeg.INITIAL), ) )
#    soft_coll_structure_qx = sub.CollStructure(
#        substructures=(soft_structure, ),
#        legs=(sub.SubtractionLeg(11,  -1, sub.SubtractionLeg.INITIAL), ) )

#    is_zero = True

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_q,)),
        }),
#        sub.Current({
#            'resolve_mother_spin_and_color'     : True,
#            'n_loops'                           : 0,
#            'squared_orders'                    : {'QCD': 2},
#            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_qx,)),
#        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=[sub.SoftStructure(
                               legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )],
                legs=(sub.SubtractionLeg(11, +1, sub.SubtractionLeg.INITIAL),)
            ),)),
            'mapping'               : torino_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_IFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        x_soft  = all_steps_info[0]['variables'][0]['xs'][1]
        x_oth = all_steps_info[0]['variables'][0]['xs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        prefactor = 1./s_rs/x_oth #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 3)
        # include the soft_collinear counterterm here, as in the torino paper
        # (see the definition of 'hard-collinear' splitting function there)
        soft_col = EpsilonExpansion({0: self.CF * 2 * x_oth / x_soft, 1: 0.})
        complete_weight = soft_col * prefactor
        evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}

# Ez
        #logger1.info("QCD_TRN_CS_IqFg.kernel")
        #logger1.info("prefactor: "+str(prefactor)+"  x_oth: "+str(x_oth)+"  x_soft: "+str(x_soft))
        #logger1.info("complete_weight[0]: "+str(complete_weight[0]))

        return evaluation



class QCD_TRN_CS_IgFg(general_current.GeneralCurrent):
    """ FI C(S(g),g)
    NLO tree-level (final) soft-collinear currents. """

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(11, 21, sub.SubtractionLeg.INITIAL), ) )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_q,)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=[sub.SoftStructure(
                    legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )],
                legs=(sub.SubtractionLeg(11, 21, sub.SubtractionLeg.INITIAL),)
            ),)),
            'mapping'               : torino_config.initial_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_IFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_initial_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]


    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        x_soft  = all_steps_info[0]['variables'][0]['xs'][1]
        x_oth = all_steps_info[0]['variables'][0]['xs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        prefactor = 1./s_rs/x_oth #
        prefactor *= compensate_sector_wgt(all_steps_info[0]['higher_PS_point'], global_variables, 3)
        # include the soft_collinear counterterm here, as in the torino paper
        # (see the definition of 'hard-collinear' splitting function there)
        soft_col = EpsilonExpansion({0: self.CA * 2 * x_oth / x_soft, 1: 0.})
        complete_weight = soft_col * prefactor
        evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

