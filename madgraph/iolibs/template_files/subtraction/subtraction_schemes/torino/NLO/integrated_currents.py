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
"""Implementation of NLO type of integrated currents for Torino subtraction"""

import os
import sys
import math

import madgraph.core.subtraction as sub
import madgraph.integrator.phase_space_generators as PS_utils
from madgraph.core.base_objects import EpsilonExpansion
import madgraph.various.misc as misc

import commons.utils as utils
import commons.QCD_local_currents as currents
import commons.factors_and_cuts as factors_and_cuts

from bidict import bidict

import commons.general_current as general_current

import torino_config
import variables as kernel_variables

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError

# a general comment: wrt the torino paper https://arxiv.org/pdf/1806.09570.pdf,
# where the following prefactor is employed:
#     (mu^2/s)^eps ,
# MadNkLO uses a different pre-factor:
#     (mu^2/Qes^2)^eps .
# In all the currents we will keep the torino-paper convention, and only
# at the very end convert to the MadNkLO one with the torino_to_madnk_epsexp function


def torino_to_madnk_epsexp(torinopoles, mu2os):
    """a function to convert an epsilon-expansion with the conventions
    of the Torino paper to the conventions used in this code
    Note that, besides the log(mu^2/s), there is also a term coming from
    the alpha_s renormalisation convention leading to an extra factor
    exp(epsilon EluerGamma) / Gamma(1-eps)
    """
    try:
        double = torinopoles[-2]
    except KeyError:
        double = 0.

    try:
        single = torinopoles[-1]
    except KeyError:
        single = 0.

    try:
        finite = torinopoles[0]
    except KeyError:
        finite = 0.

    lmu2os = math.log(mu2os)

    return EpsilonExpansion({-2: double,
                             -1: single + double * lmu2os,
                              0: finite + single * lmu2os +
                                  double * (lmu2os**2 / 2. + math.pi**2 / 12.)})


def variables_for_integrated_currents_FF(cls, reduced_process, all_steps, global_variables):
    """ the variables function for the integrated counterterms.
    It returns the recoiler number, needed in the integrated counterterms.
    It uses the get_final_state_recoilers function of torino_config
    """

    recoilers = torino_config.get_final_state_recoilers(
        reduced_process,
        excluded = global_variables['overall_parents'],
        global_variables = global_variables)

    return {'recoiler' : recoilers[0].n}


class QCD_integrated_TRN_C_FgFg(general_current.GeneralCurrent):
    """This CT corresponds to the hard-collinear integrated CT for the gg splitting
     (Eq. 2.50 in the torino paper, with p=g, keeping only the CA-proportional term),
       plus the soft-collinear contribution
     (the first line of eq 2.54 in the torino paper, with CA as prefactor)"""

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # variables should be updated
    variables = variables_for_integrated_currents_FF

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
    mapping_rules = []


    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided.
        """

        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
            all_steps_info[0]['lower_PS_point'][2])
        mu_r = global_variables['mu_r']
        mu2os = mu_r**2 / shat

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
        etapr = (p_p + p_r).square() / shat

        # hard collinear part
        overall = -self.CA / 6.
        double = 0.
        single = overall
        finite = overall * (8. / 3. - math.log(etapr))

        # soft collinear part
        overall = self.CA
        zeta2 = math.pi ** 2 / 6.
        double += overall
        single += overall * 2
        finite += overall * (6 - 7. / 2. * zeta2)

        evaluation['values'][(0, 0, 0, 0)] = torino_to_madnk_epsexp(
                                            EpsilonExpansion({-2: double,
                                                              -1: single,
                                                               0: finite}) * self.SEpsilon * (1. / (16 * math.pi**2)),
                                                                  mu2os)

        return evaluation


class QCD_integrated_TRN_C_FqFqx(general_current.GeneralCurrent):
    """this ct corresponds to eq. 2.50 in the torino paper, with p=g,
        keeping only the TR-proportional term. no Nf is included"""

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # variables should be updated
    variables = variables_for_integrated_currents_FF

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
    mapping_rules = []


    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided.
        """

        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
            all_steps_info[0]['lower_PS_point'][2])
        mu_r = global_variables['mu_r']
        mu2os = mu_r**2 / shat

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
        etapr = (p_p + p_r).square() / shat

        overall = -self.TR * 2 / 3.
        single = overall
        finite = overall * (8./3. - math.log(etapr))

        evaluation['values'][(0, 0, 0, 0)] = torino_to_madnk_epsexp(
            EpsilonExpansion({-1: single,
                              0: finite}) * self.SEpsilon * (1. / (16 * math.pi**2)), mu2os )

        return evaluation



class QCD_integrated_TRN_C_FgFq(general_current.GeneralCurrent):
    """This CT corresponds to the hard-collinear integrated CT for the qqx splitting
    (Eq. 2.50 in the torino paper, with p=q/qx)
       plus the soft-collinear contribution
     (the first line of eq 2.54 in the torino paper, with CF as prefactor)"""


    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # variables should be updated
    variables = variables_for_integrated_currents_FF

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
        sub.IntegratedCurrent({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(coll_structure_q,)),
        }),
        sub.IntegratedCurrent({
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
    mapping_rules = []


    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided.
        """

        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
            all_steps_info[0]['lower_PS_point'][2])
        mu_r = global_variables['mu_r']
        mu2os = mu_r**2 / shat

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
        etapr = (p_p + p_r).square() / shat

        # hard collinear part
        overall = -self.CF / 2.
        double = 0.
        single = overall
        finite = overall * (2 - math.log(etapr))

        # soft collinear part
        overall = self.CF
        zeta2 = math.pi**2 / 6.
        double += overall
        single += overall * 2
        finite += overall * (6 - 7./2. * zeta2)

        evaluation['values'][(0, 0, 0, 0)] = torino_to_madnk_epsexp(
            EpsilonExpansion({-2: double,
                              -1: single,
                               0: finite}) * self.SEpsilon * (1. / (16 * math.pi**2)), mu2os)

        return evaluation


#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_integrated_TRN_S_g(general_current.GeneralCurrent):
    """ The second line of eq 2.54 in the torino paper"""

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # variables should be updated
    variables = variables_for_integrated_currents_FF

    # Now define the matching singular structures
    soft_structure = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
    )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.IntegratedCurrent({
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
    mapping_rules = []

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a soft type of splitting kernel, which *does* need to know about the reduced process.
        Should be specialised by the daughter class if not dummy
        """


        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']
        higher_PS_point = all_steps_info[0]['higher_PS_point']

        # Normalization factors
        norm = -1./2.

        # All contributions from the soft counterterms are color correlated so we should remove
        # the default value None from the list of color correlations
        evaluation['color_correlations'] = []

        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
            all_steps_info[0]['lower_PS_point'][2])
        mu_r = global_variables['mu_r']
        mu2os = mu_r**2 / shat

        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        prefactor = self.SEpsilon * (1. / (16 * math.pi ** 2))

        for i, a in enumerate(colored_partons):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in colored_partons[i:]:
                # Write the eikonal for that pair
                if a != b:
                    mult_factor = 2.
                else:
                    continue
                    mult_factor = 1.

                pa = lower_PS_point[a]
                pb = lower_PS_point[b]

                etaab = (pa + pb).square() / shat

                single = math.log(etaab)
                finite = math.log(etaab) * (2 - math.log(etaab)/2.)

                kernel = torino_to_madnk_epsexp(EpsilonExpansion({-1: single, 0: finite}) * prefactor * mult_factor,
                                                                  mu2os)

                evaluation['color_correlations'].append(((a, b),))
                evaluation['values'][(0, color_correlation_index, 0, 0)] = kernel

                color_correlation_index += 1

        return evaluation


#=========================================================================================
# NLO final soft-collinear currents.
#=========================================================================================

class QCD_integrated_TRN_CS_FgFq(general_current.GeneralCurrent):
    """SC for qqx.
    We set the SC integrated CTs to zero, since we include the SC
    part in the C"""

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(11,  +1, sub.SubtractionLeg.FINAL), ) )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.IntegratedCurrent({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_q,)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # For the qg soft-collinear the color prefactor is CF
    color_factor = "CF"

    is_zero = True


class QCD_integrated_TRN_CS_FgFg(general_current.GeneralCurrent):
    """SC for gg.
    We set the SC integrated CTs to zero, since we include the SC
    part in the C"""

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL), ) )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.IntegratedCurrent({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_q,)),
        }),
    ]

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # For the gg soft-collinear the color prefactor is CA
    color_factor = "CA"

    is_zero = True

