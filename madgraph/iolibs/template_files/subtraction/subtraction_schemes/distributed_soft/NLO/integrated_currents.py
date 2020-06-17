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
"""Implementation of NLO type of integrated currents for distibuted softs."""

import os
import sys
import math

import madgraph.core.subtraction as subtraction
import madgraph.integrator.phase_space_generators as PS_utils
from madgraph.core.base_objects import EpsilonExpansion
import distributed_soft_config as config
import madgraph.various.misc as misc

import madgraph.various.math_tools.mpl as MPL

import commons.utils as utils
#import factors_and_cuts as factors_and_cuts # nofactors and cuts in distributed softs all in distributed_softs_config

import commons.general_current as general_current
import madgraph.core.subtraction as sub

#import madgraph.various.math_tools.mpl as MPL

from mpmath import polylog

import logging
logger = logging.getLogger("madgraph.integrated_currents")
# Change this to DEBUG to see the debug info
# They are quite verbose so we override the default mode
logger.setLevel(logging.INFO)

log = math.log
pi = math.pi

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# Integrated soft counterterm, which recoils against initial state
# This therefore yields a correlated convolution against initial-state beams
#=========================================================================================

class QCD_integrated_S_Fg(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = False

    is_zero = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = config.divide_by_jacobian_integrated

    # We should not need global variables for this current
    variables = None

    # Now define all the matching singular structures.
    soft_structure_g = sub.SingularStructure(
        substructures=(sub.SoftStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
            )
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        #'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        #'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... soft g(final)
    current_properties['singular_structure'] = soft_structure_g
    template_currents.append(sub.IntegratedCurrent(current_properties)) # (sub.IntegratedBeamCurrent( current_properties ))

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

#=========================================================================================
# Integrated collinear FF counterterm, which recoils against final state,
#=========================================================================================

class QCD_integrated_C_FqFq(general_current.GeneralCurrent):
    """ Integrated collinear F(q) F(q~) current. """

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = config.divide_by_jacobian_integrated

    # We should not need global variables for this current
    variables = None

    # Now define all the matching singular structures.
    coll_structure_qq = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),)
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
    }
    # Now add endpoint IntegratedBeamCurrent... # What is integrated beam current ? from parton evolution
    # ... g(initial) > q(final) qx(final)
    current_properties['singular_structure'] = coll_structure_qq
    template_currents.append(sub.IntegratedCurrent(current_properties)) # (sub.IntegratedBeamCurrent( current_properties ))

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        xi = global_variables['xis'][0]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with ' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        # This is a NLO integrated counterterm
        # it inherits from real countertems so it has
        # an `all_steps_info` with only one step
        # and the lower and higher PS points are the same since
        # there's no mapping applied. the '-1' and 'lower_PS_point' are therefore
        # not meaningful choices at this point.
        PS_point = all_steps_info[-1]['lower_PS_point']

        # Q is handled a little weirdly at the moment, we will therefore reconstruct it from PS_point,
        # which is **exactly what goes into the Matrix Element of the term at hand**.
        # We are using Gabor/colorful notation and therefore Q refers to the total initial momentum
        # of the real matrix element which was regulated by the local version of this counterterm

        n_initial_legs = global_variables['n_initial_legs']
        Q = sum(PS_point.to_list()[:n_initial_legs]) # Sum of initial state momenta :: Use global variable Q??
        Q_square = Q.square()

        color_factor = self.TR

        m_tilde = list() # math.sqrt(PS_point.to_list()[n_initial_legs:].square())
        for a in (PS_point.to_list()[n_initial_legs:]):
            if a.square() < 0.0 :
                m_tilde.append((0.0))
                # print ("mapped mass squared is below zero m_tilde^2 = " + str(a.square()) + " therefore mass is set to zero")
            else:
                m_tilde.append(math.sqrt(a.square()))
            # print(m_tilde)
        s_0 = (math.sqrt(Q_square) - sum(m_tilde))**2 #is s_0 tilde

        factor = s_0 / (mu_r ** 2)

        prefactor = EpsilonExpansion({
            0: 1.0,
            1: - log(factor),
            2: (log(factor)**2)/2.0,
        })
        prefactor *= (-4.0)/(8.0 * pi)

        kernel = EpsilonExpansion({
            -1: 1.0/(12.0 * pi),
            0: (5.0)/(36.0 * pi), # (12.0 * pi),
        })

        evaluation['values'][(0, 0, 0, 0)] = (color_factor * prefactor * kernel).truncate(max_power=0)

        return evaluation

class QCD_integrated_C_FqFg(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = config.divide_by_jacobian_integrated

    # We should not need global variables for this current
    variables = None

    # Now define all the matching singular structures.
    coll_structure_qg = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),)
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
    }

    current_properties['singular_structure'] = coll_structure_qg
    template_currents.append(sub.IntegratedCurrent(current_properties))

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        xi = global_variables['xis'][0]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with ' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        # This is a NLO integrated counterterm
        # it inherits from real countertems so it has
        # an `all_steps_info` with only one step
        # and the lower and higher PS points are the same since
        # there's no mapping applied. the '-1' and 'lower_PS_point' are therefore
        # not meaningful choices at this point.
        PS_point = all_steps_info[-1]['lower_PS_point']

        # Q is handled a little weirdly at the moment, we will therefore reconstruct it from PS_point,
        # which is **exactly what goes into the Matrix Element of the term at hand**.
        # We are using Gabor/colorful notation and therefore Q refers to the total initial momentum
        # of the real matrix element which was regulated by the local version of this counterterm

        n_initial_legs = global_variables['n_initial_legs']
        Q = sum(PS_point.to_list()[:n_initial_legs])
        Q_square = Q.square()

        color_factor = self.CF
        m_tilde = list()
        for a in (PS_point.to_list()[n_initial_legs:]):
            if a.square() < 0.0:
                m_tilde.append((0.0))
                # print ("mapped mass squared is below zero m_tilde^2 = " + str(a.square()) + " therefore mass is set to zero")
            else:
                m_tilde.append(math.sqrt(a.square()))
        s_0 = (math.sqrt(Q_square) - sum(m_tilde)) ** 2  # is s_0 tilde

        factor = s_0 / (mu_r ** 2)

        prefactor = EpsilonExpansion({
            0: 1.0,
            1: - log(factor),
            2: (log(factor) ** 2) / 2.0,
        })
        prefactor *= (-1.0) / (16.0 * pi)

        kernel = EpsilonExpansion({
            -1: 1.0 / (2.0 * pi),
            0: (1.0)/(2.0*pi)
        })

        evaluation['values'][(0, 0, 0, 0)] = (color_factor * prefactor * kernel).truncate(max_power=0)

        return evaluation

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        mu_r = global_variables['mu_r']

        PS_point = all_steps_info[-1]['lower_PS_point']
        n_initial_legs = global_variables['n_initial_legs']
        Q = sum(PS_point.to_list()[:n_initial_legs])

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']

        colinear_parent_leg = global_variables['overall_parents'][0]  # leg number of the parent of the two colinears

        Q_square = Q.square()

        m_tilde = list()
        for a in (PS_point.to_list()[n_initial_legs:]):
            if a.square() < 0.0:
                m_tilde.append((0.0))
                # print ("mapped mass squared is below zero m_tilde^2 = " + str(a.square()) + " therefore mass is set to zero")
            else:
                m_tilde.append(math.sqrt(a.square()))
        s_0 = (math.sqrt(Q_square) - sum(m_tilde)) ** 2  # is s_0 tilde

        factor = 1 / (mu_r ** 2)

        prefactor = EpsilonExpansion({
            0: 1.0,
            1: - log(factor),
            2: (log(factor) ** 2) / 2.0,
        })

        prefactor *= (1.0) / (16. * pi**2)

        norm = -1.0 # -1 from norm in the local currents
        prefactor *= norm

        color_correlation_index = 1

        #Loop over all final state hard partons exept the parent of the colinears

        for i, a in enumerate(colored_partons): # i is unused
            if a == colinear_parent_leg:
                continue
            pa = lower_PS_point[a]
            s_tilde = 2.0 * pa.dot(lower_PS_point[colinear_parent_leg])
            y_0 = s_0 / s_tilde

            int_eik = EpsilonExpansion({
                -2: 1.0,
                -1: 2.0 - log(s_tilde),
                0: (24.0 - pi ** 2 + 3.0 * (-4.0 + log(s_tilde)) * log(s_tilde) + 12.0 * (1.0 + y_0) * log(
                    1.0 + 1.0 / y_0) + 12.0 * polylog(2, -1.0 / y_0)) / 6.0,
            })

            evaluation['color_correlations'].append(((colinear_parent_leg, a),))
            result = (prefactor * int_eik).truncate(max_power=0)
            evaluation['values'][(0, color_correlation_index, 0, 0)] = result
            color_correlation_index += 1

        return evaluation


class QCD_integrated_C_FgFg(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = config.divide_by_jacobian_integrated

    # We should not need global variables for this current
    variables = None

    # Now define all the matching singular structures.
    coll_structure_gg = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),)
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
    }

    current_properties['singular_structure'] = coll_structure_gg
    template_currents.append(sub.IntegratedCurrent(current_properties)) # (sub.IntegratedBeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        xi = global_variables['xis'][0]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with ' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        # This is a NLO integrated counterterm
        # it inherits from real countertems so it has
        # an `all_steps_info` with only one step
        # and the lower and higher PS points are the same since
        # there's no mapping applied. the '-1' and 'lower_PS_point' are therefore
        # not meaningful choices at this point.
        PS_point = all_steps_info[-1]['lower_PS_point']

        # Q is handled a little weirdly at the moment, we will therefore reconstruct it from PS_point,
        # which is **exactly what goes into the Matrix Element of the term at hand**.
        # We are using Gabor/colorful notation and therefore Q refers to the total initial momentum
        # of the real matrix element which was regulated by the local version of this counterterm
        #
        # i.e. :
        #
        # * delta: M(pa_tilde+pb_tilde-> p1_tilde, ...)*delta(alpha) where
        # pa_tilde = (1-alpha) pa, pb_tilde = (1-alpha) pb. i.e. Q = sum(PS_point[initial_legs])
        #
        # * bulk:  M(pa_tilde+pb_tilde-> p1_tilde, ...) for any alpha where
        # pa_tilde = (1-alpha) pa, pb_tilde = (1-alpha) pb. i.e. Q = sum(PS_point[initial_legs])/(1-alpha)
        #
        # * counterterm: M(pa_tilde+pb_tilde-> p1_tilde, ...)|_{alpha=0} where
        # pa_tilde = (1-alpha) pa, pb_tilde = (1-alpha) pb. i.e. Q = sum(PS_point[initial_legs])

        n_initial_legs = global_variables['n_initial_legs']
        Q = sum(PS_point.to_list()[:n_initial_legs])
        Q_square = Q.square()

        color_factor = self.CA
        m_tilde = list()
        for a in (PS_point.to_list()[n_initial_legs:]):
            if a.square() < 0.0:
                m_tilde.append((0.0))
                # print ("mapped mass squared is below zero m_tilde^2 = " + str(a.square()) + " therefore mass is set to zero")
            else:
                m_tilde.append(math.sqrt(a.square()))
        s_0 = (math.sqrt(Q_square) - sum(m_tilde)) ** 2  # should be correct for all masses zero :: is s_0 tilde

        factor = s_0 / (mu_r ** 2)

        prefactor = EpsilonExpansion({
            0: 1.0,
            1: - log(factor),
            2: (log(factor) ** 2) / 2.0,
        })
        prefactor *= (-2.0) / (8.0 * pi)

        kernel = EpsilonExpansion({
            -1: 1.0 / (12.0 * pi),
            0: (5.0) / (36.0 * pi),
        })

        evaluation['values'][(0, 0, 0, 0)] = (color_factor * prefactor * kernel).truncate(max_power=0)

        return evaluation

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        mu_r = global_variables['mu_r']

        PS_point = all_steps_info[-1]['lower_PS_point']
        n_initial_legs = global_variables['n_initial_legs']
        Q = sum(PS_point.to_list()[:n_initial_legs])

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']
        colinear_parent_leg = global_variables['overall_parents'][0]  # leg number of the parent of the two colinears

        Q_square = Q.square()

        m_tilde = list()
        for a in (PS_point.to_list()[n_initial_legs:]):
            if a.square() < 0.0:
                m_tilde.append((0.0))
                # print ("mapped mass squared is below zero m_tilde^2 = " + str(a.square()) + " therefore mass is set to zero")
            else:
                m_tilde.append(math.sqrt(a.square()))
        s_0 = (math.sqrt(Q_square) - sum(m_tilde)) ** 2  # should be correct for all masses zero :: is s_0 tilde

        factor = 1.0 / (mu_r ** 2)

        prefactor = EpsilonExpansion({
            0: 1.0,
            1: - log(factor),
            2: (log(factor) ** 2) / 2.0,
        })

        prefactor *= (1.0) / (16.0 * pi**2)

        prefactor *= 2.0 # Two terms for the gluon gluon going colinear

        norm = -1.0 # -1 from norm in the local currents
        prefactor *= norm

        color_correlation_index = 1

        #Loop over all final state hard partons exept the parent of the colinears

        for i, a in enumerate(colored_partons): # i is unused
            if a == colinear_parent_leg:
                continue
            pa = lower_PS_point[a]
            s_tilde = 2.0 * pa.dot(lower_PS_point[colinear_parent_leg])
            y_0 = s_0/s_tilde
            int_eik = EpsilonExpansion({
                -2 : 1.0,
                -1 : 2.0 - log(s_tilde),
                0 : (24.0 - pi ** 2 + 3.0 * (-4.0 + log(s_tilde)) * log(s_tilde) + 12.0 * (1.0 + y_0) * log(1.0 + 1.0/y_0) + 12.0 * polylog(2, -1.0/y_0))/6.0,
            })

            evaluation['color_correlations'].append(((colinear_parent_leg, a),))
            result = (prefactor * int_eik).truncate(max_power=0)
            evaluation['values'][(0, color_correlation_index, 0, 0)] = result
            color_correlation_index += 1

        return evaluation


class QCD_integrated_CS_FF(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = False

    is_zero = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = config.divide_by_jacobian_integrated

    # We should not need global variables for this current
    variables = None

    # Now define all the matching singular structures.
    soft_structure_g = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
    )
    soft_coll_structure_gg = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=(soft_structure_g,),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),)
        ),)
    )
    soft_coll_structure_qg = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=(soft_structure_g,),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),)
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... g(initial) > g(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.IntegratedCurrent(current_properties)) # (sub.IntegratedBeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.IntegratedCurrent(current_properties)) # (sub.IntegratedBeamCurrent( current_properties ))

    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.IntegratedCurrent(current_properties)) # (sub.BeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.IntegratedCurrent(current_properties)) # (sub.BeamCurrent( current_properties ))

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]











