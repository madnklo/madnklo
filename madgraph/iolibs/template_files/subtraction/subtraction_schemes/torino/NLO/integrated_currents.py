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
import factors_and_cuts as fc


from bidict import bidict

import commons.general_current as general_current

import torino_config
import variables as kernel_variables

import logging
logger = logging.getLogger('madgraph')


pjoin = os.path.join
log = math.log
pi = math.pi

CurrentImplementationError = utils.CurrentImplementationError

logger.info("Damping factors = " + str(fc.damping_factors))


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


#######################################################################################

class QCD_F0(general_current.GeneralCurrent):
    """Implements the NLO QCD PDF counterterm of type F(xi)"""

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    beam_structure_q = sub.SingularStructure(
        substructures=(sub.BeamStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
            )
        ),)
    )
    beam_structure_g = sub.SingularStructure(
        substructures=(sub.BeamStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
            )
        ),),
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : 'proton',
        'beam_PDGs' : torino_config.beam_PDGs_supported
    }
    # Now add endpoint IntegratedBeamCurrent...
    current_properties['singular_structure'] = beam_structure_q
    # ... for an initial-state quark
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))
    # ... for an initial-state gluon
    current_properties['singular_structure'] = beam_structure_g
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... for an initial-state quark
    current_properties['singular_structure'] = beam_structure_q
    template_currents.append(sub.BeamCurrent( current_properties ))
    # ... for an initial-state gluon
    current_properties['singular_structure'] = beam_structure_g
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('PDF F0')

        # Retrieve variables
        beam_number = global_variables['overall_children'][0][0]
        if beam_number not in [1,2]:
            raise CurrentImplementationError(
                "The current %s must always be called for current for which the external" % self.__class__.__name__ +
                " leg number of the beam structure is 1 or 2, not %d."%beam_number)
        # Offset by one so as to have python list access conventions
        beam_number -= 1

        xi = global_variables['xis'][beam_number]
        mu_f = global_variables['mu_fs'][beam_number]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        distribution_type = self.currents_properties[0]['distribution_type']
        NF = self.currents_properties[0]['NF']
        logger.info('NF : ' + str(NF))
        beam_PDGs = self.currents_properties[0]['beam_PDGs']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            xis = list(global_variables['xis'])
            xis[beam_number] = 1.0
            evaluation['Bjorken_rescalings'] = [tuple(xis),]

#        # Only the order epsilon of the scales pre-factor matters here.
#        prefactor = EpsilonExpansion({
#            0: 1.,
#            1: log(mu_r ** 2 / mu_f ** 2)
#        })
#        prefactor *= EpsilonExpansion({-1: 1.}) * self.SEpsilon * (1. / (16 * pi**2) )

        # Only the order epsilon of the scales pre-factor matters here.
        mu2os = mu_r ** 2 / mu_f ** 2
        prefactor = EpsilonExpansion({-1: 1.}) * self.SEpsilon * (1. / (16 * pi**2) )


        # Assign a fake xi for now if the distribution type is 'endpoint'
        # This is not optimal, eventually we should put each of these three pieces in
        # separate currents, this will be optimal however in the low-level implementation
        if distribution_type == 'endpoint':
            xi = 0.5

        logger.info('distribution_type : ' + str(distribution_type))

        # Define the NLO QCD PDF counterterms kernels
        kernel_gg = {
            'bulk': prefactor * (
                    2. * self.CA * (1. / (1. - xi) + (1. - xi) / xi - 1. + xi * (1 - xi))
            ),
            'counterterm': prefactor * (2. * self.CA / (1. - xi)),
            'endpoint': prefactor * (11. / 6. * self.CA - 2. / 3. * NF * self.TR)
        }

        kernel_gq = {
            'bulk': prefactor * (self.CF * (1. + (1. - xi) ** 2) / xi),
            'counterterm': None,
            'endpoint': None
        }

        kernel_qg = {
            'bulk': prefactor * (self.TR * (xi ** 2 + (1. - xi) ** 2)),
            'counterterm': None,
            'endpoint': None
        }

        kernel_qq = {
            'bulk': prefactor * (self.CF * ((1. + xi ** 2) / (1. - xi))),
            'counterterm': prefactor * (self.CF * ((1. + xi ** 2) / (1. - xi))),
            'endpoint': None
        }

        if kernel_qq[distribution_type] is not None:
            kernel_qq[distribution_type] = torino_to_madnk_epsexp(kernel_qq[distribution_type], mu2os)
        elif kernel_qg[distribution_type] is not None:
            kernel_qg[distribution_type] = torino_to_madnk_epsexp(kernel_qg[distribution_type], mu2os)
        elif kernel_gq[distribution_type] is not None:
            kernel_gq[distribution_type] = torino_to_madnk_epsexp(kernel_gq[distribution_type], mu2os)
        elif kernel_gg[distribution_type] is not None:
            kernel_gg[distribution_type] = torino_to_madnk_epsexp(kernel_gg[distribution_type], mu2os)

        logger.info('kernel_gg[distribution_type] :' + str(kernel_gg[distribution_type]))
        logger.info('kernel_gq[distribution_type] :' + str(kernel_gq[distribution_type]))
        logger.info('kernel_qg[distribution_type] :' + str(kernel_qg[distribution_type]))
        logger.info('kernel_qq[distribution_type] :' + str(kernel_qq[distribution_type]))

        active_quark_PDGs = self.currents_properties[0]['active_fermions']

        # Build the NLO flavor matrix
        flavor_matrix = {}
        for reduced_flavor in beam_PDGs:
            # Gluon backward evolution
            if reduced_flavor == 21:
                gluon_dict = {}
                if kernel_gg[distribution_type] is not None:
                    gluon_dict[(21,)] = kernel_gg[distribution_type]
                if active_quark_PDGs and kernel_gq[distribution_type] is not None:
                    gluon_dict[active_quark_PDGs] = kernel_gq[distribution_type]
                if gluon_dict:
                    flavor_matrix[21] = gluon_dict

            # Quark backward evolution
            if reduced_flavor in active_quark_PDGs:
                quark_dict = {}
                if kernel_qg[distribution_type] is not None:
                    quark_dict[(21,)] = kernel_qg[distribution_type]
                if kernel_qq[distribution_type] is not None:
                    quark_dict[(reduced_flavor,)] = kernel_qq[distribution_type]
                if quark_dict:
                    flavor_matrix[reduced_flavor] = quark_dict


        # Truncate all entries of the flavor matrix so as to remove irrelevant O(\eps) terms
        for flav_in, flav_outs in flavor_matrix.items():
            for flav_out, eps_expansion in flav_outs.items():
                eps_expansion.truncate(max_power=0)

        # Now assign the flavor matrix to the result
        evaluation['values'][(0,0,0,0)] = flavor_matrix

        # Promote the format of the flavor matrix to be that of a generalised two-sided convolution
        evaluation.promote_to_two_sided_convolution(beam_number+1)

        return evaluation

class QCD_F0_lepton(general_current.GeneralCurrent):
    """Implements the dummy PDF counterterm of type F(xi) for when doing lepton collision with colorful pp."""

    # Enable the flag below to debug this current
    DEBUG = True

    # Define which leptons we want to allow
#    lepton_abs_PDGs = [11,12,13]

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    beam_structure_lepton = sub.SingularStructure(
        substructures=(sub.BeamStructure(
            substructures=tuple([]),
            legs=(
                # Set the leg PDG to +11 which belongs to the lepton equivalency set defined below
                # in the function build_equivalency_sets.
                sub.SubtractionLeg(10, 11, sub.SubtractionLeg.INITIAL),
            )
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : 'lepton',
        'beam_PDGs' : torino_config.lepton_abs_PDGs
#            [ tuple(sorted([pdg, -pdg])) for pdg in lepton_abs_PDGs ] +
#            [ tuple([pdg, ]) for pdg in lepton_abs_PDGs ] +
#            [ tuple([-pdg, ]) for pdg in lepton_abs_PDGs ]
    }

    # Now add endpoint IntegratedBeamCurrent...
    current_properties['singular_structure'] = beam_structure_lepton
    # ... for an initial-state lepton
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... for an initial-state lepton
    current_properties['singular_structure'] = beam_structure_lepton
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    # Remove this contribution at the generation level already by setting its "is_zero" flag to True.
    is_zero = False

#    @classmethod
#    def build_equivalency_sets(cls, model):
#        """ Force the PDG +1 of the template """
#        return [(set([pdg for pdg in cls.lepton_abs_PDGs])|set([-pdg for pdg in cls.lepton_abs_PDGs])),]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('PDF F0_lepton')


        # Set this contribution to zero as we are considering here QCD corrections
        evaluation['values'] = self.subtraction_current_evaluation_class.zero()['values']

        return evaluation


#=========================================================================================
# Integrated IF collinear counterterm
#=========================================================================================

class QCD_integrated_TRN_C_IqFg(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

    # Now define all the matching singular structures. We can combine in this class *all* integrated
    # initial-final counterterm thanks to the fact that *allowed_backward_evolved_flavors* is provided at
    # runtime. But it would be cleaner (and possible already with this structure) to implement each
    # kernel (qg, qq, etc...) in a separate class with its matching singular structure every time.
    coll_structure_qg = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),)
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_qg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_qg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('C_IqFg')

        # Retrieve variables
        beam_number = global_variables['overall_children'][0][0]
        if beam_number not in [1,2]:
            raise CurrentImplementationError(
                "The current %s must always be called for current for which the external" % self.__class__.__name__ +
                " leg number of the beam structure is 1 or 2, not %d."%beam_number)
        # Offset by one so as to have python list access conventions
        beam_number -= 1

        xi = global_variables['xis'][beam_number]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors[(beam_number+1)%2] != 'ALL':
            raise CurrentImplementationError('The current %s must always be called with' % self.__class__.__name__ +
                "allowed_backward_evolved_flavors[%d (i.e. inactive beam number)] = 'ALL' not %s" % ((beam_number+1)%2,str(
                allowed_backward_evolved_flavors[beam_number])))

        distribution_type = self.currents_properties[0]['distribution_type']
        beam_PDGs = self.currents_properties[0]['beam_PDGs']
        recoiler = global_variables['recoiler']

        logger.info('distribution_type : ' + str(distribution_type))

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            xis = list(global_variables['xis'])
            xis[beam_number] = 1.0
            evaluation['Bjorken_rescalings'] = [tuple(xis),]

        Q = global_variables['Q']

        # Obtain Q_square.
        Q_square = Q.square()

#        # Only up to the order epsilon^2 of the scales prefactor matters here.
#        logMuQ = log(mu_r ** 2 / Q_square)
#        prefactor = EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})
#        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )


#        logger.info('all_step_info :' + str(all_steps_info))

#        logger.info('all_step_info[] :' + str(2*all_steps_info[0]['lower_PS_point'][1]*all_steps_info[0]['lower_PS_point'][6]))
#        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
#            all_steps_info[0]['lower_PS_point'][2])
#        mu_r = global_variables['mu_r']
#        mu2os = mu_r**2 / Q_square


        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
#        etapr = (p_p + p_r).square() / Q
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))
        mu2os = mu_r**2 / (p_p + p_r).square()

        # Only up to the order epsilon^2 of the scales prefactor matters here.
#        logMuQ = log(mu_r ** 2 / (p_p + p_r).square())
#        prefactor = EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})
#        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )


        # The additional 1/x part of the prefactor is included later during the PDF
        # convolution of the event (using its 'Bjorken rescaling' attribute) because
        # we must make sure that the plus distribution hits on it.
        # Also, the same 1/x appears in the PDF counterterms as a result of the change
        # of variable necessary to bring them in the form where the plus distribution
        # only acts on the PDF. So it makes sense to keep it completely factorised.

        # Assign a fake x for now if the distribution type is 'endpoint'
        # TODO: this is not optimal, eventually we should put each of these three pieces in
        # separate currents
        if distribution_type == 'endpoint':
            x = 0.5
        else:
            x = xi

        # In MadNkLO, we use the change of variable xb' = xb*xi so that the factor
        # (Q^2)^\eps in Eq. 5.21 of https://arxiv.org/pdf/0903.1218.pdf actually reads
        # (Q^2/(xi1*xi2))^\eps and the '+' distributions also act on it, which we realize
        # by simply multiplying the Q^2 provided by the xi factor that must be set to one.
#        logMuQ_plus = log(mu_r ** 2 / (Q_square * x))
#        prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
#        prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )

        logger.info('recoiler :' + str(recoiler))

        color_factor = self.CF
        overall = 1.
        if recoiler > 2:
            kernel_qq = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: - (1 + x**2) / (1. - x) ,
                    0: (1. - x) * (log(1. - x) - fc.A1(fc.beta_IF))
                }) ,
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -1: - (1 + x**2) / (1. - x) ,
                    0: 0.
                }) ,
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -1: - 1./ 2. ,
                    0: 0.
                })
            }

        elif recoiler <=2:
            kernel_qq = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: - (1 + x**2) / (1. - x) ,
                    0: (1. - x) * (2. * log(1. - x) - fc.A1(fc.beta_II))
                }) ,
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -1: - (1 + x**2) / (1. - x) ,
                    0: 0.
                }) ,
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -1: - 1./ 2. ,
                    0: 0.
                })
            }


        kernel_gg = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
        kernel_gq = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
        kernel_qg = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
 
        if kernel_qq[distribution_type] is not None and distribution_type == 'counterterm':
            kernel_qq[distribution_type] = torino_to_madnk_epsexp(kernel_qq[distribution_type] * prefactor, mu2os / x)
        elif kernel_qq[distribution_type] is not None:
            kernel_qq[distribution_type] = torino_to_madnk_epsexp(kernel_qq[distribution_type] * prefactor, mu2os)


        active_quark_PDGs = self.currents_properties[0]['active_fermions']

        # Build the NLO flavor matrix
        flavor_matrix = {}
        for reduced_flavor in beam_PDGs:
            # Gluon backward evolution
            if reduced_flavor == 21:
                gluon_dict = {}
                if kernel_gg[distribution_type] is not None:
                    gluon_dict[(21,)] = kernel_gg[distribution_type]
                if active_quark_PDGs and kernel_gq[distribution_type] is not None:
                    gluon_dict[active_quark_PDGs] = kernel_gq[distribution_type]
                if gluon_dict:
                    flavor_matrix[21] = gluon_dict

            # Quark backward evolution
            if reduced_flavor in active_quark_PDGs:
                quark_dict = {}
                if kernel_qg[distribution_type] is not None:
                    quark_dict[(21,)] = kernel_qg[distribution_type]
                if kernel_qq[distribution_type] is not None:
                    quark_dict[(reduced_flavor,)] = kernel_qq[distribution_type]
                if quark_dict:
                    flavor_matrix[reduced_flavor] = quark_dict


        # Truncate all entries of the flavor matrix so as to remove irrelevant O(\eps) terms
        for flav_in, flav_outs in flavor_matrix.items():
            for flav_out, eps_expansion in flav_outs.items():
                eps_expansion.truncate(max_power=0)

        # Now apply the mask 'allowed_backward_evolved_flavors' if not set to 'ALL'
        filtered_flavor_matrix = self.apply_flavor_mask(flavor_matrix, allowed_backward_evolved_flavors[beam_number])

        # Now assign the flavor matrix to the result
        evaluation['values'][(0,0,0,0)] = filtered_flavor_matrix

        # Promote the format of the flavor matrix to be that of a generalised two-sided convolution
        evaluation.promote_to_two_sided_convolution(beam_number+1)

        return evaluation


class QCD_integrated_TRN_C_IgFq(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

    # Now define all the matching singular structures. We can combine in this class *all* integrated
    # initial-final counterterm thanks to the fact that *allowed_backward_evolved_flavors* is provided at
    # runtime. But it would be cleaner (and possible already with this structure) to implement each
    # kernel (qg, qq, etc...) in a separate class with its matching singular structure every time.

    coll_structure_gq = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
                sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),)
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > g(initial_after_emission) q(final)
    current_properties['singular_structure'] = coll_structure_gq
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > g(initial_after_emission) q(final)
    current_properties['singular_structure'] = coll_structure_gq
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('C_IgFq')

        # Retrieve variables
        beam_number = global_variables['overall_children'][0][0]
        if beam_number not in [1,2]:
            raise CurrentImplementationError(
                "The current %s must always be called for current for which the external" % self.__class__.__name__ +
                " leg number of the beam structure is 1 or 2, not %d."%beam_number)
        # Offset by one so as to have python list access conventions
        beam_number -= 1

        xi = global_variables['xis'][beam_number]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors[(beam_number+1)%2] != 'ALL':
            raise CurrentImplementationError('The current %s must always be called with' % self.__class__.__name__ +
                "allowed_backward_evolved_flavors[%d (i.e. inactive beam number)] = 'ALL' not %s" % ((beam_number+1)%2,str(
                allowed_backward_evolved_flavors[beam_number])))

        distribution_type = self.currents_properties[0]['distribution_type']
        beam_PDGs = self.currents_properties[0]['beam_PDGs']
        recoiler = global_variables['recoiler']


        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            xis = list(global_variables['xis'])
            xis[beam_number] = 1.0
            evaluation['Bjorken_rescalings'] = [tuple(xis),]

        Q = global_variables['Q']

        # Obtain Q_square.
        Q_square = Q.square()

#        # Only up to the order epsilon^2 of the scales prefactor matters here.
#        logMuQ = log(mu_r ** 2 / Q_square)
#        prefactor = EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})
#        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
#        etapr = (p_p + p_r).square() / Q
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))
        mu2os = mu_r**2 / (p_p + p_r).square()

        logger.info('distribution_type :' + str(distribution_type))

        # The additional 1/x part of the prefactor is included later during the PDF
        # convolution of the event (using its 'Bjorken rescaling' attribute) because
        # we must make sure that the plus distribution hits on it.
        # Also, the same 1/x appears in the PDF counterterms as a result of the change
        # of variable necessary to bring them in the form where the plus distribution
        # only acts on the PDF. So it makes sense to keep it completely factorised.


        # Assign a fake x for now if the distribution type is 'endpoint'
        # TODO: this is not optimal, eventually we should put each of these three pieces in
        # separate currents
        if distribution_type == 'endpoint':
            x = 0.5
        else:
            x = xi

        # In MadNkLO, we use the change of variable xb' = xb*xi so that the factor
        # (Q^2)^\eps in Eq. 5.21 of https://arxiv.org/pdf/0903.1218.pdf actually reads
        # (Q^2/(xi1*xi2))^\eps and the '+' distributions also act on it, which we realize
        # by simply multiplying the Q^2 provided by the xi factor that must be set to one.
#        logMuQ_plus = log(mu_r ** 2 / (Q_square * x))
#        prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
#        prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )


        color_factor = self.TR
        overall = 1. 
        if recoiler > 2:
            kernel_qg = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: - (x**2 + (1. - x)**2),
                    0: (x**2 + (1. - x)**2) *( log(1. - x) - fc.A1(fc.beta_IF) )
                }) ,
                'counterterm': None ,
                'endpoint': None
            }

        elif recoiler <=2:
            kernel_qg = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: - (x**2 + (1. - x)**2),
                    0: (x**2 + (1. - x)**2) * ( 2. * log(1. - x) - fc.A1(fc.beta_II) )
                }),
                'counterterm': None ,
                'endpoint': None
            }

        kernel_gg = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
        kernel_gq = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
        kernel_qq = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
 
        if kernel_qg[distribution_type] is not None and distribution_type == 'counterterm':
            kernel_qg[distribution_type] = torino_to_madnk_epsexp(kernel_qg[distribution_type] * prefactor, mu2os / x)
        elif kernel_qg[distribution_type] is not None:
            kernel_qg[distribution_type] = torino_to_madnk_epsexp(kernel_qg[distribution_type] * prefactor, mu2os)


        active_quark_PDGs = self.currents_properties[0]['active_fermions']

        # Build the NLO flavor matrix
        flavor_matrix = {}
        for reduced_flavor in beam_PDGs:
            # Gluon backward evolution
            if reduced_flavor == 21:
                gluon_dict = {}
                if kernel_gg[distribution_type] is not None:
                    gluon_dict[(21,)] = kernel_gg[distribution_type]
                if active_quark_PDGs and kernel_gq[distribution_type] is not None:
                    gluon_dict[active_quark_PDGs] = kernel_gq[distribution_type]
                if gluon_dict:
                    flavor_matrix[21] = gluon_dict

            # Quark backward evolution
            if reduced_flavor in active_quark_PDGs:
                quark_dict = {}
                if kernel_qg[distribution_type] is not None:
                    quark_dict[(21,)] = kernel_qg[distribution_type]
                if kernel_qq[distribution_type] is not None:
                    quark_dict[(reduced_flavor,)] = kernel_qq[distribution_type]
                if quark_dict:
                    flavor_matrix[reduced_flavor] = quark_dict

        # Truncate all entries of the flavor matrix so as to remove irrelevant O(\eps) terms
        for flav_in, flav_outs in flavor_matrix.items():
            for flav_out, eps_expansion in flav_outs.items():
                eps_expansion.truncate(max_power=0)

        # Now apply the mask 'allowed_backward_evolved_flavors' if not set to 'ALL'
        filtered_flavor_matrix = self.apply_flavor_mask(flavor_matrix, allowed_backward_evolved_flavors[beam_number])

        # Now assign the flavor matrix to the result
        evaluation['values'][(0,0,0,0)] = filtered_flavor_matrix

        # Promote the format of the flavor matrix to be that of a generalised two-sided convolution
        evaluation.promote_to_two_sided_convolution(beam_number+1)

        return evaluation

class QCD_integrated_TRN_C_IqFq(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

    # Now define all the matching singular structures. We can combine in this class *all* integrated
    # initial-final counterterm thanks to the fact that *allowed_backward_evolved_flavors* is provided at
    # runtime. But it would be cleaner (and possible already with this structure) to implement each
    # kernel (qg, qq, etc...) in a separate class with its matching singular structure every time.
    coll_structure_qq = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),)
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... g(initial) > q(initial_after_emission) qx(final)
    current_properties['singular_structure'] = coll_structure_qq
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... g(initial) > g(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_qq
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('C_IqFq')

        # Retrieve variables
        beam_number = global_variables['overall_children'][0][0]
        if beam_number not in [1,2]:
            raise CurrentImplementationError(
                "The current %s must always be called for current for which the external" % self.__class__.__name__ +
                " leg number of the beam structure is 1 or 2, not %d."%beam_number)
        # Offset by one so as to have python list access conventions
        beam_number -= 1

        xi = global_variables['xis'][beam_number]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors[(beam_number+1)%2] != 'ALL':
            raise CurrentImplementationError('The current %s must always be called with' % self.__class__.__name__ +
                "allowed_backward_evolved_flavors[%d (i.e. inactive beam number)] = 'ALL' not %s" % ((beam_number+1)%2,str(
                allowed_backward_evolved_flavors[beam_number])))

        distribution_type = self.currents_properties[0]['distribution_type']
        beam_PDGs = self.currents_properties[0]['beam_PDGs']
        recoiler = global_variables['recoiler']

#        logger.info('distribution : ' + str(distribution_type))
#        logger.info('beam_PDGs : ' + str(beam_PDGs))
#        logger.info('global_variables : ' + str(global_variables))
#        logger.info('recoiler : ' + str(recoiler))


        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            xis = list(global_variables['xis'])
            xis[beam_number] = 1.0
            evaluation['Bjorken_rescalings'] = [tuple(xis),]

        Q = global_variables['Q']

        # Obtain Q_square.
        Q_square = Q.square()

        # Only up to the order epsilon^2 of the scales prefactor matters here.
#        logMuQ = log(mu_r ** 2 / Q_square)
#        prefactor = EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})
#        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )

        # The additional 1/x part of the prefactor is included later during the PDF
        # convolution of the event (using its 'Bjorken rescaling' attribute) because
        # we must make sure that the plus distribution hits on it.
        # Also, the same 1/x appears in the PDF counterterms as a result of the change
        # of variable necessary to bring them in the form where the plus distribution
        # only acts on the PDF. So it makes sense to keep it completely factorised.

        # Assign a fake x for now if the distribution type is 'endpoint'
        # TODO: this is not optimal, eventually we should put each of these three pieces in
        # separate currents
        if distribution_type == 'endpoint':
            x = 0.5
        else:
            x = xi

#        # In MadNkLO, we use the change of variable xb' = xb*xi so that the factor
#        # (Q^2)^\eps in Eq. 5.21 of https://arxiv.org/pdf/0903.1218.pdf actually reads
#        # (Q^2/(xi1*xi2))^\eps and the '+' distributions also act on it, which we realize
#        # by simply multiplying the Q^2 provided by the xi factor that must be set to one.
#        logMuQ_plus = log(mu_r ** 2 / (Q_square * x))
#        prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
#        prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )


        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
#        etapr = (p_p + p_r).square() / Q
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))
        mu2os = mu_r**2 / (p_p + p_r).square()

        logger.info('parent :' + str(global_variables['overall_children'][0]))
        logger.info('recoiler :' + str(recoiler))
#        logger.info('distribution_type :' + str(distribution_type))
#        logger.info('p_r * p_p : :' + str((p_p + p_r).square()))


#da aggiornare
        color_factor = self.CF
        overall = 1.
        if recoiler > 2:
            kernel_gq = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: - ( 1. + ( 1. - x)**2 ) / x ,
                    0: (( 1. + ( 1. - x)**2 ) / x) * (log(1. - x)  - fc.A1(fc.beta_IF) )
                }) ,
                'counterterm': None ,
                'endpoint': None
            }

        elif recoiler <= 2:
            kernel_gq = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: - (( 1. + ( 1. - x)**2 ) / x),
                    0: (( 1. + ( 1. - x)**2 ) / x) * ( 2. * log(1. - x) - fc.A1(fc.beta_II) )
                }) ,
                'counterterm': None ,
                'endpoint': None
            }


        kernel_gg = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
        kernel_qg = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
        kernel_qq = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
 
        if kernel_gq[distribution_type] is not None and distribution_type == 'counterterm':
            kernel_gq[distribution_type] = torino_to_madnk_epsexp(kernel_gq[distribution_type] * prefactor, mu2os / x)
        elif kernel_gq[distribution_type] is not None:
            kernel_gq[distribution_type] = torino_to_madnk_epsexp(kernel_gq[distribution_type] * prefactor, mu2os)

        logger.info('kernel_gq[distribution_type] :' + str(kernel_gq[distribution_type]))

        active_quark_PDGs = self.currents_properties[0]['active_fermions']

        # Build the NLO flavor matrix
        flavor_matrix = {}
        for reduced_flavor in beam_PDGs:
            # Gluon backward evolution
            if reduced_flavor == 21:
                gluon_dict = {}
                if kernel_gg[distribution_type] is not None:
                    gluon_dict[(21,)] = kernel_gg[distribution_type]
                if active_quark_PDGs and kernel_gq[distribution_type] is not None:
                    gluon_dict[active_quark_PDGs] = kernel_gq[distribution_type]
                if gluon_dict:
                    flavor_matrix[21] = gluon_dict

            # Quark backward evolution
            if reduced_flavor in active_quark_PDGs:
                quark_dict = {}
                if kernel_qg[distribution_type] is not None:
                    quark_dict[(21,)] = kernel_qg[distribution_type]
                if kernel_qq[distribution_type] is not None:
                    quark_dict[(reduced_flavor,)] = kernel_qq[distribution_type]
                if quark_dict:
                    flavor_matrix[reduced_flavor] = quark_dict

        # Truncate all entries of the flavor matrix so as to remove irrelevant O(\eps) terms
        for flav_in, flav_outs in flavor_matrix.items():
            for flav_out, eps_expansion in flav_outs.items():
                eps_expansion.truncate(max_power=0)

        # Now apply the mask 'allowed_backward_evolved_flavors' if not set to 'ALL'
        filtered_flavor_matrix = self.apply_flavor_mask(flavor_matrix, allowed_backward_evolved_flavors[beam_number])

        # Now assign the flavor matrix to the result
        evaluation['values'][(0,0,0,0)] = filtered_flavor_matrix

        # Promote the format of the flavor matrix to be that of a generalised two-sided convolution
        evaluation.promote_to_two_sided_convolution(beam_number+1)

        return evaluation

class QCD_integrated_TRN_C_IgFg(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

    # Now define all the matching singular structures. We can combine in this class *all* integrated
    # initial-final counterterm thanks to the fact that *allowed_backward_evolved_flavors* is provided at
    # runtime. But it would be cleaner (and possible already with this structure) to implement each
    # kernel (qg, qq, etc...) in a separate class with its matching singular structure every time.
    coll_structure_gg = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
                sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),)
        ),)
    )


    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... g(initial) > g(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_gg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_gg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('C_IgFg')

        # Retrieve variables
        beam_number = global_variables['overall_children'][0][0]
        if beam_number not in [1,2]:
            raise CurrentImplementationError(
                "The current %s must always be called for current for which the external" % self.__class__.__name__ +
                " leg number of the beam structure is 1 or 2, not %d."%beam_number)
        # Offset by one so as to have python list access conventions
        beam_number -= 1

        xi = global_variables['xis'][beam_number]
        mu_r = global_variables['mu_r']
        NF = self.currents_properties[0]['NF']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors[(beam_number+1)%2] != 'ALL':
            raise CurrentImplementationError('The current %s must always be called with' % self.__class__.__name__ +
                "allowed_backward_evolved_flavors[%d (i.e. inactive beam number)] = 'ALL' not %s" % ((beam_number+1)%2,str(
                allowed_backward_evolved_flavors[beam_number])))

        distribution_type = self.currents_properties[0]['distribution_type']
        beam_PDGs = self.currents_properties[0]['beam_PDGs']
        recoiler = global_variables['recoiler']

#        logger.info('distribution : ' + str(distribution_type))
#        logger.info('beam_PDGs : ' + str(beam_PDGs))
#        logger.info('global_variables : ' + str(global_variables))
#        logger.info('recoiler : ' + str(recoiler))


        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            xis = list(global_variables['xis'])
            xis[beam_number] = 1.0
            evaluation['Bjorken_rescalings'] = [tuple(xis),]

        Q = global_variables['Q']

        # Obtain Q_square.
        Q_square = Q.square()

#        # Only up to the order epsilon^2 of the scales prefactor matters here.
#        logMuQ = log(mu_r ** 2 / Q_square)
#        prefactor = EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})
#        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
        mu2os = mu_r**2 / (p_p + p_r).square()
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))
        logger.info('p_p * p_r : ' + str((p_p + p_r).square()))


        # The additional 1/x part of the prefactor is included later during the PDF
        # convolution of the event (using its 'Bjorken rescaling' attribute) because
        # we must make sure that the plus distribution hits on it.
        # Also, the same 1/x appears in the PDF counterterms as a result of the change
        # of variable necessary to bring them in the form where the plus distribution
        # only acts on the PDF. So it makes sense to keep it completely factorised.

        # Input variables
        #y_0 = factors_and_cuts.y_0_prime
        #logy0 = log(y_0)
        # Assign a fake x for now if the distribution type is 'endpoint'
        # TODO: this is not optimal, eventually we should put each of these three pieces in
        # separate currents
        if distribution_type == 'endpoint':
            x = 0.5
        else:
            x = xi

#        # In MadNkLO, we use the change of variable xb' = xb*xi so that the factor
#        # (Q^2)^\eps in Eq. 5.21 of https://arxiv.org/pdf/0903.1218.pdf actually reads
#        # (Q^2/(xi1*xi2))^\eps and the '+' distributions also act on it, which we realize
#        # by simply multiplying the Q^2 provided by the xi factor that must be set to one.
#        logMuQ_plus = log(mu_r ** 2 / (Q_square * x))
#        prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
#        prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )


        # Define the NLO QCD integrate initial-state single collinear counterterms kernels
        logger.info('parent :' + str(global_variables['overall_children'][0]))
        logger.info('recoiler :' + str(recoiler))
#        logger.info('distribution_type :' + str(distribution_type))


        color_factor = self.CA
        overall = 2.
        if recoiler > 2:
            kernel_gg = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: - ( (1. / (1. - x) ) + (1. - x)/x - 1. + x*(1. - x) ),
                    0: ((1. / (1. -x) ) + (1. -x) * x) * (log(1. - x) - fc.A1(fc.beta_IF))
                }) ,
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -1: - (1. / (1. - x) ),
                    0: 0.
                }) ,
                'endpoint': EpsilonExpansion({
#                    -1: - (2. * self.CA + NF * self.TR * 4. / 6.) ,
                    -1: - (2. * self.CA ) ,
                    0: 0.
                })
            }

        elif recoiler <= 2:
            kernel_gg = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: - ( (1. / (1. -x) ) + (1. -x)/x - 1. + x*(1. - x) ),
                    0: ((1. / (1. -x) ) + (1. -x) * x) * (2. * log(1. - x) - fc.A1(fc.beta_II))
                }) ,
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -1: - (1. / (1. -x) ),
                    0: 0.
                }) ,
                'endpoint': EpsilonExpansion({
#                    -1: - (2. * self.CA + NF * self.TR * 4. / 6.) ,
                    -1: - (2. * self.CA ) ,
                    0: 0.
                })
            }


        kernel_qg = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
        kernel_gq = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
        kernel_qq = {
            'bulk': None ,
            'counterterm': None ,
            'endpoint': None
        }
 
        if kernel_gg[distribution_type] is not None and distribution_type == 'counterterm':
            kernel_gg[distribution_type] = torino_to_madnk_epsexp(kernel_gg[distribution_type] * prefactor, mu2os / x)
        elif kernel_gg[distribution_type] is not None:
            kernel_gg[distribution_type] = torino_to_madnk_epsexp(kernel_gg[distribution_type] * prefactor, mu2os)

        logger.info('kernel_gg[distribution_type] :' + str(kernel_gg[distribution_type]))

        active_quark_PDGs = self.currents_properties[0]['active_fermions']

        # Build the NLO flavor matrix
        flavor_matrix = {}
        for reduced_flavor in beam_PDGs:
            # Gluon backward evolution
            if reduced_flavor == 21:
                gluon_dict = {}
                if kernel_gg[distribution_type] is not None:
                    gluon_dict[(21,)] = kernel_gg[distribution_type]
                if active_quark_PDGs and kernel_gq[distribution_type] is not None:
                    gluon_dict[active_quark_PDGs] = kernel_gq[distribution_type]
                if gluon_dict:
                    flavor_matrix[21] = gluon_dict

            # Quark backward evolution
            if reduced_flavor in active_quark_PDGs:
                quark_dict = {}
                if kernel_qg[distribution_type] is not None:
                    quark_dict[(21,)] = kernel_qg[distribution_type]
                if kernel_qq[distribution_type] is not None:
                    quark_dict[(reduced_flavor,)] = kernel_qq[distribution_type]
                if quark_dict:
                    flavor_matrix[reduced_flavor] = quark_dict

        # Truncate all entries of the flavor matrix so as to remove irrelevant O(\eps) terms
        for flav_in, flav_outs in flavor_matrix.items():
            for flav_out, eps_expansion in flav_outs.items():
                eps_expansion.truncate(max_power=0)

        # Now apply the mask 'allowed_backward_evolved_flavors' if not set to 'ALL'
        filtered_flavor_matrix = self.apply_flavor_mask(flavor_matrix, allowed_backward_evolved_flavors[beam_number])

        # Now assign the flavor matrix to the result
        evaluation['values'][(0,0,0,0)] = filtered_flavor_matrix

        # Promote the format of the flavor matrix to be that of a generalised two-sided convolution
        evaluation.promote_to_two_sided_convolution(beam_number+1)

        return evaluation

#=========================================================================================
# Integrated IF soft-collinear counterterm
#=========================================================================================

class QCD_integrated_TRN_CS_IqFg(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

    # Now define all the matching singular structures. We can combine in this class *all* integrated
    # initial-final counterterm thanks to the fact that *allowed_backward_evolved_flavors* is provided at
    # runtime. But it would be cleaner (and possible already with this structure) to implement each
    # kernel (qg, qq, etc...) in a separate class with its matching singular structure every time.
    soft_structure_g = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
    )
    soft_coll_structure_qg = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=(soft_structure_g,),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),)
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('CS_IqFg')

        # Retrieve variables
        beam_number = global_variables['overall_children'][0][0]
        if beam_number not in [1,2]:
            raise CurrentImplementationError(
                "The current %s must always be called for current for which the external" % self.__class__.__name__ +
                " leg number of the beam structure is 1 or 2, not %d."%beam_number)
        # Offset by one so as to have python list access conventions
        beam_number -= 1

        xi = global_variables['xis'][beam_number]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors[(beam_number+1)%2] != 'ALL':
            raise CurrentImplementationError('The current %s must always be called with' % self.__class__.__name__ +
                "allowed_backward_evolved_flavors[%d (i.e. inactive beam number)] = 'ALL' not %s" % ((beam_number+1)%2,str(
                allowed_backward_evolved_flavors[beam_number])))

        distribution_type = self.currents_properties[0]['distribution_type']
        beam_PDGs = self.currents_properties[0]['beam_PDGs']
        recoiler = global_variables['recoiler']


        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
#        if distribution_type == 'counterterm':
#            xis = list(global_variables['xis'])
#            xis[beam_number] = 1.0
#            evaluation['Bjorken_rescalings'] = [tuple(xis),]
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]


        Q = global_variables['Q']

        # Obtain Q_square.
        Q_square = Q.square()

#        # Only up to the order epsilon^2 of the scales prefactor matters here.
#        logMuQ = log(mu_r ** 2 / Q_square)
#        prefactor = EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})
#        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))
        mu2os = mu_r**2 / (p_p + p_r).square()

        # The additional 1/x part of the prefactor is included later during the PDF
        # convolution of the event (using its 'Bjorken rescaling' attribute) because
        # we must make sure that the plus distribution hits on it.
        # Also, the same 1/x appears in the PDF counterterms as a result of the change
        # of variable necessary to bring them in the form where the plus distribution
        # only acts on the PDF. So it makes sense to keep it completely factorised.

        # Input variables
        #y_0 = factors_and_cuts.y_0_prime
        #logy0 = log(y_0)
        # Assign a fake x for now if the distribution type is 'endpoint'
        # TODO: this is not optimal, eventually we should put each of these three pieces in
        # separate currents
        if distribution_type == 'endpoint':
            x = 0.5
        else:
            x = xi

        # In MadNkLO, we use the change of variable xb' = xb*xi so that the factor
        # (Q^2)^\eps in Eq. 5.21 of https://arxiv.org/pdf/0903.1218.pdf actually reads
        # (Q^2/(xi1*xi2))^\eps and the '+' distributions also act on it, which we realize
        # by simply multiplying the Q^2 provided by the xi factor that must be set to one.
#        logMuQ_plus = log(mu_r ** 2 / (Q_square * x))
#        prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
#        prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )


        # Define the NLO QCD integrate initial-state single collinear counterterms kernels


        color_factor = self.CF
        overall = - 2.
        if recoiler > 2:
            kernel = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: x**(1. + fc.alpha) / (1. - x)  ,
                    0: - fc.A1(fc.beta_IF) * (  x / (1. - x) - x**(1. + fc.alpha) / (1. - x) ) + ( x * log(1. - x)) / (1. - x) - ( x**(1. + fc.alpha) * log(1. - x)) / (1. - x)
                }) ,
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -1: x**(1. + fc.alpha) / (1. - x) ,
                    0: - fc.A1(fc.beta_IF) * (  x / (1. - x) - x**(1. + fc.alpha) / (1. - x) ) + ( x * log(1. - x)) / (1. - x) - ( x**(1. + fc.alpha) * log(1. - x)) / (1. - x)
                }) ,
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -1: 1. - fc.A2(fc.alpha) ,
                    0: 1. - (pi**2 / 12.) + fc.A1(fc.beta_IF) * (1. - fc.A2(fc.alpha)) - (fc.A2(fc.alpha))**2/ 2. + fc.polygamma(fc.alpha) / 2.
                })
            }

        elif recoiler <= 2:
            kernel = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: x**(1. + fc.alpha) / (1. - x) ,
                    0: - fc.A1(fc.beta_II) * (  x / (1. - x) - x**(1. + fc.alpha) / (1. - x) ) + 2. * ( x * log(1. - x)) / (1. - x) - 2. * ( x**(1. + fc.alpha) * log(1. - x)) / (1. - x)
                }) ,
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -1: x**(1. + fc.alpha) / (1. - x) ,
                    0: - fc.A1(fc.beta_II) * (  x / (1. - x) - x**(1. + fc.alpha) / (1. - x) ) + ( x * log(1. - x)) / (1. - x) - ( x**(1. + fc.alpha) * log(1. - x)) / (1. - x)
                }) ,
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -1: 1. - fc.A2(fc.alpha) ,
                    0: 2. - (pi**2 / 6.) + fc.A1(fc.beta_II) * (1. - fc.A2(fc.alpha)) - (fc.A2(fc.alpha))**2 + fc.polygamma(fc.alpha)
                })
            }

        if distribution_type == 'counterterm':
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor, mu2os / xi)
        else:
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor, mu2os)

        logger.info('kernel_eval :' + str(kernel_eval))

        evaluation['values'][(0, 0, 0, 0)] = (kernel_eval).truncate(max_power=0)

        return evaluation



class QCD_integrated_TRN_CS_IgFg(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

    # Now define all the matching singular structures. We can combine in this class *all* integrated
    # initial-final counterterm thanks to the fact that *allowed_backward_evolved_flavors* is provided at
    # runtime. But it would be cleaner (and possible already with this structure) to implement each
    # kernel (qg, qq, etc...) in a separate class with its matching singular structure every time.
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
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),)
        ),)
    )


    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... g(initial) > g(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('CS_IgFg')


        # Retrieve variables
        beam_number = global_variables['overall_children'][0][0]
        if beam_number not in [1,2]:
            raise CurrentImplementationError(
                "The current %s must always be called for current for which the external" % self.__class__.__name__ +
                " leg number of the beam structure is 1 or 2, not %d."%beam_number)
        # Offset by one so as to have python list access conventions
        beam_number -= 1

        xi = global_variables['xis'][beam_number]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors[(beam_number+1)%2] != 'ALL':
            raise CurrentImplementationError('The current %s must always be called with' % self.__class__.__name__ +
                "allowed_backward_evolved_flavors[%d (i.e. inactive beam number)] = 'ALL' not %s" % ((beam_number+1)%2,str(
                allowed_backward_evolved_flavors[beam_number])))

        distribution_type = self.currents_properties[0]['distribution_type']
        beam_PDGs = self.currents_properties[0]['beam_PDGs']
        recoiler = global_variables['recoiler']


        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
#        if distribution_type == 'counterterm':
#            xis = list(global_variables['xis'])
#            xis[beam_number] = 1.0
#            evaluation['Bjorken_rescalings'] = [tuple(xis),]
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]


        Q = global_variables['Q']

        # Obtain Q_square.
        Q_square = Q.square()

#        # Only up to the order epsilon^2 of the scales prefactor matters here.
#        logMuQ = log(mu_r ** 2 / Q_square)
#        prefactor = EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})
#        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )

        # The additional 1/x part of the prefactor is included later during the PDF
        # convolution of the event (using its 'Bjorken rescaling' attribute) because
        # we must make sure that the plus distribution hits on it.
        # Also, the same 1/x appears in the PDF counterterms as a result of the change
        # of variable necessary to bring them in the form where the plus distribution
        # only acts on the PDF. So it makes sense to keep it completely factorised.

        # Input variables
        #y_0 = factors_and_cuts.y_0_prime
        #logy0 = log(y_0)
        # Assign a fake x for now if the distribution type is 'endpoint'
        # TODO: this is not optimal, eventually we should put each of these three pieces in
        # separate currents
        if distribution_type == 'endpoint':
            x = 0.5
        else:
            x = xi

#        # In MadNkLO, we use the change of variable xb' = xb*xi so that the factor
#        # (Q^2)^\eps in Eq. 5.21 of https://arxiv.org/pdf/0903.1218.pdf actually reads
#        # (Q^2/(xi1*xi2))^\eps and the '+' distributions also act on it, which we realize
#        # by simply multiplying the Q^2 provided by the xi factor that must be set to one.
#        logMuQ_plus = log(mu_r ** 2 / (Q_square * x))
#        prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
#        prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))
        mu2os = mu_r**2 / (p_p + p_r).square()

        # The additional 1/x part of the prefactor is included later during the PDF
        # convolution of the event (using its 'Bjorken rescaling' attribute) because
        # we must make sure that the plus distribution hits on it.
        # Also, the same 1/x appears in the PDF counterterms as a result of the change
        # of variable necessary to bring them in the form where the plus distribution
        # only acts on the PDF. So it makes sense to keep it completely factorised.

        # Input variables
        #y_0 = factors_and_cuts.y_0_prime
        #logy0 = log(y_0)
        # Assign a fake x for now if the distribution type is 'endpoint'
        # TODO: this is not optimal, eventually we should put each of these three pieces in
        # separate currents
        if distribution_type == 'endpoint':
            x = 0.5
        else:
            x = xi

        # In MadNkLO, we use the change of variable xb' = xb*xi so that the factor
        # (Q^2)^\eps in Eq. 5.21 of https://arxiv.org/pdf/0903.1218.pdf actually reads
        # (Q^2/(xi1*xi2))^\eps and the '+' distributions also act on it, which we realize
        # by simply multiplying the Q^2 provided by the xi factor that must be set to one.
#        logMuQ_plus = log(mu_r ** 2 / (Q_square * x))
#        prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
#        prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )


        # Define the NLO QCD integrate initial-state single collinear counterterms kernels


        color_factor = self.CA
        overall = - 2.
        if recoiler > 2:
            kernel = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: x**(1. + fc.alpha) / (1. - x) ,
                    0: - fc.A1(fc.beta_IF) * (  x / (1. - x) - x**(1. + fc.alpha) / (1. - x) ) + ( x * log(1. - x)) / (1. - x) - ( x**(1. + fc.alpha) * log(1. - x)) / (1. - x)
                }) ,
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -1: x**(1. + fc.alpha) / (1. - x)  ,
                    0: - fc.A1(fc.beta_IF) * (  x / (1. - x) - x**(1. + fc.alpha) / (1. - x) ) + ( x * log(1. - x)) / (1. - x) - ( x**(1. + fc.alpha) * log(1. - x)) / (1. - x)
                }) ,
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -1: 1. - fc.A2(fc.alpha) ,
                    0: 1. - (pi**2 / 12.) + fc.A1(fc.beta_IF) * (1. - fc.A2(fc.alpha)) - (fc.A2(fc.alpha))**2/ 2. + fc.polygamma(fc.alpha) / 2.
                })
            }

        elif recoiler <= 2:
            kernel = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -1: x**(1. + fc.alpha) / (1. - x) ,
                    0: - fc.A1(fc.beta_II) * (  x / (1. - x) - x**(1. + fc.alpha) / (1. - x) ) + 2. * ( x * log(1. - x)) / (1. - x) - 2. * ( x**(1. + fc.alpha) * log(1. - x)) / (1. - x)
                }) ,
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -1: x**(1. + fc.alpha) / (1. - x)  ,
                    0: - fc.A1(fc.beta_II) * (  x / (1. - x) - x**(1. + fc.alpha) / (1. - x) ) + ( x * log(1. - x)) / (1. - x) - ( x**(1. + fc.alpha) * log(1. - x)) / (1. - x)
                }) ,
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -1: 1. - fc.A2(fc.alpha) ,
                    0: 2. - (pi**2 / 6.) + fc.A1(fc.beta_II) * (1. - fc.A2(fc.alpha)) - (fc.A2(fc.alpha))**2 + fc.polygamma(fc.alpha)
                })
            }

        if distribution_type == 'counterterm':
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor, mu2os / xi)
        else:
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor, mu2os)

        logger.info('kernel_eval :' + str(kernel_eval))

        evaluation['values'][(0, 0, 0, 0)] = (kernel_eval).truncate(max_power=0)


        return evaluation



#=========================================================================================
# Integrated soft counterterm
#=========================================================================================


class QCD_integrated_TRN_S_g(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

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
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... soft g(final)
    current_properties['singular_structure'] = soft_structure_g
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... soft g(final)
    current_properties['singular_structure'] = soft_structure_g
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('S_g')


        xi = global_variables['xis'][0]
        xi2 = global_variables['xis'][1]
        mu_r = global_variables['mu_r']

        logger.info('xis : ' + str(global_variables['xis']))

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        # Retrieve the reduced kinematics as well as the soft momentum
        lower_PS_point = all_steps_info[0]['lower_PS_point']
        higher_PS_point = all_steps_info[0]['higher_PS_point']
#        logger.info('all_steps_info :' + str(all_steps_info[0]))
#        logger.info('all_steps_info :' + str(all_steps_info[-1]))

        # Normalization factors
        norm = -1./2.

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with ' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        distribution_type = self.currents_properties[0]['distribution_type']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]

        if distribution_type == 'endpoint':
            xi = 0.5

#        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
#            all_steps_info[0]['lower_PS_point'][2])
#        mu2os = mu_r**2 / shat

#        logger.info('shat : ' + str(shat))

        Q = global_variables['Q']

        # Obtain Q_square.
        Q_square = Q.square()

        # Only up to the order epsilon^2 of the scales prefactor matters here.
#        logMuQ = log(mu_r ** 2 / Q_square)
#        prefactor = EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})
#        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )


        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
#        prefactor = self.SEpsilon * (1. / (16 * math.pi ** 2))

        # All contributions from the soft counterterms are color correlated so we should remove
        # the default value None from the list of color correlations
        evaluation['color_correlations'] = []

        color_correlation_index = 0
#        logger.info('enumerate(colored_partons) :' + str(colored_partons))


        for i, a in enumerate(colored_partons):
            logger.info('i :' + str(i))
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in colored_partons[i:]:
                # Write the eikonal for that pair
                if a != b:
                    mult_factor = 2.
                else:
                    continue
                    mult_factor = 1.

                logger.info('a :' + str(a))
                logger.info('b :' + str(b))
                logger.info('mult_factor :' + str(mult_factor))


                pa = lower_PS_point[a]
                pb = lower_PS_point[b]
#                shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
#                    all_steps_info[0]['lower_PS_point'][2])
#                mu2os = mu_r**2 / shat
                etapr = (pa + pb).square() / shat
                mu2os = mu_r**2 / (pa + pb).square()
                #logger.info('colored partons : ' + str(all_steps_info))


                #logMuQ_plus = log(1. / etapr)
                #prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
                #prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )
                prefactor = self.SEpsilon * (1. / (16 * pi**2) )

#                logger.info('p_a * p_b : ' + str((pa + pb).square()))
                logger.info('distribution_type : ' + str(distribution_type))


                x = xi

                if a > 2 and b > 2:
                    logger.info('soft caso A')
                    kernel = {
                        'bulk': EpsilonExpansion({
                            -2: 0. ,
                            -1: 0. ,
                            0: 0. 
                        }) ,
                        'counterterm': EpsilonExpansion({
                            -2: 0. ,
                            -1: 0. ,
                            0: 0. 
                        }) ,
                        'endpoint': prefactor * EpsilonExpansion({
                            -2: - 1. ,
                            -1: - 2. * fc.A2(fc.alpha) ,
                            0: - (math.pi**2 / 12.) - 2. * (fc.A2(fc.alpha))**2 + 4. * fc.polygamma(fc.alpha) 
                        })
                    }

                if a > 2 and b <= 2:
                    logger.info('soft caso B')
                    kernel = {
                        'bulk': prefactor * EpsilonExpansion({
                            -2: 0. ,
                            -1: x**(1. + fc.alpha) / (1. - x)  ,
                            0: ( x **(1. + fc.alpha) / (1. - x) ) * fc.A2(fc.alpha) - ( log(1.-x) * x **(1. + fc.alpha) ) / (1. - x)
                        }) ,
                        'counterterm': prefactor * EpsilonExpansion({
                            -2: 0. ,
                            -1: x**(1. + fc.alpha) / (1. - x) ,
                            0: ( x **(1. + fc.alpha) / (1. - x) ) * fc.A2(fc.alpha) - ( log(1.-x) * x **(1. + fc.alpha) ) / (1. - x) 
                        }) ,
                        'endpoint': prefactor * EpsilonExpansion({
                            -2: - 1. ,
                            -1: - 2. * fc.A2(fc.alpha) ,
                            0: - (math.pi**2 / 12.) - 2. * (fc.A2(fc.alpha))**2 + 2. * fc.polygamma(fc.alpha) 
                        })                    
                    }

                if a <= 2 and b > 2:
                    logger.info('soft caso C')
                    kernel = {
                        'bulk': prefactor * EpsilonExpansion({
                            -2: 0. ,
                            -1: x **(1. + fc.alpha) / (1. - x)  ,
                            0: ( x **(1. + fc.alpha) / (1. - x) ) * fc.A2(fc.alpha) - ( log(1.-x) * x **(1. + fc.alpha) ) / (1. - x)
                        }) ,
                        'counterterm': prefactor * EpsilonExpansion({
                            -2: 0. ,
                            -1: x **(1. + fc.alpha) / (1. - x) ,
                            0: ( x **(1. + fc.alpha) / (1. - x) ) * fc.A2(fc.alpha) - ( log(1.-x) * x **(1. + fc.alpha) ) / (1. - x) 
                        }) ,
                        'endpoint': prefactor * EpsilonExpansion({
                            -2: - 1. ,
                            -1: - 2. * fc.A2(fc.alpha) ,
                            0: - (math.pi**2 / 12.) - 2. * (fc.A2(fc.alpha))**2 + 2. * fc.polygamma(fc.alpha) 
                        })                    
                    }

                if a <= 2 and b <= 2:
                    logger.info('soft caso D')
                #manca parte finita con integrazione su entrambe le x
                    kernel = {
                        'bulk': prefactor * EpsilonExpansion({
                            -2: 0. ,
                            -1: 2. * x**(1. + fc.alpha) / (1. - x) ,
                            0: 2. * fc.A2(fc.alpha) * (x **(1. + fc.alpha) / (1. - x)) - 2. * ( ( log(1.-x) * x **(1. + fc.alpha) ) / (1. - x)) - ( log(x) * x **(1. + fc.alpha) ) / (1. - x)
                        }) ,
                        'counterterm': prefactor * EpsilonExpansion({
                            -2: 0. ,
                            -1: 2. * x**(1. + fc.alpha) / (1. - x) ,
                            0: 2. * fc.A2(fc.alpha) * (x **(1. + fc.alpha) / (1. - x)) - 2. * ( ( log(1.-x) * x **(1. + fc.alpha) ) / (1. - x))
                        }) ,
                        'endpoint': prefactor * EpsilonExpansion({
                            -2: - 1. ,
                            -1: - 2. * fc.A2(fc.alpha) ,
                            0: - (math.pi**2 / 12.) - 2. * fc.A2(fc.alpha) + fc.polygamma(fc.alpha)
                        })                    
                    }


                if kernel[distribution_type] is not None and distribution_type == 'counterterm':
                    kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * mult_factor, mu2os / x)
                elif kernel[distribution_type] is not None:
                    kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * mult_factor, mu2os)
                elif kernel[distribution_type] is None:
                    kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * mult_factor, mu2os)

                logger.info('kernel_eval :' + str(kernel_eval))

                evaluation['color_correlations'].append(((a, b),))
                evaluation['values'][(0, color_correlation_index, 0, 0)] = (kernel_eval).truncate(max_power=0)
                color_correlation_index += 1

        return evaluation


#=========================================================================================
# Integrated FF collinear counterterm
#=========================================================================================

class QCD_integrated_TRN_C_FqFg(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

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
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = coll_structure_qg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = coll_structure_qg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('C_FqFg')

        xi = global_variables['xis'][0]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with ' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        distribution_type = self.currents_properties[0]['distribution_type']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]

        if distribution_type == 'endpoint':
            xi = 0.5

#        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
#            all_steps_info[0]['lower_PS_point'][2])


        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
#        etapr = (p_p + p_r).square() / shat
        mu2os = mu_r**2 / (p_p + p_r).square()

        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))

#        logMuQ_plus = log(1. / etapr)
#        prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
#        prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )


#        logger.info('Q_square : ' + str(global_variables['Q'].square()))
#        logger.info('shat : ' + str(shat))

        recoiler = global_variables['recoiler']

#        logger.info('distribution : ' + str(distribution_type))
#        logger.info('global_variables : ' + str(global_variables))
        logger.info('parent :' + str(global_variables['overall_children'][0]))
        logger.info('recoiler : ' + str(recoiler))

        color_factor = self.CF
        overall = 1. / 2.

        if recoiler > 2:

            if distribution_type == "endpoint":
                kernel = color_factor * overall * EpsilonExpansion({
                    -2: 0.,
                    -1: - 1.,
                    0: - (1. + fc.A2(fc.beta_FF))
                })
            elif distribution_type == "bulk":
                kernel = 0.
            elif distribution_type == "counterterm":
                kernel = 0. 
            else:
                raise CurrentImplementationError("Distribution type '%s' not supported." % distribution_type)


        elif recoiler <= 2:

            if distribution_type == "endpoint":
                kernel = color_factor * overall * EpsilonExpansion({
                    -2: 0.,
                    -1: - 1. ,
                    0: - (1. + fc.A2(fc.beta_FI))
                })
            elif distribution_type == "bulk":
                kernel = color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: xi **(1. + fc.beta_FI) / (1. - xi) 
                })
            elif distribution_type == "counterterm":
                kernel = color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: xi **(1. + fc.beta_FI) / (1. - xi)
                }) 
            else:
                raise CurrentImplementationError("Distribution type '%s' not supported." % distribution_type)


        if distribution_type == 'counterterm':
            kernel_eval = torino_to_madnk_epsexp(kernel * prefactor, mu2os / xi)
        else:
            kernel_eval = torino_to_madnk_epsexp(kernel * prefactor, mu2os)

        logger.info('kernel_eval :' + str(kernel_eval))

        evaluation['values'][(0, 0, 0, 0)] = (kernel_eval).truncate(max_power=0)

        return evaluation


class QCD_integrated_TRN_C_FqFqx(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

    # Now define all the matching singular structures.
    coll_structure_qqx = sub.SingularStructure(
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
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = coll_structure_qqx
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = coll_structure_qqx
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('C_FqFqx')


        xi = global_variables['xis'][0]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with ' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        distribution_type = self.currents_properties[0]['distribution_type']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]

#        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
#            all_steps_info[0]['lower_PS_point'][2])
#        mu_r = global_variables['mu_r']
#        mu2os = mu_r**2 / shat

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
        mu2os = mu_r**2 / (p_p + p_r).square() 
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))

        recoiler = global_variables['recoiler']

#        logger.info('distribution : ' + str(distribution_type))
#        logger.info('global_variables : ' + str(global_variables))
        logger.info('recoiler : ' + str(recoiler))

#        logger.info('NF : ' + str(NF))
        color_factor = self.TR
        overall = 2. / 3.

        if recoiler > 2:
            kernel = {
                'bulk': 0. ,
                'counterterm': 0. ,
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -2: 0.,
                    -1: - 1.,
                    0: - (5./3. + fc.A2(fc.beta_FF))
                })
            }

        elif recoiler <=2:
            kernel = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: x**(1. + fc.beta_FI) / (1. - x)
                }),
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: x**(1. + fc.beta_FI) / (1. - x)
                }),
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -2: 0.,
                    -1: - 1.,
                    0: - (5./3. + fc.A2(fc.beta_FI)) 
                })            
            }

        if distribution_type == 'counterterm':
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor , mu2os / xi)
        else:
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor , mu2os)

        logger.info('kernel_eval :' + str(kernel_eval))

        evaluation['values'][(0, 0, 0, 0)] = (kernel_eval).truncate(max_power=0)

        return evaluation

class QCD_integrated_TRN_C_FgFg(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

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
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = coll_structure_gg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = coll_structure_gg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('C_FgFg')


        xi = global_variables['xis'][0]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with ' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        distribution_type = self.currents_properties[0]['distribution_type']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]

#        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
#            all_steps_info[0]['lower_PS_point'][2])
#        mu2os = mu_r**2 / shat

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
#        etapr = (p_p + p_r).square() / shat
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))
        mu2os = mu_r**2 / (p_p + p_r).square()

        recoiler = global_variables['recoiler']

#        logger.info('distribution : ' + str(distribution_type))
#        logger.info('p_p * p_r : ' + str((p_p + p_r).square()))
        logger.info('parent :' + str(global_variables['overall_children'][0]))
        logger.info('recoiler : ' + str(recoiler))

        color_factor = self.CA
        overall = 1./6.

        if recoiler > 2:
            kernel = {
                'bulk': 0. ,
                'counterterm': 0. ,
                'endpoint': color_factor * EpsilonExpansion({
                    -2: 0. ,
                    -1: - overall ,
                    0: - overall * (5./3. + fc.A2(fc.beta_FF))
                })
            }


        elif recoiler <=2:
            kernel = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: x**(1. + fc.beta_FI) / (1. - x) 
                }),
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: x**(1. + fc.beta_FI) / (1. - x)
                }),
                'endpoint':  color_factor * EpsilonExpansion({
                    -2: 0. ,
                    -1: - overall ,
                    0: - overall * (5./3. + fc.A2(fc.beta_FF) )
                })
            }

        if distribution_type == 'counterterm':
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor * 2., mu2os / xi)
        else:
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor * 2. , mu2os)

        evaluation['values'][(0, 0, 0, 0)] = (kernel_eval).truncate(max_power=0)

        return evaluation


#=========================================================================================
# Integrated FF soft-collinear counterterm
#=========================================================================================

class QCD_integrated_TRN_CS_FqFg(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

    # Now define all the matching singular structures.
    soft_structure_g = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
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
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        logger.info('CS_FqFg')

        xi = global_variables['xis'][0]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with ' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        distribution_type = self.currents_properties[0]['distribution_type']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)

        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]

        if distribution_type == 'endpoint':
            xi = 0.5



#        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
#            all_steps_info[0]['lower_PS_point'][2])
#        mu_r = global_variables['mu_r']
#        mu2os = mu_r**2 / shat

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
        mu2os = mu_r**2 / (p_p + p_r).square()
#        etapr = (p_p + p_r).square() / shat
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))

        recoiler = global_variables['recoiler']

#        logger.info('shat : ' + str(shat))
#        logger.info('global_variables : ' + str(global_variables))
#        logger.info('recoiler : ' + str(recoiler))
        logger.info('distribution_type :' + str(distribution_type))

        color_factor = self.CF
        overall = - 2.

        if recoiler > 2:
            kernel = {
                'bulk': 0. ,
                'counterterm': 0. ,
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -2: 0.,
                    -1: 1. - fc.A2(fc.alpha) ,
                    0: 2. - (pi**2 / 4.) + fc.A2(fc.beta_FF) * (1. - fc.A2(fc.alpha)) - (fc.A2(fc.alpha))**2 / 2. + (3./2.) * fc.polygamma(fc.alpha)
                })
            }


        elif recoiler <=2:
            kernel = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: (xi **(1. + fc.beta_FI) / (1. - xi)) * (- 1. + fc.A2(fc.alpha)) 
                }),
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: (xi **(1. + fc.beta_FI) / (1. - xi)) * (- 1. + fc.A2(fc.alpha))
                }),
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -2: 0.,
                    -1: 1. - fc.A2(fc.alpha) ,
                    0: 2. - (pi**2 / 4.) + fc.A2(fc.beta_FI) * (1. - fc.A2(fc.alpha)) - (fc.A2(fc.alpha))**2 / 2. + (3./2.) * fc.polygamma(fc.alpha)
                })            
            }

        if distribution_type == 'counterterm':
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor, mu2os / xi)
        else:
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor, mu2os)

        logger.info('kernel_eval :' + str(kernel_eval))

        evaluation['values'][(0, 0, 0, 0)] = (kernel_eval).truncate(max_power=0)

        return evaluation


class QCD_integrated_TRN_CS_FgFg(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = variables_for_integrated_currents_FF

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

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(final) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):

        logger.info('CS_FgFg')

        """ Evaluate this counterterm given the variables provided. """

        xi = global_variables['xis'][0]
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with ' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        distribution_type = self.currents_properties[0]['distribution_type']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]

#        shat = 2 * all_steps_info[0]['lower_PS_point'][1].dot( \
#            all_steps_info[0]['lower_PS_point'][2])
#        mu2os = mu_r**2 / shat

        p_p = all_steps_info[0]['lower_PS_point'][global_variables['overall_parents'][0]]
        p_r = all_steps_info[0]['lower_PS_point'][global_variables['recoiler']]
#        etapr = (p_p + p_r).square() / shat
        mu2os = mu_r**2 / (p_p + p_r).square()
        prefactor = self.SEpsilon * (1. / (16 * pi ** 2))

        recoiler = global_variables['recoiler']

        logger.info('distribution : ' + str(distribution_type))
#        logger.info('global_variables : ' + str(global_variables))
#        logger.info('recoiler : ' + str(recoiler))

        color_factor = self.CA
        overall = - 2.

        if recoiler > 2:
            kernel = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: 0. 
                }) ,
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: 0. 
                }) ,
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -2: 0.,
                    -1: 1. - fc.A2(fc.alpha) ,
                    0: 2. - (pi**2 / 4.) + fc.A2(fc.beta_FF) * (1. - fc.A2(fc.alpha)) - (fc.A2(fc.alpha))**2 / 2. + (3./2.) * fc.polygamma(fc.alpha) 
                })
            }

        elif recoiler <=2:
            kernel = {
                'bulk': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: (x**(1. + fc.beta_FI) / (1. - x)) * (- 1. + fc.A2(fc.alpha)) 
                }),
                'counterterm': color_factor * overall * EpsilonExpansion({
                    -2: 0. ,
                    -1: 0. ,
                    0: (x**(1. + fc.beta_FI) / (1. - x)) * (- 1. + fc.A2(fc.alpha))
                }),
                'endpoint': color_factor * overall * EpsilonExpansion({
                    -2: 0.,
                    -1: 1. - fc.A2(fc.alpha) ,
                    0: 2. - (pi**2 / 4.) + fc.A2(fc.beta_FI) * (1. - fc.A2(fc.alpha)) - (fc.A2(fc.alpha))**2 / 2. + (3./2.) * fc.polygamma(fc.alpha)
                })            
            }


        if distribution_type == 'counterterm':
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor, mu2os / xi)
        else:
            kernel_eval = torino_to_madnk_epsexp(kernel[distribution_type] * prefactor, mu2os)

        logger.info('kernel_eval :' + str(kernel_eval))

        evaluation['values'][(0, 0, 0, 0)] = (kernel_eval).truncate(max_power=0)
        return evaluation


