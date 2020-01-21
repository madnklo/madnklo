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
"""Implementation of NLO type of integrated currents for colorful_pp."""

import os
import sys
import math

import madgraph.core.subtraction as subtraction
import madgraph.integrator.phase_space_generators as PS_utils
from madgraph.core.base_objects import EpsilonExpansion
import colorful_pp_config
import madgraph.various.misc as misc

import madgraph.various.math_tools.mpl as MPL

import commons.utils as utils
import commons.QCD_local_currents as currents
import commons.factors_and_cuts as factors_and_cuts

from integrated_current_expressions import HE
import commons.general_current as general_current
import madgraph.core.subtraction as sub

log = math.log
pi = math.pi

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError

# In principle we support them all, but these are the meaningful ones to consider
beam_PDGs_supported = [
    tuple(sorted([1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 21])),
    tuple(sorted([1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 21])),
    tuple(sorted([1, -1, 2, -2, 3, -3, 4, -4, 21])),
    tuple(sorted([1, -1, 2, -2, 3, -3, 21])),
    tuple(sorted([1, -1, 2, -2, 21])),
    tuple(sorted([1, -1, 21]))
]

#=========================================================================================
# NLO PDF Counterterm
#=========================================================================================

class QCD_F0(general_current.GeneralCurrent):
    """Implements the NLO QCD PDF counterterm of type F(xi)"""

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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
        'beam_PDGs' : beam_PDGs_supported
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
        beam_PDGs = self.currents_properties[0]['beam_PDGs']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            xis = list(global_variables['xis'])
            xis[beam_number] = 1.0
            evaluation['Bjorken_rescalings'] = [tuple(xis),]

        # Only the order epsilon of the scales pre-factor matters here.
        prefactor = EpsilonExpansion({
            0: 1.,
            1: log(mu_r ** 2 / mu_f ** 2)
        })
        prefactor *= EpsilonExpansion({-1: 1.}) * self.SEpsilon * (1. / (16 * pi**2) )

        # Assign a fake xi for now if the distribution type is 'endpoint'
        # This is not optimal, eventually we should put each of these three pieces in
        # separate currents, this will be optimal however in the low-level implementation
        if distribution_type == 'endpoint':
            xi = 0.5

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

#=========================================================================================
# Integrated IF collinear counterterm
#=========================================================================================

class QCD_integrated_C_IF(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define all the matching singular structures. We can combine in this class *all* integrated
    # initial-final counterterm thanks to the fact that *allowed_backward_evolved_flavors* is provided at
    # runtime. But it would be cleaner (and possible already with this structure) to implement each
    # kernel (qg, qq, etc...) in a separate class with its matching singular structure every time.
    coll_structure_qq = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
                sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),)
        ),)
    )
    coll_structure_gg = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
                sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),)
        ),)
    )
    coll_structure_qg = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
                sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),)
        ),)
    )
    coll_structure_gq = sub.SingularStructure(
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
    # ... g(initial) > g(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_gg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_qg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))
    # ... q(initial) > g(initial_after_emission) q(final)
    current_properties['singular_structure'] = coll_structure_gq
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... g(initial) > g(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_qq
    template_currents.append(sub.BeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_gg
    template_currents.append(sub.BeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = coll_structure_qg
    template_currents.append(sub.BeamCurrent( current_properties ))
    # ... q(initial) > g(initial_after_emission) q(final)
    current_properties['singular_structure'] = coll_structure_gq
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]


    # Helper function to implement the "allowed_backward_evolved_flavors" masking
    @staticmethod
    def apply_flavor_mask(flavor_matrix,allowed_backward_evolved_flavors):
        """
        Given a flavor matrix and a list of permitted flavors that can be the end point of the backward evolution
        (the flavor mask), generate a filtered flavor matrix.

        :param flavor_matrix: sparse matrix implemented as a dict of dict
                (see madgraph.integrator.ME7_integrands.ME7Event#convolve_flavors)
        :param allowed_backward_evolved_flavors: tuple of PDG ids
        :return: filtered_flavor_matrix, with the same structure as flavor_matrix
        """
        # If no mask do nothing
        if allowed_backward_evolved_flavors == 'ALL':
            return flavor_matrix
        else:
            # We will loop over a matrix M[i][j] = wgt_ij where i is the starting PDG and j a tuple of ending PDGs.
            # We filter the elements in j.
            filtered_flavor_matrix = {}
            # Loop over matrix lines i
            for reduced_flavor in flavor_matrix:
                # Each column is a dict {(PDG1.1, PDG1.2, PDG1.3) : wgt1,...} where
                # the tuple j=(PDG1.1, PDG1.2, PDG1.3) is the label of a column  and wgt the entry ij of the matrix
                for end_flavors, wgt in  flavor_matrix[reduced_flavor].items():
                    #Apply the filter on the elements of the tuple
                    allowed_end_flavors = tuple([fl for fl in end_flavors if fl in allowed_backward_evolved_flavors])
                    #If a column was entirely filtered out, do not include it
                    if allowed_end_flavors:
                        filtered_flavor_matrix[reduced_flavor] = {allowed_end_flavors:wgt}

            return filtered_flavor_matrix

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

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
        logMuQ = log(mu_r ** 2 / Q_square)
        prefactor = EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})
        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )

        # The additional 1/x part of the prefactor is included later during the PDF
        # convolution of the event (using its 'Bjorken rescaling' attribute) because
        # we must make sure that the plus distribution hits on it.
        # Also, the same 1/x appears in the PDF counterterms as a result of the change
        # of variable necessary to bring them in the form where the plus distribution
        # only acts on the PDF. So it makes sense to keep it completely factorised.

        # Input variables
        y_0 = factors_and_cuts.y_0_prime
        logy0 = log(y_0)
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
        logMuQ_plus = log(mu_r ** 2 / (Q_square * x))
        prefactor_plus = EpsilonExpansion({0: 1., 1: logMuQ_plus, 2: 0.5 * logMuQ_plus ** 2})
        prefactor_plus *= self.SEpsilon * (1. / (16 * pi**2) )

        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # Define the NLO QCD integrate initial-state single collinear counterterms kernels
        color_factor = self.CA
        kernel_gg = {
            'bulk': prefactor * color_factor * (EpsilonExpansion({
                -1: -2. * (1. / (1. - x) + (1. - x) / x - 1 + x * (1 - x)),
                0: (2. * log1mx / (1. - x)) * (1. + theta_x_1my0) + (2. * logy0 / (1. - x)) * theta_1my0_x
                   + 2. * (((1. - x) / x) - 1. + x * (1. - x)) * (log1mx * (1. + theta_x_1my0) + logy0 * theta_1my0_x)
            })),
            'counterterm': prefactor_plus * color_factor * (EpsilonExpansion({
                -1: -2. * (1. / (1. - x)),
                0: (2. * log1mx / (1. - x)) * (1. + theta_x_1my0),
            })),
            'endpoint': prefactor * color_factor * (EpsilonExpansion({
                -2: 1.,
                -1: 0.,
                0: -(math.pi ** 2 / 6.) + logy0 ** 2
            }))
        }

        color_factor = self.CA
        kernel_gq = {
            'bulk': prefactor * color_factor * (EpsilonExpansion({
                -1: -(self.CF / self.CA) * (1. + (1. - x) ** 2) / x,
                0: (self.CF / self.CA) * (
                            ((1. + (1. - x) ** 2) / x) * (log1mx * (1. + theta_x_1my0) + logy0 * theta_1my0_x) + x)
            })),
            'counterterm': None,
            'endpoint': None
        }

        color_factor = self.CF
        kernel_qg = {
            'bulk': prefactor * color_factor * (EpsilonExpansion({
                -1: -(self.TR / self.CF) * (x ** 2 + (1. - x) ** 2),
                0: (self.TR / self.CF) * ((x ** 2 + (1. - x) ** 2) * (
                            log1mx * (1. + theta_x_1my0) + logy0 * theta_1my0_x) + 2. * x * (1. - x))
            })),
            'counterterm': None,
            'endpoint': None
        }

        color_factor = self.CF
        kernel_qq = {
            'bulk': prefactor * color_factor * (EpsilonExpansion({
                -1: -((1. + x ** 2) / (1. - x)),
                0: (2. * log1mx / (1. - x)) * (1. + theta_x_1my0) + (2. * logy0 / (1. - x)) * theta_1my0_x
                   - ((1. + x) * (log1mx * (1. + theta_x_1my0) + logy0 * theta_1my0_x) - 1. + x)
            })),
            'counterterm': prefactor_plus * color_factor * (EpsilonExpansion({
                -1: -((1. + x ** 2) / (1. - x)),
                0: (2. * log1mx / (1. - x)) * (1. + theta_x_1my0),
            })),
            'endpoint': prefactor * color_factor * (EpsilonExpansion({
                -2: 1.,
                -1: 3. / 2.,
                0: -(math.pi ** 2 / 6.) + logy0 ** 2
            }))
        }

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
# Integrated soft counterterm, which recoils against initial state
# This therefore yields a correlated convolution against initial-state beams
#=========================================================================================

class QCD_integrated_S_Fg(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm of type C(g,x)(x_i)"""

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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


    ############################################################################################################
    # Beam-soft integrated counterterms
    # These are (Eikonal - its collinear limits) integrated over the unresolved phase space of Vittorio's
    # soft-recoiling-against-initial-state mapping. This is defined in sec 5.1.2 of his notes
    # Each of these should multiply a color-correlated amplitude and when summed over colors reconstruct the S+CS counterterms
    # Detailed information in Nicolas Deutschmann's handwritten notes from 03.10.2018
    # calculation in ndeutsch/workdir/research/Active/Subtractions/ReverseUnitarityCT/\
    # Vittorio_Soft/ExpandedSbars
    ############################################################################################################
    @classmethod
    def integrated_bs_endpoint_pole(cls, dipole_invariant):
        """Pole of the endpoint contribution of the collinear-subtracted eikonal integrated over the SoftVsInitial mapping unresolved PS (BS)"""
        return 2.*math.log(dipole_invariant)

    @classmethod
    def integrated_bs_endpoint_finite(cls, dipole_invariant):
        """Finite part of the endpoint contribution of the collinear-subtracted eikonal integrated over the SoftVsInitial mapping unresolved PS (BS)"""
        return  - math.log(16.)*math.log(dipole_invariant) - math.log(dipole_invariant)**2 + 2.*(MPL.G([1,0],dipole_invariant)-math.pi**2/6.)

    @classmethod
    def integrated_bs_bulk_finite(cls, dipole_invariant,xi):
        """Finite part of the bulk contribution of the collinear-subtracted eikonal integrated over the SoftVsInitial mapping unresolved PS (BS)"""
        return (8. * xi * math.log(dipole_invariant)) / ((-1. + xi) * (1. + xi))

    @classmethod
    def integrated_bs_counterterm_finite(cls, dipole_invariant,xi):
        """Finite part of the counterterm contribution of the collinear-subtracted eikonal integrated over the SoftVsInitial mapping unresolved PS (BS)"""
        return (4. * math.log(dipole_invariant) ) / (-1. + xi)

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        xi = global_variables['xis'][0]
        mu_r = global_variables['mu_r']

        PS_point = all_steps_info[-1]['lower_PS_point']

        # Represent the colored parton info dictionary as a sorted list
        colored_partons = sorted([
            (leg_number, info) for leg_number, info in colored_partons.items() ], key=lambda el: el[0]
        )

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

        Q = global_variables['Q']

        # Obtain Q_square.
        Q_square = Q.square()

        # Only up to the order epsilon^2 of the scales prefactor matters here.
        logMuQ = log(mu_r**2/Q_square)
        # Correction for the counterterm: in BS (bulk+counterterm), the variable Q_square corresponds to that
        # of the real event. However the counterterm corresponds to the residue of the bulk at xi=1.
        # This is effectively obtained by multiplying by xi: Q_residue = Q_real * xi.
        # Note for future dumb-me: log(mu_r**2/(Q_square*xi**2)) = logMuQ - log(xi**2)
        if distribution_type == 'counterterm':
            logMuQ-=log(xi**2)
        prefactor = EpsilonExpansion({ 0 : 1., 1 : logMuQ, 2 : 0.5*logMuQ**2 })
        prefactor *= self.SEpsilon * (1. / (16 * pi**2) )

        # All contributions from the soft counterterms are color correlated so we should remove
        # the default value None from the list of color correlations
        evaluation['color_correlations'] = []

        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        for i, (a, leg_info_a) in enumerate(colored_partons):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for (b, leg_info_b) in colored_partons[i+1:]:
                # Write the integrated eikonal for that pair

                evaluation['color_correlations'].append( ((a, b), ) )

                pa = PS_point[a]
                pb = PS_point[b]

                # We can only handle massless particles
                try:
                    assert pa.square()/Q_square < 1.e-09
                except AssertionError:
                    misc.sprint("No massive particles in soft currents for now")
                    raise

                # Assign the type of dipole
                # dipole_type = [bool_a, bool_b] tags the dipole tag, with True indicating initial state and False final
                dipole_type = [ leg_info_a[1] == sub.SubtractionLeg.INITIAL,
                                leg_info_b[1] == sub.SubtractionLeg.INITIAL ]

                if all(dipole_type): # Initial-initial
                    # Initial-initial: S+CS = 0
                    if distribution_type == 'bulk':
                        kernel = EpsilonExpansion({0:0})
                    elif distribution_type == 'counterterm':
                        kernel = EpsilonExpansion({0:0})
                    elif distribution_type == 'endpoint':
                        kernel = EpsilonExpansion({0:0})
                    else:
                        raise CurrentImplementationError("Distribution type '%s' not supported."%distribution_type)
                else: # At least one leg final
                    # The integrated counterterms are evaluated in terms of
                    # dipole_invariant = 1-cos(angle between the dipole momenta)
                    dipole_invariant = 0.5*pa.dot(pb)*Q.square()/(pa.dot(Q)*pb.dot(Q))
                    if distribution_type == 'bulk':
                        #The factor xi^2 below corrects the flux factor used in the bulk BS which has a 1/xi^2 too many
                        #A more permanent change is warranted after testing.
                        #See github issue #9 for reference
                        kernel = EpsilonExpansion({0:xi**2*self.integrated_bs_bulk_finite(dipole_invariant,xi)})
                    elif distribution_type == 'counterterm':
                        kernel = EpsilonExpansion({0:self.integrated_bs_counterterm_finite(dipole_invariant,xi)})
                    elif distribution_type == 'endpoint':
                        kernel = EpsilonExpansion({-1:self.integrated_bs_endpoint_pole(dipole_invariant),
                                                    0:self.integrated_bs_endpoint_finite(dipole_invariant)})
                    else:
                        raise CurrentImplementationError("Distribution type '%s' not supported."
                                                                        %distribution_type)

                    # Former implementation of the II soft+SC. Commented by Nicolas
                    # While no longer useful, this is kept for now to remember how non-zero integrated soft shoud be implemented
                    # if distribution_type == 'bulk':
                    #     kernel = EpsilonExpansion({
                    #                 0 : - 16.*xi * log(1.-xi**2) /(1.-xi**2),
                    #                 -1 : 8. * xi / (1.-xi**2),
                    #                 -2 : 0.
                    #     })
                    # elif distribution_type == 'counterterm':
                    #     kernel = EpsilonExpansion({
                    #                 0 : -8.*log(2.*(1.-xi))/(1.-xi),
                    #                 -1 : 4./(1.-xi),
                    #                 -2 : 0.
                    #     })
                    # elif distribution_type == 'endpoint':
                    #     kernel = EpsilonExpansion({
                    #                 0 : pi**2./3.-4.*log(2.)**2,
                    #                 -1 : 4.*log(2.),
                    #                 -2 : -2.
                    #     })
                    # else:
                    #     raise CurrentImplementationError("Distribution type '%s' not supported."
                    #                                                     %distribution_type)

                evaluation['values'][(0, color_correlation_index, 0, 0)] = kernel*prefactor
                color_correlation_index += 1

        return evaluation

#=========================================================================================
# Integrated collinear FF counterterm, which recoils against initial state, but using a
# different mapping as for the soft. These CTs therefore also yield a correlated
# convolution against initial-state beams.
#=========================================================================================

# @Nicolas: Specify the lengthy analytical expression of the kernels as class methods
# of the corresponding classes, like it is done for the class implementing the integrated
# soft.

class QCD_integrated_C_FqFq(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... g(initial) > q(final) qx(final)
    current_properties['singular_structure'] = coll_structure_qq
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... g(initial) > q(final) qx(final)
    current_properties['singular_structure'] = coll_structure_qq
    template_currents.append(sub.BeamCurrent( current_properties ))

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
                                             "allowed_backward_evolved_flavors=('ALL','ALL') not %s" % str(
                allowed_backward_evolved_flavors))

        distribution_type = self.currents_properties[0]['distribution_type']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]

        Q = global_variables['Q']

        PS_point = all_steps_info[-1]['lower_PS_point']

        # TODO

class QCD_integrated_C_FqFg(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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

        Q = global_variables['Q']

        PS_point = all_steps_info[-1]['lower_PS_point']

        # TODO

class QCD_integrated_C_FgFg(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... g(initial) > g(final) g(final)
    current_properties['singular_structure'] = coll_structure_gg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... g(initial) > g(final) g(final)
    current_properties['singular_structure'] = coll_structure_gg
    template_currents.append(sub.BeamCurrent( current_properties ))

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

        distribution_type = self.currents_properties[0]['distribution_type']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [(1.0, 1.0),]

        Q = global_variables['Q']

        PS_point = all_steps_info[-1]['lower_PS_point']

        # TODO

#=========================================================================================
# Integrated IF and FF soft-collinear counterterm
#=========================================================================================

class QCD_integrated_CS_IF(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single soft-collinear integrated counterterm of type C(S(g),x)(x_i)
    This is zero in this scheme as it was included directly together with the soft contribution.
    """

    # Enable the flag below to debug this current
    DEBUG = False

    # As mentioned in the docstring, this integrated counterterm is zero since it has already been accounted for
    # in the integrated soft CT.
    is_zero = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),)
        ),)
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
    # ... g(initial) > g(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.BeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

class QCD_integrated_CS_FF(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single soft-collinear integrated counterterm of type C(S(g),x)(x_i)
    This is zero in this scheme as it was included directly together with the soft contribution.
    """

    # Enable the flag below to debug this current
    DEBUG = False

    # As mentioned in the docstring, this integrated counterterm is zero since it has already been accounted for
    # in the integrated soft CT.
    is_zero = True

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... g(initial) > g(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg
    template_currents.append(sub.BeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]
