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
import factors_and_cuts as factors_and_cuts

import commons.general_current as general_current
import madgraph.core.subtraction as sub

import logging
logger = logging.getLogger("madgraph.integrated_currents")
# Change this to DEBUG to see the debug info
# They are quite verbose so we override the default mode
logger.setLevel(logging.INFO)

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

# Collinear integrated kernels

class TODO_QCD_integrated_C_IqFqpFqpx(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state single collinear integrated counterterm
    of type q(initial) > q(initial) q'(final) q'(final)"""

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, -2, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color': True,
                'n_loops': 0,
                'squared_orders': {'QCD': 4},
                'singular_structure': structure
            }),
        ),
    ]

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > qx(initial_after_emission) q'(final) q'x(final)
    current_properties['singular_structure'] = structure
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > qx(initial_after_emission) q'(final) q'x(final)
    current_properties['singular_structure'] = structure
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def get_bulk_kernel(self, x, y_0):
        """" Evaluates the bulk kernel. """

        logy0 = log(y_0)
        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # TODO code up the bulk kernel
        return EpsilonExpansion({
                0   :   1.
        })

    def get_counterterm_kernel(self, x, y_0):
        """" Evaluates the bulk kernel. """

        logy0 = log(y_0)
        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # TODO code up the counterterm kernel
        return EpsilonExpansion({
            0: 1.
        })

    def get_endpoint_kernel(self, x, y_0):
        """" Evaluates the bulk kernel. """

        logy0 = log(y_0)
        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # TODO code up the endpoint kernel
        return EpsilonExpansion({
            0: 1.
        })

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
        # Assign a fake x for now if the distribution type is 'endpoint'
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

        kernels = {
            'bulk': prefactor * self.get_bulk_kernel(x,y_0),
            'counterterm': prefactor_plus * self.get_counterterm_kernel(x,y_0),
            'endpoint': prefactor * self.get_endpoint_kernel(x,y_0)
        }

        active_quark_PDGs = self.currents_properties[0]['active_fermions']

        # Build the NLO flavor matrix
        flavor_matrix = {}

        for reduced_flavor in active_quark_PDGs:
            flavor_matrix[reduced_flavor] = {(reduced_flavor,) : kernels[distribution_type]}

        # Truncate all entries of the flavor matrix so as to remove irrelevant O(\eps) terms
        for flav_in, flav_outs in flavor_matrix.items():
            for flav_out, eps_expansion in flav_outs.items():
                eps_expansion.truncate(max_power=0)

        # Now on could apply the mask 'allowed_backward_evolved_flavors' if not set to 'ALL'
        # But that should not be necessary for this current because it refers to one specific splitting.
        #filtered_flavor_matrix = self.apply_flavor_mask(flavor_matrix, allowed_backward_evolved_flavors[beam_number])

        # Now assign the flavor matrix to the result
        evaluation['values'][(0,0,0,0)] = flavor_matrix

        # Promote the format of the flavor matrix to be that of a generalised two-sided convolution
        evaluation.promote_to_two_sided_convolution(beam_number+1)

        return evaluation

class TODO_QCD_integrated_C_IqFqFqx(TODO_QCD_integrated_C_IqFqpFqpx):
    """Implements the NLO QCD initial-state single collinear integrated counterterm
    of type q(initial) > q(initial) q(final) q(final)"""

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    coll_structure = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(12, -1, sub.SubtractionLeg.FINAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure,))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color': True,
                'n_loops': 0,
                'squared_orders': {'QCD': 4},
                'singular_structure': structure
            }),
        ),
    ]

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > qx(initial_after_emission) q(final) qx(final)
    current_properties['singular_structure'] = structure
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > qx(initial_after_emission) q(final) qx(final)
    current_properties['singular_structure'] = structure
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def get_bulk_kernel(self, x, y_0):
        """" Evaluates the bulk kernel. """

        logy0 = log(y_0)
        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # TODO code up the bulk kernel

        return EpsilonExpansion({
                0   :   1.
        })

    def get_counterterm_kernel(self, x, y_0):
        """" Evaluates the bulk kernel. """

        logy0 = log(y_0)
        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # TODO code up the counterterm kernel

        return EpsilonExpansion({
            0: 1.
        })

    def get_endpoint_kernel(self, x, y_0):
        """" Evaluates the bulk kernel. """

        logy0 = log(y_0)
        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # TODO code up the endpoint kernel

        return EpsilonExpansion({
            0: 1.
        })

# Nested collinear integrated kernels

class TODO_QCD_integrated_C_FqFqx_C_IqpFqFqx(general_current.GeneralCurrent):

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

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
    # Match both the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(sub_coll_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color': True,
                'n_loops': 0,
                'squared_orders': {'QCD': 4},
                'singular_structure': structure
            }),
        ),
    ]

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > qx(initial_after_emission) q'(final) q'x(final)
    current_properties['singular_structure'] = structure
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > qx(initial_after_emission) q'(final) q'x(final)
    current_properties['singular_structure'] = structure
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def get_bulk_kernel(self, x, y_0):
        """" Evaluates the bulk kernel. """

        logy0 = log(y_0)
        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # TODO code up the bulk kernel

        return EpsilonExpansion({
                0   :   1.
        })

    def get_counterterm_kernel(self, x, y_0):
        """" Evaluates the bulk kernel. """

        logy0 = log(y_0)
        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # TODO code up the counterterm kernel

        return EpsilonExpansion({
            0: 1.
        })

    def get_endpoint_kernel(self, x, y_0):
        """" Evaluates the bulk kernel. """

        logy0 = log(y_0)
        log1mx = log(1. - x)

        # Heaviside
        theta_x_1my0 = 1. if (x - (1 - y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1 - y_0) - x) >= 0. else 0.

        # TODO code up the endpoint kernel

        return EpsilonExpansion({
            0: 1.
        })


    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        # Retrieve variables
        initial_state_leg_number = global_variables['overall_children'][0][0]
        if initial_state_leg_number not in [1,2]:
            raise CurrentImplementationError(
                "The current %s must always be called for current for which the external" % self.__class__.__name__ +
                " leg number is 1 or 2, not %d."%beam_number)


        xi = global_variables['xis'][0]
        assert( xi == global_variables['xis'][1])

        mu_r = global_variables['mu_r']

        distribution_type = self.currents_properties[0]['distribution_type']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            evaluation['Bjorken_rescalings'] = [tuple(1.0, 1.0),]

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
        # Assign a fake x for now if the distribution type is 'endpoint'
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

        if distribution_type == 'bulk':
            kernel = prefactor * self.get_bulk_kernel(x,y_0)
        elif distribution_type == 'counterterm':
            kernel = prefactor * self.get_counterterm_kernel(x,y_0)
        elif distribution_type == 'endpoint':
            kernel = prefactor * self.get_endpoint_kernel(x,y_0)
        else:
            raise CurrentImplementationError(
                "Unsupported distirbution type '%s' in current %s" %( distribution_type, self.__class__.__name__))

        # Now assign the flavor matrix to the result
        evaluation['values'][(0,0,0,0)] = kernel.truncate(max_power=0)

        return evaluation

class TODO_QCD_integrated_C_FqFqx_C_IqFqFqx(TODO_QCD_integrated_C_FqFqx_C_IqpFqFqx):
    """ Excpet for the matching singular structure, this current should be identical to its mother class,
    i.e. TODO_QCD_integrated_C_FqFqx_C_IqpFqFqx """

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

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
    # Match both the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(sub_coll_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),))

    # This counterterm will be used if any of the currents_block of the list below matches
    defining_currents = [
        (
            sub.Current({
                'resolve_mother_spin_and_color': True,
                'n_loops': 0,
                'squared_orders': {'QCD': 4},
                'singular_structure': structure
            }),
        ),
    ]

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > qx(initial_after_emission) q(final) qx(final)
    current_properties['singular_structure'] = structure
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > qx(initial_after_emission) q(final) qx(final)
    current_properties['singular_structure'] = structure
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

# Soft integrated kernels

class TODO_QCD_integrated_S_FqFqx(general_current.GeneralCurrent):
    """Implements the NLO QCD double soft quark counterterm."""

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define all the matching singular structures.
    soft_structure_qqx = sub.SingularStructure(
        substructures=(sub.SoftStructure(
            substructures=tuple([]),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... soft g(final)
    current_properties['singular_structure'] = soft_structure_qqx
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))


    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... soft g(final)
    current_properties['singular_structure'] = soft_structure_qqx
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]


    ############################################################################################################
    # Beam-soft integrated counterterms
    # These are double quark (Eikonal - its collinear limits) integrated over the unresolved phase space of the
    # soft-recoiling-against-initial-state mapping. This is defined in sec 5.1.2 of the notes.
    # Each of these should multiply a color-correlated amplitude and when summed over colors reconstruct the S+CS counterterms
    ############################################################################################################
    def get_bulk_kernel(self, dipole_invariant, xi):

        # TODO code up the bulk kernel
        return EpsilonExpansion({
                0   :   1.
        })

    def get_counterterm_kernel(self, dipole_invariant, xi):

        # TODO code up the counterterm kernel
        return EpsilonExpansion({
                0   :   1.
        })

    def get_endpoint_kernel(self, dipole_invariant, xi):

        # TODO code up the endpoint kernel
        return EpsilonExpansion({
                0   :   1.
        })

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
                    if dipole_invariant > 1. +1.e-5:
                        raise ValueError(
                            "dipole_invariant = 0.5*pa.dot(pb)*Q.square()/(pa.dot(Q)*pb.dot(Q)) "
                            "= {dipole_invariant} > 1. something is wrong".format(
                                dipole_invariant=dipole_invariant))
                    elif dipole_invariant > 1.:
                        dipole_invariant = 1.
                    if distribution_type == 'bulk':
                        #The factor xi^2 below corrects the flux factor used in the bulk BS which has a 1/xi^2 too many
                        #A more permanent change is warranted after testing.
                        #See github issue #9 for reference
                        kernel = self.get_bulk_kernel(dipole_invariant,xi)*(xi**2)
                    elif distribution_type == 'counterterm':
                        kernel = self.get_counterterm_kernel(dipole_invariant,xi)
                    elif distribution_type == 'endpoint':
                        kernel = self.get_endpoint_kernel(dipole_invariant,xi)
                    else:
                        raise CurrentImplementationError("Distribution type '%s' not supported."%distribution_type)

                evaluation['values'][(0, color_correlation_index, 0, 0)] = (kernel*prefactor).truncate(max_power=0)
                color_correlation_index += 1

        return evaluation

# Soft-collinear integrated kernels

class QCD_integrated_S_FqFqx_C_IgFqFqx(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state soft-collinear integrated counterterms which we set to zero since
    we will include them together with the soft."""

    # As mentioned in the docstring, this integrated counterterm is zero since it has already been accounted for
    # in the integrated soft CT.
    is_zero = True

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
        )
    )
    structure = sub.SingularStructure(substructures=(coll_structure_g,))

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    current_properties['singular_structure'] = structure
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    current_properties['singular_structure'] = structure
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

class QCD_integrated_S_FqFqx_C_IqFqFqx(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state soft-collinear integrated counterterms which we set to zero since
    we will include them together with the soft."""

    # As mentioned in the docstring, this integrated counterterm is zero since it has already been accounted for
    # in the integrated soft CT.
    is_zero = True

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

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
    coll_structure_q = sub.SingularStructure(substructures=(sub.CollStructure(
        substructures=tuple([soft_structure,]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
        )
    ),))
    coll_structure_qp = sub.SingularStructure(substructures=(sub.CollStructure(
        substructures=tuple([soft_structure,]),
        legs=(
            sub.SubtractionLeg(10, +2, sub.SubtractionLeg.INITIAL),
        )
    ),))

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    # ... q(initial) > C( qx(initial_after_emission) S(q(final) qx(final)) )
    current_properties['singular_structure'] = coll_structure_q
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))
    # ... q(initial) > C( qx(initial_after_emission) S(q'(final) q'x(final)) )
    current_properties['singular_structure'] = coll_structure_qp
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > C( qx(initial_after_emission) S(q(final) qx(final)) )
    current_properties['singular_structure'] = coll_structure_q
    template_currents.append(sub.BeamCurrent( current_properties ))
    # ... q(initial) > C( qx(initial_after_emission) S(q'(final) q'x(final)) )
    current_properties['singular_structure'] = coll_structure_qp
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

class QCD_integrated_S_FqFqx_C_FqFqx(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state soft-collinear integrated counterterms which we set to zero since
    we will include them together with the soft."""

    # As mentioned in the docstring, this integrated counterterm is zero since it has already been accounted for
    # in the integrated soft CT.
    is_zero = True

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

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
    # Match both the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.SoftStructure(
            substructures=(sub_coll_structure,),
            legs=tuple([])
        ),))

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    current_properties['singular_structure'] = structure
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    current_properties['singular_structure'] = structure
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

class QCD_integrated_S_FqFqx_C_FqFqx_C_IqpFqFqx(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state soft-collinear integrated counterterms which we set to zero since
    we will include them together with the soft."""

    # As mentioned in the docstring, this integrated counterterm is zero since it has already been accounted for
    # in the integrated soft CT.
    is_zero = True

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

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
    soft_structure = sub.SoftStructure(
            substructures=(sub_coll_structure,),
            legs=tuple([])
    )
    # This counterterm will be used if any of the structures of the list below matches
    # Match the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(soft_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),))

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    current_properties['singular_structure'] = structure
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    current_properties['singular_structure'] = structure
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

class QCD_integrated_S_FqFqx_C_FqFqx_C_IqFqFqx(general_current.GeneralCurrent):
    """Implements the NLO QCD initial-state soft-collinear integrated counterterms which we set to zero since
    we will include them together with the soft."""

    # As mentioned in the docstring, this integrated counterterm is zero since it has already been accounted for
    # in the integrated soft CT.
    is_zero = True

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.SubtractionCurrentEvaluation

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
    soft_structure = sub.SoftStructure(
            substructures=(sub_coll_structure,),
            legs=tuple([])
    )
    # This counterterm will be used if any of the structures of the list below matches
    # Match the case of the initial state being a quark and/or antiquark
    structure = sub.SingularStructure(substructures=(sub.CollStructure(
            substructures=(soft_structure,),
            legs=(sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),)
        ),))

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    current_properties['singular_structure'] = structure
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    current_properties['singular_structure'] = structure
    template_currents.append(sub.BeamCurrent( current_properties ))

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

# Compound kernels

class QCD_disable_single_block_compound_PDF_counterterms(utils.VirtualCurrentImplementation):
    """ Disable integrated counterterms of the form of currents block with a *single* IntegratedBeamCurrent
    which is composed of two pieces including one PDF counterterms, for example:
       > ([((F((2, 1, i)),),(S((4, 21, f)),),) @ 0 loops]) (type: (IntegratedBeamCurrent) )
    because in colorful_pp, those are already accounted for in these pieces:
       > ([((S((4, 21, f)),),) @ 0 loops], :(F((1, 1, i)),) @ 0 loops:)
    """

    is_zero = True

    @classmethod
    def does_implement_these_currents(cls, currents_block, model):
        """ Overload the mother function so as to easily capture the pattern of all such
        single_block_compound_PDF_counterterms to disable."""

        # A single current block
        if (len(currents_block) == 1 and
            # Composed of a single integrated counterterm of type IntegratedBeamCurrent
            type(currents_block[0]) is sub.IntegratedBeamCurrent and
            # Composed of a top-level singular structure featuring more that one substructure
            len(currents_block[0]['singular_structure'].substructures)>1 and
            # with at least one that *is* a PDF counterterm *and* one that *is not* a PDF counterterm
            any( (len(ss.substructures)>0 and any(sub_ss.name()=='F' for sub_ss in ss.substructures))
                for ss in currents_block[0]['singular_structure'].substructures) and
            any((len(ss.substructures)>0 and any(sub_ss.name() != 'F' for sub_ss in ss.substructures))
                for ss in currents_block[0]['singular_structure'].substructures)
            ):
            # This will disable the counterterm since is_zero is True
            return {}

        return None

class TOTEST_QCD_integrated_Fx_Fx(general_current.GeneralCurrent):
    """Implements the endpoint QCD PDF counterterm for both the left and right beam.
    In this case it can be implemented simply as a direct orthogonal product of the
    flavour endpoint convolution matrices of each of the two PDF counterterms.
    """

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    structure_qA = sub.SingularStructure(substructures=(sub.BeamStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),
        )
    ),))
    structure_gA = sub.SingularStructure(substructures=(sub.BeamStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),
        )
    ),))
    structure_qB = sub.SingularStructure(substructures=(sub.BeamStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.INITIAL),
        )
    ),))
    structure_gB = sub.SingularStructure(substructures=(sub.BeamStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.INITIAL),
        )
    ),))

    # Now define the matching singular structures
    beam_structure_qq = sub.SingularStructure(
        substructures=(structure_qA,structure_qB)
    )
    beam_structure_qg = sub.SingularStructure(
        substructures=(structure_qA,structure_gB),
    )
    beam_structure_gq = sub.SingularStructure(
        substructures=(structure_gA,structure_qB)
    )
    beam_structure_gg = sub.SingularStructure(
        substructures=(structure_gA,structure_gB),
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 4},
        'beam_type' : 'proton',
        'beam_PDGs' : colorful_pp_config.beam_PDGs_supported
    }
    # Now add endpoint IntegratedBeamCurrent...
    current_properties['singular_structure'] = beam_structure_qq
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))
    current_properties['singular_structure'] = beam_structure_qg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))
    current_properties['singular_structure'] = beam_structure_gq
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))
    current_properties['singular_structure'] = beam_structure_gg
    template_currents.append(sub.IntegratedBeamCurrent( current_properties ))

    # The factorisation of the left and right PDF is manifest for the counterterm and bulk
    # pieces and does therefore not need to be specified here

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents, ]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        xis = global_variables['xis']
        mu_fs = global_variables['mu_fs']
        mu_r = global_variables['mu_r']

        allowed_backward_evolved_flavors = global_variables['allowed_backward_evolved_flavors']
        if allowed_backward_evolved_flavors != ('ALL','ALL'):
            raise CurrentImplementationError('The current %s must always be called with' % self.__class__.__name__ +
                                             "allowed_backward_evolved_flavors=('ALL','ALL'), not %s" % str(
                allowed_backward_evolved_flavors))

        NF = self.currents_properties[0]['NF']

        # Only the order epsilon of the scales pre-factor matters here.
        prefactorA = EpsilonExpansion({
            0: 1.,
            1: log(mu_r ** 2 / mu_fs[0] ** 2)
        })
        prefactorB = EpsilonExpansion({
            0: 1.,
            1: log(mu_r ** 2 / mu_fs[1] ** 2)
        })
        prefactor = prefactorA * prefactorB * EpsilonExpansion({-2: 1.}) * (self.SEpsilon * (1. / (16 * pi**2) ))**2

        # Define the NLO QCD PDF counterterms kernels
        kernel_gg = prefactor * (11. / 6. * self.CA - 2. / 3. * NF * self.TR)**2

        # Now assign the flavor matrix to the result
        evaluation['values'][(0,0,0,0)] = {
            (21, 21) : { ( (21, 21 ), )  : kernel_gg.truncate(max_power=0)}
        }

        return evaluation

class QCD_integrated_CS_X_F(general_current.GeneralCurrent):
    """Implements the NLO QCD soft-collinear integrated counterterms together with a PDF counterterm.
    These are all set to zero as they are included in the soft counterterms
    """

    # As mentioned before these are all zero as they are already included in the soft counterterm.
    is_zero = True

    # Enable the flag below to debug this current
    DEBUG = False

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
    soft_coll_structure_gg_final = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=(soft_structure_g,),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),)
        ),)
    )
    soft_coll_structure_qg_final = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=(soft_structure_g,),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),)
        ),)
    )
    soft_coll_structure_gg_initial = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=(soft_structure_g,),
            legs=(
                sub.SubtractionLeg(10, 21, sub.SubtractionLeg.INITIAL),)
        ),)
    )
    soft_coll_structure_qg_initial = sub.SingularStructure(
        substructures=(sub.CollStructure(
            substructures=(soft_structure_g,),
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.INITIAL),)
        ),)
    )

    # This counterterm will be used if any of the current of the list below matches
    template_currents_soft_collinear = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }
    # Now add endpoint IntegratedBeamCurrent...
    # ... g(initial) > g(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg_initial
    template_currents_soft_collinear.append(sub.IntegratedBeamCurrent( current_properties ))
    current_properties['singular_structure'] = soft_coll_structure_gg_final
    template_currents_soft_collinear.append(sub.IntegratedBeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg_initial
    template_currents_soft_collinear.append(sub.IntegratedBeamCurrent( current_properties ))
    current_properties['singular_structure'] = soft_coll_structure_qg_final
    template_currents_soft_collinear.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_gg_initial
    template_currents_soft_collinear.append(sub.BeamCurrent( current_properties ))
    current_properties['singular_structure'] = soft_coll_structure_gg_final
    template_currents_soft_collinear.append(sub.BeamCurrent( current_properties ))
    # ... q(initial) > q(initial_after_emission) g(final)
    current_properties['singular_structure'] = soft_coll_structure_qg_initial
    template_currents_soft_collinear.append(sub.BeamCurrent( current_properties ))
    current_properties['singular_structure'] = soft_coll_structure_qg_final
    template_currents_soft_collinear.append(sub.BeamCurrent( current_properties ))

    # Now add the PDF counterterm part of this compound current
    # The second piece of the compound current which is the PDF counterterm piece will be used if any of the current
    # of the list below matches
    template_currents_PDF = []

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

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : 'proton',
        'beam_PDGs' : colorful_pp_config.beam_PDGs_supported
    }

    # This compound current will always contain only the bulk PDF counterterm
    current_properties['distribution_type'] = ['bulk', ]
    # ... for an initial-state quark
    current_properties['singular_structure'] = beam_structure_q
    template_currents_PDF.append(sub.BeamCurrent( current_properties ))
    # ... for an initial-state gluon
    current_properties['singular_structure'] = beam_structure_g
    template_currents_PDF.append(sub.BeamCurrent( current_properties ))

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents_soft_collinear, template_currents_PDF]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]

class TODO_QCD_integrated_S_Fg_F(general_current.GeneralCurrent):
    """Implements the NLO QCD soft integrated counterterms together with a PDF counterterm.
    Notice that this soft counterterm also includes the integrated soft-collinear.
    """

    # Enable the flag below to debug this current
    DEBUG = False

    # Store the result for beam factorisation currents in a container that supports flavor matrices.
    subtraction_current_evaluation_class = utils.BeamFactorizationCurrentEvaluation

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
    structure = sub.SingularStructure(substructures=(soft_structure,))

    # This first piece of the compound current will be used if any of the current of the list below matches
    template_currents_soft = []

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : None, # These integrated counterterm do no refer to particular beam type or PDGs.
        'beam_PDGs' : None
    }

    # Now add endpoint IntegratedBeamCurrent...
    current_properties['singular_structure'] = structure
    template_currents_soft.append(sub.IntegratedBeamCurrent( current_properties ))

    # Now add plus distributions as BeamCurrent...
    # Specifying the distribution type as both 'bulk' and 'counterterm' means that both values
    # are acceptable when matching a target current_blocks element to this template current.
    current_properties['distribution_type'] = ['bulk', 'counterterm']
    current_properties['singular_structure'] = structure
    template_currents_soft.append(sub.BeamCurrent( current_properties ))

    # Now add the PDF counterterm part of this compound current
    # The second piece of the compound current which is the PDF counterterm piece will be used if any of the current
    # of the list below matches
    template_currents_PDF = []

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

    current_properties = {
        'resolve_mother_spin_and_color': True,
        'n_loops': 0,
        'squared_orders': {'QCD': 2},
        'beam_type' : 'proton',
        'beam_PDGs' : colorful_pp_config.beam_PDGs_supported
    }

    # This compound current will always contain only the bulk PDF counterterm
    current_properties['distribution_type'] = ['bulk', ]
    # ... for an initial-state quark
    current_properties['singular_structure'] = beam_structure_q
    template_currents_PDF.append(sub.BeamCurrent( current_properties ))
    # ... for an initial-state gluon
    current_properties['singular_structure'] = beam_structure_g
    template_currents_PDF.append(sub.BeamCurrent( current_properties ))

    # The defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [template_currents_soft, template_currents_PDF]

    # An now the mapping rules, which are not necessary in this context.
    mapping_rules = [ ]


    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        # Retrieve variables
        beam_number = self.currents_properties[1]['leg_numbers_map'][10]

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

        # Extract the distribution type from the soft template piece of this compound current
        distribution_type = self.currents_properties[0]['distribution_type']
        # Then extract the beam PDGs and NF quantity from the PDF counterterm piece of this compound current
        NF = self.currents_properties[1]['NF']
        beam_PDGs = self.currents_properties[1]['beam_PDGs']

        # Apply the Dirac delta on the convolution parameters by setting their value to the solution of the delta.
        # For these simple NLO integrated current, the delta simply sets them to one if the distribution type is
        # 'counterterm' and leaves them unchanged if the distribution type is 'bulk' (or 'endpoint' in which case
        # the convolution parameter is None anyway.)
        if distribution_type == 'counterterm':
            xis = list(global_variables['xis'])
            #TODO This compound current requires a change of variable which should be specified here!
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
                1.0 #TODO
            ),
            'counterterm': prefactor * (
                1.0 #TODO
            ),
            'endpoint': prefactor * (
                1.0 #TODO
            )
        }

        kernel_gq = {
            'bulk': prefactor * (
                1.0  # TODO
            ),
            'counterterm': None,
            'endpoint': None
        }

        kernel_qg = {
            'bulk': prefactor * (
                1.0 # TODO
            ),
            'counterterm': None,
            'endpoint': None
        }

        kernel_qq = {
            'bulk': prefactor * (
                1.0  # TODO
            ),
            'counterterm': prefactor * (
                1.0  # TODO
            ),
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