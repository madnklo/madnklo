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

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings

import commons.utils as utils
import commons.QCD_local_currents as currents
import colorful_pp_config

from madgraph.core.base_objects import EpsilonExpansion

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# Static function for accessing collinear variables with s_ij invariants
#=========================================================================================
def Q_initial_coll_variables_with_ss(PS_point, parent_momentum, children, **opts):
    na = parent_momentum
    nb = opts['Q']
    kin_variables = dict()

    # The lone initial state child is always placed first thanks to the implementation
    # of the function get_sorted_children() in the current.
    mappings.InitialCollinearVariables.get(
        PS_point, children[1:], children[0], na, nb, kin_variables)
    zs = tuple(kin_variables['z%d' % i] for i in children)
    kTs = tuple(kin_variables['kt%d' % i] for i in children)
    # This is suited for massless particles. I prefer to use explicitly the squares
    # as it would be numerically less stable
    ss = (
        2.0*PS_point[children[0]].dot(PS_point[children[1]]),
        2.0*PS_point[children[0]].dot(PS_point[children[2]]),
        2.0*PS_point[children[1]].dot(PS_point[children[2]]),
    )
    return zs, kTs, ss

#=========================================================================================
# Class listing Altarelli-Parisi Kernels
#=========================================================================================

class AltarelliParisiKernels:
    """ Implementation of AltarelliParisiKernels. Notice that the the first argument color_factors is always
    any object with the float attributes CF, CA, etc... to access color factors."""

    @staticmethod
    def P_qqpqp(color_factors, z_i, z_r, z_s, kT_i, kT_r, kT_s, s_ir, s_is, s_rs):
        """ Kernel for the q -> q qp' qp' splitting. The return value is not a float but a list of tuples:
                ( spin_correlation_vectors_with_parent, weight )
            where spin_correlation_vector_with_parent can be None if None is required.
        """

        # Compute handy term and variables
        s_irs = s_ir+s_is+s_rs
        t_rs_i = 2.*(z_r*s_is-z_s*s_ir)/(z_r+z_s) + ((z_r-z_s)/(z_r+z_s))*s_rs
        dimensional_term = z_r + z_s - s_rs/s_irs

        # Overall prefactor
        prefactor = (1./2.)*color_factors.CF*color_factors.TR*(s_irs/s_rs)

        return [
            ( None,
                EpsilonExpansion({
                    0 : -(t_rs_i**2/(s_rs*s_irs)) + (4.*z_i + (z_r - z_s)**2)/(z_r+z_s) + dimensional_term,
                    1 : -2.*dimensional_term
                })*prefactor
            )
        ]

#=========================================================================================
# NLO initial-collinear currents
#=========================================================================================

class QCD_initial_collinear_0_XXX(currents.QCDLocalCollinearCurrent):
    """Triple collinear initial-final-final tree-level current."""

    squared_orders = {'QCD': 4}
    n_loops = 0

    is_cut = staticmethod(colorful_pp_config.cut_initial_coll)
    factor = staticmethod(colorful_pp_config.factor_initial_coll)
    get_recoilers = staticmethod(colorful_pp_config.get_recoilers)
    mapping = colorful_pp_config.initial_coll_mapping
    variables = staticmethod(Q_initial_coll_variables_with_ss)
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

class QCD_initial_collinear_0_qqpqp(QCD_initial_collinear_0_XXX):
    """qg collinear ISR tree-level current. q(initial) > q(initial_after_emission) q'(final) qbar' (final) """

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, -2, sub.SubtractionLeg.FINAL),
        )),
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, -1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, -2, sub.SubtractionLeg.FINAL),
        )),
    ]

    def kernel(self, evaluation, parent, xs, ss, kTs):

        # Retrieve the collinear variable x
        x_a, x_r, x_s = xs
        s_ar, s_as, s_rs = ss
        kT_a, kT_r, kT_s = kTs

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorization formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = -1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= 1.

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqpqp(self,
                1./x_a, -x_r/x_a, -x_s/x_a, s_ar, s_as, s_rs, kT_a, kT_r, kT_s
            ):
            complete_weight = (weight*initial_state_crossing_factor).to_human_readable_dict()
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = complete_weight
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0)] = complete_weight

        return evaluation