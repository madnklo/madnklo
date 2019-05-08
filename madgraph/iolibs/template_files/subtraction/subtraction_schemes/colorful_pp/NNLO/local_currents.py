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

import commons.utils as utils
from commons.universal_kernels import AltarelliParisiKernels, SoftKernels

import commons.QCD_local_currents as currents
import colorful_pp_config

from madgraph.core.base_objects import EpsilonExpansion

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================================
# Variable functions for accessing initial and final (soft)-collinear variables with s_(a,r,s) invariants
#=========================================================================================================

def colorful_pp_initial_coll_variables(PS_point, parent_momentum, children, **opts):
    """ Initial state collinear variables as specified in Gabors Colorful ISR NNLO note."""
    
    p_a = PS_point[children[0]]
    p_fs = [PS_point[child] for child in children[1:]]
    Q = opts['Q']
    pa_dot_pb = Q.square()/2.

    # Loop over all final state momenta to obtain the corresponding variables
    xs = []
    kTs = []
    ss = {}
    for i_fs in range(len(p_fs)):
        p_r = p_fs[i_fs]
        p_s = sum(p_fs[i_fs_prime] for i_fs_prime in range(len(p_fs)) if i_fs_prime!=i_fs)
        xs.append(
            (p_r.dot(Q)/pa_dot_pb) - (p_r.dot(p_s)/pa_dot_pb)*(p_a.dot(p_r)/p_a.dot(p_r+p_s))
        )
        kTs.append(
            p_r - ((p_r.dot(Q)-2.*p_a.dot(p_r))/pa_dot_pb)*p_a - (p_a.dot(p_r)/pa_dot_pb)*Q
        )
        ss[(0,i_fs+1)] = 2.*p_a.dot(p_r)

    # Finally add the x and kT corresponding to the initial state leg
    #xs.insert(0, 1.0-sum(xs))
    xs.insert(0, 1.0 - (((p_fs[0]+p_fs[1]).dot(Q)-p_fs[0].dot(p_fs[1]))/pa_dot_pb) )
    kTs.insert(0, -sum(kTs))

    # Alternative way of getting these variables using Simone's mappings function
#    na = parent_momentum
#    nb = opts['Q']
#    kin_variables = dict()
    # The lone initial state child is always placed first thanks to the implementation
    # of the function get_sorted_children() in the current.
#    mappings.InitialCollinearVariables.get(
#        PS_point, children[1:], children[0], na, nb, kin_variables)
#    xs = tuple(kin_variables['z%d' % i] for i in children)
#    kTs = tuple(kin_variables['kt%d' % i] for i in children)

    # Add additional ss's
    for i_fs in range(len(p_fs)):
        for j_fs in range(i_fs+1,len(p_fs)):
            ss[(i_fs+1,j_fs+1)] = 2.*p_fs[i_fs].dot(p_fs[j_fs])

    return tuple(xs), tuple(kTs), ss


def colorful_pp_initial_soft_coll_variables(PS_point, parent_momentum, children, **opts):
    """ Initial state soft-collinear variables as specified in Gabors Colorful ISR NNLO note."""

    p_a_tilde = parent_momentum
    p_fs = [PS_point[child] for child in children[1:]]
    Q = opts['Q']

    # Loop over all final state momenta to obtain the corresponding variables
    xs = []
    kTs = []
    ss = {}
    for i_fs in range(len(p_fs)):
        p_r = p_fs[i_fs]
        xs.append( p_r.dot(Q) / p_a_tilde.dot(Q) )
        # KTs not needed for now, will be filled in if needed.
        kTs.append(None)
        ss[(0, i_fs + 1)] = 2. * p_a_tilde.dot(p_r)

    # Finally add the x and kT corresponding to the initial state leg
    # xs.insert(0, 1.0-sum(xs))
    xs.insert(0, 1.0-sum(xs))
    #kTs.insert(0, -sum(kTs))
    kTs.insert(0, None)

    # Alternative way of getting these variables using Simone's mappings function
    #    na = parent_momentum
    #    nb = opts['Q']
    #    kin_variables = dict()
    # The lone initial state child is always placed first thanks to the implementation
    # of the function get_sorted_children() in the current.
    #    mappings.InitialCollinearVariables.get(
    #        PS_point, children[1:], children[0], na, nb, kin_variables)
    #    xs = tuple(kin_variables['z%d' % i] for i in children)
    #    kTs = tuple(kin_variables['kt%d' % i] for i in children)

    # Add additional ss's
    for i_fs in range(len(p_fs)):
        for j_fs in range(i_fs + 1, len(p_fs)):
            ss[(i_fs + 1, j_fs + 1)] = 2. * p_fs[i_fs].dot(p_fs[j_fs])

    return tuple(xs), tuple(kTs), ss

def colorful_pp_final_coll_variables(PS_point, parent_momentum, children, **opts):
    """ For now this is identical as Q_final_coll_variables"""

    na = parent_momentum
    nb = opts['Q']
    kin_variables = dict()
    mappings.FinalCollinearVariables.get(
        PS_point, children, na, nb, kin_variables)
    zs  = tuple(kin_variables['z%d'  % i] for i in children)
    kTs = tuple(kin_variables['kt%d' % i] for i in children)

    p_fs = [PS_point[child] for child in children]
    ss = {}
    # Add additional ss's
    for i_fs in range(len(p_fs)):
        for j_fs in range(i_fs+1,len(p_fs)):
            ss[(i_fs+1,j_fs+1)] = 2.*p_fs[i_fs].dot(p_fs[j_fs])

    return zs, kTs, ss

#=========================================================================================
# NNLO initial-collinear currents
#=========================================================================================

class QCD_initial_collinear_0_XXX(currents.QCDLocalCollinearCurrent):
    """Triple collinear initial-final-final tree-level current."""

    squared_orders = {'QCD': 4}
    n_loops = 0

    is_cut = staticmethod(colorful_pp_config.cut_initial_coll)
    factor = staticmethod(colorful_pp_config.factor_initial_coll)
    get_recoilers = staticmethod(colorful_pp_config.get_recoilers)
    mapping = colorful_pp_config.initial_coll_mapping
    variables = staticmethod(colorful_pp_initial_coll_variables)
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
    ]

    def kernel(self, evaluation, parent, xs, kTs, ss):

        # Retrieve the collinear variable x
        x_a, x_r, x_s = xs
        kT_a, kT_r, kT_s = kTs
        s_ar, s_as, s_rs = ss[(0,1)], ss[(0,2)], ss[(1,2)]

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorization formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= 1.

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqpqp(self,
                1./x_a, -x_r/x_a, -x_s/x_a, -s_ar, -s_as, s_rs, kT_a, kT_r, kT_s
            ):
            complete_weight = weight*initial_state_crossing_factor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

#=========================================================================================
# NNLO final-collinear currents
#=========================================================================================

class QCD_final_collinear_0_XXX(currents.QCDLocalCollinearCurrent):
    """Triple collinear initial-final-final tree-level current."""

    squared_orders = {'QCD': 4}
    n_loops = 0

    is_cut = staticmethod(colorful_pp_config.cut_final_coll)
    factor = staticmethod(colorful_pp_config.factor_final_coll)
    get_recoilers = staticmethod(colorful_pp_config.get_recoilers)
    mapping = colorful_pp_config.final_coll_mapping
    variables = staticmethod(colorful_pp_final_coll_variables)
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

class QCD_final_collinear_0_qqpqp(QCD_final_collinear_0_XXX):
    """qg collinear FSR tree-level current. q(final) > q(final) q'(final) qbar' (final) """

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(1, +2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, -2, sub.SubtractionLeg.FINAL),
        )),
    ]

    def kernel(self, evaluation, parent, zs, kTs, ss):

        # Retrieve the collinear variable x
        z_i, z_r, z_s = zs
        kT_i, kT_r, kT_s = kTs
        s_ir, s_is, s_rs = ss[(1,2)], ss[(1,3)], ss[(2,3)]

        for spin_correlation_vector, weight in AltarelliParisiKernels.P_qqpqp(self,
                z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s
            ):
            complete_weight = weight
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append((parent, spin_correlation_vector))
                evaluation['values'][(len(evaluation['spin_correlations'])-1, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

# =========================================================================================
# NNLO soft currents
# =========================================================================================
class QCD_soft_0_kX(currents.QCDLocalCurrent):
    """ NNLO soft currents."""

    is_cut = staticmethod(colorful_pp_config.cut_soft)
    factor = staticmethod(colorful_pp_config.factor_soft)
    get_recoilers = staticmethod(colorful_pp_config.get_recoilers)

    squared_orders = {'QCD': 4}
    n_loops = 0

    mapping = colorful_pp_config.soft_mapping
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    # Must be specified by daughter classes
    soft_kernel = None

    def evaluate_subtraction_current(
            self, current,
            higher_PS_point=None, momenta_dict=None, reduced_process=None,
            hel_config=None, Q=None, **opts):

        if higher_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping.")
        if momenta_dict is None:
            raise CurrentImplementationError(
                self.name() + " requires a momentum routing dictionary.")
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment.")
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q.")

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Now find all colored leg numbers in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color') == 1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))
        soft_momenta = [ higher_PS_point[self.leg_numbers_map[soft_leg_number]] for soft_leg_number in
                                                                             self.leg_numbers_map if soft_leg_number>9 ]
        pS = sum(soft_momenta)

        # Perform mapping
        this_mapping_singular_structure = self.mapping_singular_structure.get_copy()
        this_mapping_singular_structure.legs = self.get_recoilers(reduced_process)
        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, this_mapping_singular_structure, momenta_dict,
            compute_jacobian=self.divide_by_jacobian)
        reduced_kinematics = (None, lower_PS_point)

        # Include the counterterm only in a part of the phase space
        if self.is_cut(Q=Q, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config, reduced_kinematics=('IS_CUT', lower_PS_point))

        jacobian = mapping_vars.get('jacobian', 1.)

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [],
            'reduced_kinematics': [reduced_kinematics, ],
            'values': {}
        })

        # Normalization factors
        norm = (8. * math.pi * alpha_s)**(len(soft_momenta))*(1./pS.square()**2)
        norm *= self.factor(Q=Q, pS=pS)
        if self.divide_by_jacobian:
            norm /= jacobian

        colored_partons_momenta = vectors.LorentzVectorDict()
        for colored_parton_number in all_colored_parton_numbers:
            # We want to used the reduced kinematics for our soft current
            colored_partons_momenta[colored_parton_number] = lower_PS_point[colored_parton_number]
            # Alternatively, the expression below would have given us the resolved one
            #colored_partons_momenta[colored_parton_number] = sum(higher_PS_point[child] for child in momenta_dict[colored_parton_number])

        color_correlation_index = 0
        for color_correlator, weight in self.soft_kernel(
                self, colored_partons_momenta, soft_momenta, all_colored_parton_numbers):
            evaluation['color_correlations'].append(color_correlator)
            complete_weight = weight*norm
            evaluation['values'][(0, color_correlation_index, 0)] = {'finite': complete_weight[0]}
            color_correlation_index += 1

        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

class QCD_soft_0_qqp(QCD_soft_0_kX):
    """Soft quark-anti quark pair.
    """

    structure = sub.SingularStructure(substructures=(sub.SoftStructure(
            legs=(
                sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
            )
        ),)
    )

    soft_kernel = staticmethod(SoftKernels.qqx)

# =========================================================================================
# NNLO soft-collinear currents
# =========================================================================================

class QCD_initial_soft_collinear_0_kX(currents.QCDLocalCurrent):
    """ triple soft-collinear with initial state """

    squared_orders = {'QCD': 4}
    n_loops = 0

    is_cut = staticmethod(colorful_pp_config.cut_soft)
    factor = staticmethod(colorful_pp_config.factor_soft)
    get_recoilers = staticmethod(colorful_pp_config.get_recoilers)
    mapping = colorful_pp_config.initial_soft_collinear_mapping
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    variables = staticmethod(colorful_pp_initial_soft_coll_variables)

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, Q=None, **opts ):

        if higher_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before mapping." )
        if momenta_dict is None:
            raise CurrentImplementationError(
                self.name() + " requires a momentum routing dictionary." )
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )


        """Important note about the IF CS:
        - in this scheme we want to use "tilded" momenta for the dipole legs in eikonals. This is explicitly implemented in the soft local current
        - this implies that the correct form for the local C(ir)S(r) taken as the collinear limit of the eikonals is 
        1/ (p_r + p_i_tilde)^2 (1-z_r)/z_r where z_r = p_r.Q/(p_r+p_i_tilde).Q
        - Specializing to the case where the collinear partner of the soft particle is an initial state particle (i = a ), we have 
        p_a_tilde = xi p_a and 2p_a.Q = Q^2 so that the soft-collinear takes the form
        1/(p_r+xi p_a)^2 * xi/y_rQ where y_rQ is the usual Hungarian variable
        this simplifies to 
        1/(p_r+p_a)^2 * 1/y_rQ which is exactly the soft collinear as computed *without* tilded variables (i.e. exactly eq.5.29 of 0903.1218)

        As a result we use exactly the same way of evaluating the counterterms as in honest-to-god colorful.
        """

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        children = tuple(self.leg_numbers_map[i]
                         for i in sorted(self.leg_numbers_map.keys()))
        parent = momenta_dict.inv[frozenset(children)]

        # Perform mapping
        this_mapping_singular_structure = self.mapping_singular_structure.get_copy()
        this_mapping_singular_structure.legs = self.get_recoilers(reduced_process, excluded=(parent, ))
        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, this_mapping_singular_structure, momenta_dict,
            compute_jacobian=self.divide_by_jacobian )
        reduced_kinematics = (None, lower_PS_point)

        # Include the counterterm only in a part of the phase space
        # children are the the set of particles that are going unresolved.
        # Here we have C and S going collinear with S soft.
        # The parent is the mapped C with a soft mapping, usually refered to as Ctilde.
        # S is just removed in a soft mapping.
        # Here S is a single particle but we obtain it as a list soft_children\
        # to illustrate how multiple softs would be obtained
        pCtilde = lower_PS_point[parent]
        soft_children = [ self.leg_numbers_map[soft_leg_number] for soft_leg_number in self.leg_numbers_map if soft_leg_number>9 ]
        pS = sum(higher_PS_point[child] for child in soft_children)
        collinear_final_children = [ self.leg_numbers_map[soft_leg_number] for soft_leg_number in self.leg_numbers_map if
                                                                                                0 < soft_leg_number <= 9 ]
        if len(collinear_final_children)>0:
            pCfinal = sum(higher_PS_point[child] for child in collinear_final_children)
        else:
            pCfinal = vectors.LorentzVector()
        pCinitial = higher_PS_point[self.leg_numbers_map[0]]
        pCmother = pCinitial - pCfinal - pS
        if self.is_cut(Q=Q, pC=pCmother, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config, reduced_kinematics=('IS_CUT', lower_PS_point))

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [None],
            'reduced_kinematics': [reduced_kinematics],
            'values': { }
        })

        # Evaluate kernel
        kernel_arguments = self.variables(higher_PS_point, pCtilde, children, Q=Q)

        # There is no need for the ratio of color-averaging factor between the real ME
        # initial state flavor and the one of the reduced Born ME as they are either both
        # gluons or both quarks
        evaluation = self.kernel(evaluation, parent, *kernel_arguments)

        # Add the normalization factors
        norm = (8. * math.pi * alpha_s)**(len(soft_children)+len(collinear_final_children)) / ((pCtilde+pS).square()*pS.square())
        norm *= self.factor(Q=Q, pC=pCmother, pS=pS)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result


class QCD_initial_soft_collinear_0_qqpqp(QCD_initial_soft_collinear_0_kX):
    """q' qbar' soft collinear ISR tree-level current. q(initial) > q(initial_after_emission) Soft(q'(final) qbar' (final)) """


    soft_structure = sub.SoftStructure(
        legs=(
            sub.SubtractionLeg(10, +1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(0,  1, sub.SubtractionLeg.INITIAL), ) )
    soft_coll_structure_g = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.INITIAL), ) )
    structure_q = sub.SingularStructure(substructures=(soft_coll_structure_q, ))
    structure_g = sub.SingularStructure(substructures=(soft_coll_structure_g, ))

    structure = [structure_q, structure_g]

    expected_init_opts = tuple(list(QCD_initial_soft_collinear_0_kX.expected_init_opts)+['color_charge',])

    def __init__(self, *args, **opts):
        """ Specialise constructor so as to assign value to the specified color charge. """
        super(QCD_initial_soft_collinear_0_qqpqp, self).__init__(*args, **opts)
        # Turn 'CA' into the actual numerical value self.CA for instance
        self.color_charge = getattr(self, self.color_charge)

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Overwrite this function so as to be able to specify the color-factor to consider."""

        res = super(QCD_initial_soft_collinear_0_qqpqp, cls).does_implement_this_current(current, model)
        if res is None:
            return None

        color_charge = 'CF'
        leg_numbers_map = cls.structure_q.map_leg_numbers(
            current.get('singular_structure'), [range(1, model.get_nflav()+1)])
        if leg_numbers_map is None:
            color_charge = 'CA'
            leg_numbers_map = cls.structure_g.map_leg_numbers(
                current.get('singular_structure'), [range(1, model.get_nflav()+1)])
            if leg_numbers_map is None:
                return None
        res['color_charge'] = color_charge
        return res

    def kernel(self, evaluation, parent, xs, kTs, ss):

        # Retrieve the collinear variable x
        x_a, x_r, x_s = xs
        kT_a, kT_r, kT_s = kTs
        s_ar, s_as, s_rs = ss[(0,1)], ss[(0,2)], ss[(1,2)]
        s_a_rs = s_ar + s_as

        kernel = 2. * ( 1. / (x_r + x_s) - (((s_ar * x_s - s_as * x_r) ** 2) / (s_a_rs * s_rs * ((x_r + x_s) ** 2))) )

        evaluation['values'][(0, 0, 0)] = {'finite': self.color_charge * self.TR * kernel }
        return evaluation