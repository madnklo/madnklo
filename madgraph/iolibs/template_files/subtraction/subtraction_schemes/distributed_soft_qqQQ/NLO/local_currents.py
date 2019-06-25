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
"""Implementation of NLO distributed soft currents."""

import math

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings

# from commons.universal_kernels import AltarelliParisiKernels, SoftKernels #TODO DEV
from madgraph.iolibs.template_files.subtraction.commons.universal_kernels import AltarelliParisiKernels, SoftKernels #TODO DEV
# import commons.utils as utils TODO DEV
import madgraph.iolibs.template_files.subtraction.commons.utils as utils #TODO DEV
# import commons.QCD_local_currents as currents TODO DEV
import madgraph.iolibs.template_files.subtraction.commons.QCD_local_currents as currents #TODO DEV
# import commons.factors_and_cuts as factors_and_cuts

from bidict import bidict

import madgraph.iolibs.template_files.subtraction.subtraction_schemes.distributed_soft_qqQQ.distributed_soft_qqQQ_config as scheme_config # TODO DEV
# import distributed_soft_qqQQ_config as scheme_config TODO DEV

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# Choice of recoilers
#=========================================================================================

#=========================================================================================
# NLO final-collinear currents
#=========================================================================================


class QCD_final_collinear_0_XX_soft_distrib(currents.GeneralQCDLocalCurrent):
    """Two-collinear tree-level current."""

    squared_orders = {'QCD': 2}
    n_loops = 0

    is_cut = scheme_config.final_coll_cut
    factor = scheme_config.final_coll_factor
    mapping = scheme_config.final_coll_mapping
    # variables = scheme_config.final_coll_variable
    divide_by_jacobian = scheme_config.divide_by_jacobian
    get_recoilers = scheme_config.get_recoilers

    mapping_rules = [
        {
            'singular_structure': None,
            'mapping': mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict': bidict({-1: frozenset((0, 1))}),
            'variables': currents.CompoundVariables(scheme_config.final_coll_variables),
            'is_cut': is_cut,
            'reduced_recoilers': get_recoilers,
            'additional_recoilers': sub.SubtractionLegSet([]),
        }
    ]


class QCD_final_collinear_0_gq_soft_distrib(QCD_final_collinear_0_XX_soft_distrib):
    """g q collinear tree-level current.
    This contains
    - the hard-collinear kernel for g(0)q(1) being collinear
    - the pieces of the eikonal terms for g(0) soft emitted from q(1)X(j) which also diverge
    in the g(0)q(1) collinear limit
    """

    structure = [sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL), ))]

    mapping_rules = QCD_final_collinear_0_XX_soft_distrib.mapping_rules
    mapping_rules[0]['singular_structure']=structure[0]


    def kernel(self, evaluation, all_steps_info, global_variables):

        # Retrieve the collinear variables
        z = all_steps_info[0]['variables'][0]['zs'][0]
        s = all_steps_info[0]['variables'][0]['ss'][0]
        # Compute the kernel using
        # f9d0839fc58905d67367e3e67efabee05ee390f9:madgraph/iolibs/template_files/OLD_subtraction/cataniseymour/NLO/local_currents.py:146
        evaluation['values'][(0, 0, 0)] = {'finite':self.CF*z/s}
        return evaluation


    def soft_kernel(self, evaluation, colored_parton_numbers, all_steps_info, global_variables):
        # Retrieve the gluon momentum pg
        # there is only one step and one bundle
        gluon_ID = all_steps_info[0]['bundles_info'][0]['final_state_children'][0]
        pg = all_steps_info[0]['higher_PS_point'][gluon_ID]

        # Retrieve the quark momentum pq
        # there is only one step and one bundle
        quark_ID = all_steps_info[0]['bundles_info'][0]['final_state_children'][1]
        pq = all_steps_info[0]['higher_PS_point'][quark_ID]

        parent_ID = all_steps_info[0]['bundles_info'][0]['parent']
        # Loop over the dipoles (quark,j). Only the gluon can be soft
        color_correlation_index = 1
        for j in colored_parton_numbers:
            if j == parent_ID:
                # The emission and reabsorption by the parent quark does not contribute in the massless quark
                continue
            else:
                pj=all_steps_info[0]['lower_PS_point'][j]
                # We add the color correlation (q,j) to the list of color correlations
                evaluation['color_correlations'].append(((parent_ID, j),))
                # color_correlation_index now points to that color correlation in the list
                evaluation['values'][(0, color_correlation_index, 0)] = {'finite': SoftKernels.partial_fractionned_eikonal_g(pg,pq,pj)[0]}
                # color_correlation_index now points to the next possible color correlation
                color_correlation_index += 1

        return evaluation


class QCD_final_collinear_0_qqx(QCD_final_collinear_0_XX_soft_distrib):
    """q q~ collinear tree-level current."""
    structure = [sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL), ))]
    mapping_rules = QCD_final_collinear_0_XX_soft_distrib.mapping_rules
    mapping_rules[0]['singular_structure']=structure[0]

    def kernel(self, evaluation, all_steps_info, global_variables):
        # Retrieve the collinear variables
        z = all_steps_info[0]['variables'][0]['zs'][0]
        s = all_steps_info[0]['variables'][0]['ss'][0]
        kT = all_steps_info[0]['variables'][0]['kT'][0]
        # Compute the kernel using
        # f9d0839fc58905d67367e3e67efabee05ee390f9:madgraph/iolibs/template_files/OLD_subtraction/cataniseymour/NLO/local_currents.py:146
        evaluation['values'][(0, 0, 0)] = {'finite':self.TR/s}
        evaluation['values'][(1, 0, 0)] = {'finite':4 * self.TR * z * (1-z) / kT.square()/s}
        return evaluation


class QCD_final_collinear_0_gg(QCD_final_collinear_with_soft):
    """g g collinear tree-level current."""

    structure = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), ))

    def kernel(self, evaluation, all_steps_info, global_variables):
        # Retrieve the collinear variables
        z = all_steps_info[0]['variables'][0]['zs'][0]
        s = all_steps_info[0]['variables'][0]['ss'][0]
        kT = all_steps_info[0]['variables'][0]['kT'][0]
        # Compute the kernel
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        evaluation['values'][(0, 0, 0)] = {'finite:': 0.}
        evaluation['values'][(1, 0, 0)] = { 'finite':
            -2. * self.CA * 2. * z * (1. - z) / kT.square()/s}
        return evaluation

    def soft_kernel(self, evaluation, colored_parton_numbers, all_steps_info, global_variables):
        # Retrieve the first gluon momentum pg
        # there is only one step and one bundle
        g1_ID = all_steps_info[0]['bundles_info'][0]['final_state_children'][0]
        p1 = all_steps_info[0]['higher_PS_point'][gluon_ID]

        # Retrieve the seconde gluon momentum pg
        # there is only one step and one bundle
        g2_ID = all_steps_info[0]['bundles_info'][0]['final_state_children'][1]
        p2 = all_steps_info[0]['higher_PS_point'][quark_ID]

        parent_ID = all_steps_info[0]['bundles_info'][0]['parent']
        # Loop over the dipoles (quark,j). Only the gluon can be soft
        color_correlation_index = 1
        for j in colored_parton_numbers:
            if j == parent_ID:
                # The emission and reabsorption by the parent quark does not contribute in the massless quark
                continue
            else:
                pj=all_steps_info[0]['lower_PS_point'][j]
                # We add the color correlation (q,j) to the list of color correlations
                evaluation['color_correlations'].append(((parent_ID, j),))
                # color_correlation_index now points to that color correlation in the list
                evaluation['values'][(0, color_correlation_index, 0)] = {'finite': # Sum over p1/p2 soft from (2/1,j)
                                                        SoftKernels.partial_fractionned_eikonal_g(p1,p2,pj)[0]
                                                      + SoftKernels.partial_fractionned_eikonal_g(p2,p1,pj)[0]
                                                                         }
                # color_correlation_index now points to the next possible color correlation
                color_correlation_index += 1

        return evaluation

##TODO OLD
# class QCD_final_collinear_with_soft(QCD_final_collinear_0_XX):
#     """Generic class subtracting all the NLO counterterms
#      that diverge when a pair of particles become unresolved"""
#
#     @staticmethod
#     def partial_fractionned_eikonal(pi, pj, pr):
#         """Partial fractionned eikonal factor for soft particle with momentum pr
#         emitted from the dipole with momenta pi and pj.
#         This is obtained starting from the eikonal and:
#         - ignoring 1 / sir, which is already included in the normalisation factor;
#         - multiplying by the partial fraction sjr / (sir + sjr) to regulate for sjr -> 0.
#         """
#
#         pipj = pi.dot(pj)
#         pijpr = pr.dot(pi + pj)
#         return 2 * pipj / pijpr
#
#     def evaluate_subtraction_current(
#         self, current,
#         higher_PS_point=None, momenta_dict=None, reduced_process=None,
#         hel_config=None, **opts):
#
#         if higher_PS_point is None:
#             raise CurrentImplementationError(
#                 self.name() + " needs the higher phase-space point.")
#         if momenta_dict is None:
#             raise CurrentImplementationError(
#                 self.name() + " requires a momenta dictionary.")
#         if reduced_process is None:
#             raise CurrentImplementationError(
#                 self.name() + " requires a reduced_process.")
#         if not hel_config is None:
#             raise CurrentImplementationError(
#                 self.name() + " does not support helicity assignment.")
#
#         # Retrieve alpha_s and mu_r
#         model_param_dict = self.model.get('parameter_dict')
#         alpha_s = model_param_dict['aS']
#         mu_r = model_param_dict['MU_R']
#         # Include the counterterm only in a part of the phase space
#         # Retrieve leg numbers
#         children = tuple(self.leg_numbers_map[i]
#                          for i in sorted(self.leg_numbers_map.keys()))
#         parent = momenta_dict.inv[frozenset(children)]
#
#         # Perform mapping
#         self.mapping_singular_structure.legs = self.get_recoilers(
#             reduced_process, excluded=(parent, ) )
#         lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
#             higher_PS_point, self.mapping_singular_structure, momenta_dict,
#             compute_jacobian=self.divide_by_jacobian)
#
#         # Retrieve kinematics
#         Q = mapping_vars['Q']
#         pC = sum(higher_PS_point[child] for child in children)
#         qC = lower_PS_point[parent]
#         if self.is_cut(Q=Q, pC=pC):
#             return utils.SubtractionCurrentResult.zero(
#                 current=current, hel_config=hel_config)
#         reduced_kinematics = (None, lower_PS_point)
#
#         # Evaluate collinear subtracted kernel
#         zs, kTs = self.variables(higher_PS_point, qC, children, Q=Q)
#         evaluation = self.kernel(zs, kTs, parent, reduced_kinematics)
#
#         # Start handling the soft counterterms
#         # First check which particles can go soft
#
#         gluons = [self.is_gluon(l,self.model) for l in self.structure.get_all_legs()]
#         # If any can go soft, we need to add eikonal pieces
#         if any(gluons):
#             # Find all colored leg numbers ?except for the parent? in the reduced process
#             all_colored_parton_numbers = []
#             for leg in reduced_process.get('legs'):
#                 if self.model.get_particle(leg.get('id')).get('color') == 1:
#                     continue
#                 all_colored_parton_numbers.append(leg.get('number'))
#
#             color_correlation_index = 1
#             p0 = higher_PS_point[children[0]]
#             p1 = higher_PS_point[children[1]]
#
#             # Now loop over the colored parton number pairs (parent, j)
#             # and add the corresponding contributions to this current
#             # NB: we are using a (0,1) collinear mapping so in the S(0) limit, parent=p1
#             # and in the S(1) limit, parent = p0.
#             # As a result, S(0) emitted from the dipole (1,j) has the unresolved color correlation (parent,j)
#             # as well as S(1) emitted from the dipole (0,j). As a result, for each j, we have a single color
#             # correlation for the two eikonal pieces pj.p1/(pj.p0)/(pj.p0+p1.p0) and pj.p0/(pj.p1)/(pj.p1+p0.p1)
#
#             for j in all_colored_parton_numbers:
#                 # Write the eikonal for that pair
#                 if j == parent:
#                     # j is already the other leg of the dipole,
#                     # we skip as we don't handle massive emitters
#                     continue
#                 pj = higher_PS_point[j]
#                 evaluation['color_correlations'].append(((parent, j),))
#                 eiks = 0
#                 if gluons[0]:
#                     # Soft 0 emitted from (1,j) with C(0,j) screened
#                     eiks -= self.partial_fractionned_eikonal(p1, pj, p0)
#                 if gluons[1]:
#                     # Soft 1 emitted from (0,j) with C(1,j) screened
#                     eiks -= self.partial_fractionned_eikonal(p0, pj, p1)
#                 evaluation['values'][(0, color_correlation_index,0)] = {'finite': eiks}
#                 color_correlation_index += 1
#
#         # Add the normalization factors
#         pC2 = pC.square()
#         norm = 8. * math.pi * alpha_s / pC2
#         norm *= self.factor(Q=Q, pC=pC, qC=qC)
#         for k in evaluation['values']:
#             evaluation['values'][k]['finite'] *= norm
#
#         # Construct and return result
#         result = utils.SubtractionCurrentResult()
#         result.add_result(
#             evaluation,
#             hel_config=hel_config,
#             squared_orders=tuple(sorted(current.get('squared_orders').items())))
#         return result
#
#
# class QCD_final_collinear_0_qqx(QCD_final_collinear_with_soft):
#     """q q~ collinear tree-level current."""
#
#     structure = sub.SingularStructure(sub.CollStructure(
#         sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
#         sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL), ))
#
#     def kernel(self, zs, kTs, parent, reduced_kinematics):
#
#         # Retrieve the collinear variables
#         z = zs[0]
#         kT = kTs[0]
#         # Instantiate the structure of the result
#         evaluation = utils.SubtractionCurrentEvaluation({
#             'spin_correlations': [None, ((parent, (kT,)),), ],
#             'color_correlations': [None],
#             'reduced_kinematics': [reduced_kinematics],
#             'values': {(0, 0, 0): {'finite': None},
#                        (1, 0, 0): {'finite': None}, }
#         })
#         # Compute the kernel
#         # The line below implements the g_{\mu\nu} part of the splitting kernel.
#         # Notice that the extra longitudinal terms included in the spin-correlation 'None'
#         # from the relation:
#         #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
#         #    = g^{\mu\nu} + longitudinal terms
#         # are irrelevant because Ward identities evaluate them to zero anyway.
#         evaluation['values'][(0, 0, 0)]['finite'] = self.TR
#         evaluation['values'][(1, 0, 0)]['finite'] = 4 * self.TR * z * (1-z) / kT.square()
#         return evaluation
#
#
# class QCD_final_collinear_0_gq(QCD_final_collinear_with_soft):
#     """g q collinear tree-level current.
#     This contains
#     - the hard-collinear kernel for g(0)q(1) being collinear
#     - the pieces of the eikonal terms for g(0) soft emitted from q(1)X(j) which also diverge
#     in the g(0)q(1) collinear limit
#     """
#
#     structure = sub.SingularStructure(sub.CollStructure(
#         sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL),
#         sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL), ))
#
#     def kernel(self, zs, kTs, parent, reduced_kinematics):
#
#         # Retrieve the collinear variables
#         z = zs[0]
#         # Instantiate the structure of the result
#         evaluation = utils.SubtractionCurrentEvaluation({
#             'spin_correlations': [None],
#             'color_correlations': [None],
#             'reduced_kinematics': [reduced_kinematics],
#             'values': {(0, 0, 0): {'finite': None}}
#         })
#         # Compute the kernel using
#         # f9d0839fc58905d67367e3e67efabee05ee390f9:madgraph/iolibs/template_files/OLD_subtraction/cataniseymour/NLO/local_currents.py:146
#         evaluation['values'][(0, 0, 0)]['finite'] = \
#             self.CF*z
#         return evaluation
#
#
# class QCD_final_collinear_0_gg(QCD_final_collinear_with_soft):
#     """g g collinear tree-level current."""
#
#     structure = sub.SingularStructure(sub.CollStructure(
#         sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL),
#         sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), ))
#
#     def kernel(self, zs, kTs, parent, reduced_kinematics):
#
#         # Retrieve the collinear variables
#         z = zs[0]
#         kT = kTs[0]
#         # Instantiate the structure of the result
#         evaluation = utils.SubtractionCurrentEvaluation({
#             'spin_correlations': [None, ((parent, (kT,)),), ],
#             'color_correlations': [None],
#             'reduced_kinematics': [reduced_kinematics],
#             'values': {(0, 0, 0): {'finite': None},
#                        (1, 0, 0): {'finite': None}, }
#         })
#         # Compute the kernel
#         # The line below implements the g_{\mu\nu} part of the splitting kernel.
#         # Notice that the extra longitudinal terms included in the spin-correlation 'None'
#         # from the relation:
#         #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
#         #    = g^{\mu\nu} + longitudinal terms
#         # are irrelevant because Ward identities evaluate them to zero anyway.
#         evaluation['values'][(0, 0, 0)]['finite'] = \
#             0.
#         evaluation['values'][(1, 0, 0)]['finite'] = \
#             -2. * self.CA * 2. * z * (1. - z) / kT.square()
#         return evaluation
#
# #=========================================================================================
# NLO soft current
#=========================================================================================

class NoSoft(currents.QCDLocalCurrent):
    """There is no explicit soft counterterm in this scheme. It is distributed into collinear counterterms"""

    is_zero = True
    squared_orders = {'QCD': 2}
    n_loops = 0
    structure = sub.SingularStructure(
        sub.SoftStructure(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL)) )

    
    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, **opts ):

        # Just return 0
        result = utils.SubtractionCurrentResult()
        result.add_result(
            utils.SubtractionCurrentEvaluation.zero(),
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())) )
        return result

#=========================================================================================
# NLO final soft-collinear currents
#=========================================================================================

class NoSoftCollinear(currents.QCDLocalCurrent):
    """There is no explicit soft counterterm in this scheme. It is distributed into collinear counterterms"""

    is_zero=True
    squared_orders = {'QCD': 2}
    n_loops = 0

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(1,  1, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_g = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), ) )
    structure_q = sub.SingularStructure(substructures=(soft_coll_structure_q, ))
    structure_g = sub.SingularStructure(substructures=(soft_coll_structure_g, ))

    @classmethod
    def does_implement_this_current(cls, current, model):

        if not cls.check_current_properties(current): return None

        color_charge = 'CF'
        leg_numbers_map = cls.structure_q.map_leg_numbers(
            current.get('singular_structure'), [range(1, model.get_nflav()+1)])
        if leg_numbers_map is None:
            color_charge = 'CA'
            leg_numbers_map = cls.structure_g.map_leg_numbers(
                current.get('singular_structure'), [range(1, model.get_nflav()+1)])
            if leg_numbers_map is None:
                return None
        mapping_singular_structure = current.get('singular_structure').get_copy()
        return {
            'leg_numbers_map': leg_numbers_map,
            'color_charge': color_charge,
            'mapping_singular_structure': mapping_singular_structure
        }

    def __init__(self, *args, **opts):

        for opt_name in ['color_charge', ]:
            try:
                setattr(self, opt_name, opts.pop(opt_name))
            except KeyError:
                raise CurrentImplementationError(
                    "__init__ of " + self.__class__.__name__ + " requires " + opt_name)

        super(NoSoftCollinear, self).__init__(*args, **opts)
        # At this state color_charge is the string of the group factor ('CA' or 'CF');
        # now that the super constructor has been called,
        # the group factors have been initialized and we can retrieve them.
        self.color_charge = getattr(self, self.color_charge)
    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, **opts ):

        # Just return 0
        result = utils.SubtractionCurrentResult()
        result.add_result(
            utils.SubtractionCurrentEvaluation.zero(),
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())) )
        return result