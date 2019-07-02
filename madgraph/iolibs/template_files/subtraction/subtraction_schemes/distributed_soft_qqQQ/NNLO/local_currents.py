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

from commons.universal_kernels import AltarelliParisiKernels, SoftKernels

import commons.utils as utils

import commons.QCD_local_currents as currents

import commons.factors_and_cuts as factors_and_cuts

from bidict import bidict

import distributed_soft_qqQQ_config as scheme_config

import NNLO.QQqq_kernels as kernels

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#==================================================================
# NNLO final collinear mother class
#==================================================================

class QCD_final_collinear_0_XXX_soft_distrib(currents.GeneralQCDLocalCurrent):
    """Two-collinear tree-level current."""

    squared_orders = {'QCD': 4}
    n_loops = 0

    is_cut = scheme_config.final_coll_cut
    factor = scheme_config.final_coll_factor
    mapping = scheme_config.final_coll_mapping
    variables = staticmethod(scheme_config.global_final_coll_variables)
    divide_by_jacobian = scheme_config.divide_by_jacobian
    get_recoilers = scheme_config.get_recoilers

    mapping_rules = [
        {
            'singular_structure': None,
            'mapping': mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict': bidict({1001: frozenset((0, 1))}),
            'variables': currents.CompoundVariables(scheme_config.final_coll_variables),
            'is_cut': is_cut,
            'reduced_recoilers': get_recoilers,
            'additional_recoilers': sub.SubtractionLegSet([sub.SubtractionLeg(2, +2, sub.SubtractionLeg.FINAL)]),
        },
        {
            'singular_structure': None,
            'mapping': mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict': bidict({-1: frozenset((1001, 2))}),
            'variables': currents.CompoundVariables(scheme_config.final_coll_variables),
            'is_cut': is_cut,
            'reduced_recoilers': get_recoilers,
            'additional_recoilers': sub.SubtractionLegSet([]),
        }
    ]





class QCD_final_triple_collinear_qqxQ(QCD_final_collinear_0_XXX_soft_distrib):
    structure = [sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(2, +2, sub.SubtractionLeg.FINAL), ))]

    mapping_rules = QCD_final_collinear_0_XXX_soft_distrib.mapping_rules[:]
    mapping_rules[0]['singular_structure'] = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),),)
    mapping_rules[1]['singular_structure'] = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(2, +2, sub.SubtractionLeg.FINAL),),)

    def kernel(self, evaluation, all_steps, global_variables):
        """Combining the whole qbar-q-Q forest
        The implementation is not memory-cautious at all but aims at readability
        """
        # Variables of the triple splitting (Q->Q q qbar)

        (z1,z2,z3)=global_variables['zs']
        (s12,s13,s23)=global_variables['ss']
        # TODO DEV Variables used in the nested currents not used as long as we keep things separated
        # (k1perp, k2perp, k3perp) = global_variables['kTs']
        #
        # # Variables of the nested double splitting (g->q qbar)
        # (z1_qq,z2_qq) = all_steps[0]['variables'][0]['zs']
        # (k1perp_qq,k2perp_qq) = all_steps[0]['variables'][0]['kTs']
        # (s12_qq,) = all_steps[0]['variables'][0]['ss'] #This is the same as s12 above
        #
        # # Variables of the mapped double splitting (Q->Qg)
        # (z3hat,z12hat) = all_steps[1]['variables'][0]['zs']
        # (k3perphat,k12perphat) = all_steps[1]['variables'][0]['kTs']
        # (s_12hat_3hat,) = all_steps[1]['variables'][0]['ss']

        # Triple-collinear kernel
        result = kernels.C123_minus_C123S12_Qqqx_kernel(z1, z2, z3, s12, s13, s23)

        # The color factor should be better handled
        result *= self.CF * self.TR
        evaluation['values'][(0,0,0)] = {'finite':0.5*result}
        return evaluation

    def soft_kernel(self, evaluation, colored_parton_numbers, all_steps_info, global_variables):
        # in Q -> Q q qx, the q qx pair (0,1) can go soft
        # The parent is the double unresolved mapped Q, obtained after
        # the full mapping (all_steps_info[1])
        parent_ID = all_steps_info[1]['bundles_info'][0]['parent']
        # Obtain the momenta of the Quark, quark and anti-quark
        q_ID = global_variables['overall_children'][0][0]
        qbar_ID = global_variables['overall_children'][0][1]
        Q_ID = global_variables['overall_children'][0][2]
        pq = all_steps_info[0]['higher_PS_point'][q_ID]
        pqbar = all_steps_info[0]['higher_PS_point'][qbar_ID]
        pQ = all_steps_info[0]['higher_PS_point'][Q_ID]

        # Building the variables
        # Final state invariants
        (sqqbar,sqQ,sqbarQ) = global_variables['ss']
        # Variables of the g->q qx splitting
        z1_qq = all_steps_info[0]['variables'][0]['zs'][0]
        k1perp = all_steps_info[0]['variables'][0]['kTs'][0]
        # Variables of the Q->Qg splitting)
        (p3hat,p12hat) = [all_steps_info[1]['higher_PS_point'][i] for i in all_steps_info[1]['bundles_info'][0]['final_state_children']]


        # We loop over dipoles (k,parent) in the reduced process
        # The dipoles are labelled by color_correlation_index
        color_correlation_index = 1
        for k in colored_parton_numbers:
            if k == parent_ID:
                # No color-diagonal emission
                continue
            else:
                # Add the label (parent_ID,k) to the ordered list of
                # emitting dipoles
                evaluation['color_correlations'].append(((parent_ID, k),))
                # Evaluate the eikonal factor
                pk = all_steps_info[0]['higher_PS_point'][k]
                sQk = 2*pk.dot(pQ)
                sqk = 2*pk.dot(pq)
                sqbark = 2*pk.dot(pqbar)
                weight = 0
                weight += kernels.S12_qqx_PF_kernel(sqQ, sqbarQ, sQk, sqk, sqbark, sqqbar)
                weight -= kernels.S12C12_Qqqx_PF_kernel(z1_qq, sqqbar, p12hat, p3hat, pk, k1perp)
                weight *= self.TR
                evaluation['values'][(0, color_correlation_index,0)] = {'finite': weight}
                color_correlation_index += 1


        return evaluation



class NoFinalNestedCollinear(QCD_final_collinear_0_XXX_soft_distrib):
    structure = [sub.SingularStructure(sub.CollStructure(
            sub.CollStructure(
                sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),),
            sub.SubtractionLeg(2, +2, sub.SubtractionLeg.FINAL), ))]

    # mapping_rules = QCD_final_collinear_0_XXX_soft_distrib.mapping_rules[:]
    # mapping_rules[0]['singular_structure'] = structure[0]

    mapping_rules = QCD_final_collinear_0_XXX_soft_distrib.mapping_rules[:]
    mapping_rules[0]['singular_structure'] = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),),)
    mapping_rules[1]['singular_structure'] = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(2, +2, sub.SubtractionLeg.FINAL),),)

    def kernel(self, evaluation, all_steps, global_variables):
        """TODO DEV COMBINE ABOVE AFTER DEBUG Combining the whole qbar-q-Q forest
        The implementation is not memory-cautious at all but aims at readability
        """
        # Variables of the triple splitting (Q->Q q qbar)
        (z1,z2,z3)=global_variables['zs']
        (s12,s13,s23)=global_variables['ss']
        (k1perp, k2perp, k3perp) = global_variables['kTs']

        # Variables of the nested double splitting (g->q qbar)
        (z1_qq,z2_qq) = all_steps[0]['variables'][0]['zs']
        (k1perp_qq,k2perp_qq) = all_steps[0]['variables'][0]['kTs']
        (s12_qq,) = all_steps[0]['variables'][0]['ss'] #This is the same as s12 above

        # Variables of the mapped double splitting (Q->Qg)
        (z3hat,z12hat) = all_steps[1]['variables'][0]['zs']
        (k3perphat,k12perphat) = all_steps[1]['variables'][0]['kTs']
        (s_12hat_3hat,) = all_steps[1]['variables'][0]['ss']

        # Triple-collinear kernel
        result = 0.25*kernels.C123C12_Qqqx_kernel(z1_qq, z12hat, k1perp_qq, k12perphat, s12, s_12hat_3hat,s12+s13+s23) # TODO TIDY UP THESE FACTORS OF 2
        result *=  self.CF * self.TR
        evaluation['values'][(0,0,0)] = {'finite':result}
        return evaluation


class NoFinalDoubleSoftinTripleCollinear(QCD_final_collinear_0_XXX_soft_distrib):
    structure = [sub.SingularStructure(sub.CollStructure(
            sub.SoftStructure(
                sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),),
            sub.SubtractionLeg(2, +2, sub.SubtractionLeg.FINAL), ))]

    # mapping_rules = QCD_final_collinear_0_XXX_soft_distrib.mapping_rules[:]
    # mapping_rules[0]['singular_structure'] = structure[0]

    is_zero = True

class NoFinalDoubleSoftCollinearinTripleCollinear(QCD_final_collinear_0_XXX_soft_distrib):
    structure = [sub.SingularStructure(sub.CollStructure(
            sub.SoftStructure(sub.CollStructure(
                sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),),),
            sub.SubtractionLeg(2, +2, sub.SubtractionLeg.FINAL), ))]

    # mapping_rules = QCD_final_collinear_0_XXX_soft_distrib.mapping_rules[:]
    # mapping_rules[0]['singular_structure'] = structure[0]
    # --------------------------- TODO DEV
    # TMP
    # --------------------------- TODO DEV
    mapping_rules = QCD_final_collinear_0_XXX_soft_distrib.mapping_rules[:]
    mapping_rules[0]['singular_structure'] = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),),)
    mapping_rules[1]['singular_structure'] = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(1001, 21, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(2, +2, sub.SubtractionLeg.FINAL),),)

    def kernel(self, evaluation, all_steps, global_variables):
        """TODO DEV COMBINE ABOVE AFTER DEBUG Combining the whole qbar-q-Q forest
        The implementation is not memory-cautious at all but aims at readability
        """
        # Variables of the triple splitting (Q->Q q qbar)
        (z1,z2,z3)=global_variables['zs']
        (s12,s13,s23)=global_variables['ss']
        (k1perp, k2perp, k3perp) = global_variables['kTs']

        # Variables of the nested double splitting (g->q qbar)
        (z1_qq,z2_qq) = all_steps[0]['variables'][0]['zs']
        (k1perp_qq,k2perp_qq) = all_steps[0]['variables'][0]['kTs']
        (s12_qq,) = all_steps[0]['variables'][0]['ss'] #This is the same as s12 above

        # Variables of the mapped double splitting (Q->Qg)
        (z3hat,z12hat) = all_steps[1]['variables'][0]['zs']
        (k3perphat,k12perphat) = all_steps[1]['variables'][0]['kTs']
        (s_12hat_3hat,) = all_steps[1]['variables'][0]['ss']
        (p3hat,p12hat) = [all_steps[1]['higher_PS_point'][i] for i in all_steps[1]['bundles_info'][0]['final_state_children']]
        # Triple-collinear kernel
        result = 0.25*kernels.C123S12C12_Qqqx_PF_kernel(z1_qq,z12hat,s12,p12hat,p3hat,k1perp_qq,k12perphat)
        result *=  self.CF * self.TR
        evaluation['values'][(0,0,0)] = {'finite':result}
        return evaluation

    # --------------------------- TODO DEV



    #is_zero = True TODO DEV

class NoFinalDoubleSoft(QCD_final_collinear_0_XXX_soft_distrib):
    structure = [sub.SingularStructure(
            sub.SoftStructure(
                sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),
            ),)]

    # mapping_rules = QCD_final_collinear_0_XXX_soft_distrib.mapping_rules[:]
    # mapping_rules[0]['singular_structure'] = structure[0]

    is_zero = True

class NoFinalDoubleSoftCollinear(QCD_final_collinear_0_XXX_soft_distrib):
    structure = [sub.SingularStructure(
            sub.SoftStructure(sub.CollStructure(
                sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
                sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),),),)]

    # mapping_rules = QCD_final_collinear_0_XXX_soft_distrib.mapping_rules[:]
    # mapping_rules[0]['singular_structure'] = structure[0]

    is_zero = True
