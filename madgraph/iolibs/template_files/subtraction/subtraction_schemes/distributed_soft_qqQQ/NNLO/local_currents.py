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

#from commons.universal_kernels import AltarelliParisiKernels, SoftKernels #TODO DEV
from madgraph.iolibs.template_files.subtraction.commons.universal_kernels import AltarelliParisiKernels, SoftKernels #TODO DEV

#import commons.utils as utils# TODO DEV
import madgraph.iolibs.template_files.subtraction.commons.utils as utils #TODO DEV

#import commons.QCD_local_currents as currents #TODO DEV
import madgraph.iolibs.template_files.subtraction.commons.QCD_local_currents as currents #TODO DEV

#import commons.factors_and_cuts as factors_and_cuts #TODO DEV
import madgraph.iolibs.template_files.subtraction.commons.factors_and_cuts as factors_and_cuts #TODO DEV

from bidict import bidict

import madgraph.iolibs.template_files.subtraction.subtraction_schemes.distributed_soft_qqQQ.distributed_soft_qqQQ_config as scheme_config # TODO DEV
#import distributed_soft_qqQQ_config as scheme_config #TODO DEV

import madgraph.iolibs.template_files.subtraction.subtraction_schemes.distributed_soft_qqQQ.NNLO.QQqq_kernels as kernels # TODO DEV
#import distributed_soft_qqQQ.NNLO.QQqq_kernels as kernels # TODO DEV

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
        result = kernels.C123_Qqqx_kernel(z1,z2,z3,s12,s13,s23)
        result *=  self.CF * self.TR
        evaluation['values'][(0,0,0)] = {'finite':result}
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

    is_zero = True

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