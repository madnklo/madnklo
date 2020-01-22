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

import commons.utils as utils
import commons.QCD_local_currents as currents
# import commons.factors_and_cuts as factors_and_cuts

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError
#=========================================================================================
# Choice of recoilers
#=========================================================================================

def get_recoilers(reduced_process, excluded=()):

    model = reduced_process.get('model')
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            model.get_particle(leg['id']).get('mass').upper() == 'ZERO',
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded
        ])
    ])


#=========================================================================================
# Variables, mappings, jacobians, factors and cuts
#=========================================================================================
# Note that variables, factors and cuts will be class members by design
# so they can easily be overridden by subclasses.
# They will be taken from the following variables
# so we can quickly switch them coherently across the entire subtraction scheme.

variables = currents.n_final_coll_variables
mapping = mappings.FinalGroupingMapping
# soft_mapping = mappings.SoftVsFinalPureRescalingMapping
# soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
#    soft_mapping=soft_mapping, collinear_mapping=coll_mapping)
divide_by_jacobian = True
# factor_coll = factors_and_cuts.factor_coll
# factor_soft = factors_and_cuts.factor_soft
# is_cut_coll = factors_and_cuts.cut_coll
# is_cut_soft = factors_and_cuts.cut_soft
no_cut = currents.no_cut
no_factor = currents.no_factor

#=========================================================================================
# NLO final-collinear currents
#=========================================================================================

class QCD_final_collinear_0_XX(currents.QCDLocalCollinearCurrent):
    """Two-collinear tree-level current."""
    squared_orders = {'QCD': 4}
    n_loops = 0

    is_cut = staticmethod(no_cut)
    factor = staticmethod(no_factor)
    mapping = mapping
    variables = staticmethod(variables)
    divide_by_jacobian = divide_by_jacobian
    get_recoilers = staticmethod(get_recoilers)

class QCD_final_collinear_0_QQxq(QCD_final_collinear_0_XX):
    """q Q Q~ triple-collinear tree-level current."""

    structure = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(2, +2, sub.SubtractionLeg.FINAL), ))

    def __init__(self,*args,**kwargs):
        # TODO THIS SHOULD NOT BE USED YET
        raise NotImplementedError

    def C123_Qqqx_kernel(self, z1, z2, z3, s12, s13, s23, s123):

        t123 = tijk(z1, z2, s12, s13, s23)
        sqrbrk  = -(t123 ** 2)/(s12*s123)
        sqrbrk += (4*z3 + (z1-z2)**2) / (z1+z2)
        sqrbrk += z1 + z2 - s12/s123
        return sqrbrk / (s12*s123)


    # def kernel(self, zs, kTs, invariants, parent, reduced_kinematics):
    #     #raise NotImplementedError TODO DEV
    #     # Retrieve the collinear variables
    #     z = zs[0]
    #     kT = kTs[0]
    #     # Instantiate the structure of the result
    #     evaluation = utils.SubtractionCurrentEvaluation({
    #         'spin_correlations': [None, ((parent, (kT,)),), ],
    #         'color_correlations': [None],
    #         'reduced_kinematics': [reduced_kinematics],
    #         'values': {(0, 0, 0): {'finite': None},
    #                    (1, 0, 0): {'finite': None}, }
    #     })
    #     # Compute the kernel
    #     # The line below implements the g_{\mu\nu} part of the splitting kernel.
    #     # Notice that the extra longitudinal terms included in the spin-correlation 'None'
    #     # from the relation:
    #     #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
    #     #    = g^{\mu\nu} + longitudinal terms
    #     # are irrelevant because Ward identities evaluate them to zero anyway.
    #     evaluation['values'][(0, 0, 0)]['finite'] = self.TR
    #     evaluation['values'][(1, 0, 0)]['finite'] = 4 * self.TR * z * (1-z) / kT.square()
    #     return evaluation