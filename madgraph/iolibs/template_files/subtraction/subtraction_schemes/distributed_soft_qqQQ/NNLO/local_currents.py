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
# Helper functions
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

def kTdiff(i, j, zs, kTs):
    """Compute the difference vector (kTi/zi - kTj/zj)."""

    return kTs[i-1]/zs[i-1] - kTs[j-1]/zs[j-1]

def sij(i, j, zs, kTs, kTdiffij=None):
    """Compute the invariant mass sij from zs and kTs."""

    if kTdiffij: kTdiffij2 = kTdiffij.square()
    else: kTdiffij2 = kTdiff(i, j, zs, kTs).square()
    return -zs[i-1]*zs[j-1]*kTdiffij2

def tijk(zi, zj, sij, sik, sjk):

    return (2*(zi*sjk-zj*sik)+(zi-zj)*sij)/(zi+zj)

def t(i, j, k, zs, kTs, s_ij=None, s_ik=None, s_jk=None):
    """Compute the invariant t_{ij,k} from zs and kTs."""

    if not s_ij: s_ij = sij(i, j, zs, kTs)
    if not s_ik: s_ik = sij(i, k, zs, kTs)
    if not s_jk: s_jk = sij(j, k, zs, kTs)
    return tijk(zs[i-1],zs[j-1],s_ij,s_ik,s_jk)

def get_intermediate_PS_point(self, higher_PS_point, children):

    recoilers = tuple(
        i for i in higher_PS_point.keys()
        if i not in (1, 2, children[0], children[1]) )
    def dummy_leg(i):
        return subtraction.SubtractionLeg(i, 21, subtraction.SubtractionLeg.FINAL)
    structure = subtraction.SingularStructure(
        substructures=[subtraction.CollStructure(
            legs=(dummy_leg(children[0]), dummy_leg(children[1])) ), ],
        legs=(dummy_leg(r) for r in recoilers) )
    leg_numbers_map = subtraction.bidict({
        i: frozenset({i,})
        for i in higher_PS_point.keys()
        if i not in (children[0], children[1])})
    leg_numbers_map[1000] = frozenset({children[0], children[1]})
    mapping = mappings.FinalRescalingOneMapping
    intermediate_ps_point, mapping_vars = mapping.map_to_lower_multiplicity(
        higher_PS_point, structure, leg_numbers_map, compute_jacobian=True)
    return intermediate_ps_point, mapping_vars

def n_final_coll_variables_triple(higher_PS_point, qC, children,**opts):
    """Obtain the n_final_coll_variables definition of z and kT and also provide the pair-wise invariants from the higher_PS_point.

    This is redundant information but screens possible issues when using the freedom to redefine kT in a way that would break the relationship used in the function sij above.

    In the spirit of future low-level hardcoded structures, we hardcode the expectation that there are three relevant invariants
    """
    return currents.n_final_coll_variables(higher_PS_point, qC, children,**opts)+((
        2. * higher_PS_point[children[0]].dot(higher_PS_point[children[1]]),
        2. * higher_PS_point[children[0]].dot(higher_PS_point[children[2]]),
        2. * higher_PS_point[children[1]].dot(higher_PS_point[children[2]])
    ),)

#=========================================================================================
# Variables, mappings, jacobians, factors and cuts
#=========================================================================================
# Note that variables, factors and cuts will be class members by design
# so they can easily be overridden by subclasses.
# They will be taken from the following variables
# so we can quickly switch them coherently across the entire subtraction scheme.

variables = n_final_coll_variables_triple
#currents.n_final_coll_variables
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
# NNLO final-collinear currents
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

#    def __init__(self,*args,**kwargs):
        # TODO THIS SHOULD NOT BE USED YET
#        raise NotImplementedError

    def C123_Qqqx_kernel(self, z1, z2, z3, s12, s13, s23, s123):
        """Triple collinear splitting kernel"""
        t123 = tijk(z1, z2, s12, s13, s23)
        sqrbrk  = -(t123 ** 2)/(s12*s123)
        sqrbrk += (4*z3 + (z1-z2)**2) / (z1+z2)
        sqrbrk += z1 + z2 - s12/s123
        return 1./2.*self.CF*self.TR*s123*sqrbrk / (s12)

    def C123C12_kernel(self, higher_PS_point, parent_momentum, children, **opts):

        # Rebuild the two-stage mapping
        interm_PS_point, interm_mapping_vars = self.get_intermediate_PS_point(
            higher_PS_point, children )
        final__PS_point, final__mapping_vars = self.get_final_PS_point(
            interm_PS_point, children )

        # Retrieve momenta
        p1 = higher_PS_point[children[0]]
        p2 = higher_PS_point[children[1]]
        p3 = higher_PS_point[children[2]]
        p12 = p1 + p2
        p123 = p12 + p3
        p12hat = interm_PS_point[1000]
        p3hat = interm_PS_point[children[2]]
        p123hat = p12hat + p3hat
        p123tilde = final__PS_point[2000]
        Q = interm_mapping_vars['Q']

        # Compute momentum fractions and transverse momenta
        zs_interm, kTs_interm = self.variables(
            higher_PS_point, p12hat, children[:2], Q=Q)
        zs_final_, kTs_final_ = self.variables(
            interm_PS_point, p123tilde, (1000, children[2]), Q=Q)
        z1 = zs_interm[0]
        k1perp = kTs_interm[0]
        z12 = zs_final_[0]
        k12perp = kTs_final_[0]

        # Build scalar products
        s12 = p12.square()
        s12hat_3hat = 2*p3hat.dot(p12hat)
        k1perp2 = k1perp.square()
        k12perp2 = k12perp.square()
        kperpSP = 2*k1perp.dot(k12perp)

        # Construct the iterated current C(C(1,2),3)
        perpterm = ((1-z12)/z12) * (kperpSP**2)/(k1perp2*k12perp2)
        pqg = (1+(1-z12)**2) / z12
        brk = z12 + perpterm
        C123C12_current = 4*(pqg - 2*z1*(1-z1)*brk) / (s12hat_3hat*s12)

        # If jacobians are active, correct the one-step jacobian to the two-step one
        # and correct current factors altogether
        try:
            jacobian = opts['jacobian']
            jacobian /= (final__mapping_vars['jacobian']*interm_mapping_vars['jacobian'])
        except KeyError:
            jacobian = 1
        # Correct current factors if they are active
        factor_interm = self.factor(Q=Q, pC=p12, qC=p12hat)
        factor_final_ = self.factor(Q=Q, pC=p123hat, qC=p123tilde)
        factor_direct = self.factor(Q=Q, pC=p123, qC=p123tilde)
        factor = factor_interm * factor_final_ / factor_direct

        return factor*C123C12_current*jacobian


    def kernel(self, evaluation, parent, zs, kTs , sijs):
        ker = self.C123_Qqqx_kernel(zs[0],zs[1],zs[2], sijs[0], sijs[1], sijs[2], sum(sijs))
        evaluation['values'][(0, 0, 0)] = {'finite' : ker}
        return evaluation
