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
"""Implementation of NNLO type of currents."""

import os
import math

import madgraph.core.subtraction as subtraction
import madgraph.integrator.mappings as mappings
import madgraph.various.misc as misc

try:
    # First try to import this in the context of the exported currents
    import SubtractionCurrents.subtraction_current_implementations_utils as utils
    import SubtractionCurrents.QCD_local_currents as currents
except ImportError:
    # If not working, then it must be within MG5_aMC context:
    import madgraph.iolibs.template_files.\
                   subtraction.subtraction_current_implementations_utils as utils
    import madgraph.iolibs.template_files.\
                   subtraction.QCD_local_currents as currents

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# Auxiliary functions
#=========================================================================================

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

#=========================================================================================
# NNLO final-collinear currents, containing the soft limits
#=========================================================================================

class QCD_final_collinear_0_QQxq(currents.QCDLocalCollinearCurrent):
    """Q Q~ q collinear tree-level current."""

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NNLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 4, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that there are 3 massless final state quarks or antiquarks
        if len(ss.legs) != 3: return None
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
            if not cls.is_quark(leg, model): return None
        # Look for a quark/antiquark pair
        pair = None
        for i in range(len(ss.legs)):
            for j in range(i+1, len(ss.legs)):
                if cls.are_antiparticles(ss.legs[i], ss.legs[j]):
                    pair = (ss.legs[i], ss.legs[j])
                    continue
            if pair is not None: continue
        if pair is None: return None
        # Identify the remaining quark
        other_quarks = [leg for leg in ss.legs if leg not in pair]
        # Since leg numbers have been discarded, equal legs will not appear here
        # Thus if the quark species were the same, other_quarks = []
        if len(other_quarks) != 1: return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        if cls.are_antiparticles(legs[0], legs[1]):
            return (legs[0].n, legs[1].n, legs[2].n)
        elif cls.are_antiparticles(legs[0], legs[2]):
            return (legs[0].n, legs[2].n, legs[1].n)
        else:
            return (legs[1].n, legs[2].n, legs[0].n)

    def S12_kernel(self, sir, sis, sik, skr, sks, srs):

        sk_rs = skr + sks
        si_rs = sir + sis
        skrs = sk_rs + srs
        sirs = si_rs + srs
        frac = skrs / (sirs + skrs)
        return frac * 4 * self.TR * (
            (sir*sks + skr*sis - sik*srs) / (sirs*sk_rs)
            - sir * sis / (sirs ** 2)
            - skr * sks / (sk_rs ** 2)
        ) / (srs ** 2)

    def C123_kernel(self, z1, z2, z3, s12, s13, s23, s123):

        t123 = tijk(z1, z2, s12, s13, s23)
        sqrbrk  = -(t123 ** 2)/(s12*s123)
        sqrbrk += (4*z3 + (z1-z2)**2) / (z1+z2)
        sqrbrk += z1 + z2 - s12/s123
        return sqrbrk / (s12*s123)

    def C123_minus_C123S12_kernel(self, z1, z2, z3, s12, s13, s23, s123):

        term1 = -s123**(-2)
        term2 = (z1**2+z2**2)/(z1+z2)/(s12*s123)
        return 2*(term1+term2)

    def C123S12_kernel(self, z1, z2, z3, s12, s13, s23, s123):

        c123 = self.C123_kernel(z1, z2, z3, s12, s13, s23, s123)
        c123_m_c123s12 = self.C123_minus_C123S12_kernel(z1, z2, z3, s12, s13, s23, s123)
        return c123 - c123_m_c123s12

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
        mapping = mappings.FinalGroupingMapping
        intermediate_ps_point, mapping_vars = mapping.map_to_lower_multiplicity(
            higher_PS_point, structure, leg_numbers_map)
        return intermediate_ps_point, mapping_vars

    def C123C12_kernel(self, higher_PS_point, parent_momentum, children, **opts):

        p1 = higher_PS_point[children[0]]
        p2 = higher_PS_point[children[1]]
        p12 = p1+p2
        s12 = p12.square()
        intermediate_PS_point, mapping_vars = self.get_intermediate_PS_point(
            higher_PS_point, children )
        p12hat = intermediate_PS_point[1000]
        p3hat = intermediate_PS_point[children[2]]
        zs, kTs = self.variables(
            higher_PS_point, p12hat, children[:2], Q=mapping_vars['Q'])
        zs2, kTs2 = self.variables(
            intermediate_PS_point, parent_momentum,
            (1000, children[2]), Q=mapping_vars['Q'])
        s12hat_3hat = 2*p3hat.dot(p12hat)
        z1 = zs[0]
        k1perp = kTs[0]
        z12 = zs2[0]
        k12perp = kTs2[0]
        # misc.sprint(z1, k1perp, z12, k12perp)
        k1perp2 = k1perp.square()
        # sperp = 2*p3hat.dot(k1perp)
        k12perp2 = k12perp.square()
        kperpSP = 2*k1perp.dot(k12perp)
        # perpterm = sperp**2/(k1perp2*s12hat_3hat)
        perptermalt = ((1-z12)/z12) * (kperpSP**2)/(k1perp2*k12perp2)
        # perpratio = perpterm / perptermalt
        # n = mappings.LorentzVector([1,0,0,1])
        # Q = mapping_vars['Q']
        # misc.sprint(p12.dot(k12perp), n.dot(k12perp), Q.dot(k12perp))
        # misc.sprint(perpterm, perptermalt, perpratio)
        pqg = (1+(1-z12)**2) / z12
        brk = z12 + perptermalt
        # brk = z12 - sperp**2/(k1perp2*s12hat_3hat)
        return 4*(pqg - 2*z1*(1-z1) *brk) / (s12hat_3hat*s12)

    def S12C12_kernel(self, higher_PS_point, children, emitter, spectator, **opts):

        p1 = higher_PS_point[children[0]]
        p2 = higher_PS_point[children[1]]
        p3 = higher_PS_point[emitter]
        p4 = higher_PS_point[spectator]
        p12  = p1+p2
        p34  = p3+p4
        p123 = p1+p2+p3
        p124 = p1+p2+p4
        s12  = p12.square()
        s34  = p34.square()
        s123 = p123.square()
        s124 = p124.square()
        intermediate_PS_point, mapping_vars = self.get_intermediate_PS_point(
            higher_PS_point, children )
        p12hat = intermediate_PS_point[1000]
        p3hat = intermediate_PS_point[emitter]
        p4hat = intermediate_PS_point[spectator]
        zs, kTs = self.variables(
            higher_PS_point, p12hat, children[:2], Q=mapping_vars['Q'])
        # zs, kTs = self.variables(
        #     higher_PS_point, p12hat, children[:2], Q=opts['Q'])
        s12hat_3hat = 2*p3hat.dot(p12hat)
        s12hat_4hat = 2*p4hat.dot(p12hat)
        s3hat_4hat  = 2*p3hat.dot(p4hat)
        z = zs[0]
        kperp = kTs[0]
        eik_num = -s3hat_4hat
        # eik_num = -2*s34
        sperp3hat = 2*kperp.dot(p3hat)
        sperp4hat = 2*kperp.dot(p4hat)
        kperp2 = kperp.square()
        kperp_part = sperp3hat*sperp4hat/kperp2
        x = (sperp3hat*s12hat_4hat)/(sperp4hat*s12hat_3hat)
        offdiag = 1. - 0.5*(x + 1./x)
        # kperp_part = k1perp.dot(p3)*k1perp.dot(p4)/k1perp.square()
        cor_num = 2*z*(1-z)*kperp_part
        den = s12*s12hat_3hat*s12hat_4hat
        # den = s12*s123*s124
        frac = 1
        frac = s12hat_4hat / (s12hat_3hat+s12hat_4hat)
        # frac = s124 / (s123+s124)
        return 2*self.TR*frac*(eik_num+cor_num*offdiag)/den

    def C123S12C12_kernel(self, higher_PS_point, parent_momentum, children, **opts):

        p1 = higher_PS_point[children[0]]
        p2 = higher_PS_point[children[1]]
        p12 = p1+p2
        s12 = p12.square()
        intermediate_PS_point, mapping_vars = self.get_intermediate_PS_point(
            higher_PS_point, children )
        p12hat = intermediate_PS_point[1000]
        p3hat = intermediate_PS_point[children[2]]
        zs, kTs = self.variables(
            higher_PS_point, p12hat, children[:2], Q=mapping_vars['Q'])
        zhats, kThats = self.variables(
            intermediate_PS_point, parent_momentum,
            (1000, children[2]), Q=mapping_vars['Q'])
        s12hat_3hat = 2*p3hat.dot(p12hat)
        z = zs[0]
        kT = kTs[0]
        zhat = zhats[0]
        kThat = kThats[0]
        sTThat = 2*kT.dot(kThat)
        kT2 = kT.square()
        kThat2 = kThat.square()
        sqrbrk = 1. - z*(1-z)*sTThat*sTThat/(kT2*kThat2)
        return 8*(1-zhat)/zhat/(s12hat_3hat*s12)*sqrbrk

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables and compute basic quantities
        z1, z2, z3 = zs
        s12 = sij(1, 2, zs, kTs)
        s13 = sij(1, 3, zs, kTs)
        s23 = sij(2, 3, zs, kTs)
        misc.sprint(s12,s13,s23)
        s123 = s12 + s13 + s23
        # Assemble kernel
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None],
            'color_correlations' : [None],
            'values'             : {}
        })
        ker = 0
        # ker += 2*C123(z1, z2, z3, s12, s13, s23, s123)
        evaluation['values'][(0, 0)] = {'finite': 0.5*self.CF*self.TR*ker}
        return evaluation

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, lower_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):

        if higher_PS_point is None or lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping." )
        if leg_numbers_map is None:
            raise CurrentImplementationError(
                self.name() + " requires a leg numbers map, i.e. a momentum dictionary." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Include the counterterm only in a part of the phase space
        children = self.get_sorted_children(current, self.model)
        parent = leg_numbers_map.inv[frozenset(children)]
        pC = sum(higher_PS_point[child] for child in children)
        qC = lower_PS_point[parent]
        if self.is_cut(Q=Q, pC=pC):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Evaluate collinear subtracted kernel
        zs, kTs = self.variables(higher_PS_point, qC, children, Q=Q)
        pr = higher_PS_point[children[0]]
        ps = higher_PS_point[children[1]]
        pi = higher_PS_point[children[2]]
        sir = 2*pi.dot(pr)
        sis = 2*pi.dot(ps)
        srs = 2*pr.dot(ps)
        evaluation = self.evaluate_kernel(zs, kTs, parent)
        ker = 0
        misc.sprint(srs,sir,sis)
        ker += 2*self.C123_kernel(zs[0], zs[1], zs[2], srs, sir, sis, sir+sis+srs)
        ker -= 2*self.C123S12_kernel(zs[0], zs[1], zs[2], srs, sir, sis, sir+sis+srs)
        ker -= self.C123C12_kernel(
            higher_PS_point, lower_PS_point[parent], children, Q=Q)
        ker += self.C123S12C12_kernel(
            higher_PS_point, lower_PS_point[parent], children, Q=Q)
        evaluation['values'][(0, 0)]['finite'] += 0.5*self.CF*self.TR*ker

        # Find all colored leg numbers except for the parent in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color') == 1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        color_correlation_index = 1

        # Now loop over the colored parton number pairs (parent, k)
        # and add the corresponding contributions to this current
        emitter = children[2]
        for k in all_colored_parton_numbers:
            spectator = k
            if k == parent:
                spectator = children[2]
            evaluation['color_correlations'].append(((parent, k),))
            weight = 0
            if k != parent:
                pk = higher_PS_point[k]
                sik = 2*pi.dot(pk)
                skr = 2*pk.dot(pr)
                sks = 2*pk.dot(ps)
                weight += self.S12_kernel(sir, sis, sik, skr, sks, srs)
                weight -= 2*self.S12C12_kernel(
                    higher_PS_point, children, emitter, spectator, Q=Q)
            evaluation['values'][(0, color_correlation_index)] = {'finite': weight}
            color_correlation_index += 1

        # Add the normalization factors
        norm = (8. * math.pi * alpha_s) ** 2
        norm *= self.factor(Q=Q, pC=pC, qC=qC)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result
