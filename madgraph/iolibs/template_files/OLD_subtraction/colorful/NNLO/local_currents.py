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

import madgraph.various.misc as misc

import commons.utils as utils
import commons.QCD_local_currents as currents

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

def tijk(i, j, k, zs, kTs, s_ij=None, s_ik=None, s_jk=None):
    """Compute the invariant t_{ij,k} from zs and kTs."""

    if not s_ij: s_ij = sij(i, j, zs, kTs)
    if not s_ik: s_ik = sij(i, k, zs, kTs)
    if not s_jk: s_jk = sij(j, k, zs, kTs)
    return (2*(zs[i-1]*s_jk-zs[j-1]*s_ik)+(zs[i-1]-zs[j-1])*s_ij)/(zs[i-1]+zs[j-1])

#=========================================================================================
# NNLO final-collinear currents
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

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables and compute basic quantities
        z1, z2, z3 = zs
        s12 = sij(1, 2, zs, kTs)
        s13 = sij(1, 3, zs, kTs)
        s23 = sij(2, 3, zs, kTs)
        s123 = s12 + s13 + s23
        t123 = tijk(1, 2, 3, zs, kTs, s_ij=s12, s_ik=s13, s_jk=s23)
        # Assemble kernel
        sqrbrk  = -(t123 ** 2)/(s12*s123)
        sqrbrk += (4*z3 + (z1-z2)**2) / (z1+z2)
        sqrbrk += z1 + z2 - s12/s123
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None}}
        })
        evaluation['values'][(0, 0)]['finite'] = self.CF*self.TR * s123 / (2*s12) * sqrbrk
        return evaluation



class QCD_final_collinear_0_qqxq(currents.QCDLocalCollinearCurrent):
    """q q~ q (same flavour) collinear tree-level current."""

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NNLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 4, 0)  # gs^4 and 0 loops
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get(
            'singular_structure')  #: C((0,21,F),(0,4,F)) -> F: fijal state , 21: gluon, 4: c  going collinear
        # Check that there are 3 massless final state quarks or antiquarks
        if len(ss.legs) != 3: return None  # better have 3 legs
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
            if not cls.is_quark(leg, model): return None  # better be a (anti)quark
        # Look for a quark/antiquark pair
        pair = None
        for i in range(len(ss.legs)):
            for j in range(i + 1, len(ss.legs)):
                if cls.are_antiparticles(ss.legs[i], ss.legs[j]):
                    pair = (ss.legs[i], ss.legs[j])
                    continue
            if pair is not None: continue
        if pair is None: return None
        # Identify the remaining quark
        # remaining quark only appears in other_quarks if it has a different flavor than the q-qbar pair
        other_quarks = [leg for leg in ss.legs if leg not in pair]
        # we want same flavour, so we want the other_quarks to be empty
        if len(other_quarks) != 0: return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):
        # sorts out particles within current so, for example in q q' Q', q is always number 1
        legs = current.get('singular_structure').legs
        # we have a q q qx current, or a q qx qx current and we want to spot the "special" qx, or q respectively
        if cls.are_antiparticles(legs[0], legs[1]):  # 0 and 1 are anti
            if cls.are_antiparticles(legs[0], legs[2]):  # 0 and 2 are anti, so 0 is special
                return (legs[0].n, legs[1].n, legs[2].n)
            else:  # 0 and 2 are same, so 1 is special
                return (legs[1].n, legs[2].n, legs[0].n)
        else:  # 0 and 1 are same , so 2 is special
            return (legs[2].n, legs[0].n, legs[1].n)

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables and compute basic quantities
        z1, z2, z3 = zs
        k1, k2, k3 = kTs
        s12 = sij(1, 2, zs, kTs)
        s13 = sij(1, 3, zs, kTs)
        s23 = sij(2, 3, zs, kTs)
        s123 = s12 + s13 + s23

        # Assemble kernel
        kernel = self.evaluate_unsymmetrized_kernel((z1, z2, z3), (s12, s13, s23, s123), kTs)
        kernel += self.evaluate_unsymmetrized_kernel((z1, z3, z2), (s13, s12, s23, s123), (k1, k3, k2))

        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [None],
            'values': {(0, 0): {'finite': None}}
        })
        evaluation['values'][(0, 0)]['finite'] = kernel
        return evaluation

    def evaluate_unsymmetrized_kernel(self, zs, sijs, k_ts):
        z1, z2, z3 = zs
        s12, s13, s23, s123 = sijs
        t123 = tijk(1, 2, 3, zs, k_ts, s_ij=s12, s_ik=s13, s_jk=s23)
        # Assemble unsymmetrized kernel - piece identical to the kernel q->q q' qbar'
        sqrbrk = -(t123 ** 2) / (s12 * s123)
        sqrbrk += (4 * z3 + (z1 - z2) ** 2) / (z1 + z2)
        sqrbrk += z1 + z2 - s12 / s123
        # Assemble unsymmetrized kernel - piece special to same flavour case
        res = 2.0 * s23 / s12
        res += s123 / s12 * ((1 + z1 ** 2) / (1 - z2) - 2. * z2 / (1 - z3))
        res += - s123 ** 2 / s12 / s13 * z1 / 2. * (1 + z1 ** 2) / (1 - z2) / (1 - z3)
        return 0.5 * self.CF * self.TR * s123 / s12 * sqrbrk + self.CF * (self.CF - self.CA / 2.) * res


class QCD_final_collinear_0_qgg(currents.QCDLocalCollinearCurrent):
    """q g g  collinear tree-level current."""

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NNLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 4, 0)  # gs^4 and 0 loops
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get(
            'singular_structure')  #: C((0,21,F),(0,4,F)) -> F: fijal state , 21: gluon, 4: c  going collinear
        # Check that there are 3 massless final state partons
        if len(ss.legs) != 3: return None  # better have 3 legs
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
        # Check that there are two gluons and one (anti)quark
        number_of_quarks = 0
        number_of_gluons = 0
        for leg in ss.legs:
            if cls.is_quark(leg,model):number_of_quarks += 1
            if cls.is_gluon(leg,model):number_of_gluons += 1
        if number_of_quarks != 1: return None
        if number_of_gluons != 2: return None

        # The current is valid if we've made it till here
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):
        # sorts out particles within current so, for example in q q' Q', q is always number 1
        legs = current.get('singular_structure').legs
        # we have a q g g current  and we want to spot the "special" q (and place it at position 3)
        if cls.is_quark(legs[0], model):  # 0 is quark
            return (legs[1].n, legs[2].n, legs[0].n)
        elif cls.is_quark(legs[1], model): # 1 is quark
            return (legs[2].n, legs[0].n, legs[1].n)
        else:  # it must be that 2 is a quark
            return (legs[0].n, legs[1].n, legs[2].n)

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables and compute basic quantities
        z1, z2, z3 = zs
        k1, k2, k3 = kTs
        s12 = sij(1, 2, zs, kTs)
        s13 = sij(1, 3, zs, kTs)
        s23 = sij(2, 3, zs, kTs)
        s123 = s12 + s13 + s23

        # Assemble kernel
        kernel = self.evaluate_unsymmetrized_kernel((z1, z2, z3), (s12, s13, s23, s123), kTs)
        kernel += self.evaluate_unsymmetrized_kernel((z2, z1, z3), (s12, s23, s13, s123), (k2, k1, k3))

        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [None],
            'values': {(0, 0): {'finite': None}}
        })
        evaluation['values'][(0, 0)]['finite'] = kernel
        return evaluation

    def evaluate_unsymmetrized_kernel(self, zs, sijs, k_ts):
        z1, z2, z3 = zs
        s12, s13, s23, s123 = sijs
        t123 = tijk(1, 2, 3, zs, k_ts, s_ij=s12, s_ik=s13, s_jk=s23)
        # Assemble unsymmetrized kernel - Abelian piece
        abelian = s123**2 / (2 * s13 * s23) * z3 * (1+z3**2) / (z1 * z2)
        abelian += s123/s13 * (z3*(1-z1)+(1-z2)**3) / (z1 * z2)
        abelian += - s23/s13
        # Assemble unsymmetrized kernel - Non-abelian piece
        nonabelian = t123**2 / (4 * s12**2) + 1. / 4. + s123**2 / (2. * s12 * s13) * (  ((1-z3)**2 + 2 * z3) / z2 + (z2**2 + 2 * (1. - z2)) / (1-z3))
        nonabelian += - s123**2 / (4 * s13 * s23) * z3 * ((1-z3)**2 + 2 * z3 ) / (z1 * z2)
        nonabelian += s123 / (2 * s12) * ( (z1 * (2 - 2 * z1 + z1**2) - z2 * (6 -6*z2 + z2**2)) / (z2 * (1-z3))   )
        nonabelian += s123 / (2 * s13) * (  ((1-z2)**3 + z3**2 - z2) / (z2 * (1-z3))     - ( z3 * (1-z1) + (1-z2)**3) / (z1 * z2)  )
        return  self.CF**2 * abelian + self.CF *  self.CA  * nonabelian


class QCD_final_collinear_0_gqqx(currents.QCDLocalCollinearCurrent):
    """g -> g q q~  collinear tree-level current."""

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NNLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 4, 0)  # gs^4 and 0 loops
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get(
            'singular_structure')  #: C((0,21,F),(0,4,F)) -> F: fijal state , 21: gluon, 4: c  going collinear
        # Check that there are 3 massless final state partons
        if len(ss.legs) != 3: return None  # better have 3 legs
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
        # Check that there are two gluons and one (anti)quark
        number_of_quarks = 0
        number_of_gluons = 0
        for leg in ss.legs:
            if cls.is_quark(leg,model):number_of_quarks += 1
            if cls.is_gluon(leg,model):number_of_gluons += 1
        if number_of_quarks != 2: return None
        if number_of_gluons != 1: return None
        # Look for a quark/antiquark pair
        pair = None
        for i in range(len(ss.legs)):
            for j in range(i + 1, len(ss.legs)):
                if cls.are_antiparticles(ss.legs[i], ss.legs[j]):
                    pair = (ss.legs[i], ss.legs[j])
                    continue
            if pair is not None: continue
        if pair is None: return None
        # The current is valid if we've made it till here
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):
        # sorts out particles within current so, for example in g q Q, g is always number 1
        #we don;t care about the order of q and Q
        legs = current.get('singular_structure').legs
        # we have a g q Q current  and we want to spot the "special" g (and place it at position 1)
        if cls.is_gluon(legs[0], model):  # 0 is gluon
            return (legs[0].n, legs[1].n, legs[2].n)
        elif cls.is_gluon(legs[1], model): # 1 is gluon
            return (legs[1].n, legs[2].n, legs[0].n)
        else:  # it must be that 2 is a gluon
            return (legs[2].n, legs[0].n, legs[1].n)

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables and compute basic quantities
        z1, z2, z3 = zs
        k1, k2, k3 = kTs
        s12 = sij(1, 2, zs, kTs)
        s13 = sij(1, 3, zs, kTs)
        s23 = sij(2, 3, zs, kTs)
        s123 = s12 + s13 + s23
        CF = self.CF
        CA = self.CA
        TR = self.TR
        t231 = tijk(2, 3, 1, zs, kTs, s_ij=s23, s_ik=s12, s_jk=s13)
        t321 = tijk(3, 2, 1, zs, kTs, s_ij=s23, s_ik=s13, s_jk=s12)

        k23 = k2 + k3
        u23 = k2/z2 - k3/z3

        # Assemble kernel
        kernel_gmunu = (CF*TR*(2 - s23**2/(s12*s13) - s123**2/(s12*s13))
                        + CA*TR*(0.5 + s123**2/(s12*s13) + t231**2/(4.*s23**2)
                                 + t321**2/(4.*s23**2) + s123/(s23*(1 - z1))
                                 - s123/(2.*s12*(1 - z1)*z1) - s123/(2.*s13*(1 - z1)*z1)
                                 - s123/(s23*(1 - z1)*z1) - (2*s123*z1)/(s23*(1 - z1))
                                 + (s123**2*z2)/(s12*s23*(1 - z1))
                                 + (s123*z2)/(2.*s13*(1 - z1)*z1)
                                 - (s123**2*z2)/(2.*s12*s23*(1 - z1)*z1)
                                 + (s123**2*z3)/(s13*s23*(1 - z1))
                                 + (s123*z3)/(2.*s12*(1 - z1)*z1)
                                 - (s123**2*z3)/(2.*s13*s23*(1 - z1)*z1))
                        )
        kernel_k1mu_k1nu = ((2*CA*TR*s123)/(s12*s13) - (4*CF*TR*s123)/(s12*s13))
        kernel_k2mu_k2nu = ((-4*CF*TR*s123)/(s12*s13)
                            + CA*TR*(((4*s123)/(s12*s13)
                                      + (s123*((-16*z3**2)/((1 - z1)*z1)
                                               - 4*(1 + (2*(-z1 + z2)*z3)/((1 - z1)*z1))))/(s13*s23))/4.
                                     + ((4*s123)/(s12*s13)
                                        + (s123*(8 - 4*(1 + (2*z2*(-z1 + z3))/((1 - z1)*z1))))/(s12*s23))/4.)
                            )
        kernel_k3mu_k3nu = ((-4*CF*TR*s123)/(s12*s13)
                            + CA*TR*(((4*s123)/(s12*s13)
                                      + (s123*(8 - 4*(1 + (2*(-z1 + z2)*z3)/((1 - z1)*z1))))/(s13*s23))/4.
                                     + ((4*s123)/(s12*s13)
                                        + (s123*((-16*z2**2)/((1 - z1)*z1)
                                                 - 4*(1 + (2*z2*(-z1 + z3))/((1 - z1)*z1))))/(s12*s23))/4.)
                            )
        kernel_k23mu_k23nu = ((4*CF*TR*s123)/(s12*s13)
                              + CA*TR*(((-4*s123)/(s12*s13)
                                        + (4*s123*(1 + (2*(-z1 + z2)*z3)/((1 - z1)*z1)))/(s13*s23))/4.
                                       + ((-4*s123)/(s12*s13)
                                          + (4*s123*(1 + (2*z2*(-z1 + z3))/((1 - z1)*z1)))/(s12*s23))/4.)
                              )

        kernel_u23mu_u23nu = ((-8*CA*TR*s123*z2**2*z3**2)/(s23**2*(1 - z1)*z1))

        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None, ((parent,( k1, )), ),
                                        ((parent,( k2, )), ),
                                        ((parent,( k3, )), ),
                                        ((parent,( k23, )), ),
                                        ((parent,( u23, )), ), ],
            'color_correlations': [None],
            'values': {(0, 0): {'finite': -kernel_gmunu},  # for the minus sign: Madgraph multiplies this by -g_mu_nu
                       (1, 0): {'finite': kernel_k1mu_k1nu},
                       (2, 0): {'finite': kernel_k2mu_k2nu},
                       (3, 0): {'finite': kernel_k3mu_k3nu},
                       (4, 0): {'finite': kernel_k23mu_k23nu},
                       (5, 0): {'finite': kernel_u23mu_u23nu},}
        })

        return evaluation


class QCD_final_collinear_0_ggg(currents.QCDLocalCollinearCurrent):
    """g -> g g g  collinear tree-level current."""

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NNLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 4, 0)  # gs^4 and 0 loops
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get(
            'singular_structure')  #: C((0,21,F),(0,4,F)) -> F: fijal state , 21: gluon, 4: c  going collinear
        # Check that there are 3 massless final state partons
        if len(ss.legs) != 3: return None  # better have 3 legs
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
        # Check that there are no (anti)quarks and that there are 3 gluons
        number_of_gluons = 0
        for leg in ss.legs:
            if cls.is_quark(leg, model): return None
            if cls.is_gluon(leg, model): number_of_gluons = number_of_gluons + 1
        if number_of_gluons != 3: return None

        # The current is valid if we've made it till here
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):
        # sorts out particles within current
        # we don't care about ordering here - the g g g vertex is symmetric
        legs =  current.get('singular_structure').legs
        return (legs[0].n, legs[1].n, legs[2].n)

    def ttijk(self,i, j, k, zs, s_ij, s_ik, s_jk):
        """Compute the invariant t_{ij,k} from zs and sij."""
        return (2 * (zs[i - 1] * s_jk - zs[j - 1] * s_ik) + (zs[i - 1] - zs[j - 1]) * s_ij) / (zs[i - 1] + zs[j - 1])

    def evaluate_kernel(self, zs, kTs, parent):
        #misc.sprint("HELLO")

        # Retrieve the collinear variables and compute basic quantities
        z1, z2, z3 = zs
        k1, k2, k3 = kTs
        s12 = sij(1, 2, zs, kTs)
        s13 = sij(1, 3, zs, kTs)
        s23 = sij(2, 3, zs, kTs)
        s123 = s12 + s13 + s23
        CA = self.CA

        k12 = k1 + k2
        k23 = k2 + k3
        k13 = k1 + k3
        u12 = k1 / z1 - k2 / z2
        u23 = k2 / z2 - k3 / z3
        u31 = k1 / z1 - k3 / z3

        # Assemble kernel
        kernel_gmunu = self.evaluate_permuted_gmunu(z1,z2,z3,s12,s13,s23)

        kernel_u12mu_u12nu = (8*CA**2*s123*z1**2*z2**2)/(s12**2*(1 - z3)*z3)

        kernel_u23mu_u23nu = (8*CA**2*s123*z2**2*z3**2)/(s23**2*(1 - z1)*z1)

        kernel_u13mu_u13nu = (8*CA**2*s123*z1**2*z3**2)/(s13**2*(1 - z2)*z2)

        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None, ((parent, (k1,)),),
                                  ((parent, (k2,)),),
                                  ((parent, (k3,)),),
                                  ((parent, (u12,)),),
                                  ((parent, (u23,)),),
                                  ((parent, (u31,)),),
                                  ],
            'color_correlations': [None],
            'values': {(0, 0): {'finite': -kernel_gmunu},  # for the minus sign: Madgraph multiplies this by -g_mu_nu
                       (1, 0): {'finite': self.evaluatecoefk1k1(z1,z2,z3,s12,s13,s23)},
                       (2, 0): {'finite': self.evaluatecoefk2k2(z1,z2,z3,s12,s13,s23)},
                       (3, 0): {'finite': self.evaluatecoefk3k3(z1,z2,z3,s12,s13,s23)},
                       (4, 0): {'finite': kernel_u12mu_u12nu},
                       (5, 0): {'finite': kernel_u23mu_u23nu},
                       (6, 0): {'finite': kernel_u13mu_u13nu},
                       }
        })

        return evaluation

    def evaluate_permuted_gmunu(self,z1,z2,z3,s12,s13,s23):
        s123 = s12 + s13 + s23
        zs = (z1,z2,z3)
        tt123 = self.ttijk(1, 2, 3, zs, s12, s13, s23)
        tt312 = self.ttijk(3, 1, 2, zs, s13, s23, s12)
        tt231 = self.ttijk(2, 3, 1, zs, s23, s12, s13)
        tt213 = - tt123
        tt132 = - tt312
        tt321 = - tt231
        return (          self.evaluate_gmunu(z1, z2, z3, s12, s13, s23, s123, tt123)
                        + self.evaluate_gmunu(z2, z3, z1, s23, s12, s13, s123, tt231)
                        + self.evaluate_gmunu(z3, z1, z2, s13, s23, s12, s123, tt312)
                        + self.evaluate_gmunu(z3, z2, z1, s23, s13, s12, s123, tt321)
                        + self.evaluate_gmunu(z1, z3, z2, s13, s12, s23, s123, tt132)
                        + self.evaluate_gmunu(z2, z1, z3, s12, s23, s13, s123, tt213)
                        )

    def evaluate_gmunu(self, z1, z2, z3, s12, s13, s23,s123,tt123):
        return self.CA**2*(-0.75 - tt123**2/(4.*s12**2)
                           + (2*s123)/(s12*(1 - z1)*z1)
                           - (2*s123)/(s12*(1 - z3))
                           - s123**2/(2.*s12*s13*(1 - z2)*(1 - z3))
                           + (s123**2*z1)/(s12*s13*(1 - z2)*(1 - z3))
                           - (s123**2*z1**2)/(s12*s13*(1 - z2)*(1 - z3))
                           - s123/(s12*(1 - z1)*z1*z3)
                           - s123**2/(2.*s12*s13*z2*z3)
                           + (s123**2*z1)/(s12*s13*z2*z3)
                           - (s123**2*z1**2)/(s12*s13*z2*z3)
                           + (2*s123)/(s12*(1 - z3)*z3)
                           - (2*s123*z3)/(s12*(1 - z1)*z1)
                           + (4*s123*z3)/(s12*(1 - z3))
                           + (2*s123**2*z2*z3)/(s12*s13*(1 - z2)*(1 - z3))
                           )


    def evaluatecoefk1k1(self,z1,z2,z3,s12,s13,s23):
        s123 = s12 + s23 + s13
        return self.CA**2*s123*((-3 + (2*(1 - z2)*z2)/((1 - z3)*z3))/(s12*s13)
                                + (3 - (2*(1 - z1)*z1)/((1 - z3)*z3) + (2*z2*(1 - 2*z3))/((1 - z3)*z3))/(s12*s23)
                                + (3 - (2*(1 - z1)*z1)/((1 - z2)*z2) + (2*(1 - 2*z2)*z3)/((1 - z2)*z2))/(s13*s23)
                                + (3 - (2*(1 - z2)*z2)/((1 - z1)*z1) + (2*(1 - 2*z2)*z3)/((1 - z2)*z2))/(s13*s23)
                                + (3 + (2*z2*(1 - 2*z3))/((1 - z3)*z3) - (2*(1 - z3)*z3)/((1 - z1)*z1))/(s12*s23)
                                + (-3 + (2*(1 - z3)*z3)/((1 - z2)*z2))/(s12*s13)
                                )

    def evaluatecoefk2k2(self,z1,z2,z3,s12,s13,s23):
        s123 = s12 + s23 + s13
        return self.CA**2*s123*((-3 + (2*(1 - z1)*z1)/((1 - z3)*z3))/(s12*s23)
                                + (3 - (2*(1 - z2)*z2)/((1 - z3)*z3) + (2*z1*(1 - 2*z3))/((1 - z3)*z3))/(s12*s13)
                                + (3 - (2*(1 - z1)*z1)/((1 - z2)*z2) + (2*(1 - 2*z1)*z3)/((1 - z1)*z1))/(s13*s23)
                                + (3 - (2*(1 - z2)*z2)/((1 - z1)*z1) + (2*(1 - 2*z1)*z3)/((1 - z1)*z1))/(s13*s23)
                                + (-3 + (2*(1 - z3)*z3)/((1 - z1)*z1))/(s12*s23)
                                + (3 + (2*z1*(1 - 2*z3))/((1 - z3)*z3)
                                   - (2*(1 - z3)*z3)/((1 - z2)*z2))/(s12*s13)
                                )

    def evaluatecoefk3k3(self,z1,z2,z3,s12,s13,s23):
        s123 = s12 + s23 + s13
        return self.CA**2*s123*((-3 + (2*(1 - z1)*z1)/((1 - z2)*z2))/(s13*s23)
                                + (-3 + (2*(1 - z2)*z2)/((1 - z1)*z1))/(s13*s23)
                                + (3 + (2*(1 - 2*z1)*z2)/((1 - z1)*z1) - (2*(1 - z1)*z1)/((1 - z3)*z3))/(s12*s23)
                                + (3 + (2*z1*(1 - 2*z2))/((1 - z2)*z2) - (2*(1 - z2)*z2)/((1 - z3)*z3))/(s12*s13)
                                + (3 + (2*(1 - 2*z1)*z2)/((1 - z1)*z1) - (2*(1 - z3)*z3)/((1 - z1)*z1))/(s12*s23)
                                + (3 + (2*z1*(1 - 2*z2))/((1 - z2)*z2) - (2*(1 - z3)*z3)/((1 - z2)*z2))/(s12*s13)
                                )




#=========================================================================================
# NNLO soft currents
#=========================================================================================

class QCD_final_soft_0_gg(currents.QCDLocalSoftCurrent):
    """Double soft tree-level current."""
    
    #TODO Implement meaningful cut and factors for the double soft
    # The single-soft versions below are not expected to work directly
    # THIS IS PROBABLY CRUCIAL IN ORDER TO HAVE THE EXPECTED CANCELLATION
    #    S(5) + S(5)S(6) + S(5,6) = 0 in the limit S(6)
    ###is_cut = staticmethod(currents.SomogyiChoices.cut_soft)
    ###factor = staticmethod(currents.SomogyiChoices.factor_soft)
    
    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NNLO QCD soft tree-level current
        init_vars = cls.common_does_implement_this_current(current, QCD_squared_order=4, n_loops=0)
        if init_vars is None: return None
        
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that there are 2 massless gluons
        if len(ss.legs) != 2: return None
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
            if not cls.is_gluon(leg, model): return None

        # The current is valid
        return init_vars
   
    @staticmethod
    def eikonal(pi, pj, pr):
        """Eikonal factor for soft particle with number 'r'
        emitted from 'i' and reconnecting to 'j'.
        """

        pipj = pi.dot(pj)
        pipr = pi.dot(pr)
        pjpr = pj.dot(pr)
        return pipj/(pipr*pjpr)

    @staticmethod
    def strongly_ordered_double_soft_kernel(pi, pj, pr1, pr2):
        """Strongly ordered piece of the double soft current
        with first soft particle with momentum 'pr1' and second with momentum 'pr2'
        with emitting dipole of momenta 'pi' and 'pj'.
        See Eq.(109-111) of hep-ph/9908523v1 (Catani-Grazzini).
        """

        pipj = pi.dot(pj)
        pr1pr2 = pr1.dot(pr2)
        pipr1 = pi.dot(pr1)
        pipr2 = pi.dot(pr2)
        pjpr1 = pj.dot(pr1)
        pjpr2 = pj.dot(pr2)
        
        return ( 
            (pipj/pr1pr2) * ( 1./(pipr1*pjpr2) + 1./(pjpr1*pipr2) )
            -  (pipj)**2 / (pipr1*pjpr2*pipr2*pjpr1)
        )

    @staticmethod
    def non_abelian_eikonal(pi, pj, pr1, pr2):
        """Non abelian piece of the double soft current
        with first soft particle with momentum 'pr1' and second with momentum 'pr2'
        with emitting dipole of momenta 'pi' and 'pj'.
        See Eq.(109-111) of hep-ph/9908523v1 (Catani-Grazzini).
        """

        pipj = pi.dot(pj)
        pr1pr2 = pr1.dot(pr2)
        pipr1 = pi.dot(pr1)
        pipr2 = pi.dot(pr2)
        pjpr1 = pj.dot(pr1)
        pjpr2 = pj.dot(pr2)

        so_double_soft = QCD_final_soft_0_gg.strongly_ordered_double_soft_kernel(
            pi, pj, pr1, pr2)
        
        pure_double_soft_non_abelian = (
            so_double_soft + ((pipr1*pjpr2+pipr2*pjpr1)/((pipr1+pipr2)*(pjpr1+pjpr2)))*(
                (1./(pr1pr2**2)) - 0.5*so_double_soft
            ) - ((2.*pipj)/(pr1pr2*(pipr1+pipr2)*(pjpr1+pjpr2)))
        )
        
        # We want to absorb here also the non-abelian strongly ordered part of the 
        # two-iterated soft counterterms (i.e. S(r) \otimes S(s) and S(s) \otimes S(r)) 
        # which are each '-so_double_soft' when using the current normalisation of 
        # this non-abelian eikonal term.
        return pure_double_soft_non_abelian - 2.*so_double_soft
   
    def create_CataniGrazzini_correlator(self, (i,j),(k,l)):
        """ Returns the correlator of Catani-Grazzini (Eq.113 of hep-ph/9908523v1)
                <M| ( T^-1_i \dot T^-1_j ) * ( T^-1_k \dot T^-1_l ) | M > 
                
            converted into MadGraph's conventions:
            
              ( (a,-1,a),(b,-2,b) ) , ( (c,-1,c),(d,-2,d) ) --> T^-1_a T^-2_b T^-1_c T^-2_d
            """

        # It is important to never commute two color operators acting on the same index, so we must chose carefully which
        # index to pick to carry the gluon index '-2' of the first connection. This can be either 'k' or 'l'.
        if j!=k and j!=l:
            # If all indices are different, we can pick either k or l, it is irrelevant
            index1, index2, index3, index4 = i, k, j, l
        elif j==k and j!=l:
            # If j is equal to k, we must pick l
            index1, index2, index3, index4 = i, l, j, k
        elif j==l and j!=k:
            # If j is equal to l, we must pick k
            index1, index2, index3, index4 = i, k, j, l
        elif j==l and j==k:
            # If j is equal to both l and k, then agin it doesn't matter and we can pick k
            index1, index2, index3, index4 = i, k, j, l

        # The sorting according to the first index of each tuple of each of the two convention is to match
        # Madgraph's convention for sorting color connection in the color correlators definition
        return (
            tuple(sorted([ (index1,-1,index1), (index2,-2,index2) ], key = lambda el: el[0], reverse=True)),
            tuple(sorted([ (index3,-1,index3), (index4,-2,index4) ], key = lambda el: el[0], reverse=True))
        )

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
        mu_r    = model_param_dict['MU_R']

        # Now find all colored leg numbers in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color')==1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        # Identify the soft leg numbers and momenta
        soft_leg_number_A = current.get('singular_structure').legs[0].n
        soft_leg_number_B = current.get('singular_structure').legs[1].n
        pA = higher_PS_point[soft_leg_number_A]
        pB = higher_PS_point[soft_leg_number_B]
        pS = pA + pB

        # Include the counterterm only in a part of the phase space
        if self.is_cut(Q=Q, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [],
            'values'              : {}
        })
        
        # Normalization factors
        couplings_factors = 4. * math.pi * alpha_s
        norm = couplings_factors**2 * self.factor(Q=Q, pS=pS)

        # Keep track of the color correlators added
        color_correlators_added = {}
        color_correlation_max_index = 0
        # Now loop over the colored parton numbers to implement the squared double-soft
        # current flow. Note that significant improvement can be obtained by taking advantage
        # of the symmetries, as well as writing fully expanded hard-coded expressions with
        # all dot-products cached.
        # For efficiency, we choose here to define coefficients of each single product of
        # correlators. So, typically in the abelian piece, the color correlator reads:
        #       ( '{}' denotes an anti-commutator below and '.' denotes a dot-product )
        #    { T_i . T_j, T_k . T_l }
        # And we will register *two* entries:
        #    evaluation['color_correlations'].append(
        #      ( 
        #         ( ( (i,-1,i), (k,-2,k) ), ( (j,-1,j), (l,-2,l) ) ),
        #      )
        #    )
        # And:
        #    evaluation['color_correlations'].append(
        #      (     
        #         ( ( (i,-2,i), (k,-1,k) ), ( (j,-2,j), (l,-1,l) ) ),
        #      )
        #    )
        # As opposed to directly defining their sum:
        #    evaluation['color_correlations'].append(
        #            ( 
        #                    ( ( (i,-1,i), (k,-2,k) ), ( (j,-1,j), (l,-2,l) ) ),
        #                    ( ( (i,-2,i), (k,-1,k) ), ( (j,-2,j), (l,-1,l) ) )
        #            )
        #    )
        for i in all_colored_parton_numbers:
            for j in all_colored_parton_numbers:

                # Compute the non-abelian eikonal
                pi = sum(higher_PS_point[child] for child in leg_numbers_map[i])
                pj = sum(higher_PS_point[child] for child in leg_numbers_map[j])
                # pi = lower_PS_point[i]
                # pj = lower_PS_point[j]
                non_abelian_eikonal = self.non_abelian_eikonal(pi, pj, pA, pB)
                non_abelian_kernel = - self.CA * norm * non_abelian_eikonal
                
                # Implement the non-abelian piece
                non_abelian_correlator = ( 
                    ( ( (i,-1,i), ), ( (j,-1,j), ) ), 
                )
                if non_abelian_correlator in color_correlators_added:
                    color_correlation_index = color_correlators_added[non_abelian_correlator]
                    evaluation['values'][(0, color_correlation_index)]['finite'] += non_abelian_kernel
                else:
                    evaluation['color_correlations'].append(non_abelian_correlator)
                    color_correlation_index = color_correlation_max_index
                    color_correlators_added[non_abelian_correlator] = color_correlation_max_index
                    color_correlation_max_index += 1
                    evaluation['values'][(0, color_correlation_index)] = {
                                                             'finite': non_abelian_kernel }
    
                for k in all_colored_parton_numbers:
                    for l in all_colored_parton_numbers:

                        # Compute the abelian eikonal
                        pk = sum(higher_PS_point[child] for child in leg_numbers_map[k])
                        pl = sum(higher_PS_point[child] for child in leg_numbers_map[l])
                        # pk = lower_PS_point[k]
                        # pl = lower_PS_point[l]
                        eik_ij = self.eikonal(pi, pj, pA)
                        eik_kl = self.eikonal(pk, pl, pB)
                        abelian_kernel = 0.5 * norm * eik_ij * eik_kl
                        
                        # Implement the abelian piece
                        abelian_correlator_A = ( self.create_CataniGrazzini_correlator((i,j),(k,l)), )
                        abelian_correlator_B = ( self.create_CataniGrazzini_correlator((k,l),(i,j)), )

                        for correlator in [abelian_correlator_A, abelian_correlator_B]:
                            if correlator in color_correlators_added:
                                color_correlation_index = color_correlators_added[correlator]
                                #misc.sprint('Adding %f ((%d,%d,%d)->%f, (%d,%d,%d)->%f) to CC: %d, %s'%\
                                #            (abelian_kernel,
                                #             i,j,soft_leg_number_A,self.eikonal(PS_point, i, j, soft_leg_number_A),
                                #             k,l,soft_leg_number_B,self.eikonal(PS_point, k, l, soft_leg_number_B),
                                #             color_correlation_index, str(correlator)
                                #            ))
                                evaluation['values'][(0, color_correlation_index)]['finite'] += abelian_kernel
                            else:
                                evaluation['color_correlations'].append(correlator)
                                color_correlation_index = color_correlation_max_index
                                color_correlators_added[correlator] = color_correlation_max_index
                                color_correlation_max_index += 1
                                evaluation['values'][(0, color_correlation_index)] = {
                                                                 'finite': abelian_kernel }
                                #misc.sprint('Adding %f ((%d,%d,%d)->%f, (%d,%d,%d)->%f) to CC: %d, %s'%\
                                #            (abelian_kernel,
                                #             i,j,soft_leg_number_A,self.eikonal(PS_point, i, j, soft_leg_number_A),
                                #             k,l,soft_leg_number_B,self.eikonal(PS_point, k, l, soft_leg_number_B),
                                #             color_correlation_index, str(correlator)
                                #            ))
        
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())) )
        
        #misc.sprint('==BELOW LIST CC==')
        #for i, cc in enumerate(evaluation['color_correlations']):
        #    misc.sprint('Color correlator: %d : %s = %f'%(i, cc, evaluation['values'][(0,i)]['finite']))
        #misc.sprint('==ABOVE LIST CC==')

        return result
