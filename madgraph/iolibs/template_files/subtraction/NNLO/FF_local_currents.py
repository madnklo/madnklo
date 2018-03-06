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
