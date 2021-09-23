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

import os
import math

import madgraph.integrator.vectors as vectors
import madgraph.various.misc as misc
import madgraph.integrator.mappings as mappings

from madgraph.integrator.vectors import Vector, LorentzVector

import commons.utils as utils

import logging



logger = logging.getLogger('madgraph')


CurrentImplementationError = utils.CurrentImplementationError


def get_sudakov_decomp_FF(ki, kj, kr):
    """returns the sudakov decomposition (x and kt) of ki and kj, with kr being the other
    reference momentum.
    Follows eq 2.24-2.27 of the torino paper"""

    sij = ki.dot(kj)
    sir = ki.dot(kr)
    sjr = kj.dot(kr)

    k = ki + kj
    kbar = k - sij / (sir + sjr) * kr

    xi = sir / (sir + sjr)
    #xj = sjr / (sir + sjr)
    xj = 1. - xi 

    kitil = ki - k * xi  - kr * (k.dot(ki) / k.square() - xi) * k.square() / k.dot(kr)
    kjtil = kj - k * xj  - kr * (k.dot(kj) / k.square() - xj) * k.square() / k.dot(kr)

    return xi, xj, kitil, kjtil

def TRN_FFn_variables(higher_PS_point, qC, children, **opts):
    """
    returns variables for the torino subtraction scheme, in particular x's and
    ktil as returned from get_sudakov_decomp

    """

    all_p_fs = [higher_PS_point[child] for child in children]
    Q = opts['Q']
    # the recoiler
    p_rec = higher_PS_point[opts['ids']['c']]

    # Add additional ss's
    ss_i_j = {}
    for i_fs in range(len(all_p_fs)):
        for j_fs in range(i_fs+1,len(all_p_fs)):
            ss_i_j[(i_fs,j_fs)] = 2.*all_p_fs[i_fs].dot(all_p_fs[j_fs])
    ss_i_j_tot = sum(ss_i_j.values())
    for k,v in ss_i_j.items():
        ss_i_j[(k[1],k[0])] = v

    # now x's and kts
    xi, xj, kitil, kjtil = get_sudakov_decomp_FF(all_p_fs[0], all_p_fs[1], p_rec)

    xis = (xi, xj)
    kTs = (kitil, kjtil)

    return [{'xs':xis, 'kTs':kTs, 'ss':ss_i_j},]

def get_sudakov_decomp_IF(ki, kj, kr):
    """returns the sudakov decomposition (x and kt) of ki and kj, with kr being the other
    reference momentum for the IF case. ki=initial, kj=final
    Follows eq 1.18-1.2 of the nlo_v4 note"""

    sij = ki.dot(kj)
    sir = ki.dot(kr)
    sjr = kj.dot(kr)

    #xi_old = (sir - sjr)/sir
    #xj = 1. - xi


    na, nb = mappings.InitialCollinearVariables.collinear_and_reference(ki)

    pa = ki
    # Compute the sum of momenta
    pA = LorentzVector(pa)
    pA -= kj
    # Pre-compute variables
    napA = na.dot(pA)
    nbpA = nb.dot(pA)
    nanb = na.dot(nb)
    nbpa = nb.dot(pa)
    zA = nbpA / nbpa
    ktA = pA - (nbpA * na + napA * nb) / nanb
    # ptA = pA - (nbpA * na + napA * nb) / nanb
    # ktA = ptA / zA
    # Initialize variables for sum rules check
    zsum = 0
    ktsum = LorentzVector()
    ktabssum = LorentzVector()
    # Fill in A data, using child number improperly
    zsum += zA
    ktsum += ktA
    for j in range(len(ktA)):
        ktabssum[j] += abs(ktA[j])
    # Compute all kinematic variables
    pi = kj
    napi = na.dot(pi)
    nbpi = nb.dot(pi)
    zi = nbpi / nbpa
    kti = pi - (nbpi*na+napi*nb) / nanb
    # pti = pi - (nbpi*na+napi*nb) / nanb
    # kti = pti - zi * ktA
    zsum += zi
    ktsum += kti
    for j in range(len(kti)):
        ktabssum[j] += abs(kti[j])

    #logger.info("kti_new = " + str(ktA))
    #logger.info("ktj_new = " + str(kti))

    #kjtil_old = kj - ki * (1. - xi_old)  - kr * (sij / sir)

    #logger.info("kjtil_old = " + str(kjtil_old))
    #logger.info("xi_old = " + str(xi_old))
    #logger.info("xi_new = " + str(zA))

    xi = zA
    xj = zi
    kitil = ktA
    kjtil = kti

    return xi, xj, kitil, kjtil

def TRN_IFn_variables(higher_PS_point, qC, children, **opts):
    """
    returns variables for the torino subtraction scheme, in particular x's and
    ktil as returned from get_sudakov_decomp_IF

    """
#    logger1.info('children: '+str(children))

    all_p_fs = [higher_PS_point[child] for child in children]
    Q = opts['Q']
    # the recoiler
    p_rec = higher_PS_point[opts['ids']['c']]


    # Add additional ss's
    ss_i_j = {}
    for i_fs in range(len(all_p_fs)):
        for j_fs in range(i_fs+1,len(all_p_fs)):
            ss_i_j[(i_fs,j_fs)] = 2.*all_p_fs[i_fs].dot(all_p_fs[j_fs])
    ss_i_j_tot = sum(ss_i_j.values())
    for k,v in ss_i_j.items():
        ss_i_j[(k[1],k[0])] = v

    # now x's and kts where all_p_fs[0] = IN , all_p_fs[1] = F
    xi, xj, kitil, kjtil = get_sudakov_decomp_IF(all_p_fs[0], all_p_fs[1], p_rec)

    xis = (xi, xj)
    kTs = (kitil, kjtil)


    return [{'xs':xis, 'kTs':kTs, 'ss':ss_i_j,},]
