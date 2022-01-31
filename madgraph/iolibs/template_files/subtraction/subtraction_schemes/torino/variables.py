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
import factors_and_cuts as fc

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

    y = sij / (sij + sir + sjr)    #for final kr
    x = 1. - sij / (sir + sjr)     #for initial kr

    return xi, xj, kitil, kjtil, y, x

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
    xi, xj, kitil, kjtil, y, x = get_sudakov_decomp_FF(all_p_fs[0], all_p_fs[1], p_rec)

    xis = (xi, xj)
    kTs = (kitil, kjtil)

    return [{'xs':xis, 'kTs':kTs, 'ss':ss_i_j, 'y':y, 'x':x},]

def get_sudakov_decomp_IF(kj, ki, kr, kr_leg_number):
    """returns the sudakov decomposition (x and kt) of ki and kj, with kr being the other
    reference momentum for the IF case. kj=initial, ki=final 
    Follows eq 1.18-1.2 of the nlo_v4 note"""

    sij = ki.dot(kj)
    sir = ki.dot(kr)
    sjr = kj.dot(kr)

    na, nb = mappings.InitialCollinearVariables.collinear_and_reference(kj)
    pa = kj
        # Compute the sum of momenta
    pA = LorentzVector(pa)
    pA -= ki
    # Pre-compute variables
    napA = na.dot(pA)
    nbpA = nb.dot(pA)
    nanb = na.dot(nb)
    nbpa = nb.dot(pa)
    # zA = nbpA / nbpa
    ktA = pA - (nbpA * na + napA * nb) / nanb
    # Compute all kinematic variables
    pi = ki
    napi = na.dot(pi)
    nbpi = nb.dot(pi)
    # zi = nbpi / nbpa
    kti = pi - (nbpi*na+napi*nb) / nanb
    # kts
    kjtil = ktA
    kitil = kti

    if kr_leg_number > 2:
        xj = 1. - sir / (sij+sjr)
        xi = 1. - xj
    else:
        xj = 1. - ((sij + sir) / sjr)
        xi = 1. - xj

    # print('VAR - xj : ' + str(xj))
    # print('VAR - xi : ' + str(xi))

    # Test on kt
    # print('VAR - kjtil     : ' + str(kjtil))
    # print('VAR - test kt   : ' + str(ki-kj*(sir/sjr)-kr*(sij/sjr)))
    # print('VAR - test kt 2 : ' + str(kjtil.dot(kr)))
    # print('VAR - test kt 3 : ' + str(kjtil.dot(kj)))


    # Variable for initial splitting - final rec
    z = sij / (sij + sjr)
    # Variables for initial splitting - initial rec
    v = sij / (sij + sir)

    # x_r_final = 1. - (sjr / (sij + sir))
    # x_r_initial = 1. - ((sij + sjr) / sir)
    # #x1 = 1. - (sjr / (sir - sij))    #for soft initial case
    # x1 = 1. - ((sij + sjr) / sir)    #for soft initial case

    
    return xj, xi, kjtil, kitil, z, v

    # return xj, xi, kjtil, kitil, x_r_final, x_r_initial, z, v, x1

def TRN_IFn_variables(higher_PS_point, qC, children, **opts):
    """
    returns variables for the torino subtraction scheme, in particular x's and
    ktil as returned from get_sudakov_decomp_IF

    """

    all_p_fs = [higher_PS_point[child] for child in children]
    Q = opts['Q']
    # the recoiler
    p_rec = higher_PS_point[opts['ids']['c']]
    p_rec_leg_number = opts['ids']['c']

    # print('VAR : all_p_fs : ' + str(all_p_fs))
    # print('VAR : Q        : ' + str(Q))
    # print('VAR : p_rec    : ' + str(p_rec))

    # Add additional ss's
    ss_i_j = {}
    for i_fs in range(len(all_p_fs)):
        for j_fs in range(i_fs+1,len(all_p_fs)):
            ss_i_j[(i_fs,j_fs)] = 2.*all_p_fs[i_fs].dot(all_p_fs[j_fs])
    ss_i_j_tot = sum(ss_i_j.values())
    for k,vv in ss_i_j.items():
        ss_i_j[(k[1],k[0])] = vv

    # now x's and kts where all_p_fs[0] = IN , all_p_fs[1] = F
    # xi, xj, kitil, kjtil, x_r_final, x_r_initial, z, v, x1 = get_sudakov_decomp_IF(all_p_fs[0], all_p_fs[1], p_rec, p_rec_leg_number)

    # xis = (xi, xj)
    # kTs = (kitil, kjtil)
    # x_mapping = (x_r_final, x_r_initial)
    # return [{'xs':xis, 'kTs':kTs, 'ss':ss_i_j, 'x_mapping':x_mapping, 'z':z, 'v':v, 'x1':x1},]

    xj, xi, kjtil, kitil, z, v = get_sudakov_decomp_IF(all_p_fs[0], all_p_fs[1], p_rec, p_rec_leg_number)

    xis = (xj, xi)
    kTs = (kjtil, kitil)
 
    return [{'xs':xis, 'kTs':kTs, 'ss':ss_i_j, 'z':z, 'v':v},]
