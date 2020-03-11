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

import commons.utils as utils

CurrentImplementationError = utils.CurrentImplementationError


def get_sudakov_decomp(ki, kj, kr):
    """returns the sudakov decomposition (x and kt) of ki and kj, with kr being the other
    reference momentum.
    Follows eq 2.24-2.27 of the torino paper"""

    sij = ki.dot(kj)
    sir = ki.dot(kr)
    sjr = kj.dot(kr)

    k = ki + kj
    kbar = k - sij / (sir + sjr) * kr

    xi = sir / (sir + sjr)
    xj = sjr / (sir + sjr)

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
    xi, xj, kitil, kjtil = get_sudakov_decomp(all_p_fs[0], all_p_fs[1], p_rec)

    xis = (xi, xj)
    kTs = (kitil, kjtil)

    return [{'xs':xis, 'kTs':kTs, 'ss':ss_i_j},]
