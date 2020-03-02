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

    assert(ki != kj)
    assert(kr != kj)
    assert(ki != kr)

    sij = ki.dot(kj)
    sir = ki.dot(kr)
    sjr = kj.dot(kr)

    k = ki + kj
    kbar = k - sij / (sir + sjr) * kr

    xi = sir / (sir + sjr)
    xj = sjr / (sir + sjr)

    assert(xi+xj == 1)

    kitil = ki - k * xi  - kr * (k.dot(ki) / k.square() - xi) * k.square() / k.dot(kr)
    kjtil = kj - k * xj  - kr * (k.dot(kj) / k.square() - xj) * k.square() / k.dot(kr)

    assert(kitil == -kjtil)
    assert(kitil.dot(kbar) == kitil.dot(kr) == 0.)
    assert(kjtil.dot(kbar) == kjtil.dot(kr) == 0.)

    return xi, xj, kitil, kjtil

def TRN_FFn_variables(higher_PS_point, qC, children, **opts):
    """MZ: a copy of olorful_pp_FFn_variables, where also some other
    variables (E.g. the x's of eq 2.25)
    are added Variables for 'n' *pure final state* collinear
    recoiling exclusively against the initial state."""

    all_p_fs = [higher_PS_point[child] for child in children]
    Q = opts['Q']


    alpha =0.1 ## MZ

    # Loop over all final state momenta to obtain the corresponding variables
    zs = []
    for p_fs in all_p_fs:
        zs.append(p_fs.dot(Q))
    normalisation = sum(zs)
    for i in range(len(zs)):
        zs[i] /= normalisation

    # Add additional ss's
    ss_i_j = {}
    for i_fs in range(len(all_p_fs)):
        for j_fs in range(i_fs+1,len(all_p_fs)):
            ss_i_j[(i_fs,j_fs)] = 2.*all_p_fs[i_fs].dot(all_p_fs[j_fs])
    ss_i_j_tot = sum(ss_i_j.values())
    for k,v in ss_i_j.items():
        ss_i_j[(k[1],k[0])] = v

    ss_i_others = []
    for i_fs in range(len(all_p_fs)):
        p_other = vectors.LorentzVector()
        for j_fs in range(len(all_p_fs)):
            if j_fs==i_fs: continue
            p_other += all_p_fs[j_fs]
        ss_i_others.append(2.*all_p_fs[i_fs].dot(p_other))

    xis = [(z-ss_i_others[i]/(alpha*(2.*sum(all_p_fs).dot(Q)))) for i, z in enumerate(zs)]

    if len(all_p_fs)==3:
        # WARNING function below is *not* symmetric in r and s.
        def bigZ(i,r,s):
            return (ss_i_j[(i,r)] - ss_i_j[(r,s)] - 2.*zs[r]*ss_i_j_tot)/(alpha*2.*qC.dot(Q))
    elif len(all_p_fs)==2:
        def bigZ(i,r):
            return (ss_i_j[(i,r)]*(zs[r]-zs[i]))/(alpha*2.*qC.dot(Q))
    else:
        raise NotImplementedError

    kTs = {}
    for i in range(len(all_p_fs)):
        kT = vectors.LorentzVector()
        other_indices = [k for k in range(len(all_p_fs)) if k!=i]
        for j in other_indices:
            kT += xis[i]*all_p_fs[j]
            kT -= xis[j]*all_p_fs[i]

        kT += bigZ(i,*other_indices)*qC
        kTs[(i,tuple(other_indices))] = kT


    # now the fks variables
    # fist, assert that we are in the partonic c.o.m frame
    assert (higher_PS_point[1][0] == higher_PS_point[2][0] and
            higher_PS_point[1][3] == - higher_PS_point[2][3])

    #MZ: maybe it is enough to compute these variables only for the children, as for
    # the previous ones

    # then, the rescaled energies (fks_xi_i)
    sqrts = math.sqrt(2 * higher_PS_point[1].dot(higher_PS_point[2]))
    fks_xi_i = [2 * p[0] / sqrts  for p in all_p_fs]

    fks_y_ij = {}
    # we fill assuming i>j. Note that i and j here may NOT
    #  coincide with the i/j in the FKS convention

    for j in range(len(all_p_fs)-1):
        for i in range(len(all_p_fs))[j+1:]:
            fks_y_ij[(i, j)] = vectors.LorentzVector.cos(all_p_fs[i], all_p_fs[j])
            fks_y_ij[(j, i)] = fks_y_ij[(i, j)]

    fks_z_ij = {} # madfks paper eq 5.9
    # we fill assuming i>j. Note that i and j here may NOT
    #  coincide with the i/j in the
    #  FKS convention
    for j in range(len(all_p_fs)-1):
        for i in range(len(all_p_fs))[j+1:]:
            fks_z_ij[(i, j)] = all_p_fs[i][0] / \
                    (all_p_fs[i][0] + all_p_fs[j][0])
            fks_z_ij[(j, i)] = 1. - fks_z_ij[(i, j)]

    return [{'zs':tuple(zs), 'kTs':kTs, 'ss':ss_i_j, 'ss_i_others':ss_i_others,
        'fks_xi': tuple(fks_xi_i), 'fks_y': fks_y_ij, 'fks_z': fks_z_ij},]
