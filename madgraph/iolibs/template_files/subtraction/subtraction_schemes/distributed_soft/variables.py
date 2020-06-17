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
import distributed_soft_config as config

import commons.utils as utils

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# Variables as in Gabor's note, suited for a mapping of final-state collinears
# recoiling against the INITIAL states. See Gabor's notes.
#=========================================================================================
def distributed_soft_FFn_variables(higher_PS_point, qC, children, **opts):
    """ Variables for 'n' *pure final state* collinear recoiling exclusively against the initial state.""" #final state recoilers problem??
    # import ipdb
    # ipdb.set_trace()

    all_p_fs = [higher_PS_point[child] for child in children] # read out momenta of the children
    Q = opts['Q']
    # TODO SET ALPHA CORRECTLY :: What is alpha ?
    alpha = 1.0
    # alpha = 2.0
    # there is not opts['alpha']

    # print("alpha is "+str(alpha))
    
    # Loop over all final state momenta to obtain the corresponding variables
    zs = [] # momentum fractions of the children of the singular parent with respect to the light like reference vector n
    n = config.ref_n # now with massless reference vector :: before massive Q
    for p_fs in all_p_fs:
        zs.append(p_fs.dot(n))
    normalisation = sum(zs)
    for i in range(len(zs)):
        zs[i] /= normalisation

    # Add additional ss's
    # 2 * the product of the momenta of the children no doubles and masses
    # s_{ij}
    ss_i_j = {}
    for i_fs in range(len(all_p_fs)):
        for j_fs in range(i_fs+1,len(all_p_fs)):
            ss_i_j[(i_fs,j_fs)] = 2.*all_p_fs[i_fs].dot(all_p_fs[j_fs])
    ss_i_j_tot = sum(ss_i_j.values())
    for k,v in ss_i_j.items():
        ss_i_j[(k[1],k[0])] = v

    # product of the momenta of the children with the sum of the moment of all other children times 2
    # ss_i_others[i] = S_{i(1, ..., i-1, i+1, ..., n)}
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
        other_indices = [k for k in range(len(all_p_fs)) if k!=i] # all colinear momenta exept i
        for j in other_indices:
            kT += xis[i]*all_p_fs[j]
            kT -= xis[j]*all_p_fs[i]

        kT += bigZ(i,*other_indices)*qC
        kTs[(i,tuple(other_indices))] = kT

    return [{'zs':tuple(zs), 'kTs':kTs, 'ss':ss_i_j, 'ss_i_others':ss_i_others },]






























