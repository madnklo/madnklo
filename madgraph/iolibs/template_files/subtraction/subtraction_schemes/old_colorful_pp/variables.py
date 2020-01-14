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

#=========================================================================================
# Variables as in Gabor's note, suited for a mapping of final-state collinears
# recoiling against the INITIAL states. See Gabor's notes.
#=========================================================================================
def colorful_pp_FFn_variables(higher_PS_point, qC, children, **opts):
    """ Variables for 'n' *pure final state* collinear recoiling exclusively against the initial state."""

    all_p_fs = [higher_PS_point[child] for child in children]
    Q = opts['Q']
    alpha = opts['alpha']

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

    return [{'zs':tuple(zs), 'kTs':kTs, 'ss':ss_i_j, 'ss_i_others':ss_i_others },]

def colorful_pp_IFn_variables(PS_point, parent_momentum, children, **opts):
    """ Variables for 'n' initial-state collinear recoiling exclusively against the initial state."""

    p_a = PS_point[children[0]]
    p_fs = [PS_point[child] for child in children[1:]]
    Q = opts['Q']
    pa_dot_pb = Q.square()/2.

    # Loop over all final state momenta to obtain the corresponding variables
    xs = []
    kTs = []
    ss = {}
    for i_fs in range(len(p_fs)):
        p_r = p_fs[i_fs]
        if len(p_fs)>1:
            p_s = sum(p_fs[i_fs_prime] for i_fs_prime in range(len(p_fs)) if i_fs_prime!=i_fs)
        else:
            p_s = vectors.LorentzVector()
        xs.append(
            (p_r.dot(Q)/pa_dot_pb) - (p_r.dot(p_s)/pa_dot_pb)*(p_a.dot(p_r)/p_a.dot(p_r+p_s))
        )
        kTs.append(
            p_r - ((p_r.dot(Q)-2.*p_a.dot(p_r))/pa_dot_pb)*p_a - (p_a.dot(p_r)/pa_dot_pb)*Q
        )
        ss[(0,i_fs+1)] = 2.*p_a.dot(p_r)

    # Finally add the x and kT corresponding to the initial state leg
    xs.insert(0, 1.0-sum(xs))
    #xs.insert(0, 1.0 - (((p_fs[0]+p_fs[1]).dot(Q)-p_fs[0].dot(p_fs[1]))/pa_dot_pb) )
    kTs.insert(0, -sum(kTs))

    # Alternative way of getting these variables using Simone's mappings function
#    na = parent_momentum
#    nb = opts['Q']
#    kin_variables = dict()
    # The lone initial state child is always placed first thanks to the implementation
    # of the function get_sorted_children() in the current.
#    mappings.InitialCollinearVariables.get(
#        PS_point, children[1:], children[0], na, nb, kin_variables)
#    xs = tuple(kin_variables['z%d' % i] for i in children)
#    kTs = tuple(kin_variables['kt%d' % i] for i in children)

    # Add additional ss's
    for i_fs in range(len(p_fs)):
        for j_fs in range(i_fs+1,len(p_fs)):
            ss[(i_fs+1,j_fs+1)] = 2.*p_fs[i_fs].dot(p_fs[j_fs])

    return [{'xs':tuple(xs), 'kTs':tuple(kTs), 'ss':ss },]

def colorful_pp_FFF_softFF_variables(higher_PS_point, qC, children, **opts):
    """ Variables for the *pure final state* collinear recoiling exclusively against the initial state,
    with a nest soft-structure
    """
    p_i_tilde = qC
    # Only include here the momenta from the final_states *not* involved in the soft structures
    # which *by convention* is the first child
    p_fs = [higher_PS_point[child] for child in children[1:]]
    Q = opts['Q']

    # Loop over all final state momenta to obtain the corresponding variables
    zs = []
    kTs = []
    ss = {}
    for i_fs in range(len(p_fs)):
        p_r = p_fs[i_fs]
        zs.append( p_r.dot(Q) / p_i_tilde.dot(Q) )
        # KTs not needed for now, will be filled in if needed.
        kTs.append(None)
        ss[(0, i_fs + 1)] = 2. * p_i_tilde.dot(p_r)

    # Finally add the x and kT corresponding to the initial state leg
    # xs.insert(0, 1.0-sum(xs))
    zs.insert(0, 1.0-sum(zs))
    #kTs.insert(0, -sum(kTs))
    kTs.insert(0, None)

    # Add additional ss's
    for i_fs in range(len(p_fs)):
        for j_fs in range(i_fs + 1, len(p_fs)):
            ss[(i_fs + 1, j_fs + 1)] = 2. * p_fs[i_fs].dot(p_fs[j_fs])

    return [{'xs':tuple(zs), 'kTs':tuple(kTs), 'ss':ss},]

def colorful_pp_IFF_softFF_variables(PS_point, parent_momentum, children, **opts):
    """ Variables for the initial-state collinear recoiling exclusively against the initial state,
    with a nested soft-structure
    """

    p_a_tilde = parent_momentum
    p_fs = [PS_point[child] for child in children[1:]]
    Q = opts['Q']

    # Loop over all final state momenta to obtain the corresponding variables
    xs = []
    kTs = []
    ss = {}
    for i_fs in range(len(p_fs)):
        p_r = p_fs[i_fs]
        xs.append( p_r.dot(Q) / p_a_tilde.dot(Q) )
        # KTs not needed for now, will be filled in if needed.
        kTs.append(None)
        ss[(0, i_fs + 1)] = 2. * p_a_tilde.dot(p_r)

    # Finally add the x and kT corresponding to the initial state leg
    # xs.insert(0, 1.0-sum(xs))
    xs.insert(0, 1.0-sum(xs))
    #kTs.insert(0, -sum(kTs))
    kTs.insert(0, None)

    # Alternative way of getting these variables using Simone's mappings function
    #    na = parent_momentum
    #    nb = opts['Q']
    #    kin_variables = dict()
    # The lone initial state child is always placed first thanks to the implementation
    # of the function get_sorted_children() in the current.
    #    mappings.InitialCollinearVariables.get(
    #        PS_point, children[1:], children[0], na, nb, kin_variables)
    #    xs = tuple(kin_variables['z%d' % i] for i in children)
    #    kTs = tuple(kin_variables['kt%d' % i] for i in children)

    # Add additional ss's
    for i_fs in range(len(p_fs)):
        for j_fs in range(i_fs + 1, len(p_fs)):
            ss[(i_fs + 1, j_fs + 1)] = 2. * p_fs[i_fs].dot(p_fs[j_fs])

    return [{'xs':tuple(xs), 'kTs':tuple(kTs), 'ss':ss},]