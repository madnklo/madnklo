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
"""Common factors and cuts for subtraction schemes."""
# TODO Ideally these functions should be part of colorful and colorful_pp

import madgraph.integrator.mappings as mappings

alpha_0     = 0.5
y_0         = 0.5
y_0_prime   = 0.5
d_0         = 1
d_0_prime   = 2

def cut_coll(**opts):

    try:
        alpha = opts['alpha']
    except KeyError:
        pC    = opts['pC']
        Q     = opts['Q']
        alpha = mappings.FinalRescalingOneMapping.alpha(pC, Q)
    # Include the counterterm only up to alpha_0
    return alpha > alpha_0

def cut_initial_coll(**opts):

    pA    = opts['pA']
    pR    = opts['pR']
    Q     = opts['Q']
    y_0p  = (2.*pA.dot(pR))/Q.square()
    # Include the counterterm only up to y_0_prime
    return y_0p > y_0_prime

def cut_soft(**opts):

    try:
        y  = opts['y']
    except KeyError:
        pS = opts['pS']
        Q  = opts['Q']
        y = mappings.SoftVsFinalPureRescalingMapping.y(pS, Q)
    # Include the counterterm only up to y_0
    return y > y_0

def factor_coll(**opts):

    try:
        alpha = opts['alpha']
    except KeyError:
        pC    = opts['pC']
        Q     = opts['Q']
        alpha = mappings.FinalRescalingOneMapping.alpha(pC, Q)
    norm = (1 - alpha) ** (2 * (d_0 - 1))
    return norm

def factor_soft(**opts):

    try:
        y  = opts['y']
    except KeyError:
        pS = opts['pS']
        Q  = opts['Q']
        y = mappings.SoftVsFinalPureRescalingMapping.y(pS, Q)
    norm = (1 - y) ** (d_0_prime - 2)
    return norm
