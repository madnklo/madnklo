###########################################################
#
# colorful_pp subtraction scheme -- configuration variables
#
###########################################################

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings
import commons.utils as utils
import commons.QCD_local_currents as currents

#=========================================================================================
# Variables, mappings, jacobians, factors and cuts
#=========================================================================================
# Note that variables, factors and cuts will be class members by design
# so they can easily be overridden by subclasses.
# They will be taken from the following variables
# so we can quickly switch them coherently across the entire subtraction scheme.

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

def no_cut(**opts):
    return True

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

def no_factor(**opts):
    return 1.0

divide_by_jacobian = True

# Initial-collinear configuration
initial_coll_variables = currents.Q_initial_coll_variables
factor_initial_coll = factor_coll
cut_initial_coll = cut_coll

# Soft configuration
factor_soft = no_factor
cut_soft = no_cut

initial_coll_mapping = mappings.InitialLorentzOneMapping
soft_mapping = mappings.SoftVsInitialMapping
soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=initial_coll_mapping)