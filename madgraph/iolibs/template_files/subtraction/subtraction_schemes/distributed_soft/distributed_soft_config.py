###########################################################
#
# distributed soft subtraction scheme -- configuration variables
#
###########################################################

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings
import commons.utils as utils
import madgraph.various.misc as misc
import madgraph.integrator.vectors as vec

import math

import numpy as np

# from scipy.optimize import newton

CurrentImplementationError = utils.CurrentImplementationError

# In principle we support them all, but these are the meaningful ones to consider
beam_PDGs_supported= [
    tuple(sorted([1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 21])),
    tuple(sorted([1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 21])),
    tuple(sorted([1, -1, 2, -2, 3, -3, 4, -4, 21])),
    tuple(sorted([1, -1, 2, -2, 3, -3, 21])),
    tuple(sorted([1, -1, 2, -2, 21])),
    tuple(sorted([1, -1, 21]))
]

#=========================================================================================
# mappings, jacobians, factors and cuts
#=========================================================================================
# Note that factors and cuts will be class members by design
# so they can easily be overridden by subclasses.
# They will be taken from the entities below
# so we can quickly switch them coherently across the entire subtraction scheme.

import madgraph.integrator.mappings as mappings

# We consider only the two initial state of the reduced process as recoilers,
# assuming that the initial state numbers of the lower multiplicity process
# are identical to those of the final state numbers.

# No initial state recoilers
# def get_initial_state_recoilers(reduced_process, excluded=(), **opts):
#
#     model = reduced_process.get('model')
#     return sub.SubtractionLegSet([
#         leg for leg in reduced_process.get('legs') if all([
#             leg['state'] == leg.INITIAL,
#             leg['number'] not in excluded
#         ])
#     ])

def get_final_state_recoilers(reduced_process, excluded=(), **opts):

    model = reduced_process.get('model')
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded
        ])
    ])

def no_factor(**opts): #same as no factor and no cut in factors and cuts for colourful
    return 1.0

def no_cut(**opts):
    return False

divide_by_jacobian_local = False # What does this do?
divide_by_jacobian_integrated = False

# Initial collinear configuration
# initial_coll_factor = no_factor
# initial_coll_cut = no_cut # = False # Different from colourful
# initial_coll_mapping = mappings.InitialLorentzOneMapping

# Soft configuration
# soft_factor = no_factor
# # The soft counterterm must not have cuts as we did not implement the soft integrated counterterm with cuts
# soft_cut = no_cut # = False
# soft_mapping = mappings.SoftVsInitialMapping

# Final collinear configuration
# WARNING: This is *not* the same final-collinear mapping as in colorful, where one has 'FinalRescalingOneMapping' instead.
final_coll_mapping = mappings.FinalGroupingMapping #Different from colourful
final_coll_factor = no_factor # = 1.0
# A cut on the final-final would make little sense since we integrate numerically over the symmetric rescaling of the initial
# momenta, so this is the same story as for the soft.
final_coll_cut = no_cut # = False

# # Final soft-collinear configuration (not strictly speaking necessary)
# final_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
#     soft_mapping=soft_mapping, collinear_mapping=final_coll_mapping)
# initial_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
#     soft_mapping=soft_mapping, collinear_mapping=initial_coll_mapping)

# light like reference vector
ref_n = vec.LorentzVector(([1., 0., 0., 1.]))

ref_m = vec.LorentzVector(([1., 1., 0., 0.]))

def generalised_cuts(cut_inputs, global_variables):
    """ Function applying the correct cut for each bundles depending on the variables passed for each which can be:
        {'pA': ..., 'pC': ....} for the initial-state momenta -pA (if any, otherwise absent) and the final-state ones pC for collines
        {'pS': ...} for the soft
    """
    Q = global_variables['Q']
    # print("generalised_cuts is called with cut inputs")
    # print(cut_inputs)

    for bundle_info in cut_inputs['bundles_info']:
        bundle_cut_info = bundle_info['cut_inputs']
        if 'pS' in bundle_cut_info:
            if soft_cut(pS=bundle_cut_info['pS'], Q=Q):
                return True
        elif 'pA' in bundle_cut_info:
            if initial_coll_cut(pA=bundle_cut_info['pA'], pR=bundle_cut_info['pC'], Q=Q):
                return True
        elif 'pC' in bundle_cut_info:
            if final_coll_cut(pC=bundle_cut_info['pC'], Q=Q):
                return True
        else:
            raise CurrentImplementationError("The generalised_cuts function is only applicable for mapping structure "+
                                                                                              " are collinear or soft.")
    return False

def Jacobian_determinant (all_steps_info, global_variables, all_mapped_masses_are_zero=True) :

    # returns the Jacobian of the generalized rescaling mapping for the give phase space point and mapped phase
    # space point.

    # print (all_steps_info[0]['bundles_info'][0]['final_state_children'])

    Q = global_variables['Q']
    Q2 = Q.square()
    # print(Q)

    n_initial_legs = global_variables['n_initial_legs']

    lower_PS_point = all_steps_info[0]['lower_PS_point']
    higher_PS_point = all_steps_info[0]['higher_PS_point']
    # print (lower_PS_point)
    # print (higher_PS_point)

    parent = all_steps_info[0]['bundles_info'][0]['parent']

    m = len(lower_PS_point) - n_initial_legs # number of resolved final state particles

    d = 4 # dimension of space time

    kappa = 0.0
    a = 1.0
    b = 0.0
    c = 0.0

    final_state = dict()
    parent_added = False # only for one parent

    for i, particle in enumerate (higher_PS_point):
        if particle <= n_initial_legs:
            continue

        if particle in all_steps_info[0]['bundles_info'][0]['final_state_children']:
            if parent_added:
                continue
            final_state.update({parent : sum(higher_PS_point[u] for u in all_steps_info[0]['bundles_info'][0]['final_state_children'])})
            parent_added = True
            continue

        final_state.update({particle : higher_PS_point[particle]})

    mapped_final_state = dict()
    for i, particle in enumerate (lower_PS_point):
        # print(i)
        # print(particle)
        if particle <= n_initial_legs:
            continue

        mapped_final_state.update({particle : lower_PS_point[particle]})

    # print (final_state)
    # print (mapped_final_state)

    for i in final_state:
        # print (a)
        a *= (mapped_final_state[i].dot(Q))/(final_state[i].dot(Q))
        # print(b)
        b += (((mapped_final_state[i].dot(Q))**2)/Q2 - mapped_final_state[i].square())/(mapped_final_state[i].dot(Q))
        # print(c)
        c += (((mapped_final_state[i].dot(Q))**2)/Q2 - mapped_final_state[i].square())/(final_state[i].dot(Q))

    factor = a * (b / c)

    if all_mapped_masses_are_zero:
        for i in final_state:
            # print(kappa)
            kappa += math.sqrt(((final_state[i].dot(Q))**2)/Q2 - final_state[i].square())
        kappa /= math.sqrt(Q2)
    else: # solve numerically :: check sign of \vec{pi}^2

        # if m == 2:
        #     kappa = 1.0
        #
        # else:
        kappa_0 = 1.0 # starting point for numerical solution
        from scipy.optimize import newton

        def func(k):
            return sum (math.sqrt((final_state[i].square() - ((final_state[i].dot(Q))**2)/Q2)/k**2 + mapped_final_state[i].square()) for i in final_state) - math.sqrt(Q2)

        def df(k):
            return sum ((((final_state[i].square() - ((final_state[i].dot(Q))**2)/Q2) * (-2.0))/k**3)/math.sqrt((final_state[i].square() - ((final_state[i].dot(Q))**2)/Q2)/k**2 + mapped_final_state[i].square()) for i in final_state)
        kappa = newton(func,x_0=kappa_0,fprime=df)
        # print("not all masses zero mapping Jacobian not implemented")

    determinant = kappa ** (d*(m-1)-(m+1)) * factor
    # print(determinant)

    return determinant

















