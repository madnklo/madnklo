###########################################################
#
# distributed soft subtraction scheme -- configuration variables
#
###########################################################

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings
import commons.utils as utils
import madgraph.various.misc as misc

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
def get_initial_state_recoilers(reduced_process, excluded=(), **opts):

    model = reduced_process.get('model')
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            leg['state'] == leg.INITIAL,
            leg['number'] not in excluded
        ])
    ])

def get_final_state_recoilers(reduced_process, excluded=(), **opts):

    model = reduced_process.get('model')
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded
        ])
    ])

def no_factor(**opts):
    return 1.0

def no_cut(**opts):
    return False

divide_by_jacobian = False

# Initial collinear configuration
initial_coll_factor = no_factor
initial_coll_cut = no_cut
initial_coll_mapping = mappings.InitialLorentzOneMapping

# Soft configuration
soft_factor = no_factor
# The soft counterterm must not have cuts as we did not implement the soft integrated counterterm with cuts
soft_cut = no_cut
soft_mapping = mappings.SoftVsInitialMapping

# Final collinear configuration
# WARNING: This is *not* the same final-collinear mapping as in colorful, where one has 'FinalRescalingOneMapping' instead.
final_coll_mapping = mappings.FinalCollinearVsInitialMapping
final_coll_factor = no_factor
# A cut on the final-final would make little sense since we integrate numerically over the symmetric rescaling of the initial
# momenta, so thhis is the same story as for the soft.
final_coll_cut = no_cut

# Final soft-collinear configuration (not strictly speaking necessary)
final_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=final_coll_mapping)
initial_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=initial_coll_mapping)

def generalised_cuts(cut_inputs, global_variables):
    """ Function applying the correct cut for each bundles depending on the variables passed for each which can be:
        {'pA': ..., 'pC': ....} for the initial-state momenta -pA (if any, otherwise absent) and the final-state ones pC for collines
        {'pS': ...} for the soft
    """
    Q = global_variables['Q']

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