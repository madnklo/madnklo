###########################################################
#
# colorful_pp subtraction scheme -- configuration variables
#
###########################################################

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings
import commons.utils as utils
import commons.QCD_local_currents as currents
import commons.factors_and_cuts as factors_and_cuts
import madgraph.various.misc as misc

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
def get_initial_state_recoilers(reduced_process, excluded=()):

    model = reduced_process.get('model')
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            leg['state'] == leg.INITIAL,
            leg['number'] not in excluded
        ])
    ])

def get_final_state_recoilers(reduced_process, excluded=()):
    model = reduced_process.get('model')
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded
        ])
    ])

divide_by_jacobian = False

# Initial collinear configuration
initial_coll_factor = factors_and_cuts.no_factor
initial_coll_cut = factors_and_cuts.cut_initial_coll
initial_coll_mapping = mappings.InitialLorentzOneMapping

# Soft configuration
soft_factor = factors_and_cuts.no_factor
soft_cut = factors_and_cuts.no_cut
soft_mapping = mappings.SoftVsInitialMapping

# Final collinear configuration
# WARNING: This is *not* the same final-collinear mapping as in colorful, where one has 'FinalRescalingOneMapping' instead.
final_coll_mapping = mappings.FinalCollinearVsInitialMapping
final_coll_factor = factors_and_cuts.no_factor
final_coll_cut = factors_and_cuts.cut_coll

# Final soft-collinear configuration (not strictly speaking necessary)
final_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=final_coll_mapping)
initial_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=initial_coll_mapping)

def generalised_cuts(bundles_cut_info,Q):
    """ Function applying the correct cut for each bundles depending on the variables passed for each which can be:
        {'pA': ..., 'pC': ....} for the initial-state momenta -pA (if any, otherwise absent) and the final-state ones pC for collines
        {'pS': ...} for the soft
    """
    for bundle_cut_info in bundles_cut_info:
        if 'pS' in bundle_cut_info:
            return soft_cut(pS=bundle_cut_info['pS'], Q=Q)
        elif 'pA' in bundle_cut_info:
            return initial_coll_cut(pA=bundle_cut_info['pA'], pR=bundle_cut_info['pC'], Q=Q)
        else:
            return final_coll_cut(pC=bundle_cut_info['pC'], Q=Q)