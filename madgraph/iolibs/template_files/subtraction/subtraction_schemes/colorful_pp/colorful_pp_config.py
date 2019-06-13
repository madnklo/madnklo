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

#=========================================================================================
# Variables, mappings, jacobians, factors and cuts
#=========================================================================================
# Note that variables, factors and cuts will be class members by design
# so they can easily be overridden by subclasses.
# They will be taken from the following variables
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

divide_by_jacobian = False

get_recoilers = get_initial_state_recoilers

# Initial collinear configuration
initial_coll_variables = currents.Q_initial_coll_variables
initial_coll_factor = factors_and_cuts.no_factor
initial_coll_cut = factors_and_cuts.cut_initial_coll
initial_coll_mapping = mappings.InitialLorentzOneMapping

# Soft configuration
soft_factor = factors_and_cuts.no_factor
soft_cut = factors_and_cuts.no_cut
soft_mapping = mappings.SoftVsInitialMapping

# Final collinear configuration
# WARNING: This is *not* the same final-collinear mapping as in colorful, where one has 'FinalRescalingOneMapping' instead.
# The __init__.py of colorful_pp will make sure to overwrite this mapping choice for the final collinear imported from
# coloful. For colorful_pp we need even final state collinears to recoil against initial states.
final_coll_mapping = mappings.FinalCollinearVsInitialMapping
final_coll_factor = factors_and_cuts.no_factor
final_coll_cut = factors_and_cuts.cut_coll
final_coll_variables = currents.Q_final_coll_mapping_recoiling_against_initial_state

# Final soft-collinear configuration (not strictly speaking necessary)
final_soft_coll_variables = currents.compute_energy_fractions
final_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=final_coll_mapping)
initial_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=initial_coll_mapping)