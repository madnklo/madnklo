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

# We consider all final state particles as recoiler, including massive and colorless ones
def get_all_final_recoilers(reduced_process, excluded=()):

    model = reduced_process.get('model')
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded
        ])
    ])

divide_by_jacobian = True

get_recoilers = get_all_final_recoilers

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
# coloful. Notice then that the final_coll quantity specified here apply only then to the *NNLO* final collinear.
# For the NLO ones, as explained above, these properties will be overwritten appropriately irrespectively of what is below.
final_coll_mapping = mappings.FinalLorentzOneMapping
final_coll_factor = factors_and_cuts.no_factor
final_coll_cut = factors_and_cuts.cut_coll

# Final soft-collinear configuration
final_soft_coll_variables = currents.compute_energy_fractions
final_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=final_coll_mapping)
initial_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=initial_coll_mapping)