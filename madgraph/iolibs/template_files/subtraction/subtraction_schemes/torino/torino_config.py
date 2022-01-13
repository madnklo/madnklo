###########################################################
#
# FKS subtraction scheme -- configuration variables
#
###########################################################

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings
import commons.utils as utils
import factors_and_cuts as factors_and_cuts
import madgraph.various.misc as misc
import logging

logger = logging.getLogger('madgraph')

CurrentImplementationError = utils.CurrentImplementationError

#gl
# In principle we support them all, but these are the meaningful ones to consider
beam_PDGs_supported= [
    tuple(sorted([1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 21])),
    tuple(sorted([1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 21])),
    tuple(sorted([1, -1, 2, -2, 3, -3, 4, -4, 21])),
    tuple(sorted([1, -1, 2, -2, 3, -3, 21])),
    tuple(sorted([1, -1, 2, -2, 21])),
    tuple(sorted([1, -1, 21]))
]

lepton_abs_PDGs = [11, -11, 12, -12, 13, -13]

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


#DEFAULT
#The following two function are the DEFAULT ones:
# - 'get_initial_state_recoilers' just choose initial-state rec;
# - 'get_final_state_recoilers' just choose final-state rec.

def pause_get_initial_state_recoilers(reduced_process, excluded=(), global_variables={}):

    model = reduced_process.get('model')
    try:
        sector_legs = global_variables['sector_info'][0].leg_numbers
    except AttributeError:
        # this is in the case of the integrated CT's
        sector_legs = ()

    recoilers = [
        leg for leg in reduced_process.get('legs') if all([
            leg['state'] == leg.INITIAL,
            leg['number'] not in excluded,
            leg['number'] not in sector_legs])
                ]

    return sub.SubtractionLegSet([recoilers[0]])


def pause_get_final_state_recoilers(reduced_process, excluded=(), global_variables={}):
    """in the Torino Subtraction scheme, a single recoiler has to be choosen
    At the moment we simply pick it among massless, final-state legs which
    do not belong to the sector leg. By convention, we pick the one with smallest
    PDG id """
    model = reduced_process.get('model')
    try:
        sector_legs = global_variables['sector_info'][0].leg_numbers
    except AttributeError:
        # this is in the case of the integrated CT's
        sector_legs = ()

    recoilers = [
        leg for leg in reduced_process.get('legs') if all([
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded,
            leg['number'] not in sector_legs,
             model.get_particle(leg['id'])['mass'].lower() == 'zero'])
                ]

    #logger.info('Recoilers set= ' + str(recoilers))

    # check that recoilers exist
    if not recoilers :
        raise CurrentImplementationError("Recoilers cannot be found for process %s" % reduced_process.nice_string())

    # sort the recoilers according to their id, and return the first one
    recoilers.sort(key = lambda l: l['id'])
    return sub.SubtractionLegSet([recoilers[0]])


#TRY 1
#The following two function are the same and
#ALWAYS PREFER initial-state rec over final-state rec when possible.
# - specific for partonic recoilers

def pause_get_initial_state_recoilers(reduced_process, excluded=(), global_variables={}):
# We prefere to choose initial state recoiler over final state recoilers when possible

    model = reduced_process.get('model')
    try:
        sector_legs = global_variables['sector_info'][0].leg_numbers
    except AttributeError:
        # this is in the case of the integrated CT's
        sector_legs = ()

    #logger.info('Sector legs= ' + str(sector_legs))

    partons_id = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 21]

    recoilers = [
        leg for leg in reduced_process.get('legs') if all([
            leg['id'] in partons_id,
            leg['number'] not in excluded,
            leg['number'] not in sector_legs])
                ]

    logger.info('Sector Legs= ' + str(sector_legs))
    if recoilers[0]['number'] > 2:
        # sort the recoilers according to their id, and return the first one
        recoilers.sort(key = lambda l: l['id'])

    return sub.SubtractionLegSet([recoilers[0]])


def pause_get_final_state_recoilers(reduced_process, excluded=(), global_variables={}):
# We prefere to choose initial state recoiler over final state recoilers when possible

    model = reduced_process.get('model')
    try:
        sector_legs = global_variables['sector_info'][0].leg_numbers
    except AttributeError:
        # this is in the case of the integrated CT's
        sector_legs = ()

    #logger.info('Sector legs= ' + str(sector_legs))

    partons_id = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 21]

    recoilers = [
        leg for leg in reduced_process.get('legs') if all([
            leg['id'] in partons_id,
            leg['number'] not in excluded,
            leg['number'] not in sector_legs])
                ]

    logger.info('Sector Legs= ' + str(sector_legs))
    if recoilers[0]['number'] > 2:
        # sort the recoilers according to their id, and return the first one
        recoilers.sort(key = lambda l: l['id'])

    return sub.SubtractionLegSet([recoilers[0]])


#CHOOSEN METHOD
#The following two function are the same:
# - for final-state collinear singularity they look for final-state rec over initial-state one when possible;
# - for initial-state collinear singularity they look for initial-state rec over final-state one when possible.
# - specific for partonic recoilers.

def get_initial_state_recoilers(reduced_process, excluded=(), global_variables={}):

    model = reduced_process.get('model')
    try:
        sector_legs = global_variables['sector_info'][0].leg_numbers
    except AttributeError:
        # this is in the case of the integrated CT's
        sector_legs = ()

    #logger.info('Sector Legs= ' + str(sector_legs))

    partons_id = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 21]

    initial_recoilers = [
        leg for leg in reduced_process.get('legs') if all([
            leg['id'] in partons_id,
            leg['state'] == leg.INITIAL,
            leg['number'] not in excluded,
            leg['number'] not in sector_legs])
                ]

    final_recoilers = [
        leg for leg in reduced_process.get('legs') if all([
            leg['id'] in partons_id,
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded,
            leg['number'] not in sector_legs])
                ]

    # for integrated currents
    if len(sector_legs) == 0:
        if global_variables['overall_children'][0][0] > 2 and global_variables['overall_children'][0][1] > 2:
            if len(final_recoilers) > 0:
                # sort the recoilers according to their id, and return the first one
                final_recoilers.sort(key = lambda l: l['id'])
                return sub.SubtractionLegSet([final_recoilers[0]])
            else: 
                return sub.SubtractionLegSet([initial_recoilers[0]])
        if global_variables['overall_children'][0][0] <= 2 or global_variables['overall_children'][0][1] <= 2:
            if len(initial_recoilers) > 0:
                return sub.SubtractionLegSet([initial_recoilers[0]])
            else: 
                # sort the recoilers according to their id, and return the first one
                final_recoilers.sort(key = lambda l: l['id'])
                return sub.SubtractionLegSet([final_recoilers[0]])


    # for local currents
    if sector_legs[0] > 2 and sector_legs[1] > 2:
        if len(final_recoilers) > 0:
            # sort the recoilers according to their id, and return the first one
            final_recoilers.sort(key = lambda l: l['id'])
            return sub.SubtractionLegSet([final_recoilers[0]])
        else:
            return sub.SubtractionLegSet([initial_recoilers[0]])
    elif sector_legs[0] <= 2 or sector_legs[1] <= 2:
        if len(initial_recoilers) > 0:
            return sub.SubtractionLegSet([initial_recoilers[0]])
        else:
            # sort the recoilers according to their id, and return the first one
            final_recoilers.sort(key = lambda l: l['id'])
            return sub.SubtractionLegSet([final_recoilers[0]])



def get_final_state_recoilers(reduced_process, excluded=(), global_variables={}):

    model = reduced_process.get('model')
    try:
        sector_legs = global_variables['sector_info'][0].leg_numbers
    except AttributeError:
        # this is in the case of the integrated CT's
        sector_legs = ()

    #logger.info('Sector Legs= ' + str(sector_legs))

    partons_id = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 21]

    initial_recoilers = [
        leg for leg in reduced_process.get('legs') if all([
            leg['id'] in partons_id,
            leg['state'] == leg.INITIAL,
            leg['number'] not in excluded,
            leg['number'] not in sector_legs])
                ]

    final_recoilers = [
        leg for leg in reduced_process.get('legs') if all([
            leg['id'] in partons_id,
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded,
            leg['number'] not in sector_legs])
                ]

    # for integrated currents
    if len(sector_legs) == 0:
        if global_variables['overall_children'][0][0] > 2 and global_variables['overall_children'][0][1] > 2:
            if len(final_recoilers) > 0:
                # sort the recoilers according to their id, and return the first one
                final_recoilers.sort(key = lambda l: l['id'])
                return sub.SubtractionLegSet([final_recoilers[0]])
            else: 
                return sub.SubtractionLegSet([initial_recoilers[0]])
        if global_variables['overall_children'][0][0] <= 2 or global_variables['overall_children'][0][1] <= 2:
            if len(initial_recoilers) > 0:
                return sub.SubtractionLegSet([initial_recoilers[0]])
            else: 
                # sort the recoilers according to their id, and return the first one
                final_recoilers.sort(key = lambda l: l['id'])
                return sub.SubtractionLegSet([final_recoilers[0]])


    # for local currents
    if sector_legs[0] > 2 and sector_legs[1] > 2:
        if len(final_recoilers) > 0:
            # sort the recoilers according to their id, and return the first one
            final_recoilers.sort(key = lambda l: l['id'])
            return sub.SubtractionLegSet([final_recoilers[0]])
        else:
            return sub.SubtractionLegSet([initial_recoilers[0]])
    elif sector_legs[0] <= 2 or sector_legs[1] <= 2:
        if len(initial_recoilers) > 0:
            return sub.SubtractionLegSet([initial_recoilers[0]])
        else:
            # sort the recoilers according to their id, and return the first one
            final_recoilers.sort(key = lambda l: l['id'])
            return sub.SubtractionLegSet([final_recoilers[0]])


def get_soft_recoilers(reduced_process, excluded=(), **opts):

    model = reduced_process.get('model')
    # logger.info('Soft_recoilers : ' + str([
    #     leg for leg in reduced_process.get('legs') if all([
    #         leg['number'] not in excluded
    #     ])
    # ]))
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            leg['number'] not in excluded
        ])
    ]) 



#TODO cuts are disabled for now because they need to be rethought for double unresolved limits
divide_by_jacobian = False

# Initial-Final collinear configuration
initial_coll_mapping = mappings.CollinearTRNMapping
initial_coll_factor = factors_and_cuts.no_factor
#initial_coll_cut = factors_and_cuts.cut_initial_coll
initial_coll_cut = factors_and_cuts.no_cut


# Soft configuration
soft_factor = factors_and_cuts.no_factor
soft_cut = factors_and_cuts.no_cut
soft_mapping = mappings.SoftTRNMapping


# Final collinear configuration
# WARNING: This is *not* the same final-collinear mapping as in FKS, where one has 'FinalRescalingOneMapping' instead.
#final_coll_mapping = mappings.FinalTRNMapping
final_coll_mapping = mappings.CollinearTRNMapping
final_coll_factor = factors_and_cuts.no_factor
#final_coll_cut = factors_and_cuts.cut_coll
final_coll_cut = factors_and_cuts.no_cut

# Final soft-collinear configuration (not strictly speaking necessary)
#final_soft_coll_mapping = mappings.FinalTRNMapping()
#initial_soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
#    soft_mapping=soft_mapping, collinear_mapping=initial_coll_mapping)

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
