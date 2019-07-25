###########################################################
#
# distributed_soft_qqQQ subtraction scheme -- configuration variables
#
###########################################################

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings
import commons.utils as utils
import commons.QCD_local_currents as currents
import commons.factors_and_cuts as factors_and_cuts

import madgraph.integrator.mappings as mappings

#=========================================================================================
# Choice of recoilers
#=========================================================================================


def get_recoilers(reduced_process, excluded=()):
    model = reduced_process.get('model')
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            model.get_particle(leg['id']).get('mass').upper() == 'ZERO',
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded
        ])
    ])


#===========================
# Triple collinear variables
#===========================
low_level_coll_variables= currents.n_final_coll_variables
def triple_collinear_global_variables(all_steps, global_info):
    PS_point=all_steps[0]['higher_PS_point']
    children=global_info['overall_children'][0]
    parent=global_info['overall_parents'][0]
    pC = all_steps[1]['lower_PS_point'][parent]
    return low_level_coll_variables(PS_point,pC,children)[0]

divide_by_jacobian = False

# Final collinear configuration
# NB: There are only collinear CTs in this scheme
final_coll_cut = currents.no_cut
final_coll_factor = currents.no_factor
final_coll_variables = low_level_coll_variables
global_final_coll_variables = triple_collinear_global_variables
final_coll_mapping = mappings.FinalGroupingMapping
