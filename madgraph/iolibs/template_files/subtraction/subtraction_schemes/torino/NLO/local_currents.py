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
"""Implementation of NLO fks local currents."""

import math
from bidict import bidict

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings

import commons.utils as utils
import commons.general_current as general_current
import dipole_current as dipole_current
import commons.factors_and_cuts as factors_and_cuts
import madgraph.integrator.vectors as vectors


import torino_config
import variables as kernel_variables

from commons.universal_kernels import AltarelliParisiKernels, SoftKernels

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# NLO final collinears
#=========================================================================================



class QCD_TRN_C_FgFq(general_current.GeneralCurrent):
    """ FF C(q_g)"""

    # Enable the flag below to debug this current
    DEBUG = True

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures, we want to match both a quark or an antiquark
    coll_structure_q = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
        )
    )
    coll_structure_qx = sub.CollStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(11, -1, sub.SubtractionLeg.FINAL),
        )
    )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(coll_structure_q,)),
        }),
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(coll_structure_qx,)),
        }),
    ]

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL),
                    sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : torino_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : general_current.CompoundVariables(kernel_variables.TRN_FFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        # inside variables we assume always i > j, but we filled (i,j), not (j,i);
        # in the mapping rules, the quark(j) comes second, the gluon first(i)
        z_ji = all_steps_info[0]['variables'][0]['fks_z'][(1,0)]

        prefactor = 1./s_rs # this comes putting together 5.16 and 5.18 of the madfks paper
        # The P_gq kernel is the one that matches the z definitions above.
        for spin_correlation_vector, weight in AltarelliParisiKernels.P_gq(self, z_ji, kT_FF):
            complete_weight = weight * prefactor
            if spin_correlation_vector is None:
                evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
            else:
                evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
                evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation

#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_TRN_S_g(dipole_current.DipoleCurrent):
    """ F S(q_g)"""

    # Enable the flag below to debug this current
    DEBUG = False

    divide_by_jacobian = torino_config.divide_by_jacobian

    # We should not need global variables for this current
    variables = None

    # Now define the matching singular structures
    soft_structure = sub.SoftStructure(
        substructures=tuple([]),
        legs=(
            sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
        )
    )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_structure,)),
        }),
    ]

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.SoftStructure(
                substructures=tuple([]),
                legs=(
                    sub.SubtractionLeg(11, 21, sub.SubtractionLeg.FINAL),
                )
            ),)),
            'mapping'               : torino_config.soft_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({}),
            'variables'             : None,
            'is_cut'                : torino_config.generalised_cuts,
            # Note that even though the soft mapping recoils against the final state, it still needs to
            # receive the list of *final state* as the recoilers, because those are the ones which will absorb
            # the transverse momentum recoil.
            'reduced_recoilers'     : torino_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a soft type of splitting kernel, which *does* need to know about the reduced process.
        Should be specialised by the daughter class if not dummy
        """

        # The colored_partons argument gives a dictionary with keys being leg numbers and values being
        # their color representation. At NLO, all we care about is what are the colored leg numbers.
        colored_partons = sorted(colored_partons.keys())

        # Normalization factors
        norm = -1./2.

        # All contributions from the soft counterterms are color correlated so we should remove
        # the default value None from the list of color correlations
        evaluation['color_correlations'] = []

        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        for i, a in enumerate(colored_partons):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in colored_partons[i:]:
                # Write the eikonal for that pair
                if a != b:
                    mult_factor = 2.
                else:
                    continue # MZ as long as we just care for massless partons, this is irrelevant.
                    # note that this has to be taken care of for the case of massive one, in
                    # particular how to identify the extra recoiler for the momentum mapping
                    mult_factor = 1.

                # Retrieve the reduced kinematics as well as the soft momentum
                lower_PS_point = all_steps_info[(a, b)][0]['lower_PS_point']
                higher_PS_point = all_steps_info[(a, b)][0]['higher_PS_point']

                ###pS = higher_PS_point[all_steps_info[(a, b)][0]['bundles_info'][0]['final_state_children'][0]]
                # the above expression cannot be used any more since we use a mapping structure of collinear
                #  type also for this current
                pS = higher_PS_point[ self.map_leg_number(
                                    global_variables['leg_numbers_map'],
                                    self.soft_structure.legs[0].n,
                                    global_variables['overall_parents'])]

                pa = higher_PS_point[a]
                pb = higher_PS_point[b]
                # some special behaviour may be needed in the case
                # either a or b are one of the particles defining
                # the sector, in particular for the SC limit.

                # At the moment, we do not implement anything special
                eikonal = SoftKernels.eikonal_dipole(pa, pb, pS)
                evaluation['color_correlations'].append(((a, b),))
                evaluation['values'][(0, color_correlation_index, 0, 0)] = {
                    'finite': norm * mult_factor * eikonal}
                evaluation['reduced_kinematics'].append(('Dip %d-%d' % (a, b), lower_PS_point))
                color_correlation_index += 1

        return evaluation

#=========================================================================================
# NLO final soft-collinear currents
#=========================================================================================


class QCD_FKS_CS_FgFq(general_current.GeneralCurrent):
    """ FF C(S(g),q)
    NLO tree-level (final) soft-collinear currents. Uses a collinear kernel,
    however wrt C_FqFg the order of q/g is exchanged"""

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(11,  +1, sub.SubtractionLeg.FINAL), ) )

    # This counterterm will be used if any of the current of the list below matches
    currents = [
        sub.Current({
            'resolve_mother_spin_and_color'     : True,
            'n_loops'                           : 0,
            'squared_orders'                    : {'QCD': 2},
            'singular_structure'                : sub.SingularStructure(substructures=(soft_coll_structure_q,)),
        }),
    ]

    # Te defining currents correspond to a currents block composed of several lists of template currents
    # for matching each currents of the target currents block (in an unordered fashion)
    defining_currents = [ currents, ]

    # An now the mapping rules
    mapping_rules = [
        {
            'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
                substructures=[sub.SoftStructure(
                               legs=(sub.SubtractionLeg(10, 21, sub.SubtractionLeg.FINAL), ) )],
                legs=(sub.SubtractionLeg(11, +1, sub.SubtractionLeg.FINAL),)
            ),)),
            'mapping'               : torino_config.final_coll_mapping,
            # Intermediate legs should be strictly superior to a 1000
            'momenta_dict'          : bidict({-1:frozenset((10,11))}),
            'variables'             : None, #general_current.CompoundVariables(kernel_variables.fks_FFn_variables),
            'is_cut'                : torino_config.generalised_cuts,
            'reduced_recoilers'     : torino_config.get_final_state_recoilers,
            'additional_recoilers'  : sub.SubtractionLegSet([]),
        },
    ]

    # For the qg soft-collinear the color prefactor is CF
    color_factor = "CF"

    def kernel(self, evaluation, all_steps_info, global_variables):
        """ Evaluate this counterterm given the variables provided. """

        kT_FF = all_steps_info[0]['variables'][0]['kTs'][(0,(1,))]
        z_FF  = all_steps_info[0]['variables'][0]['zs'][0]
        s_rs  = all_steps_info[0]['variables'][0]['ss'][(0,1)]

        # inside variables we assume always i > j, but we filled (i,j), not (j,i);
        # in the mapping rules, the gluon(i) comes first, then the quark(j)
        z_ji = all_steps_info[0]['variables'][0]['fks_z'][(1,0)]
        y_ij = all_steps_info[0]['variables'][0]['fks_y'][(1,0)]
        xi_i = all_steps_info[0]['variables'][0]['fks_xi'][1] #<-- why 1?? it should be 0
        ### CHECK!
        shat = 2 * all_steps_info[0]['higher_PS_point'][1].dot( \
            all_steps_info[0]['higher_PS_point'][2])

        prefactor = 2./shat/xi_i**2/(1-y_ij)
        # this comes from 5.19 of the madfks paper in the zij->1 limit
        # the extra (1-z) is from the definition of the barred kernels

        # The P_gq kernel is the one that matches the z definitions above.
        evaluation['values'][(0, 0, 0, 0)] = {'finite': prefactor * 2 * self.CF}
        ##for spin_correlation_vector, weight in AltarelliParisiKernels.P_gq(self, z_ji, kT_FF):
        ##    complete_weight = weight * prefactor
        ##    if spin_correlation_vector is None:
        ##        evaluation['values'][(0, 0, 0, 0)] = {'finite': complete_weight[0]}
        ##    else:
        ##        evaluation['spin_correlations'].append( ((parent, spin_correlation_vector),) )
        ##        evaluation['values'][(len(evaluation['spin_correlations']) - 1, 0, 0, 0)] = {'finite': complete_weight[0]}

        return evaluation


#=========================================================================================
# NLO initial-collinear currents
#=========================================================================================


#=========================================================================================
# NLO soft initial-collinear currents
#=========================================================================================

