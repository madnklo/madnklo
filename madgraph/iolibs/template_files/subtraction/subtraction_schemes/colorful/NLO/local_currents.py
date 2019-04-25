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
"""Implementation of NLO colorful currents."""

import math

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings

import commons.utils as utils
import commons.QCD_local_currents as currents
import commons.factors_and_cuts as factors_and_cuts

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

def get_recoilers(reduced_process, excluded=()):

    model = reduced_process.get('model')
    return sub.SubtractionLegSet([
        leg for leg in reduced_process.get('legs') if all([
            model.get_particle(leg['id']).get('mass').upper() == 'ZERO',
            leg['state'] == leg.FINAL,
            leg['number'] not in excluded
        ])
    ])

variables = staticmethod(currents.Q_final_coll_variables)

#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_soft_0_g(currents.QCDCurrent):
    """Soft gluon eikonal current at tree level, eq.4.12-4.13 of arXiv:0903.1218."""
    # TODO Ideally set recoilers for mapping_singular_structure
    # within does_implement_this_current and __init__

    structure = sub.SingularStructure(substructures=(
        sub.SoftStructure(legs=[sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), ]),
    ))
    mapping = mappings.SoftVsFinalPureRescalingMapping
    factor = staticmethod(factors_and_cuts.factor_soft)
    is_cut = staticmethod(factors_and_cuts.cut_soft)

    @classmethod
    def does_implement_this_current(cls, current, model):

        if not all([
            current.get('squared_orders') == {'QCD': 2},
            current.get('n_loops') == 0,
            current.get('resolve_mother_spin_and_color') == True,
            not isinstance(current, (sub.BeamCurrent, sub.IntegratedBeamCurrent, sub.IntegratedCurrent))
        ]):
            return None

        leg_numbers_map = cls.structure.map_leg_numbers(
            current.get('singular_structure'), [range(1, model.get_nflav()+1)])
        if leg_numbers_map is None:
            return None
        mapping_singular_structure = current.get('singular_structure').get_copy()
        return {
            'leg_numbers_map': leg_numbers_map,
            'mapping_singular_structure': mapping_singular_structure }

    def __init__(self, *args, **opts):

        for opt_name in ['leg_numbers_map', 'mapping_singular_structure']:
            try:
                setattr(self, opt_name, opts.pop(opt_name))
            except KeyError:
                raise CurrentImplementationError(
                    "__init__ of " + self.__class__.__name__ + " requires " + opt_name)
        self.supports_helicity_assignment = False
        super(QCD_soft_0_g, self).__init__(*args, **opts)

    @staticmethod
    def eikonal(pi, pj, ps):
        """Eikonal factor for soft particle with momentum ps
        emitted from the dipole with momenta pi and pj.
        """

        pipj = pi.dot(pj)
        pips = pi.dot(ps)
        pjps = pj.dot(ps)
        return pipj/(pips*pjps)
    
    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, **opts ):

        if higher_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the higher phase-space point." )
        if momenta_dict is None:
            raise CurrentImplementationError(
                self.name() + " requires a momenta dictionary." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Perform mapping
        if not self.mapping_singular_structure.legs:
            self.mapping_singular_structure.legs = get_recoilers(reduced_process)
        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, self.mapping_singular_structure, momenta_dict,
            compute_jacobian=True )
        Q = mapping_vars['Q']
        jacobian = mapping_vars['jacobian']

        soft_leg_number = self.leg_numbers_map[1]
        pS = higher_PS_point[soft_leg_number]

        # Include the counterterm only in a part of the phase space
        if self.is_cut(Q=Q, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config,
                reduced_kinematics=(None, lower_PS_point))

        # Normalization factors
        norm = -4. * math.pi * alpha_s
        norm *= self.factor(Q=Q, pS=pS)
        norm /= jacobian

        # Find all colored leg numbers in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color') == 1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        # Initialize the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [],
            'reduced_kinematics': [(None, lower_PS_point)],
            'values': {}
        })

        # Loop over colored parton number pairs (a, b)
        # and add the corresponding contributions to this current
        color_correlation_index = 0
        for i, a in enumerate(all_colored_parton_numbers):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in all_colored_parton_numbers[i:]:
                # Write the eikonal for that pair
                if a != b:
                    mult_factor = 2.
                else:
                    mult_factor = 1.
                pa = higher_PS_point[a]
                pb = higher_PS_point[b]
                eikonal = self.eikonal(pa, pb, pS)
                evaluation['color_correlations'].append( ((a, b), ) )
                evaluation['values'][(0, color_correlation_index, 0)] = {
                    'finite': norm * mult_factor * eikonal }
                color_correlation_index += 1
        
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())) )
        return result

#=========================================================================================
# NLO final soft-collinear currents
#=========================================================================================

class QCD_final_softcollinear_0_gX(currents.QCDLocalSoftCollinearCurrent):
    """NLO tree-level (final) soft-collinear currents."""
    # TODO Ideally set recoilers for mapping_singular_structure
    # within does_implement_this_current and __init__

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(2,  1, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_g = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(2, 21, sub.SubtractionLeg.FINAL), ) )
    structure_q = sub.SingularStructure(substructures=(soft_coll_structure_q, ))
    structure_g = sub.SingularStructure(substructures=(soft_coll_structure_g, ))
    mapping = mappings.SoftVsFinalPureRescalingMapping
    factor = staticmethod(factors_and_cuts.factor_soft)
    is_cut = staticmethod(factors_and_cuts.cut_soft)

    @classmethod
    def does_implement_this_current(cls, current, model):

        if not all([
            current.get('squared_orders') == {'QCD': 2},
            current.get('n_loops') == 0,
            current.get('resolve_mother_spin_and_color') == True,
            not isinstance(current, (sub.BeamCurrent, sub.IntegratedBeamCurrent, sub.IntegratedCurrent))
        ]):
            return None

        color_charge = 'CF'
        leg_numbers_map = cls.structure_q.map_leg_numbers(
            current.get('singular_structure'), [range(1, model.get_nflav()+1)])
        if leg_numbers_map is None:
            color_charge = 'CA'
            leg_numbers_map = cls.structure_g.map_leg_numbers(
                current.get('singular_structure'), [range(1, model.get_nflav()+1)])
            if leg_numbers_map is None:
                return None
        mapping_soft_structure = sub.SoftStructure(
            legs=(sub.SubtractionLeg(leg_numbers_map[1], 21, sub.SubtractionLeg.FINAL),))
        mapping_singular_structure = sub.SingularStructure(
            substructures=[mapping_soft_structure, ])
        return {
            'leg_numbers_map': leg_numbers_map,
            'color_charge': color_charge,
            'mapping_singular_structure': mapping_singular_structure
        }

    def __init__(self, *args, **opts):

        for opt_name in ['leg_numbers_map', 'color_charge', 'mapping_singular_structure']:
            try:
                setattr(self, opt_name, opts.pop(opt_name))
            except KeyError:
                raise CurrentImplementationError(
                    "__init__ of " + self.__class__.__name__ + " requires " + opt_name)

        super(QCD_final_softcollinear_0_gX, self).__init__(*args, **opts)
        # At this state color_charge is the string of the group factor ('CA' or 'CF');
        # now that the super constructor has been called,
        # the group factors have been initialized and we can retrieve them.
        self.color_charge = getattr(self, self.color_charge)

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, **opts ):
        if higher_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the higher phase-space point." )
        if momenta_dict is None:
            raise CurrentImplementationError(
                self.name() + " requires a momenta dictionary." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        soft_leg_number = self.leg_numbers_map[1]
        coll_leg_number = self.leg_numbers_map[2]
        children = (soft_leg_number, coll_leg_number, )
        parent = momenta_dict.inv[frozenset(children)]
        pC = higher_PS_point[soft_leg_number] + higher_PS_point[coll_leg_number]
        pS = higher_PS_point[soft_leg_number]

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation.zero()

        # Perform mapping
        if not self.mapping_singular_structure.legs:
            self.mapping_singular_structure.legs = get_recoilers(
                reduced_process, excluded=(parent, ))
        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, self.mapping_singular_structure, momenta_dict,
            compute_jacobian=True )
        Q = mapping_vars['Q']
        jacobian = mapping_vars['jacobian']

        # Include the counterterm only in a part of the phase space
        if self.is_cut(Q=Q, pC=pC, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Evaluate kernel
        zs, kTs = variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        z = zs[0]
        evaluation['values'][(0, 0, 0)]['finite'] = self.color_charge * 2. * (1.-z) / z

        # Add the normalization factors
        s12 = pC.square()
        norm = 8. * math.pi * alpha_s / s12
        norm *= self.factor(Q=Q, pC=pC, pS=pS)
        norm /= jacobian
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result
