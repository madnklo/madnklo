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

#=========================================================================================
# Variables, mappings, jacobians, factors and cuts
#=========================================================================================
# Note that variables, factors and cuts will be class members by design
# so they can easily be overridden by subclasses.
# They will be taken from the following variables
# so we can quickly switch them coherently across the entire subtraction scheme.

variables = currents.Q_final_coll_variables
coll_mapping = mappings.FinalRescalingOneMapping
soft_mapping = mappings.SoftVsFinalPureRescalingMapping
soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
    soft_mapping=soft_mapping, collinear_mapping=coll_mapping)
divide_by_jacobian = True
factor_coll = factors_and_cuts.factor_coll
factor_soft = factors_and_cuts.factor_soft
is_cut_coll = factors_and_cuts.cut_coll
is_cut_soft = factors_and_cuts.cut_soft

#=========================================================================================
# NLO final-collinear currents
#=========================================================================================

class QCD_final_collinear_0_XX(currents.QCDLocalCollinearCurrent):
    """Two-collinear tree-level current."""

    squared_orders = {'QCD': 2}
    n_loops = 0

    is_cut = staticmethod(is_cut_coll)
    factor = staticmethod(factor_coll)
    mapping = coll_mapping
    variables = staticmethod(variables)
    divide_by_jacobian = divide_by_jacobian
    get_recoilers = staticmethod(get_recoilers)


class QCD_final_collinear_0_qqx(QCD_final_collinear_0_XX):
    """q q~ collinear tree-level current."""

    structure = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL), ))

    def kernel(self, evaluation, parent, zs, kTs):

        # Retrieve the collinear variables
        z = zs[0]
        kT = kTs[0]
        # Instantiate the structure of the result
        evaluation['spin_correlations'] = [None, ((parent, (kT,)),), ]

        # Compute the kernel
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        evaluation['values'][(0, 0, 0)] = { 'finite' : self.TR }
        evaluation['values'][(1, 0, 0)] = { 'finite' : 4 * self.TR * z * (1-z) / kT.square() }
        return evaluation


class QCD_final_collinear_0_gq(QCD_final_collinear_0_XX):
    """g q collinear tree-level current."""

    structure = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL), ))

    def kernel(self, evaluation, parent, zs, kTs):

        # Retrieve the collinear variables
        z = zs[0]

        # Compute the kernel
        evaluation['values'][(0, 0, 0)] = { 'finite' : \
            self.CF * (1 + (1-z)**2) / z }

        return evaluation


class QCD_final_collinear_0_gg(QCD_final_collinear_0_XX):
    """g g collinear tree-level current."""

    structure = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), ))

    def kernel(self, evaluation, parent, zs, kTs):

        # Retrieve the collinear variables
        z = zs[0]
        kT = kTs[0]

        # Instantiate the structure of the result
        evaluation['spin_correlations'] = [None, ((parent, (kT,)),), ]

        # Compute the kernel
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        evaluation['values'][(0, 0, 0)] = { 'finite' : \
            +2 * self.CA * (z/(1-z) + (1-z)/z) }
        evaluation['values'][(1, 0, 0)] = { 'finite' : \
            -2 * self.CA * 2 * z * (1-z) / kT.square() }
        return evaluation

#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_soft_0_g(currents.QCDLocalCurrent):
    """Soft gluon eikonal current at tree level, eq.4.12-4.13 of arXiv:0903.1218."""

    squared_orders = {'QCD': 2}
    n_loops = 0

    structure = sub.SingularStructure(
        sub.SoftStructure(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL)) )

    is_cut = staticmethod(is_cut_soft)
    factor = staticmethod(factor_soft)
    mapping = soft_mapping
    divide_by_jacobian = divide_by_jacobian
    get_recoilers = staticmethod(get_recoilers)

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

        # Retrieve leg numbers
        soft_leg_number = self.leg_numbers_map[0]

        # Perform mapping
        self.mapping_singular_structure.legs = self.get_recoilers(reduced_process)
        lower_PS_point, mapping_vars = soft_mapping.map_to_lower_multiplicity(
            higher_PS_point, self.mapping_singular_structure, momenta_dict,
            compute_jacobian=self.divide_by_jacobian )

        # Retrieve kinematics
        Q = mapping_vars['Q']
        pS = higher_PS_point[soft_leg_number]
        jacobian = mapping_vars.get('jacobian', 1)

        # Include the counterterm only in a part of the phase space
        if self.is_cut(Q=Q, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config,
                reduced_kinematics=(None, lower_PS_point))

        # Normalization factors
        norm = -4 * math.pi * alpha_s
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

class QCD_final_softcollinear_0_gX(currents.QCDLocalCurrent):
    """NLO tree-level (final) soft-collinear currents."""

    squared_orders = {'QCD': 2}
    n_loops = 0

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(1,  1, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_g = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), ) )
    structure_q = sub.SingularStructure(substructures=(soft_coll_structure_q, ))
    structure_g = sub.SingularStructure(substructures=(soft_coll_structure_g, ))

    is_cut = staticmethod(is_cut_soft)
    factor = staticmethod(factor_soft)
    mapping = soft_coll_mapping
    variables = staticmethod(variables)
    divide_by_jacobian = divide_by_jacobian
    get_recoilers = staticmethod(get_recoilers)

    @classmethod
    def does_implement_this_current(cls, current, model):

        if not cls.check_current_properties(current): return None

        color_charge = 'CF'
        leg_numbers_map = cls.structure_q.map_leg_numbers(
            current.get('singular_structure'), [range(1, model.get_nflav()+1)])
        if leg_numbers_map is None:
            color_charge = 'CA'
            leg_numbers_map = cls.structure_g.map_leg_numbers(
                current.get('singular_structure'), [range(1, model.get_nflav()+1)])
            if leg_numbers_map is None:
                return None
        mapping_singular_structure = current.get('singular_structure').get_copy()
        return {
            'leg_numbers_map': leg_numbers_map,
            'color_charge': color_charge,
            'mapping_singular_structure': mapping_singular_structure
        }

    def __init__(self, *args, **opts):

        for opt_name in ['color_charge', ]:
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

        # Retrieve leg numbers
        soft_leg_number = self.leg_numbers_map[0]
        coll_leg_number = self.leg_numbers_map[1]
        children = (soft_leg_number, coll_leg_number, )
        parent = momenta_dict.inv[frozenset(children)]

        # Perform mapping
        self.mapping_singular_structure.legs = self.get_recoilers(
            reduced_process, excluded=(parent, ))
        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, self.mapping_singular_structure, momenta_dict,
            compute_jacobian=self.divide_by_jacobian )

        # Retrieve kinematics
        Q = mapping_vars['Q']
        pS = higher_PS_point[soft_leg_number]
        pC = pS + higher_PS_point[coll_leg_number]
        jacobian = mapping_vars.get('jacobian', 1)

        # Include the counterterm only in a part of the phase space
        if self.is_cut(Q=Q, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config,
                reduced_kinematics=(None, lower_PS_point))

        # Evaluate kernel
        zs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        z = zs[0]
        evaluation = utils.SubtractionCurrentEvaluation.zero(
            reduced_kinematics=(None, lower_PS_point))
        evaluation['values'][(0, 0, 0)]['finite'] = self.color_charge * 2 * (1-z) / z

        # Add the normalization factors
        s12 = pC.square()
        norm = 8 * math.pi * alpha_s / s12
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
