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
"""Implementation of NLO distributed soft currents."""

import math

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings

import commons.utils as utils
import commons.QCD_local_currents as currents
# import commons.factors_and_cuts as factors_and_cuts

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

variables = currents.n_final_coll_variables
mapping = mappings.FinalGroupingMapping
# soft_mapping = mappings.SoftVsFinalPureRescalingMapping
# soft_coll_mapping = mappings.SoftCollinearVsFinalMapping(
#    soft_mapping=soft_mapping, collinear_mapping=coll_mapping)
divide_by_jacobian = True
# factor_coll = factors_and_cuts.factor_coll
# factor_soft = factors_and_cuts.factor_soft
# is_cut_coll = factors_and_cuts.cut_coll
# is_cut_soft = factors_and_cuts.cut_soft
no_cut = currents.no_cut
no_factor = currents.no_factor

#=========================================================================================
# NLO final-collinear currents
#=========================================================================================

class QCD_final_collinear_0_XX(currents.QCDLocalCollinearCurrent):
    """Two-collinear tree-level current."""

    squared_orders = {'QCD': 2}
    n_loops = 0

    is_cut = staticmethod(no_cut)
    factor = staticmethod(no_factor)
    mapping = mapping
    variables = staticmethod(variables)
    divide_by_jacobian = divide_by_jacobian
    get_recoilers = staticmethod(get_recoilers)


class QCD_final_collinear_with_soft(QCD_final_collinear_0_XX):
    """Generic class subtracting all the NLO counterterms
     that diverge when a pair of particles become unresolved"""

    @staticmethod
    def partial_fractionned_eikonal(pi, pj, pr):
        """Partial fractionned eikonal factor for soft particle with momentum pr
        emitted from the dipole with momenta pi and pj.
        This is obtained starting from the eikonal and:
        - ignoring 1 / sir, which is already included in the normalisation factor;
        - multiplying by the partial fraction sjr / (sir + sjr) to regulate for sjr -> 0.
        """

        pipj = pi.dot(pj)
        pijpr = pr.dot(pi + pj)
        return 2 * pipj / pijpr

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, **opts):

        if higher_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the higher phase-space point.")
        if momenta_dict is None:
            raise CurrentImplementationError(
                self.name() + " requires a momenta dictionary.")
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment.")

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Include the counterterm only in a part of the phase space
        # Retrieve leg numbers
        children = tuple(self.leg_numbers_map[i]
                         for i in sorted(self.leg_numbers_map.keys()))
        parent = momenta_dict.inv[frozenset(children)]

        # Perform mapping
        self.mapping_singular_structure.legs = self.get_recoilers(
            reduced_process, excluded=(parent, ) )
        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, self.mapping_singular_structure, momenta_dict,
            compute_jacobian=self.divide_by_jacobian)

        # Retrieve kinematics
        Q = mapping_vars['Q']
        pC = sum(higher_PS_point[child] for child in children)
        qC = lower_PS_point[parent]
        if self.is_cut(Q=Q, pC=pC):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)
        reduced_kinematics = (None, lower_PS_point)

        # Evaluate collinear subtracted kernel
        zs, kTs = self.variables(higher_PS_point, qC, children, Q=Q)
        evaluation = self.kernel(zs, kTs, parent, reduced_kinematics)

        # Start handling the soft counterterms
        # First check which particles can go soft

        gluons = [self.is_gluon(l,self.model) for l in self.structure.get_all_legs()]
        # If any can go soft, we need to add eikonal pieces
        if any(gluons):
            # Find all colored leg numbers ?except for the parent? in the reduced process
            all_colored_parton_numbers = []
            for leg in reduced_process.get('legs'):
                if self.model.get_particle(leg.get('id')).get('color') == 1:
                    continue
                all_colored_parton_numbers.append(leg.get('number'))

            color_correlation_index = 1
            p0 = higher_PS_point[children[0]]
            p1 = higher_PS_point[children[1]]

            # Now loop over the colored parton number pairs (parent, j)
            # and add the corresponding contributions to this current
            # NB: we are using a (0,1) collinear mapping so in the S(0) limit, parent=p1
            # and in the S(1) limit, parent = p0.
            # As a result, S(0) emitted from the dipole (1,j) has the unresolved color correlation (parent,j)
            # as well as S(1) emitted from the dipole (0,j). As a result, for each j, we have a single color
            # correlation for the two eikonal pieces pj.p1/(pj.p0)/(pj.p0+p1.p0) and pj.p0/(pj.p1)/(pj.p1+p0.p1)

            for j in all_colored_parton_numbers:
                # Write the eikonal for that pair
                if j == parent:
                    # j is already the other leg of the dipole,
                    # we skip as we don't handle massive emitters
                    continue
                pj = higher_PS_point[j]
                evaluation['color_correlations'].append(((parent, j),))
                eiks = 0
                if gluons[0]:
                    # Soft 0 emitted from (1,j) with C(0,j) screened
                    eiks -= self.partial_fractionned_eikonal(p1, pj, p0)
                if gluons[1]:
                    # Soft 1 emitted from (0,j) with C(1,j) screened
                    eiks -= self.partial_fractionned_eikonal(p0, pj, p1)
                evaluation['values'][(0, color_correlation_index,0)] = {'finite': eiks}
                color_correlation_index += 1

        # Add the normalization factors
        pC2 = pC.square()
        norm = 8. * math.pi * alpha_s / pC2
        norm *= self.factor(Q=Q, pC=pC, qC=qC)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result


class QCD_final_collinear_0_qqx(QCD_final_collinear_with_soft):
    """q q~ collinear tree-level current."""

    structure = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, +1, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL), ))

    def kernel(self, zs, kTs, parent, reduced_kinematics):

        # Retrieve the collinear variables
        z = zs[0]
        kT = kTs[0]
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None, ((parent, (kT,)),), ],
            'color_correlations': [None],
            'reduced_kinematics': [reduced_kinematics],
            'values': {(0, 0, 0): {'finite': None},
                       (1, 0, 0): {'finite': None}, }
        })
        # Compute the kernel
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        evaluation['values'][(0, 0, 0)]['finite'] = self.TR
        evaluation['values'][(1, 0, 0)]['finite'] = 4 * self.TR * z * (1-z) / kT.square()
        return evaluation


class QCD_final_collinear_0_gq(QCD_final_collinear_with_soft):
    """g q collinear tree-level current.
    This contains
    - the hard-collinear kernel for g(0)q(1) being collinear
    - the pieces of the eikonal terms for g(0) soft emitted from q(1)X(j) which also diverge
    in the g(0)q(1) collinear limit
    """

    structure = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL), ))

    def kernel(self, zs, kTs, parent, reduced_kinematics):

        # Retrieve the collinear variables
        z = zs[0]
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [None],
            'reduced_kinematics': [reduced_kinematics],
            'values': {(0, 0, 0): {'finite': None}}
        })
        # Compute the kernel using
        # f9d0839fc58905d67367e3e67efabee05ee390f9:madgraph/iolibs/template_files/OLD_subtraction/cataniseymour/NLO/local_currents.py:146
        evaluation['values'][(0, 0, 0)]['finite'] = \
            self.CF*z
        return evaluation


class QCD_final_collinear_0_gg(QCD_final_collinear_with_soft):
    """g g collinear tree-level current."""

    structure = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL),
        sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), ))

    def kernel(self, zs, kTs, parent, reduced_kinematics):

        # Retrieve the collinear variables
        z = zs[0]
        kT = kTs[0]
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None, ((parent, (kT,)),), ],
            'color_correlations': [None],
            'reduced_kinematics': [reduced_kinematics],
            'values': {(0, 0, 0): {'finite': None},
                       (1, 0, 0): {'finite': None}, }
        })
        # Compute the kernel
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        evaluation['values'][(0, 0, 0)]['finite'] = \
            -4. * self.CA
        evaluation['values'][(1, 0, 0)]['finite'] = \
            -2. * self.CA * 2. * z * (1. - z) / kT.square()
        return evaluation

#=========================================================================================
# NLO soft current
#=========================================================================================

class NoSoft(currents.QCDLocalCurrent):
    """There is no explicit soft counterterm in this scheme. It is distributed into collinear counterterms"""

    is_zero = True
    squared_orders = {'QCD': 2}
    n_loops = 0
    structure = sub.SingularStructure(
        sub.SoftStructure(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL)) )

    
    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, **opts ):

        # Just return 0
        result = utils.SubtractionCurrentResult()
        result.add_result(
            utils.SubtractionCurrentEvaluation.zero(),
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())) )
        return result

#=========================================================================================
# NLO final soft-collinear currents
#=========================================================================================

class NoSoftCollinear(currents.QCDLocalCurrent):
    """There is no explicit soft counterterm in this scheme. It is distributed into collinear counterterms"""

    is_zero=True
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

        super(NoSoftCollinear, self).__init__(*args, **opts)
        # At this state color_charge is the string of the group factor ('CA' or 'CF');
        # now that the super constructor has been called,
        # the group factors have been initialized and we can retrieve them.
        self.color_charge = getattr(self, self.color_charge)
    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, **opts ):

        # Just return 0
        result = utils.SubtractionCurrentResult()
        result.add_result(
            utils.SubtractionCurrentEvaluation.zero(),
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())) )
        return result