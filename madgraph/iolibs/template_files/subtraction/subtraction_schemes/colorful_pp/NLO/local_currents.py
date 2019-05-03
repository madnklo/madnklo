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
"""Implementation of NLO colorful_pp local currents."""

import math

import madgraph.core.subtraction as sub
import madgraph.integrator.mappings as mappings

import commons.utils as utils
import commons.QCD_local_currents as currents
import colorful_pp_config

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_soft_0_g(currents.QCDLocalCurrent):
    """Soft gluon eikonal current at tree level, eq.4.12-4.13 of arXiv:0903.1218.
    This is modified with respect to the above reference because we replace all the two legs of the dipole of each
    eikonal by their mapped version.
    """

    is_cut = staticmethod(colorful_pp_config.cut_soft)
    factor = staticmethod(colorful_pp_config.factor_soft)
    get_recoilers = staticmethod(colorful_pp_config.get_recoilers)

    squared_orders = {'QCD': 2}
    n_loops = 0

    structure = sub.SingularStructure(
        sub.SoftStructure(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL)) )

    mapping = colorful_pp_config.soft_mapping
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian
   
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
        higher_PS_point = None, momenta_dict = None, reduced_process = None,
        hel_config = None, Q=None, **opts ):

        if higher_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping." )
        if momenta_dict is None:
            raise CurrentImplementationError(
                self.name() + " requires a momentum routing dictionary." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )

        """Important note about the IF CS:
        - in this scheme we want to use "tilded" momenta for the dipole legs in eikonals. This is explicitly implemented in the soft local current
        - this implies that the correct form for the local C(ir)S(r) taken as the collinear limit of the eikonals is 
        1/ (p_r + p_i_tilde)^2 (1-z_r)/z_r where z_r = p_r.Q/(p_r+p_i_tilde).Q
        - Specializing to the case where the collinear partner of the soft particle is an initial state particle (i = a ), we have 
        p_a_tilde = xi p_a and 2p_a.Q = Q^2 so that the soft-collinear takes the form
        1/(p_r+xi p_a)^2 * xi/y_rQ where y_rQ is the usual Hungarian variable
        this simplifies to 
        1/(p_r+p_a)^2 * 1/y_rQ which is exactly the soft collinear as computed *without* tilded variables (i.e. exactly eq.5.29 of 0903.1218)
        
        As a result we use exactly the same way of evaluating the counterterms as in honest-to-god colorful.
        """

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Now find all colored leg numbers in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color')==1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))
        soft_leg_number = self.leg_numbers_map[0]

        pS = higher_PS_point[soft_leg_number]

        # Perform mapping
        this_mapping_singular_structure = self.mapping_singular_structure.get_copy()
        this_mapping_singular_structure.legs = self.get_recoilers(reduced_process)
        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, this_mapping_singular_structure, momenta_dict,
            compute_jacobian=self.divide_by_jacobian )
        reduced_kinematics = (None, lower_PS_point)

        # Include the counterterm only in a part of the phase space
        if self.is_cut(Q=Q, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config, reduced_kinematics=('IS_CUT', lower_PS_point))

        # Retrieve kinematics
        pS = higher_PS_point[soft_leg_number]
        jacobian = mapping_vars.get('jacobian', 1.)

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [],
            'reduced_kinematics'  : [reduced_kinematics, ],
            'values'              : {}
        })
        
        # Normalization factors
        norm = -4. * math.pi * alpha_s
        norm *= self.factor(Q=Q, pS=pS)
        if self.divide_by_jacobian:
            norm /= jacobian

        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        for i, a in enumerate(all_colored_parton_numbers):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in all_colored_parton_numbers[i:]:
                # Write the eikonal for that pair
                if a!=b:
                    mult_factor = 2.
                else:
                    mult_factor = 1.
                #pa = sum(higher_PS_point[child] for child in momenta_dict[a])
                #pb = sum(higher_PS_point[child] for child in momenta_dict[b])
                pa = lower_PS_point[a]
                pb = lower_PS_point[b]
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
    """NLO tree-level (final) soft-collinear currents. The momenta used in this current are the mapped momenta from the soft mapping."""

    is_cut = staticmethod(colorful_pp_config.cut_soft)
    factor = staticmethod(colorful_pp_config.factor_soft)
    get_recoilers = staticmethod(colorful_pp_config.get_recoilers)
    variables = staticmethod(colorful_pp_config.final_soft_coll_variables)

    squared_orders = {'QCD': 2}
    n_loops = 0

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(0,  1, sub.SubtractionLeg.FINAL), ) )
    soft_coll_structure_g = sub.CollStructure(
        substructures=(soft_structure, ),
        legs=(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.FINAL), ) )
    structure_q = sub.SingularStructure(substructures=(soft_coll_structure_q, ))
    structure_g = sub.SingularStructure(substructures=(soft_coll_structure_g, ))

    mapping = colorful_pp_config.final_soft_collinear_mapping
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    def __init__(self, *args, **opts):

        # Make sure it is initialized with the proper set of options and remove them
        # before calling the mother constructor
        if 'color_charge' not in opts:
            raise CurrentImplementationError(
                "The current '%s' must be instantiated with " % self.__class__.__name__ +
                "a 'color_charge' option specified.")
        color_charge = opts.pop('color_charge')

        super(QCD_final_softcollinear_0_gX, self).__init__(*args, **opts)
        # At this state color_charge is the string of the group factor ('CA' or 'CF');
        # now that the mother constructor has been called,
        # the group factors have been initialized and we can retrieve them.
        self.color_charge = getattr(self, color_charge)

    @classmethod
    def does_implement_this_current(cls, current, model):

        if not cls.check_current_properties(current): return None

        leg_numbers_map = cls.structure_q.map_leg_numbers(
            current.get('singular_structure'), [range(1, model.get_nflav() + 1)])

        if leg_numbers_map is not None:
            color_charge = 'CF'
        else:
            leg_numbers_map = cls.structure_g.map_leg_numbers(
                current.get('singular_structure'), [range(1, model.get_nflav() + 1)])
            if leg_numbers_map is not None:
                color_charge = 'CA'

        if leg_numbers_map is None:
            return None

        #mapping_singular_structure = sub.SingularStructure(sub.SoftStructure(
        #    legs=(sub.SubtractionLeg(leg_numbers_map[0], 21, sub.SubtractionLeg.FINAL),)))
        mapping_singular_structure = current.get('singular_structure').get_copy()
        return {
            'leg_numbers_map': leg_numbers_map,
            'color_charge': color_charge,
            'mapping_singular_structure': mapping_singular_structure
        }

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, Q=None, **opts ):

        if higher_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before mapping." )
        if momenta_dict is None:
            raise CurrentImplementationError(
                self.name() + " requires a momentum routing dictionary." )
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        children = tuple(self.leg_numbers_map[i]
                         for i in sorted(self.leg_numbers_map.keys()))
        parent = momenta_dict.inv[frozenset(children)]

        # Perform mapping
        this_mapping_singular_structure = self.mapping_singular_structure.get_copy()
        this_mapping_singular_structure.legs = self.get_recoilers(reduced_process, excluded=(parent, ))
        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, this_mapping_singular_structure, momenta_dict,
            compute_jacobian=self.divide_by_jacobian )
        reduced_kinematics = (None, lower_PS_point)

        # Include the counterterm only in a part of the phase space
        # children are the the set of particles that are going unresolved.
        # Here we have C and S going collinear with S soft.
        # The parent is the mapped C with a soft mapping, usually refered to as Ctilde.
        # S is just removed in a soft mapping.
        # Here S is a single particle but we obtain it as a list soft_children\
        # to illustrate how multiple softs would be obtained
        pCtilde = lower_PS_point[parent]
        soft_children = [ self.leg_numbers_map[1], ]

        pS = sum(higher_PS_point[child] for child in soft_children)
        if self.is_cut(Q=Q, pC=pCtilde, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config, reduced_kinematics=('IS_CUT', lower_PS_point))

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [None],
            'reduced_kinematics': [reduced_kinematics,],
            'values': {(0, 0, 0): {'finite': None}}
        })

        # Evaluate kernel
        zs = self.variables([pS,pCtilde],Q)
        z = zs[0]
        evaluation['values'][(0, 0, 0)]['finite'] = self.color_charge * 2.*(1.-z) / z

        # Add the normalization factors
        s12 = (pCtilde+pS).square()
        norm = 8. * math.pi * alpha_s / s12
        norm *= self.factor(Q=Q, pC=pCtilde, pS=pS)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

#=========================================================================================
# NLO initial-collinear currents
#=========================================================================================

class QCD_initial_collinear_0_XX(currents.QCDLocalCollinearCurrent):
    """Two-collinear tree-level current."""

    squared_orders = {'QCD': 2}
    n_loops = 0

    is_cut = staticmethod(colorful_pp_config.cut_initial_coll)
    factor = staticmethod(colorful_pp_config.factor_initial_coll)
    get_recoilers = staticmethod(colorful_pp_config.get_recoilers)
    mapping = colorful_pp_config.initial_coll_mapping
    variables = staticmethod(colorful_pp_config.initial_coll_variables)
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

class QCD_initial_collinear_0_qg(QCD_initial_collinear_0_XX):
    """qg collinear ISR tree-level current. q(initial) > q(initial_after_emission) g(final)"""

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), )),
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, -1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), )),
    ]

    def kernel(self, evaluation, parent, xs, kTs):
        # Retrieve the collinear variable x
        x = xs[0]

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorization formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = -1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= 1.

        z = 1. / x

        norm = initial_state_crossing_factor * self.CF
        # We re-use here the Altarelli-Parisi Kernel of the P_qg final state kernel
        evaluation['values'][(0, 0, 0)] = {'finite' : norm * ((1. + z ** 2) / (1. - z))}

        return evaluation

class QCD_initial_collinear_0_gq(QCD_initial_collinear_0_XX):
    """gq collinear ISR tree-level current. q(initial) > g(initial_after_emission) q(final)"""

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, +1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL), )),
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, -1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL), )),
    ]

    def kernel(self, evaluation, parent, xs, kTs):

        # Retrieve the collinear variable x
        x = xs[0]
        kT = kTs[0]

        # Set spin correlations
        evaluation['spin_correlations'] = [None, ((parent, (kT,)),), ]

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorization formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (gluon)
        initial_state_crossing_factor *= (self.NC**2-1)/float(self.NC)
        
        z = 1./x

        norm = initial_state_crossing_factor * self.TR
        # We re-use here the Altarelli-Parisi Kernel of the P_q\bar{q} final state kernel
        evaluation['values'][(0, 0, 0)] = { 'finite' : norm }
        evaluation['values'][(1, 0, 0)] = { 'finite' : 4. * norm * z*(1.-z) / kT.square() }

        return evaluation
    
class QCD_initial_collinear_0_qq(QCD_initial_collinear_0_XX):
    """qq collinear ISR tree-level current. g(initial) > q(initial_after_emission) qx(final)"""

    # Make sure to have the initial particle with the lowest index
    structure = [
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, 21, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, +1, sub.SubtractionLeg.FINAL), )),
        sub.SingularStructure(sub.CollStructure(
            sub.SubtractionLeg(0, 21, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(1, -1, sub.SubtractionLeg.FINAL), ))
    ]

#    Below is an example on how to efficiently debug 'does_implement_this_current'
#    @classmethod
#    def does_implement_this_current(cls, current, model):
#
#        res = super(QCD_initial_collinear_0_qq, cls).does_implement_this_current(current, model)
#
#        misc.sprint("For this structure: %s"%current.get('singular_structure').__str__(print_n=True, print_pdg=True, print_state=True, sort=True))
#        misc.sprint("I get:",res)
#        import pdb
#        pdb.set_trace()

    def kernel(self, evaluation, parent, xs, kTs):

        # Retrieve the collinear variable x
        x = xs[0]

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorization formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (gluon) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= (self.NC/float(self.NC**2-1))
        
        z = 1./x

        norm = initial_state_crossing_factor * self.CF
        # We re-use here the Altarelli-Parisi Kernel of the P_gq final state kernel
        evaluation['values'][(0, 0, 0)] = {'finite' : norm * ( (1.+(1.-z)**2)/z ) }

        return evaluation

class QCD_initial_collinear_0_gg(QCD_initial_collinear_0_XX):
    """gg collinear ISR tree-level current. g(initial) > g(initial_after_emission) g(final)"""


    # Make sure to have the initial particle with the lowest index
    structure = sub.SingularStructure(sub.CollStructure(
        sub.SubtractionLeg(0, 21, sub.SubtractionLeg.INITIAL),
        sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL), )),

    def kernel(self, evaluation, parent, xs, kTs):

        # Retrieve the collinear variable x
        x = xs[0]
        kT = kTs[0]

        # Set spin correlations
        evaluation['spin_correlations'] = [None, ((parent, (kT,)),), ]

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorization formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = -1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (gluon) and the one of the reduced Born ME (gluon)
        initial_state_crossing_factor *= 1.

        z = 1. / x

        # We re-use here the Altarelli-Parisi Kernel of the P_g\bar{g} final state kernel

        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        norm = initial_state_crossing_factor * 2. * self.CA
        evaluation['values'][(0, 0, 0)] = {'finite': norm * ((z / (1. - z)) + ((1. - z) / z))}
        evaluation['values'][(1, 0, 0)] = {'finite': -norm * 2. * z * (1. - z) / kT.square()}

        return evaluation

#=========================================================================================
# NLO soft initial-collinear currents
#=========================================================================================

class QCD_initial_softcollinear_0_Xg(currents.QCDLocalCurrent):
    """NLO tree-level (initial) soft-collinear currents."""

    is_cut = staticmethod(colorful_pp_config.cut_soft)
    factor = staticmethod(colorful_pp_config.factor_soft)
    get_recoilers = staticmethod(colorful_pp_config.get_recoilers)

    squared_orders = {'QCD': 2}
    n_loops = 0

    soft_structure = sub.SoftStructure(
        legs=(sub.SubtractionLeg(1, 21, sub.SubtractionLeg.FINAL),))
    soft_coll_structure_q = sub.CollStructure(
        substructures=(soft_structure,),
        legs=(sub.SubtractionLeg(0, 1, sub.SubtractionLeg.INITIAL),))
    soft_coll_structure_g = sub.CollStructure(
        substructures=(soft_structure,),
        legs=(sub.SubtractionLeg(0, 21, sub.SubtractionLeg.INITIAL),))
    structure_q = sub.SingularStructure(substructures=(soft_coll_structure_q,))
    structure_g = sub.SingularStructure(substructures=(soft_coll_structure_g,))

    mapping = colorful_pp_config.initial_soft_collinear_mapping
    divide_by_jacobian = colorful_pp_config.divide_by_jacobian

    variables = staticmethod(colorful_pp_config.initial_coll_variables)

    def __init__(self, *args, **opts):

        # Make sure it is initialized with the proper set of options and remove them
        # before calling the mother constructor
        if 'color_charge' not in opts:
            raise CurrentImplementationError(
                "The current '%s' must be instantiated with " % self.__class__.__name__ +
                "a 'color_charge' option specified.")
        color_charge = opts.pop('color_charge')

        super(QCD_initial_softcollinear_0_Xg, self).__init__(*args, **opts)
        # At this state color_charge is the string of the group factor ('CA' or 'CF');
        # now that the mother constructor has been called,
        # the group factors have been initialized and we can retrieve them.
        self.color_charge = getattr(self, color_charge)

    @classmethod
    def does_implement_this_current(cls, current, model):

        if not cls.check_current_properties(current): return None

        leg_numbers_map = cls.structure_q.map_leg_numbers(
            current.get('singular_structure'), [range(1, model.get_nflav() + 1)])

        if leg_numbers_map is not None:
            color_charge = 'CF'
        else:
            leg_numbers_map = cls.structure_g.map_leg_numbers(
                current.get('singular_structure'), [range(1, model.get_nflav() + 1)])
            if leg_numbers_map is not None:
                color_charge = 'CA'

        if leg_numbers_map is None:
            return None

        #soft_structure = sub.SingularStructure(sub.SoftStructure(
        #    legs=(sub.SubtractionLeg(leg_numbers_map[0], 21, sub.SubtractionLeg.FINAL),)))
        mapping_singular_structure = current.get('singular_structure').get_copy()
        return {
            'leg_numbers_map': leg_numbers_map,
            'color_charge': color_charge,
            'mapping_singular_structure': mapping_singular_structure
        }

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, Q=None, **opts ):

        if higher_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before mapping." )
        if momenta_dict is None:
            raise CurrentImplementationError(
                self.name() + " requires a momentum routing dictionary." )
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Perform mapping
        this_mapping_singular_structure = self.mapping_singular_structure.get_copy()
        this_mapping_singular_structure.legs = self.get_recoilers(reduced_process)
        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, this_mapping_singular_structure, momenta_dict,
            compute_jacobian=self.divide_by_jacobian )
        reduced_kinematics = (None, lower_PS_point)

        # Include the counterterm only in a part of the phase space
        children = tuple(self.leg_numbers_map[i]
                         for i in sorted(self.leg_numbers_map.keys()))
        pC = higher_PS_point[children[0]]
        pS = higher_PS_point[children[1]]
        pC = pC + pS
        parent = momenta_dict.inv[frozenset(children)]
        if self.is_cut(Q=Q, pC=pC, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config, reduced_kinematics=('IS_CUT',lower_PS_point))

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [None],
            'reduced_kinematics': [reduced_kinematics ],
            'values': {(0, 0, 0): {'finite': None}}
        })

        # Evaluate kernel
        xs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        x = xs[0]
        # See Eq. (4.17) of NNLO compatible NLO scheme publication arXiv:0903.1218v2

        # There is no need for the ratio of color-averaging factor between the real ME
        # initial state flavor and the one of the reduced Born ME as they are either both
        # gluons or both quarks
        
        evaluation['values'][(0, 0,0)]['finite'] = self.color_charge * ( 2. / (1. - x) )

        # Add the normalization factors
        s12 = pC.square()
        norm = 8. * math.pi * alpha_s / s12
        norm *= self.factor(Q=Q, pC=pC, pS=pS)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result
