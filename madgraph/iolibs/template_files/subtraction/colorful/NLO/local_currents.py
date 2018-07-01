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

import os
import math

try:
    # First try to import this in the context of the exported currents
    import SubtractionCurrents.subtraction_current_implementations_utils as utils
    import SubtractionCurrents.QCD_local_currents as currents
except ImportError:
    # If not working, then it must be within MG5_aMC context:
    import madgraph.iolibs.template_files.\
                   subtraction.subtraction_current_implementations_utils as utils
    import madgraph.iolibs.template_files.\
                   subtraction.QCD_local_currents as currents

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# NLO final-collinear currents
#=========================================================================================

class QCD_final_collinear_0_qqx(currents.QCDLocalCollinearCurrent):
    """q q~ collinear tree-level current."""

    is_cut = staticmethod(currents.SomogyiChoices.cut_coll)
    factor = staticmethod(currents.SomogyiChoices.factor_coll)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that the particles are a massless quark and its anti-quark in final-state
        if len(ss.legs) != 2: return None
        for leg in ss.legs:
            if not cls.is_quark(leg, model): return None
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
        if not cls.are_antiparticles(ss.legs[0], ss.legs[1]): return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        return tuple(leg.n for leg in legs)

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables
        z = zs[0]
        kT = kTs[0]
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None, ((parent, (kT, )), ), ],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None},
                                    (1, 0): {'finite': None}, }
        })
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        evaluation['values'][(0, 0)]['finite'] = self.TR
        evaluation['values'][(1, 0)]['finite'] = 4. * self.TR * z*(1.-z) / kT.square()
        return evaluation

class QCD_final_collinear_0_gq(currents.QCDLocalCollinearCurrent):
    """g q collinear tree-level current."""

    is_cut = staticmethod(currents.SomogyiChoices.cut_coll)
    factor = staticmethod(currents.SomogyiChoices.factor_coll)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that all particles are massless final state
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
        # Check that there are a quark and a gluon
        if len(ss.legs) != 2: return None
        if (cls.is_gluon(ss.legs[0], model) and cls.is_quark(ss.legs[1], model)):
            pass
        elif (cls.is_quark(ss.legs[0], model) and cls.is_gluon(ss.legs[1], model)):
            pass
        else:
            return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        if cls.is_gluon(legs[0], model): return (legs[0].n, legs[1].n)
        else: return (legs[1].n, legs[0].n)

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables
        z = zs[0]
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None}}
        })
        evaluation['values'][(0, 0)]['finite'] = self.CF * (1.+(1.-z)**2)/z
        return evaluation

class QCD_final_collinear_0_gg(currents.QCDLocalCollinearCurrent):
    """g g collinear tree-level current."""

    is_cut = staticmethod(currents.SomogyiChoices.cut_coll)
    factor = staticmethod(currents.SomogyiChoices.factor_coll)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that the particles are two final-state massless gluons
        if len(ss.legs) != 2: return None
        for leg in ss.legs:
            if not cls.is_gluon(leg, model): return None
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        return tuple(leg.n for leg in legs)

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables
        z = zs[0]
        kT = kTs[0]
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None, ((parent,( kT, )), ), ],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None},
                                    (1, 0): {'finite': None}, }
        })
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        evaluation['values'][(0, 0)]['finite'] =  2.*self.CA * ( (z/(1.-z)) + ((1.-z)/z) )
        evaluation['values'][(1, 0)]['finite'] = -2.*self.CA * 2.*z*(1.-z) / kT.square()
        return evaluation

#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_soft_0_g(currents.QCDLocalSoftCurrent):
    """Soft gluon eikonal current at tree level, eq.4.12-4.13 of arXiv:0903.1218."""

    is_cut = staticmethod(currents.SomogyiChoices.cut_soft)
    factor = staticmethod(currents.SomogyiChoices.factor_soft)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD soft tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve the singular structure
        singular_structure = current.get('singular_structure')
        # It should consist in exactly one legs going soft
        if len(singular_structure.legs) != 1:
            return None
        # The leg going soft should be a massless gluon in the final state
        leg = singular_structure.legs[0]
        if not cls.is_gluon(leg, model): return None
        if not cls.is_massless(leg, model): return None
        if cls.is_initial(leg): return None
        # The current is valid
        return init_vars
   
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
        higher_PS_point=None, lower_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):
        if higher_PS_point is None or lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping." )
        if leg_numbers_map is None:
            raise CurrentImplementationError(
                self.name() + " requires a leg numbers map, i.e. a momentum dictionary." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )

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
        soft_leg_number = current.get('singular_structure').legs[0].n

        pS = higher_PS_point[soft_leg_number]

        # Include the counterterm only in a part of the phase space
        if self.is_cut(Q=Q, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [],
            'values'              : {}
        })
        
        # Normalization factors
        norm = -4. * math.pi * alpha_s
        norm *= self.factor(Q=Q, pS=pS)

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
                pa = sum(higher_PS_point[child] for child in leg_numbers_map[a])
                pb = sum(higher_PS_point[child] for child in leg_numbers_map[b])
                # pa = lower_PS_point[a]
                # pb = lower_PS_point[b]
                eikonal = self.eikonal(pa, pb, pS)
                evaluation['color_correlations'].append( ((a, b), ) )
                evaluation['values'][(0, color_correlation_index)] = {
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

    is_cut = staticmethod(currents.SomogyiChoices.cut_soft)
    variables = staticmethod(currents.Q_final_coll_variables)
    factor = staticmethod(currents.SomogyiChoices.factor_soft)

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

        # Check the general properties common to NLO QCD
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve the singular structure
        singular_structure = current.get('singular_structure')
        # It should have only one leg and one nested substructure with one soft leg
        if len(singular_structure.legs) != 1: return None
        if len(singular_structure.substructures) != 1: return None
        sub_singular_structure = singular_structure.substructures[0]
        if len(sub_singular_structure.legs) != 1: return None
        # The hard and soft legs are now identified
        hard_leg = singular_structure.legs[0]
        soft_leg = sub_singular_structure.legs[0]
        # Make sure legs are massless and final state
        if not cls.is_massless(hard_leg, model): return None
        if not cls.is_massless(soft_leg, model): return None
        if cls.is_initial(hard_leg) or cls.is_initial(soft_leg): return None
        # The hard leg must be quark or a gluon
        if not (cls.is_gluon(hard_leg, model) or cls.is_quark(hard_leg, model)):
            return None
        # The soft leg must be a gluon
        if not cls.is_gluon(soft_leg, model): return None
        # Check if hard_leg is a quark or a gluon, and set the color factor accordingly
        if   cls.is_gluon(hard_leg, model): init_vars['color_charge'] = 'CA'
        elif cls.is_quark(hard_leg, model): init_vars['color_charge'] = 'CF'
        else: return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        ss = current.get('singular_structure')
        return (ss.substructures[0].legs[0].n, ss.legs[0].n)

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, lower_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):
        if higher_PS_point is None or lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping." )
        if leg_numbers_map is None:
            raise CurrentImplementationError(
                self.name() + " requires a leg numbers map, i.e. a momentum dictionary." )
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

        # Include the counterterm only in a part of the phase space
        children = self.get_sorted_children(current, self.model)
        pC = sum(higher_PS_point[child] for child in children)
        soft_children = []
        for substructure in current.get('singular_structure').substructures:
            soft_children += [leg.n for leg in substructure.get_all_legs()]
        pS = sum(higher_PS_point[child] for child in soft_children)
        parent = leg_numbers_map.inv[frozenset(children)]
        if self.is_cut(Q=Q, pC=pC, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [None],
            'values': {(0, 0): {'finite': None}}
        })

        # Evaluate kernel
        zs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        z = zs[0]
        evaluation['values'][(0, 0)]['finite'] = self.color_charge * 2.*(1.-z) / z

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

#=========================================================================================
# NLO initial-collinear currents
#=========================================================================================

class QCD_initial_collinear_0_qg(currents.QCDLocalCollinearCurrent):
    """qg collinear ISR tree-level current. q(initial) > q(initial_after_emission) g(final)"""

    is_cut = staticmethod(currents.SomogyiChoices.cut_coll)
    factor = staticmethod(currents.SomogyiChoices.factor_coll)
    variables = staticmethod(currents.Q_initial_coll_variables)
    
    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        
        # Retrieve singular structure
        ss = current.get('singular_structure')

        # Check that the particles are a massless quark and its anti-quark in final-state
        if len(ss.legs) != 2: return None
        
        n_initial_state_quarks = 0
        for leg in ss.legs:
            if cls.is_quark(leg, model) and cls.is_initial(leg):
                n_initial_state_quarks += 1
        if n_initial_state_quarks != 1: return None
        
        n_final_state_gluon = 0
        for leg in ss.legs:
            if cls.is_gluon(leg, model) and not cls.is_initial(leg):
                n_final_state_gluon += 1
        if n_final_state_gluon != 1: return None

        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        # Always put the initial state child first
        children_numbers = [leg.n for leg in legs if leg.state == leg.INITIAL]
        # Then the final state ones
        children_numbers.extend([leg.n for leg in legs if leg.state == leg.FINAL])

        return tuple(children_numbers)

    def evaluate_kernel(self, xs, kTs, parent):

        # Retrieve the collinear variable x
        x = xs[0]
        
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None}}
        })
        
        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorisation formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = -1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= 1.
        
        z = 1./x

        norm = initial_state_crossing_factor * self.CF
        # We re-use here the Altarelli-Parisi Kernel of the P_qg final state kernel
        evaluation['values'][(0, 0)]['finite'] = norm * ( (1.+z**2)/(1.-z) )

        return evaluation

class QCD_initial_collinear_0_gq(currents.QCDLocalCollinearCurrent):
    """gq collinear ISR tree-level current. q(initial) > g(initial_after_emission) q(final)"""

    is_cut = staticmethod(currents.SomogyiChoices.cut_coll)
    factor = staticmethod(currents.SomogyiChoices.factor_coll)
    variables = staticmethod(currents.Q_initial_coll_variables)
    
    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        
        # Retrieve singular structure
        ss = current.get('singular_structure')

        # Check that the particles are a massless quark and its anti-quark in final-state
        if len(ss.legs) != 2: return None
        
        n_initial_state_quarks = 0
        for leg in ss.legs:
            if cls.is_quark(leg, model) and cls.is_initial(leg):
                n_initial_state_quarks += 1
        if n_initial_state_quarks != 1: return None
        
        n_final_state_quarks = 0
        for leg in ss.legs:
            if cls.is_quark(leg, model) and not cls.is_initial(leg):
                n_final_state_quarks += 1
        if n_final_state_quarks != 1: return None

        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        # Always put the initial state child first
        children_numbers = [leg.n for leg in legs if leg.state == leg.INITIAL]
        # Then the final state ones
        children_numbers.extend([leg.n for leg in legs if leg.state == leg.FINAL])

        return tuple(children_numbers)

    def evaluate_kernel(self, xs, kTs, parent):

        # Retrieve the collinear variable x
        x = xs[0]
        kT = kTs[0]

        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None, ((parent, (kT, )), ), ],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None},
                                    (1, 0): {'finite': None}, }
        })

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorisation formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (gluon)
        initial_state_crossing_factor *= (self.NC**2-1)/float(self.NC)
        
        z = 1./x

        norm = initial_state_crossing_factor * self.TR
        # We re-use here the Altarelli-Parisi Kernel of the P_q\bar{q} final state kernel
        evaluation['values'][(0, 0)]['finite'] = norm
        evaluation['values'][(1, 0)]['finite'] = 4. * norm * z*(1.-z) / kT.square()

        return evaluation
    
class QCD_initial_collinear_0_qq(currents.QCDLocalCollinearCurrent):
    """qq collinear ISR tree-level current. g(initial) > q(initial_after_emission) qx(final)"""

    is_cut = staticmethod(currents.SomogyiChoices.cut_coll)
    factor = staticmethod(currents.SomogyiChoices.factor_coll)
    variables = staticmethod(currents.Q_initial_coll_variables)
    
    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        
        # Retrieve singular structure
        ss = current.get('singular_structure')

        # Check that the particles are a massless quark and its anti-quark in final-state
        if len(ss.legs) != 2: return None
        
        n_initial_state_gluons = 0
        for leg in ss.legs:
            if cls.is_gluon(leg, model) and cls.is_initial(leg):
                n_initial_state_gluons += 1
        if n_initial_state_gluons != 1: return None
        
        n_final_state_quarks = 0
        for leg in ss.legs:
            if cls.is_quark(leg, model) and not cls.is_initial(leg):
                n_final_state_quarks += 1
        if n_final_state_quarks != 1: return None

        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        # Always put the initial state child first
        children_numbers = [leg.n for leg in legs if leg.state == leg.INITIAL]
        # Then the final state ones
        children_numbers.extend([leg.n for leg in legs if leg.state == leg.FINAL])

        return tuple(children_numbers)

    def evaluate_kernel(self, xs, kTs, parent):

        # Retrieve the collinear variable x
        x = xs[0]
        
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None}}
        })
        
        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorisation formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (gluon) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= (self.NC/float(self.NC**2-1))
        
        z = 1./x

        norm = initial_state_crossing_factor * self.CF
        # We re-use here the Altarelli-Parisi Kernel of the P_gq final state kernel
        evaluation['values'][(0, 0)]['finite'] = norm * ( (1.+(1.-z)**2)/z )

        return evaluation

class QCD_initial_collinear_0_gg(currents.QCDLocalCollinearCurrent):
    """gg collinear ISR tree-level current. g(initial) > g(initial_after_emission) g(final)"""

    is_cut = staticmethod(currents.SomogyiChoices.cut_coll)
    factor = staticmethod(currents.SomogyiChoices.factor_coll)
    variables = staticmethod(currents.Q_initial_coll_variables)
    
    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        
        # Retrieve singular structure
        ss = current.get('singular_structure')

        # Check that the particles are a massless quark and its anti-quark in final-state
        if len(ss.legs) != 2: return None
        
        n_initial_state_gluons = 0
        for leg in ss.legs:
            if cls.is_gluon(leg, model) and cls.is_initial(leg):
                n_initial_state_gluons += 1
        if n_initial_state_gluons != 1: return None
        
        n_final_state_gluons = 0
        for leg in ss.legs:
            if cls.is_gluon(leg, model) and not cls.is_initial(leg):
                n_final_state_gluons += 1
        if n_final_state_gluons != 1: return None

        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        # Always put the initial state child first
        children_numbers = [leg.n for leg in legs if leg.state == leg.INITIAL]
        # Then the final state ones
        children_numbers.extend([leg.n for leg in legs if leg.state == leg.FINAL])

        return tuple(children_numbers)

    def evaluate_kernel(self, xs, kTs, parent):

        # Retrieve the collinear variable x
        x = xs[0]
        kT = kTs[0]

        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None, ((parent,( kT, )), ), ],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None},
                                    (1, 0): {'finite': None}, }
        })

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorisation formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = -1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (gluon) and the one of the reduced Born ME (gluon)
        initial_state_crossing_factor *= 1.
        
        z = 1./x

        # We re-use here the Altarelli-Parisi Kernel of the P_g\bar{g} final state kernel
        
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        norm = initial_state_crossing_factor * 2. * self.CA
        evaluation['values'][(0, 0)]['finite'] =  norm * ( (z/(1.-z)) + ((1.-z)/z) )
        evaluation['values'][(1, 0)]['finite'] = -norm * 2.*z*(1.-z) / kT.square()

        return evaluation
    
#=========================================================================================
# NLO soft initial-collinear currents
#=========================================================================================

class QCD_initial_softcollinear_0_Xg(currents.QCDLocalSoftCollinearCurrent):
    """NLO tree-level (initial) soft-collinear currents."""

    is_cut = staticmethod(currents.SomogyiChoices.cut_soft)
    factor = staticmethod(currents.SomogyiChoices.factor_soft)
    variables = staticmethod(currents.Q_initial_coll_variables)

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

        # Check the general properties common to NLO QCD
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve the singular structure
        singular_structure = current.get('singular_structure')
        # It should have only one leg and one nested substructure with one soft leg
        if len(singular_structure.legs) != 1: return None
        if len(singular_structure.substructures) != 1: return None
        sub_singular_structure = singular_structure.substructures[0]
        if len(sub_singular_structure.legs) != 1: return None
        # The hard and soft legs are now identified
        hard_leg = singular_structure.legs[0]
        soft_leg = sub_singular_structure.legs[0]
        # Make sure legs are massless and the hard one final state
        if not cls.is_massless(hard_leg, model): return None
        if not cls.is_massless(soft_leg, model): return None
        if not cls.is_initial(hard_leg): return None
        # The hard leg must be quark or a gluon
        if not (cls.is_gluon(hard_leg, model) or cls.is_quark(hard_leg, model)):
            return None
        # The soft leg must be a gluon
        if not cls.is_gluon(soft_leg, model): return None
        # Check if hard_leg is a quark or a gluon, and set the color factor accordingly
        if   cls.is_gluon(hard_leg, model): init_vars['color_charge'] = 'CA'
        elif cls.is_quark(hard_leg, model): init_vars['color_charge'] = 'CF'
        else: return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        ss = current.get('singular_structure')
        return (ss.legs[0].n, ss.substructures[0].legs[0].n)

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, lower_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):
        if higher_PS_point is None or lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping." )
        if leg_numbers_map is None:
            raise CurrentImplementationError(
                self.name() + " requires a leg numbers map, i.e. a momentum dictionary." )
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

        # Include the counterterm only in a part of the phase space
        children = self.get_sorted_children(current, self.model)
        pC = higher_PS_point[children[0]]
        pS = higher_PS_point[children[1]]
        pC -= pS
        parent = leg_numbers_map.inv[frozenset(children)]
        if self.is_cut(Q=Q, pC=pC, pS=pS):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [None],
            'values': {(0, 0): {'finite': None}}
        })

        # Evaluate kernel
        xs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        x = xs[0]
    
        # See Eq. (4.17) of NNLO compatible NLO scheme publication arXiv:0903.1218v2

        # There is no need for the ratio of color-averaging factor between the real ME
        # initial state flavor and the one of the reduced Born ME as they are either both
        # gluons or both quarks
        
        evaluation['values'][(0, 0)]['finite'] = self.color_charge * ( 2. / (1. - x) )

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
