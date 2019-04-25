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


import commons.utils as utils
import commons.QCD_local_currents as currents
import commons.factors_and_cuts as factors_and_cuts

import madgraph.various.misc as misc

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# NLO soft current
#=========================================================================================

class QCD_soft_0_g(currents.QCDLocalSoftCurrent):
    """Soft gluon eikonal current at tree level, eq.4.12-4.13 of arXiv:0903.1218."""

    factor = staticmethod(factors_and_cuts.factor_soft)
    is_cut = staticmethod(factors_and_cuts.cut_soft)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD soft tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve the singular structure
        singular_structure = current.get('singular_structure').substructures[0]
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
        higher_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):
        if higher_PS_point is None:
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
        soft_leg_number = current.get('singular_structure').substructures[0].legs[0].n

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
    factor = staticmethod(factors_and_cuts.factor_soft)
    is_cut = staticmethod(factors_and_cuts.cut_soft)


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
        singular_structure = current.get('singular_structure').substructures[0]
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

        ss = current.get('singular_structure').substructures[0]
        return (ss.substructures[0].legs[0].n, ss.legs[0].n)

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):
        if higher_PS_point is None:
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
        for substructure in current.get('singular_structure').substructures[0].substructures:
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

        ######
        # TODO lower_PS_point no longer available here, mapping must be invoked!
        ######
        # Evaluate kernel
        # zs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        zs, kTs = self.variables(higher_PS_point, vectors.LorentzVector(), children, Q=Q)

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
