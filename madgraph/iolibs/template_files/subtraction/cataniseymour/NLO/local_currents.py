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
# Eikonal factor, modified by partial fractioning and without divergence
#=========================================================================================

def mod_eikonal(PS_point, i, j, r):
    """Modified eikonal factor for soft particle with number 'r'
    emitted from 'i' and reconnecting to 'j'.
    This is obtained starting from the eikonal and:
    - ignoring 1 / sir, which is already included in the normalisation factor;
    - multiplying by the partial fraction sjr / (sir + sjr) to regulate for sjr -> 0.
    """

    sij = 2*PS_point[i].dot(PS_point[j])
    sir = 2*PS_point[i].dot(PS_point[r])
    sjr = 2*PS_point[j].dot(PS_point[r])
    return 2 * sij / (sir + sjr)

#=========================================================================================
# NLO final-collinear currents, containing the soft limits
#=========================================================================================

class QCD_final_collinear_0_qqx(currents.QCDLocalCollinearCurrent):
    """q q~ collinear tree-level current."""

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
        evaluation['values'][(0, 0)]['finite'] = self.CF * ((1.-z)**2 - 1.)/z
        return evaluation

    def evaluate_subtraction_current(
        self, current, PS_point,
        reduced_process=None, hel_config=None,
        mapping_variables=None, leg_numbers_map=None):

        if not hel_config is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s"
                "does not support helicity assignment." % self.__class__.__name__)

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Include the counterterm only in a part of the phase space
        children = self.get_sorted_children(current, self.model)
        parent = leg_numbers_map.inv[frozenset(children)]
        if self.is_cut(mapping_variables, parent):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Evaluate collinear subtracted kernel
        zs, kTs = self.variables(PS_point, parent, children, mapping_variables)
        evaluation = self.evaluate_kernel(zs, kTs, parent)

        # Find all colored leg numbers except for the parent in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color') == 1:
                continue
            number = leg.get('number')
            if number == parent:
                continue
            all_colored_parton_numbers.append(number)

        color_correlation_index = 1
        # Now loop over the colored parton number pairs (parent, a)
        # and add the corresponding contributions to this current
        for a in all_colored_parton_numbers:
            evaluation['color_correlations'].append(((parent, a),))
            # Write the eikonal for that pair
            evaluation['values'][(0, color_correlation_index)] = {
                'finite': -mod_eikonal(PS_point, parent, a, children[0]) }
            color_correlation_index += 1

        # Add the normalization factors
        pC = mapping_variables['pC' + str(parent)]
        pC2 = pC.square()
        norm = 8. * math.pi * alpha_s / pC2
        norm *= self.factor(mapping_variables, parent)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

class QCD_final_collinear_0_gg(currents.QCDLocalCollinearCurrent):
    """g g collinear tree-level current."""

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
        # full_00 = (z/(1.-z)) + ((1.-z)/z)
        # limit_00 = 1./z + 1./(1.-z)
        # evaluation['values'][(0, 0)]['finite'] =  2.*self.CA * (full_00-limit_00)
        evaluation['values'][(0, 0)]['finite'] = -4.*self.CA
        evaluation['values'][(1, 0)]['finite'] = -2.*self.CA * 2.*z*(1.-z) / kT.square()
        return evaluation

    def evaluate_subtraction_current(
        self, current, PS_point,
        reduced_process=None, hel_config=None,
        mapping_variables=None, leg_numbers_map=None):

        if not hel_config is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s"
                "does not support helicity assignment." % self.__class__.__name__)

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Include the counterterm only in a part of the phase space
        children = self.get_sorted_children(current, self.model)
        parent = leg_numbers_map.inv[frozenset(children)]
        if self.is_cut(mapping_variables, parent):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Evaluate collinear subtracted kernel
        zs, kTs = self.variables(PS_point, parent, children, mapping_variables)
        evaluation = self.evaluate_kernel(zs, kTs, parent)

        # Find all colored leg numbers except for the parent in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color') == 1:
                continue
            number = leg.get('number')
            if number == parent:
                continue
            all_colored_parton_numbers.append(number)

        color_correlation_index = 1
        # Loop over the colored parton number pairs (parent, a)
        # and add the corresponding contributions to this current
        for a in all_colored_parton_numbers:
            evaluation['color_correlations'].append(((parent, a),))
            # Write the eikonal for that pair
            eik0 = -mod_eikonal(PS_point, parent, a, children[0])
            eik1 = -mod_eikonal(PS_point, parent, a, children[1])
            evaluation['values'][(0, color_correlation_index)] = {'finite': eik0 + eik1}
            color_correlation_index += 1

        # Add the normalization factors
        pC = mapping_variables['pC' + str(parent)]
        pC2 = pC.square()
        norm = 8. * math.pi * alpha_s / pC2
        norm *= self.factor(mapping_variables, parent)
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
# NLO soft current
#=========================================================================================

class NoSoftCurrent(currents.QCDCurrent):
    """Trivial current returning zero for any NLO limit containing softs."""

    is_zero = True

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD soft tree-level currents
        init_vars = currents.QCDLocalSoftCurrent.\
            common_does_implement_this_current(current, 2, 0)
        if init_vars is not None:
            return init_vars
        init_vars = currents.QCDLocalSoftCollinearCurrent.\
            common_does_implement_this_current(current, 2, 0)
        return init_vars

    def evaluate_subtraction_current(
        self, current, PS_point,
        reduced_process=None, hel_config=None,
        mapping_variables=None, leg_numbers_map=None):
        """Return 0 for this current."""

        return utils.SubtractionCurrentResult.zero(current=current, hel_config=hel_config)
