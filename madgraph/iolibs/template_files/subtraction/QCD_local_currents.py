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
"""Common functions for QCD local currents."""

import os
import math

import madgraph.integrator.mappings as mappings
import madgraph.various.misc as misc

try:
    # First try to import this in the context of the exported currents
    import SubtractionCurrents.subtraction_current_implementations_utils as utils
except ImportError:
    # If not working, then it must be within MG5_aMC context:
    import madgraph.iolibs.template_files.\
                   subtraction.subtraction_current_implementations_utils as utils

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# Cuts and factors functions
#=========================================================================================

def no_cut(mapping_variables, parent):

    # The default behavior is to include counterterms everywhere in phase space
    return False

def no_factor(mapping_variables, parent):

    # The default behavior is that basic currents include all necessary factors
    return 1

def n_final_coll_variables(PS_point, parent, children, mapping_variables):

    na, nb = mappings.FinalCollinearVariables.collinear_and_reference(PS_point[parent])
    kin_variables = dict()
    mappings.FinalCollinearVariables.get(
        PS_point, children, na, nb, kin_variables)
    zs  = tuple(kin_variables['z%d'  % i] for i in children)
    kTs = tuple(kin_variables['kt%d' % i] for i in children)
    return zs, kTs

def Q_final_coll_variables(PS_point, parent, children, mapping_variables):

    na = PS_point[parent]
    nb = mapping_variables['Q']
    kin_variables = dict()
    mappings.FinalCollinearVariables.get(
        PS_point, children, na, nb, kin_variables)
    zs  = tuple(kin_variables['z%d'  % i] for i in children)
    kTs = tuple(kin_variables['kt%d' % i] for i in children)
    return zs, kTs

def n_initial_coll_variables(PS_point, parent, children, mapping_variables):

    na, nb = mappings.InitialCollinearVariables.collinear_and_reference(PS_point[parent])
    kin_variables = dict()
    # The lone initial state child is always placed first thanks to the implementation
    # of the function get_sorted_children() in the current.
    mappings.FinalCollinearVariables.get(
        PS_point, children[1:], children[0], na, nb, kin_variables)
    zs  = tuple(kin_variables['z%d'  % i] for i in children)
    kTs = tuple(kin_variables['kt%d' % i] for i in children)
    return zs, kTs

def Q_initial_coll_variables(PS_point, parent, children, mapping_variables):
    
    na = PS_point[parent]
    nb = mapping_variables['Q']
    kin_variables = dict()
    # The lone initial state child is always placed first thanks to the implementation
    # of the function get_sorted_children() in the current.
    mappings.InitialCollinearVariables.get(
        PS_point, children[1:], children[0], na, nb, kin_variables)
    zs  = tuple(kin_variables['z%d'  % i] for i in children)
    kTs = tuple(kin_variables['kt%d' % i] for i in children)
    return zs, kTs

#=========================================================================================
# QCDCurrent
#=========================================================================================

class QCDCurrent(utils.VirtualCurrentImplementation):
    """Common functions for QCD currents."""

    def __init__(self, model, **opts):

        super(QCDCurrent, self).__init__(model, **opts)
        # Extract constants from the UFO model if present, otherwise take default values
        try:
            model_param_dict = self.model.get('parameter_dict')
        except:
            model_param_dict = dict()
        self.TR = model_param_dict.get('TR', 0.5)
        self.NC = model_param_dict.get('NC', 3.0)
        self.CF = model_param_dict.get('CF', (self.NC ** 2 - 1) / (2 * self.NC))
        self.CA = model_param_dict.get('CA', self.NC)

    @staticmethod
    def is_quark(leg, model):
        # This will return True both for a quark and anti-quark because it is not the function
        # get_color() which is called but get('color'), and the sign of the dictionary
        # value is always positive.
        return model.get_particle(leg.pdg).get('color') == 3

    @staticmethod
    def is_gluon(leg, model):

        return model.get_particle(leg.pdg).get('color') == 8

    @staticmethod
    def is_massless(leg, model):

        return model.get_particle(leg.pdg).get('mass').upper() == 'ZERO'

    @staticmethod
    def are_antiparticles(leg1, leg2):
        
        # Notice that for this function and two below, one should in principle use the
        # safer functions 'get_pdg_code()' and 'get_anti_pdg_code()' of the particle objects
        # retreived with: model.get_particle(leg.pdg)
        # But that should be irrelevant.
        return leg1.pdg == -leg2.pdg

    @staticmethod
    def are_same_particle(leg1, leg2):

        return leg1.pdg == leg2.pdg

    @staticmethod
    def both_particle_or_antiparticle(leg1, leg2):

        return leg1.pdg * leg2.pdg > 0

    @staticmethod
    def is_initial(leg):

        return leg.state == leg.INITIAL

    is_cut = staticmethod(no_cut)
    factor = staticmethod(no_factor)

    # Prefix this base function with 'common' to screen it from the lookup
    # performed by the MG5aMC current exporter.
    @classmethod
    def common_does_implement_this_current(
        cls, current, QCD_squared_order=None, n_loops=None):
        """General checks common to all QCD currents."""
        
        # Check the current is pure QCD
        squared_orders = current.get('squared_orders')
        for order in squared_orders:
            if order != 'QCD' and squared_orders[order] != 0:
                return None
        # Check the coupling order
        if not QCD_squared_order is None:
            if squared_orders['QCD'] != QCD_squared_order:
                return None
        # Check the number of loops
        if not n_loops is None:
            if current.get('n_loops') != n_loops:
                return None
        # Make sure we don't need to sum over the quantum number of the mother leg
        if not current.get('resolve_mother_spin_and_color'):
            return None
        # All checks passed
        return {}

#=========================================================================================
# QCDLocalCollinearCurrent
#=========================================================================================

class QCDLocalCollinearCurrent(QCDCurrent):
    """Common functions for QCD local collinear currents."""

    variables = staticmethod(Q_final_coll_variables)

    def __init__(self, *args, **opts):

        super(QCDLocalCollinearCurrent, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    # Prefix this base function with 'common' to screen it from the lookup
    # performed by the MG5aMC current exporter.
    @classmethod
    def common_does_implement_this_current(
        cls, current, QCD_squared_order=None, n_loops=None):
        """General checks common to all QCD collinear currents."""

        # Check the general properties common to QCD currents
        init_vars = super(
            QCDLocalCollinearCurrent, cls).common_does_implement_this_current(
            current, QCD_squared_order, n_loops)
        if init_vars is None: return None
        # Check the structure is a simple collinear
        singular_structure = current.get('singular_structure')
        if singular_structure.name() != 'C': return None
        if singular_structure.substructures: return None
        # All checks passed
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):
        """Return a tuple with the leg numbers of canonically sorted children."""

        raise NotImplemented

    def evaluate_kernel(self, zs, kTs, parent):
        """Evaluate the basic splitting kernel, return a SubtractionCurrentEvaluation."""

        raise NotImplemented

    def evaluate_subtraction_current(
        self, current, PS_point,
        reduced_process=None, hel_config=None,
        mapping_variables=None, leg_numbers_map=None ):

        if not hel_config is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s"
                "does not support helicity assignment." % self.__class__.__name__ )

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

        # Evaluate kernel
        zs, kTs = self.variables(PS_point, parent, children, mapping_variables)
        evaluation = self.evaluate_kernel(zs, kTs, parent)

        # Add the normalization factors
        pC = mapping_variables['pC' + str(parent)]
        pC2 = pC.square()
        norm = (8. * math.pi * alpha_s / pC2) ** (len(children) - 1)
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
# QCDLocalSoftCurrent
#=========================================================================================

class QCDLocalSoftCurrent(QCDCurrent):
    """Common functions for QCD local soft currents."""

    def __init__(self, *args, **opts):

        super(QCDLocalSoftCurrent, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    # Prefix this base function with 'common' to screen it from the lookup
    # performed by the MG5aMC current exporter.
    @classmethod
    def common_does_implement_this_current(
        cls, current, QCD_squared_order=None, n_loops=None):
        """General checks common to all QCD soft currents."""

        # Check the general properties common to QCD currents
        init_vars = super(QCDLocalSoftCurrent, cls).common_does_implement_this_current(
            current, QCD_squared_order, n_loops)
        if init_vars is None: return None
        # Check the structure is a simple soft
        singular_structure = current.get('singular_structure')
        if singular_structure.name() != 'S': return None
        if singular_structure.substructures: return None
        # All checks passed
        return init_vars

#=========================================================================================
# QCDLocalSoftCollinearCurrent
#=========================================================================================

class QCDLocalSoftCollinearCurrent(QCDCurrent):
    """Common functions for QCD local soft-collinear currents."""

    def __init__(self, *args, **opts):

        super(QCDLocalSoftCollinearCurrent, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    # Prefix this base function with 'common' to screen it from the lookup
    # performed by the MG5aMC current exporter.
    @classmethod
    def common_does_implement_this_current(
        cls, current, QCD_squared_order=None, n_loops=None):
        """General checks common to all QCD soft-collinear currents."""

        # Check the general properties common to QCD currents
        init_vars = super(
            QCDLocalSoftCollinearCurrent, cls).common_does_implement_this_current(
            current, QCD_squared_order, n_loops)
        if init_vars is None: return None
        # Retrieve the singular structure
        singular_structure = current.get('singular_structure')
        # The main structure should be collinear
        if singular_structure.name() != 'C': return None
        if not singular_structure.substructures: return None
        # Substructures should be simple soft structures
        for sub_singular_structure in singular_structure.substructures:
            if sub_singular_structure.name() != 'S': return None
            if sub_singular_structure.substructures: return None
        # All checks passed
        return init_vars

#=========================================================================================
# Original Somogyi choices
#=========================================================================================

class SomogyiChoices(object):
    """Original Somogyi choices."""

    alpha_0 = 0.5
    y_0 = 0.5
    divide_by_jacobian = True
    d_0 = 1
    d_0_prime = 2

    @staticmethod
    def cut_coll(mapping_variables, parent):

        try:
            alpha = mapping_variables['alpha' + str(parent)]
        except:
            pC    = mapping_variables['pC'    + str(parent)]
            Q     = mapping_variables['Q']
            alpha = mappings.FinalRescalingOneMapping.alpha(pC, Q)
        # Include the counterterm only up to alpha_0
        return alpha > SomogyiChoices.alpha_0

    @staticmethod
    def cut_soft(mapping_variables, parent):

        y = mapping_variables['y']
        # Include the counterterm only up to y_0
        return y > SomogyiChoices.y_0

    @staticmethod
    def factor_coll(mapping_variables, parent):

        try:
            alpha = mapping_variables['alpha' + str(parent)]
        except:
            pC    = mapping_variables['pC'    + str(parent)]
            Q     = mapping_variables['Q']
            alpha = mappings.FinalRescalingOneMapping.alpha(pC, Q)
        norm = (1 - alpha) ** (2 * (SomogyiChoices.d_0 - 1))
        # Jacobian power is a HACK,
        # because the jacobian is now considered part of the current
        # which leads to double counting for disjoint currents
        # with the same mapping.
        # Consider moving it to the mapping or the hike.
        jacobian_power = 1./mapping_variables.get('pow', 1)
        if SomogyiChoices.divide_by_jacobian:
            norm /= mapping_variables['jacobian'] ** jacobian_power
        return norm

    @staticmethod
    def factor_soft(mapping_variables, parent):

        y = mapping_variables['y']
        norm = (1 - y) ** (SomogyiChoices.d_0_prime - 2)
        if SomogyiChoices.divide_by_jacobian:
            norm /= mapping_variables['jacobian']
        return norm
