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
from madgraph.core.subtraction import BeamCurrent, IntegratedBeamCurrent, Counterterm, SubtractionLeg
from madgraph.core.base_objects import EpsilonExpansion

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

def no_cut(**opts):

    # The default behavior is to include counterterms everywhere in phase space
    return False

def no_factor(**opts):

    # The default behavior is that basic currents include all necessary factors
    return 1

def n_final_coll_variables(PS_point, parent_momentum, children, **opts):

    na, nb = mappings.FinalCollinearVariables.collinear_and_reference(parent_momentum)
    kin_variables = dict()
    mappings.FinalCollinearVariables.get(
        PS_point, children, na, nb, kin_variables)
    zs  = tuple(kin_variables['z%d'  % i] for i in children)
    kTs = tuple(kin_variables['kt%d' % i] for i in children)
    return zs, kTs

def Q_final_coll_variables(PS_point, parent_momentum, children, **opts):

    na = parent_momentum
    nb = opts['Q']
    kin_variables = dict()
    mappings.FinalCollinearVariables.get(
        PS_point, children, na, nb, kin_variables)
    zs  = tuple(kin_variables['z%d'  % i] for i in children)
    kTs = tuple(kin_variables['kt%d' % i] for i in children)
    return zs, kTs

def n_initial_coll_variables(PS_point, parent_momentum, children, **opts):

    na, nb = mappings.InitialCollinearVariables.collinear_and_reference(parent_momentum)
    kin_variables = dict()
    # The lone initial state child is always placed first thanks to the implementation
    # of the function get_sorted_children() in the current.
    mappings.InitialCollinearVariables.get(
        PS_point, children[1:], children[0], na, nb, kin_variables)
    zs  = tuple(kin_variables['z%d'  % i] for i in children)
    kTs = tuple(kin_variables['kt%d' % i] for i in children)
    return zs, kTs

def Q_initial_coll_variables(PS_point, parent_momentum, children, **opts):
    
    na = parent_momentum
    nb = opts['Q']
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

        self.TR = model_param_dict.get('TR', utils.Constants.TR)
        self.NC = model_param_dict.get('NC', utils.Constants.NC)
        self.CF = model_param_dict.get('CF', (self.NC ** 2 - 1) / (2 * self.NC))
        self.CA = model_param_dict.get('CA', self.NC)
        
        self.EulerGamma = utils.Constants.EulerGamma
        # S_\eps = (4 \[Pi])^\[Epsilon] E^(-\[Epsilon] EulerGamma)
        self.SEpsilon   = utils.Constants.SEpsilon
        # The SEpsilon volume factor is factorized from all virtual and integrated contributions
        # so that the poles between the two cancel even before multiplying by SEpsilon, hence
        # making the result equivalent to what one would have obtained by multiplying by SEpsilon=1.
        self.SEpsilon = EpsilonExpansion({ 0 : 1.})

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
# QCDBeamFactorizationCurrent
#=========================================================================================
class QCDBeamFactorizationCurrent(QCDCurrent):
    """ Mother class for QCD beam factorization currents (integrated ISR and PDF counter-terms
    characterized by the fact that they all have a xi-dependence (or are the end-point of
    a xi-dependent function."""

    # When chosing a pattern where a single class accepts all distribution type, then
    # keep the variable below set to ['bulk','counterterm','endpoint'], otherwise define
    # a subset.
    distribution_types_implemented_in_this_class = ['bulk','counterterm','endpoint']

    # Similarly for the beam_type
    beam_types_implemented_in_this_class = ['proton']
    
    # And for the beam_PDGs
    beam_PDGs_implemented_in_this_class = [
        tuple(sorted([1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,21])),
        tuple(sorted([1,-1,2,-2,3,-3,4,-4,5,-5,21])),
        tuple(sorted([1,-1,2,-2,3,-3,4,-4,21])),
        tuple(sorted([1,-1,2,-2,3,-3,21])),
    ]

    def __init__(self, model, **opts):
        """ Initial general attributes common to all beam factorization counterterms."""
        
        super(QCDBeamFactorizationCurrent, self).__init__(model, **opts)
        
        self.beam_type = opts.get('beam_type', 'Unknown')
        self.beam_PDGs = opts.get('beam_PDGs', tuple([]))

        self.NF = len([1 for pdg in self.beam_PDGs if
                                      model.get_particle(pdg).get('mass').upper()=='ZERO'])

        # This entry *must* be specified.
        self.distribution_type = opts['distribution_type']

    # Prefix this base function with 'common' to screen it from the lookup
    # performed by the MG5aMC current exporter. Implement general checks common to all
    # derived beam factorization currents.
    @classmethod
    def common_does_implement_this_current(
        cls, current, QCD_squared_order=None, n_loops=None):
        """General checks common to all QCD currents."""

        # Check the general properties common to QCD currents
        init_vars = super(
            QCDBeamFactorizationCurrent, cls).common_does_implement_this_current(
            current, QCD_squared_order, n_loops)
        if init_vars is None:
            return None

        if not isinstance(current, (BeamCurrent, IntegratedBeamCurrent)):
            return None

        init_vars = {}
        
        # All checks passed
        if isinstance(current, IntegratedBeamCurrent):
            init_vars['distribution_type'] = 'endpoint'
        elif isinstance(current, BeamCurrent):
            # The two possible types in this case are 'bulk' and 'counterterm'
            init_vars['distribution_type'] = current['distribution_type']
        else:
            # Then this is not an ISR
            return None
        
        if not init_vars['distribution_type'] in cls.distribution_types_implemented_in_this_class:
            return None

        if cls.beam_types_implemented_in_this_class!='ALL':
            if current['beam_type'] in cls.beam_types_implemented_in_this_class:
                init_vars['beam_type'] = current['beam_type']
            else:
                return None
        
        if cls.beam_PDGs_implemented_in_this_class!='ALL':
            if tuple(sorted(current['beam_PDGs'])) in cls.beam_PDGs_implemented_in_this_class:
                init_vars['beam_PDGs'] = tuple(sorted(current['beam_PDGs']))
            else:
                return None

        return init_vars

    def evaluate_subtraction_current(self, current, higher_PS_point=None, lower_PS_point=None, 
            reduced_process = None, xi=None, mu_r=None, mu_f=None, hel_config=None, **opts ):
        """ This implementation of the main function call in the base class preprocess
        the inputs so as to define the variable generically useful for all beam factorization
        current."""

        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if self.distribution_type != 'endpoint' and xi is None:
            raise CurrentImplementationError(
                self.name() + " requires the rescaling variable xi." )
        if mu_f is None:
            raise CurrentImplementationError(
                self.name() + " requires the factorization scale mu_f." )
        if mu_r is None:
            raise CurrentImplementationError(
                self.name() + " requires the factorization scale mu_r." )
        if lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " requires a lower PS point to be specified" )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a process instance to be specified" )

        # Retrieve alpha_s
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']

        # Compute the normalization factor
        normalization = ( alpha_s / (2. * math.pi) ) ** (current['n_loops'] + 1)

        # For beam factorization terms, this function returns an instance of
        # BeamFactorizationCurrentEvaluation which can specify color-correlations as 
        # well as reduced and resolved flavors.
        evaluation = self.evaluate_kernel(
                          lower_PS_point, reduced_process, xi, mu_r, mu_f, normalization)

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

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
        
        # Make sure this is not a beam factorization current
        if isinstance(current, (BeamCurrent, IntegratedBeamCurrent)):
            return None
        
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
        parent = leg_numbers_map.inv[frozenset(children)]
        pC = sum(higher_PS_point[child] for child in children)
        if self.is_cut(Q=Q, pC=pC):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Evaluate kernel
        zs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        evaluation = self.evaluate_kernel(zs, kTs, parent)

        # Add the normalization factors
        pC2 = pC.square()
        norm = (8. * math.pi * alpha_s / pC2) ** (len(children) - 1)
        norm *= self.factor(Q=Q, pC=pC)
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

        # Make sure this is not a beam factorization current
        if isinstance(current, (BeamCurrent, IntegratedBeamCurrent)):
            return None

        # Make sure we don't need to sum over the quantum number of the mother leg
        if not current.get('resolve_mother_spin_and_color'):
            return None
        
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
        
        # Make sure we don't need to sum over the quantum number of the mother leg
        if not current.get('resolve_mother_spin_and_color'):
            return None
        
        # Make sure this is not a beam factorization current
        if isinstance(current, (BeamCurrent, IntegratedBeamCurrent)):
            return None
        
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
    d_0 = 1
    d_0_prime = 2

    @staticmethod
    def cut_coll(**opts):

        try:
            alpha = opts['alpha']
        except KeyError:
            pC    = opts['pC']
            Q     = opts['Q']
            alpha = mappings.FinalRescalingOneMapping.alpha(pC, Q)
        # Include the counterterm only up to alpha_0
        return alpha > SomogyiChoices.alpha_0

    @staticmethod
    def cut_soft(**opts):

        try:
            y  = opts['y']
        except KeyError:
            pS = opts['pS']
            Q  = opts['Q']
            y = mappings.SoftVsFinalMapping.y(pS, Q)
        # Include the counterterm only up to y_0
        return y > SomogyiChoices.y_0

    @staticmethod
    def factor_coll(**opts):

        try:
            alpha = opts['alpha']
        except KeyError:
            pC    = opts['pC']
            Q     = opts['Q']
            alpha = mappings.FinalRescalingOneMapping.alpha(pC, Q)
        norm = (1 - alpha) ** (2 * (SomogyiChoices.d_0 - 1))
        return norm

    @staticmethod
    def factor_soft(**opts):

        try:
            y  = opts['y']
        except KeyError:
            pS = opts['pS']
            Q  = opts['Q']
            y = mappings.SoftVsFinalMapping.y(pS, Q)
        norm = (1 - y) ** (SomogyiChoices.d_0_prime - 2)
        return norm
