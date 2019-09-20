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
from bidict import bidict
from pprint import pformat

import madgraph.integrator.mappings as mappings
import madgraph.integrator.vectors as vectors
import madgraph.various.misc as misc
from madgraph.core.subtraction import BeamCurrent, IntegratedBeamCurrent, \
    IntegratedCurrent, Counterterm, SubtractionLeg
from madgraph.core.base_objects import EpsilonExpansion

import commons.utils as utils

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError
import madgraph.core.subtraction as sub

#=========================================================================================
# Cuts and factors functions
#=========================================================================================

def no_cut(**opts):

    # The default behavior is to include counterterms everywhere in phase space
    return False

def no_factor(**opts):

    # The default behavior is that basic currents include all necessary factors
    return 1

def alpha_jacobian(**opts):
    """The jacobian of the change of variables
    between virtuality and the alpha parameter of the rescaling collinear mapping.
    """

    Q = opts['Q']
    pC = opts['pC']
    qC = opts['qC']
    return Q.dot(qC)/Q.dot(pC)

#=========================================================================================
# Current variables
#=========================================================================================

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

def anti_final_coll_variables(PS_point, parent_momentum, children, **opts):

    Q = opts['Q']
    Q2 = Q.square()
    Qnorm = Q2 ** 0.5
    pvec = parent_momentum - (Q.dot(parent_momentum)/Q2) * Q
    pnorm = (-pvec.square()) ** 0.5
    n = pvec / pnorm
    t = Q / Qnorm
    na = t + n
    nb = t - n
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

def compute_energy_fractions(momenta,reference,**opts):
    """Compute the energy fractions of a set of momenta with respect to a reference

    Given a set of momenta p1...pk and a reference vector n, the energy fractions are defined as zi = pi.n/(p1.n+...+pk.n)
    :param momenta: momenta whose energy fraction we compute
    :type momenta: list of LorentzVector
    :param reference: reference vector
    :type reference: LorentzVector
    :return: list of energy fractions ordered like the list of momenta
    :rtype: list of float
    """
    energy_fractions = []
    for p in momenta:
        z = p.dot(reference)
        energy_fractions.append(z)
    normalization = sum(energy_fractions)
    return [z/normalization for z in energy_fractions]

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
        # self.SEpsilon   = utils.Constants.SEpsilon
        # The SEpsilon volume factor is factorized from all virtual and integrated contributions
        # so that the poles between the two cancel even before multiplying by SEpsilon, hence
        # making the result equivalent to what one would have obtained by multiplying by SEpsilon=1.
        self.SEpsilon = EpsilonExpansion({0 : 1.})

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

    # TODO This is mantained for backward compatibility, DEPRECATED
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
# QCDLocalCurrent
#=========================================================================================
class QCDLocalCurrent(QCDCurrent):
    """Parent class for QCD local currents."""

    expected_init_opts = ('leg_numbers_map', 'mapping_singular_structure')
    non_local_currents = (BeamCurrent, IntegratedBeamCurrent, IntegratedCurrent)
    resolve_mother_spin_and_color = True

    @classmethod
    def check_current_properties(cls, current):

        try:
            return all([
                current['squared_orders'] == cls.squared_orders,
                current['n_loops'] == cls.n_loops,
                current['resolve_mother_spin_and_color'] == cls.resolve_mother_spin_and_color,
                not isinstance(current, cls.non_local_currents)
            ])
        except KeyError as ke:
            raise CurrentImplementationError(
                "The current " + cls.__name__ + " called check_current_properties " +
                "without setting " + str(ke) + ", please review your implementation")

    @classmethod
    def build_equivalency_sets(cls, model):
        """ By default, the equivalency class is such that light quarks and antiquarks will be
        matched indiscriminatory. In the case of electroweak theory you may
        not want to put the quarks and antiquarks in the same equivalency class. For QCD
        you do."""
        return [set(range(1, model.get_nflav()+1)) | set(range(-model.get_nflav(), 0))]

    @classmethod
    def does_implement_this_current(cls, current, model):

        if not cls.check_current_properties(current):
            return None

#TZdebugger gets past this point

        all_template_structures = cls.structure
        if isinstance(all_template_structures, sub.SingularStructure):
            all_template_structures = [cls.structure, ]

        leg_numbers_map = None
        for template_structure in all_template_structures:
            try:
                leg_numbers_map = template_structure.map_leg_numbers(
                    current.get('singular_structure'), cls.build_equivalency_sets(model))
            except AttributeError:
                continue
            if leg_numbers_map is None:
                continue
            else:
                break

        if leg_numbers_map is None:
#TZ debugger stuck here
#            if len(current["singular_structure"].substructures[0].substructures)>0:
#                if len(current["singular_structure"].substructures[0].substructures[0].legs)==2:
#                    import ipdb
#                    ipdb.set_trace()
            return None

        mapping_singular_structure = current.get('singular_structure').get_copy()
        return {
            'leg_numbers_map': leg_numbers_map,
            'mapping_singular_structure': mapping_singular_structure }

    def __init__(self, *args, **opts):

        for opt_name in self.expected_init_opts:
            try:
                setattr(self, opt_name, opts.pop(opt_name))
            except KeyError:
                raise CurrentImplementationError(
                    "__init__ of " + self.__class__.__name__ + " requires " + opt_name)
        self.supports_helicity_assignment = False
        super(QCDLocalCurrent, self).__init__(*args, **opts)

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
        tuple(sorted([1, -1, 2, -2, 3, -3, 21])),
        tuple(sorted([1, -1, 2, -2, 21])),
        tuple(sorted([1, -1, 21]))
    ]

    def __init__(self, model, **opts):
        """ Initial general attributes common to all beam factorization counterterms."""
        
        super(QCDBeamFactorizationCurrent, self).__init__(model, **opts)
        
        self.beam_type = opts.get('beam_type', 'Unknown')
        self.beam_PDGs = opts.get('beam_PDGs', tuple([]))

        # Retrieve NF from the active beam_PDGs, ommitting all massless gauge vectors
        # and grouping particles and anti particles
        self.active_quarks = sorted(list(set( abs(pdg) for pdg in self.beam_PDGs if 
                                                model.get_particle(pdg).get('spin')==2  )))
        self.NF = len([1 for pdg in self.active_quarks if
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

    def evaluate_subtraction_current(self, current, lower_PS_point=None,
            reduced_process = None, xi=None, mu_r=None, mu_f=None, Q=None, hel_config=None, 
            allowed_backward_evolved_flavors = 'ALL', **opts ):
        """ This implementation of the main function call in the base class pre-process
        the inputs so as to define the variable generically useful for all beam factorization
        current."""

        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if self.distribution_type != 'endpoint' and xi is None:
            raise CurrentImplementationError(
                self.name() + " requires the rescaling variable xi." )
        if lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the lower phase-space point.")
        if mu_f is None:
            raise CurrentImplementationError(
                self.name() + " requires the factorization scale mu_f." )
        if mu_r is None:
            raise CurrentImplementationError(
                self.name() + " requires the factorization scale mu_r." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a process instance to be specified" )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total initial momentum Q." )

        # Retrieve alpha_s
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']

        # Compute the normalization factor
        normalization = self.SEpsilon * ( alpha_s / (2. * math.pi) ) ** (current['n_loops'] + 1)

        # For beam factorization terms, this function returns an instance of
        # BeamFactorizationCurrentEvaluation which can specify color-correlations as 
        # well as reduced and resolved flavors.
        evaluation = self.evaluate_kernel(
            lower_PS_point, reduced_process, xi, mu_r, mu_f, Q, normalization,
            allowed_backward_evolved_flavors = allowed_backward_evolved_flavors)

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

    @staticmethod
    def apply_flavor_mask(flavor_matrix,allowed_backward_evolved_flavors):
        """
        Given a flavor matrix and a list of permitted flavors that can be the end point of the backward evolution (the flavor mask), generate a filtered flavor matrix.

        :param flavor_matrix: sparse matrix implemented as a dict of dict (see madgraph.integrator.ME7_integrands.ME7Event#convolve_flavors)
        :param allowed_backward_evolved_flavors: tuple of PDG ids
        :return: filtered_flavor_matrix, with the same structure as flavor_matrix
        """
        # If no mask do nothing
        if allowed_backward_evolved_flavors == 'ALL':
            return flavor_matrix
        else:
            # We will loop over a matrix M[i][j] = wgt_ij where i is the starting PDG and j a tuple of ending PDGs. We filter the elements in j.
            filtered_flavor_matrix = {}
            # Loop over matrix lines i
            for reduced_flavor in flavor_matrix:
                # Each column is a dict {(PDG1.1, PDG1.2, PDG1.3) : wgt1,...} where the tuple j=(PDG1.1, PDG1.2, PDG1.3) is the label of a column  and wgt the entry ij of the matrix
                for end_flavors, wgt in  flavor_matrix[reduced_flavor].items():
                    #Apply the filter on the elements of the tuple
                    allowed_end_flavors = tuple([fl for fl in end_flavors if fl in allowed_backward_evolved_flavors])
                    #If a column was entirely filtered out, do not include it
                    if allowed_end_flavors:
                        filtered_flavor_matrix[reduced_flavor] = {allowed_end_flavors:wgt}

            return filtered_flavor_matrix

#=========================================================================================
# QCDLocalCollinearCurrent
#=========================================================================================
class QCDLocalCollinearCurrent(QCDLocalCurrent):
    """Common functions for QCD local collinear currents."""

    expected_init_opts = ('leg_numbers_map', 'mapping_singular_structure', 'has_initial_state')

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Add information about whether or not this local collinear current contains initial states."""
        res = super(QCDLocalCollinearCurrent, cls).does_implement_this_current(current, model)

        if res is not None:
            res['has_initial_state'] = current.get('singular_structure').get_all_legs().has_initial_state_leg()

        return res

    def kernel(self, zs, kTs, parent, reduced_kinematics):
        """Evaluate the basic splitting kernel.

        :param zs: momentum fractions
        :param kTs: transverse momenta
        :param parent: parent leg number
        :param reduced_kinematics: reduced kinematics associated to this kernel
            in the form (reduced_kinematics_identifier, lower_PS_point)
        :return: value of the splitting kernel
        :return type: SubtractionCurrentEvaluation
        """

        raise NotImplementedError

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, Q=None, **opts):

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
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires specification of the total incoming momentum Q.")

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Retrieve leg numbers
        children = tuple(self.leg_numbers_map[i]
                         for i in sorted(self.leg_numbers_map.keys()))
        parent = momenta_dict.inv[frozenset(children)]

        # Perform mapping
        self.mapping_singular_structure.legs = self.get_recoilers( reduced_process, excluded=(parent, ) )

        lower_PS_point, mapping_vars = self.mapping.map_to_lower_multiplicity(
            higher_PS_point, self.mapping_singular_structure, momenta_dict,
            compute_jacobian=self.divide_by_jacobian)

        # Retrieve kinematics
        # The Q variable of the mapping cannot be relied upon
        #Q = mapping_vars['Q']
        if self.has_initial_state:
            pC = higher_PS_point[children[0]]-sum(higher_PS_point[child] for child in children[1:])
        else:
            pC = sum(higher_PS_point[child] for child in children)

        qC = lower_PS_point[parent]
        jacobian = mapping_vars.get('jacobian', 1)
        reduced_kinematics = (None, lower_PS_point)

        # Include the counterterm only in a part of the phase space
        if self.has_initial_state:
            pA = higher_PS_point[children[0]]
            pR = sum(higher_PS_point[child] for child in children[1:])
            # Initial state collinear cut
            if self.is_cut(Q=Q, pA=pA, pR=pR):
                return utils.SubtractionCurrentResult.zero(
                    current=current, hel_config=hel_config,
                    reduced_kinematics=('IS_CUT', lower_PS_point))
        else:
            # Final state collinear cut
            if self.is_cut(Q=Q, pC=pC):
                return utils.SubtractionCurrentResult.zero(
                    current=current, hel_config=hel_config,
                    reduced_kinematics=('IS_CUT', lower_PS_point))

        # Evaluate kernel
        # First construct variables necessary for its evaluation
        # pass the mapping variables to the variables function, except for Q which is provided externally
        mapping_vars.pop('Q',None)
        kernel_arguments = self.variables(higher_PS_point, qC, children, Q=Q, **mapping_vars)

        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [ None, ],
            'color_correlations': [ None, ],
            'reduced_kinematics': [ reduced_kinematics ],
            'values': { }
        })
        evaluation = self.kernel(evaluation, parent, *kernel_arguments)

        # Add the normalization factors
        pC2 = pC.square()
        norm = (8. * math.pi * alpha_s / pC2) ** (len(children) - 1)
        norm *= self.factor(Q=Q, pC=pC, qC=qC)
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

#=========================================================================================
#
# Nested currents
#
# These are more advanced currents that require dedicated implementations of the cuts,
# variables and class attributes
#
#=========================================================================================

class CompoundVariables(object):

    def __init__(self, *variables_generators):
        self.variables_generators = variables_generators

    def __call__(self, higher_PS_point, lower_PS_point, bundles_info, **opts):
        """ Calls in sequence the variables for each bundles, and return the aggregated list."""

        assert(len(bundles_info)==len(self.variables_generators))

        all_variables = []
        for bundle_info, variables_generator in zip(bundles_info,self.variables_generators):
            if variables_generator is not None:
                all_variables.append(variables_generator(
                    higher_PS_point,
                    lower_PS_point[bundle_info['parent']],
                    tuple(list(bundle_info['initial_state_children'])+list(bundle_info['final_state_children'])),
                    **opts)[0]
                )
            else:
                all_variables.append({})
        return all_variables

class GeneralQCDLocalCurrent(QCDLocalCurrent):
    """Common functions for QCD local collinear currents."""

    expected_init_opts = ('leg_numbers_map', 'mapping_singular_structure', 'has_initial_state')

    # This should be defined by daughter class, notice that if a list of structures are specified
    # they can only differ according to their flavors.
    structure = NotImplemented

    # The class attribute mapping rules specifies each level of mapping in an ordered list, each
    # composed of a dictionary with the necessary component for performing this mapping:
    mapping_rules = NotImplemented

    # Decide whether one must divide by the jacobian or not. By default we don't
    divide_by_jacobian = False

    # The attirbutes below are not used in this implementation
    get_recoilers   = None
    factor          = None
    mapping         = None
    is_cut          = None

    # A global variables builder can also be supplied and will be called to build the overall variables
    # that need to know about the entire mapping history. It can be left to None if not necessary.
    variables = None

    # Example: for the nested collinear FF within IF
    # mapping_rules = [
    #     {
    #         'singular_structure'    : sub.SingularStructure(substructures=(sub.CollStructure(
    #             substructures=tuple([]),
    #             legs=(
    #                 sub.SubtractionLeg(10, +2, sub.SubtractionLeg.FINAL),
    #                 sub.SubtractionLeg(11, -2, sub.SubtractionLeg.FINAL),
    #             )
    #         ),)),
    #         'mapping'               : colorful_pp_config.final_coll_mapping,
    #         # Intermediate legs should be strictly superior to a 1000
    #         'momenta_dict'          : bidict({1001:frozenset((10,11))}),
    #         'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_FFn_variables),
    #         'is_cut'                : colorful_pp_config.generalised_cuts,
    #         'reduced_recoilers'     : colorful_pp_config.get_initial_state_recoilers,
    #         'additional_recoilers'  : sub.SubtractionLegSet([sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL)]),
    #     },
    #     {
    #         'singular_structure': sub.SingularStructure(substructures=(sub.CollStructure(
    #             substructures=tuple([]),
    #             legs=(
    #                 sub.SubtractionLeg(1, +1, sub.SubtractionLeg.INITIAL),
    #                 sub.SubtractionLeg(1001, +1, sub.SubtractionLeg.FINAL),
    #             )
    #         ),)),
    #         'mapping'               : colorful_pp_config.initial_coll_mapping,
    #         # -1 indicates that this ID should be replaced by the first overall parent connecting to the ME
    #         'momenta_dict'          : bidict({-1: frozenset((1001, 1))}),
    #         'variables'             : currents.CompoundVariables(kernel_variables.colorful_pp_IFn_variables),
    #         'is_cut'                : colorful_pp_config.generalised_cuts,
    #         'reduced_recoilers'     : colorful_pp_config.get_final_state_recoilers,
    #         'additional_recoilers'  : sub.SubtractionLegSet([]),
    #     },
    # ]

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Add information about whether or not this local collinear current contains initial states."""

        res = super(GeneralQCDLocalCurrent, cls).does_implement_this_current(current, model)
#TZdebugger, checked
#        if len(current["singular_structure"].substructures[0].substructures)>0:
#            if len(current["singular_structure"].substructures[0].substructures[0].legs)==2:
#                import ipdb
#                ipdb.set_trace()

        if res is not None:
            res['has_initial_state'] = current.get('singular_structure').get_all_legs().has_initial_state_leg()

        return res

    # Dummy soft-kernel, it can be specialised by the daughter class if it should not be dummy

    #removed "parent" argument that was just after "evaluation"
    def kernel(self, evaluation, steps_and_bundles_variables, global_variables):
        """Evaluate a collinear type of splitting kernel, which does *not* need to know about the reduced process
        Should be specialised by the daughter class if not dummy
        """

        return evaluation

    def call_soft_kernel(self, evaluation, reduced_process, *args):
        """Evaluate a collinear type of splitting kernel, which *does* need to know about the reduced process
        Should be specialised by the daughter class if more info than just the colored partons must be extracted
        by from the reduced process.
        """

        colored_parton_numbers = {}
        for leg in reduced_process.get('legs'):
            leg_color_quantum_number = self.model.get_particle(leg.get('id')).get('color')
            if leg_color_quantum_number==1:
                continue
            colored_parton_numbers[leg.get('number')] = leg_color_quantum_number

        return self.soft_kernel(evaluation,colored_parton_numbers,*args)

    def soft_kernel(self, evaluation, colored_partons, all_steps_info, global_variables):
        """Evaluate a collinear type of splitting kernel, which *does* need to know about the reduced process
        Should be specialised by the daughter class if not dummy
        """

        return evaluation

    def map_leg_number(self, leg_number, parents):
        """
        Map a leg number to its true ID in a specific process.
        Inputs:
            - leg_number: ID of the leg within the singular structure
            - parents: list of process-specific parent momentum in the top-level reduced process. They
        Rules are as follows:
            leg_number < 0 : Should be replaced by the overall parent #abs(leg_number) of the substructure considered.
            leg_number >= 0 : Should be replaced according to self.leg_numbers_map
            leg_number > 1000: Should be left as is given that this is an intermediate leg number
        """
        if leg_number > 1000:
            return leg_number
        elif leg_number < 0:
            return parents[abs(leg_number)-1]
        else:
            return self.leg_numbers_map[leg_number]

    def map_leg_numbers_in_singular_structure(self, singular_structure, parents):
        """ Recursively walk through the singular structure and map the generic subtraction leg numbers to either:
        - the explicit indices of the resovled process
        - the explicit indices of the unresolved process
        - dummy intermediate indices
        """

        new_subtraction_legs = []
        for leg in singular_structure.legs:
            new_subtraction_legs.append(sub.SubtractionLeg(self.map_leg_number(leg.n, parents),  leg.pdg, leg.state))
        singular_structure.legs = sub.SubtractionLegSet(new_subtraction_legs)
        for ss in singular_structure.substructures:
            self.map_leg_numbers_in_singular_structure(ss,parents)

    @classmethod
    def get_parent(cls, children_set, momenta_dict):
        """ retrieves parent index according to the momenta_dict supplied given set of children.
        Returns None if it has no parent. This is very similar to the functiong get_ancestor of the
        Counterm class, except that this function follows relabelling rules if present."""

        if len(children_set)==1:
            if children_set in momenta_dict.inv:
                # Follow relabelling rule(s)
                return cls.get_parent(frozenset([ momenta_dict.inv[children_set], ]),momenta_dict)
            else:
                return list(children_set)[0]

        if children_set in momenta_dict.inv:
            # Follow relabelling rule(s)
            return cls.get_parent(frozenset([ momenta_dict.inv[children_set], ]), momenta_dict)

        new_children_set = frozenset(children_set)
        for key in momenta_dict.inv.keys():
            if len(key)==1 or key.isdisjoint(new_children_set):
                continue
            if key.issubset(new_children_set):
                new_children_set = new_children_set.difference(key)
                new_children_set = new_children_set.union(frozenset([ momenta_dict.inv[key], ]))
        if new_children_set!=children_set:
            return cls.get_parent(new_children_set, momenta_dict)
        else:
            raise CurrentImplementationError("Could not determine the parents of children %s"%str(children_set)+
                                             " from the momenta_dict supplied: %s"%str(momenta_dict))


    @classmethod
    def has_parent(cls, singular_structure_bundle, n_children):
        """ Tests whether the singular structure of a bundle has a parent and raises an error if it
        should have more than one. """

        n_parents = n_children - singular_structure_bundle.count_unresolved()
        if n_parents==0:
            return False
        elif n_parents==1:
            return True
        else:
            raise CurrentImplementationError("A bundle singular structure should have one or zero parents, but"+
                                             " %s has %d."%(str(singular_structure_bundle), n_parents))

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, Q=None, **opts):

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
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires specification of the total incoming momentum Q.")

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Retrieve leg numbers
        overall_children = [] # the legs becoming unresolved in the resolved process (momenta in the real-emission ME)
        overall_parents = [] # the parent legs in the top-level reduced process (momenta in the mapped ME of this CT)
        # We take one defining structure (self.structure[0]). Its top-level substructures are the independent bundles
        # eg (C(1,2),C(4,5),S(6)).
        # Using structure[0] as the defining structure is fine for the purpose below since they should
        # all use the same leg numbers
        for bundle in self.structure[0].substructures:
            overall_children.append(tuple(self.leg_numbers_map[l.n] for l in bundle.get_all_legs()))
            if self.has_parent(bundle, len(overall_children[-1])):
                overall_parents.append(self.get_parent(frozenset(overall_children[-1]), momenta_dict))
            else:
                overall_parents.append(None)

        all_steps = [{'higher_PS_point': higher_PS_point},]
        overall_jacobian = 1.
        for i_step, mapping_information in enumerate(self.mapping_rules):
            # Now obtain recoilers
            reduced_recoilers = mapping_information['reduced_recoilers'](reduced_process, excluded=tuple(overall_parents))
            additional_recoilers = sub.SubtractionLegSet( SubtractionLeg(self.map_leg_number(l.n, overall_parents),l.pdg,l.state)
                                                                    for l in mapping_information['additional_recoilers'] )
            all_recoilers = sub.SubtractionLegSet( list(reduced_recoilers)+list(additional_recoilers))

            # Now recursively apply leg numbers mappings
            mapping_singular_structure = mapping_information['singular_structure'].get_copy()
            self.map_leg_numbers_in_singular_structure(mapping_singular_structure, overall_parents)
            # Now assign the recoilers (whose leg numbers have already been mapped)
            mapping_singular_structure.legs = all_recoilers

            # Build the momenta_dict by also substituting leg numbers
            this_momenta_dict = bidict({self.map_leg_number(k,overall_parents):frozenset([
                self.map_leg_number(n, overall_parents) for n in v]) for k,v in mapping_information['momenta_dict'].items()})

            lower_PS_point, mapping_vars = mapping_information['mapping'].map_to_lower_multiplicity(
                all_steps[-1]['higher_PS_point'], mapping_singular_structure, this_momenta_dict,
                compute_jacobian=self.divide_by_jacobian)
            if mappings.RelabellingMapping.needs_relabelling(this_momenta_dict):
                lower_PS_point = mappings.RelabellingMapping.map_to_lower_multiplicity(lower_PS_point, None, this_momenta_dict)

            all_steps[-1]['lower_PS_point'] = lower_PS_point
            # Q is provided externally
            mapping_vars.pop('Q', None)
            all_steps[-1]['mapping_vars'] = mapping_vars

            overall_jacobian *= mapping_vars.get('jacobian', 1.)

            bundles_info = []
            for bundle in mapping_singular_structure.substructures:
                bundles_info.append({})
                all_legs = bundle.get_all_legs()
                # This sorting is important so that the variables generated can be related to the legs specified
                # in the mapping singular structures of the mapping rules
                all_initial_legs = sorted([l for l in all_legs if l.state==l.INITIAL], key = lambda l: l.n)
                all_final_legs = sorted([l for l in all_legs if l.state==l.FINAL], key = lambda l: l.n)
                bundles_info[-1]['initial_state_children'] = tuple(l.n for l in all_initial_legs)
                bundles_info[-1]['final_state_children'] = tuple(l.n for l in all_final_legs)
                if self.has_parent(bundle, len(all_legs)):
                    bundles_info[-1]['parent'] = self.get_parent(frozenset(l.n for l in all_legs),this_momenta_dict)
                else:
                    bundles_info[-1]['parent'] = None

                # Retrieve kinematics
                bundles_info[-1]['cut_inputs'] = {}
                if bundle.name() == 'C':
                    if len(all_initial_legs)>0:
                        bundles_info[-1]['cut_inputs']['pA'] = -sum(all_steps[-1]['higher_PS_point'][l.n] for l in all_initial_legs)
                    bundles_info[-1]['cut_inputs']['pC'] = sum(all_steps[-1]['higher_PS_point'][l.n] for l in all_final_legs)
                elif bundle.name() == 'S':
                    bundles_info[-1]['cut_inputs']['pS'] = sum(all_steps[-1]['higher_PS_point'][l.n] for l in all_final_legs)

            all_steps[-1]['bundles_info'] = bundles_info

            # Get all variables for this level
            if mapping_information['variables'] is not None:
                all_steps[-1]['variables'] = mapping_information['variables'](
                    all_steps[-1]['higher_PS_point'], all_steps[-1]['lower_PS_point'],
                    bundles_info, Q=Q, **mapping_vars)
            else:
                all_steps[-1]['variables'] = {}

            # Add the next higher PS point for the next level if necessary
            if i_step<(len(self.mapping_rules)-1):
                all_steps.append({'higher_PS_point':lower_PS_point})

        overall_lower_PS_point = all_steps[-1]['lower_PS_point']
        reduced_kinematics = (None, overall_lower_PS_point)

        global_variables = {
                'overall_children': overall_children,
                'overall_parents' : overall_parents,
                'leg_numbers_map' : self.leg_numbers_map,
                'Q' : Q,
        }
        # Build global variables if necessary
        if self.variables is not None:
            global_variables.update(self.variables(all_steps, global_variables))

        # Apply cuts: include the counterterm only in a part of the phase-space
        for i_step, step_info in enumerate(all_steps):
            cut_inputs = dict(step_info)
            if self.mapping_rules[i_step]['is_cut'](cut_inputs, global_variables):
                return utils.SubtractionCurrentResult.zero(
                    current=current, hel_config=hel_config,
                    reduced_kinematics=('IS_CUT', overall_lower_PS_point))

#        for i_step, step_info in enumerate(all_steps):
#            misc.sprint("Higher PS point at step #%d: %s"%(i_step, str(step_info['higher_PS_point'])))
#            misc.sprint("Lower PS point at step  #%d: %s"%(i_step, str(step_info['lower_PS_point'])))

        # Evaluate kernel
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations' : [ None, ],
            'color_correlations': [ None, ],
            'reduced_kinematics': [ reduced_kinematics, ],
            'values': { }
        })

        # Apply collinear kernel (can be dummy)
        evaluation = self.kernel(evaluation, all_steps, global_variables)

        # Apply soft kernel (can be dummy), which also knows about the reduced process
        evaluation = self.call_soft_kernel(evaluation, reduced_process, all_steps, global_variables)

        # Add the normalization factors
        # WARNING! In this implementation the propagator denominators must be included in the kernel evaluation.
        norm = (8. * math.pi * alpha_s) ** (self.squared_orders['QCD']/2)
        norm /= overall_jacobian
        for k in evaluation['values']:
            for term in evaluation['values'][k]:
                evaluation['values'][k][term] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result
