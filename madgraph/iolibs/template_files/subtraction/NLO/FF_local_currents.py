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
"""Implementation of NLO type of currents."""

import os
import math

import madgraph.core.subtraction as subtraction
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
# Common functions for FF NLO QCD currents
#=========================================================================================

class NLO_FF_QCD_local_subtraction_current(utils.VirtualCurrentImplementation):
    """Template class for all Final-Final NLO QCD local subtraction current."""
    
    # Prefix this base function with 'common' to screen it from the lookup performed
    # by the MG5aMC current exporter.
    @classmethod
    def common_does_implement_this_current(cls, current, model):
        """General class of checks common to all currents inheriting from this class."""
        
        # Make sure it is a local subtraction counterterm and not an integrated one.
        if isinstance(current, subtraction.IntegratedCurrent):
            return None

        squared_orders = current.get('squared_orders')

        # First check that we indeed have a pure NLO QCD current
        if squared_orders['QCD'] != 2:
            return None
        for order in squared_orders:
            if order != 'QCD' and squared_orders[order] != 0:
                return None

        # Now check that it is tree-level
        if current.get('n_loops') != 0:
            return None

        # Make sure we don't need to sum over the quantum number of the mother leg
        if not current.get('resolve_mother_spin_and_color'):
            return None
        
        # Finally check that the singular structure and PDG matches
        singular_structure = current.get('singular_structure')
        
        # Check that all legs are final states
        for leg in singular_structure.get_all_legs():
            if leg.state != leg.FINAL:
                return None
        
        return {}

#=========================================================================================
# Common functions for NLO final-collinear currents
#=========================================================================================

class NLO_FF_QCD_collinear(NLO_FF_QCD_local_subtraction_current):
    """Common functions for NLO final-collinear currents."""

    # The currents will include a Heaviside theta function of argument (alpha_0 - alpha)
    # It _should_ be possible to override the value in daughter currents
    alpha_0 = 0.1
    divide_by_jacobian = False
    d_0 = 0

    def __init__(self, *args, **opts):

        super(NLO_FF_QCD_collinear, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    @classmethod
    def common_does_implement_this_current(cls, current, model):
        """Return None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning
        a dictionary, it specifies potential options that must passed upon instantiating
        this implementation for the current given in argument. """

        singular_structure = current.get('singular_structure')

        # Check the general properties common to FF NLO QCD
        if super(NLO_FF_QCD_collinear, cls).common_does_implement_this_current(
            current, model) is None:
            return None
        # It should be a collinear type of structure
        if singular_structure.name() != 'C':
            return None
        # It should not have nested structures
        if singular_structure.substructures:
            return None
        # It should consist in exactly two legs going collinear
        if len(singular_structure.legs) != 2:
            return None
        # We now know that this current is implemented here, so we return
        # an empty dictionary which could potentially have contained specific
        # options to specify upon instantiating this class.
        return {}

    def get_cache_and_result_key(self, current,
                                 PS_point,
                                 reduced_process=None,
                                 leg_numbers_map=None,
                                 hel_config=None,
                                 mapping_variables={},
                                 **opts
                                 ):
        """If this subtraction current depends on more than just the PS point and alpha_s,
        complement the cache key here with the additional necessary information.
        """

        cache_key, result_key = super(NLO_FF_QCD_collinear,
                                      self).get_cache_and_result_key(
            current, PS_point,
            reduced_process=reduced_process, leg_numbers_map=leg_numbers_map,
            hel_config=hel_config, mapping_variables=mapping_variables)

        # cache_key['another_call_specifier'] = 'set_it_here'

        return cache_key, result_key

    def evaluate_kernel(self, parent, zs, kTs):

        raise NotImplemented

    def evaluate_subtraction_current(
        self, current, PS_point,
        reduced_process=None, leg_numbers_map=None,
        hel_config=None, mapping_variables=None ):
        """Evaluate the current and return the corresponding instance of
        SubtractionCurrentResult. See documentation of the mother function for more details.
        """

        if not hel_config is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s"
                "does not support helicity assignment." % self.__class__.__name__ )

        if leg_numbers_map is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s"
                "requires the leg_number_map." % self.__class__.__name__ )

        ss = current.get('singular_structure')

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        parent, children, _ = mappings.get_structure_numbers(ss, leg_numbers_map)
        # Always sort the children in descending PDG order to identify variables
        children = tuple(leg.n
                         for leg in sorted(ss.legs, key=lambda x: x.pdg, reverse=True)
                         if leg.n in children)

        # Read mapping variables
        try:
            pC       = mapping_variables['pC' + str(parent)]
            Q        = mapping_variables['Q']
            alpha    = mapping_variables['alpha' + str(parent)]
            jacobian = mapping_variables['jacobian']
        except:
            raise CurrentImplementationError(
                "Subtraction current %s" % self.__class__.__name__ +
                "requires mapping variables pC, Q, alpha and jacobian" )

        # Include the counterterm only up to alpha_0
        if alpha > self.alpha_0:
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Compute kinematic variables
        na = PS_point[parent]
        nb = Q
        kin_variables = dict()
        mappings.FinalCollinearVariables.get(PS_point, children, na, nb, kin_variables)

        zs = tuple(kin_variables['z%d' % i] for i in children)
        kTs = tuple(kin_variables['kt%d' % i] for i in children)
        pC2 = pC.square()

        evaluation = self.evaluate_kernel(parent, zs, kTs)

        # Now add the normalization factors
        norm = 4. * math.pi * alpha_s * (2. / pC2)
        if self.divide_by_jacobian: norm /= jacobian
        norm *= (1-alpha)**(2*(self.d_0-1))
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

#=========================================================================================
# NLO final-collinear currents
#=========================================================================================

class NLO_FF_QCD_collinear_qqx(NLO_FF_QCD_collinear):
    """q q~ collinear tree-level current."""

    # Override alpha_0, divide_by_jacobian and d_0 here if needed

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to FF NLO QCD collinear
        if cls.common_does_implement_this_current(current, model) is None:
            return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that all legs are massless quarks
        for leg in ss.legs:
            if not (cls.is_quark(leg, model) and cls.is_massless(leg, model)):
                return None
        # The current is valid and needs no options
        return {}

    def evaluate_kernel(self, parent, zs, kTs):

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


class NLO_FF_QCD_collinear_gq(NLO_FF_QCD_collinear):
    """g q collinear tree-level current."""

    # Override alpha_0, divide_by_jacobian and d_0 here if needed

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to FF NLO QCD collinear
        if cls.common_does_implement_this_current(current, model) is None:
            return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that there are a massless gluon and a massless quark
        if not(
            (cls.is_quark(ss.legs[0], model) and cls.is_gluon(ss.legs[1], model)) or
            (cls.is_quark(ss.legs[1], model) and cls.is_gluon(ss.legs[0], model)) ):
            return None
        # Check that all particles are massless
        for leg in ss.legs:
            if not cls.is_massless(leg, model):
                return None
        # The current is valid and needs no options
        return {}

    def evaluate_kernel(self, parent, zs, kTs):

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


class NLO_FF_QCD_collinear_gg(NLO_FF_QCD_collinear):
    """g g collinear tree-level current."""

    # Override alpha_0, divide_by_jacobian and d_0 here if needed

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to FF NLO QCD collinear
        if cls.common_does_implement_this_current(current, model) is None:
            return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that the unresolved particles are massless gluons
        for leg in ss.legs:
            if not (cls.is_gluon(leg, model) and cls.is_massless(leg, model)):
                return None
        # The current is valid and needs no options
        return {}

    def evaluate_kernel(self, parent, zs, kTs):

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

class NLO_QCD_soft_gluon(NLO_FF_QCD_local_subtraction_current):
    """Soft gluon eikonal current, eq.4.12-4.13 of arXiv:0903.1218."""

    y_0 = 1.
    divide_by_jacobian = False
    d_0_prime = 0

    def __init__(self, *args, **opts):

        super(NLO_QCD_soft_gluon, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False
        
    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to FF NLO QCD
        if cls.common_does_implement_this_current(current, model) is None:
            return None
        # Retrieve the singular structure
        singular_structure = current.get('singular_structure')
        # It should be a soft type of structure
        if singular_structure.name()!='S':
            return None
        # It should not have nested structures
        if singular_structure.substructures:
            return None
        # It should consist in exactly one legs going soft
        if len(singular_structure.legs) != 1:
            return None
        # Legs should be massless gluons
        for leg in singular_structure.legs:
            if not cls.is_gluon(leg, model) or not cls.is_massless(leg, model):
                return None
        # The current is valid and needs no options
        return {}
   
    @staticmethod
    def eikonal(PS_point, i, j, r):
        """Eikonal factor for soft particle with number 'r'
        emitted from 'i' and reconnecting to 'j'.
        """

        pipj = PS_point[i].dot(PS_point[j])
        pipr = PS_point[i].dot(PS_point[r])
        pjpr = PS_point[j].dot(PS_point[r])
        return pipj/(pipr*pjpr)
    
    def evaluate_subtraction_current(self,  current, 
                                            PS_point, 
                                            reduced_process = None,
                                            leg_numbers_map = None,
                                            hel_config = None,
                                            mapping_variables = {}
                                     ):
        """ Evaluates this current and return the corresponding instance of
        SubtractionCurrentResult. See documentation of the mother function for more details."""
        
        if not hel_config is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s" % self.__class__.__name__ +
                " does not support helicity assignment." )
        if leg_numbers_map is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s" % self.__class__.__name__ +
                " requires the leg_number_map." )
        if reduced_process is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s" % self.__class__.__name__ +
                " requires a reduced_process." )
        

        ss = current.get('singular_structure')
        
        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        soft_leg_number = ss.legs[0].n
        # Use the momenta map, in case it has been remapped.
        # Although for the soft current it's typically not the case
        soft_leg_number = leg_numbers_map.inv[frozenset([soft_leg_number,])]
        
        # Now find all colored leg numbers in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color')==1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        # Read mapping variables
        try:
            y        = mapping_variables['y']
            jacobian = mapping_variables['jacobian']
        except:
            raise CurrentImplementationError(
                "Subtraction current %s" % self.__class__.__name__ +
                "requires mapping variables y and jacobian" )

        # Include the counterterm only up to y_0
        if y > self.y_0:
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [],
            'values'              : {}
        })
        
        # Normalization factors
        norm = -8.*math.pi*alpha_s
        if self.divide_by_jacobian: norm /= jacobian
        norm *= (1-y)**(self.d_0_prime-2)

        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        for i, a in enumerate(all_colored_parton_numbers):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in all_colored_parton_numbers[i+1:]:
                evaluation['color_correlations'].append( ((a, b), ) )
                # Write the eikonal for that pair
                evaluation['values'][(0, color_correlation_index)] = {
                    'finite': norm * self.eikonal(PS_point, a, b, soft_leg_number)
                }
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

class NLO_FF_QCD_softcollinear(NLO_QCD_soft_gluon):
    """NLO final soft-collinear currents."""

    def __init__(self, *args, **opts):

        # Make sure it is initialized with the proper set of options and remove them
        # before calling the mother constructor
        if 'color_charge' not in opts:
            raise CurrentImplementationError(
                "The current '%s' must be instantiated with " % self.__class__.__name__ +
                "a 'color_charge' option specified.")
        color_charge = opts.pop('color_charge')

        super(NLO_FF_QCD_softcollinear, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False
        # At this state color_charge is the string of the group factor ('CA' or 'CF');
        # now that the mother constructor has been called,
        # the group factors have been initialized and we can retrieve them.
        self.color_charge = getattr(self, color_charge)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to FF NLO QCD
        if cls.common_does_implement_this_current(current, model) is None:
            return None
        # Retrieve the singular structure
        singular_structure = current.get('singular_structure')
        # It main structure should be of collinear type
        if singular_structure.name() != 'C':
            return None
        # It should have only one leg
        if len(singular_structure.legs) != 1:
            return None
        # It should have exactly one nested structures
        if len(singular_structure.substructures) != 1:
            return None
        # Substructure identified
        sub_singular_structure = singular_structure.substructures[0]
        # Make sure the substructure is soft
        if sub_singular_structure.name() != 'S':
            return None
        # Make sure it contains a single soft leg
        if len(sub_singular_structure.legs) != 1:
            return None
        # Make sure it has no deeper structure
        if sub_singular_structure.substructures:
            return None
        # The hard and soft legs are now identified
        hard_leg = singular_structure.legs[0]
        soft_leg = sub_singular_structure.legs[0]
        # Make sure legs are massless
        if not (cls.is_massless(hard_leg, model) and cls.is_massless(soft_leg, model)):
            return None
        # The leg not soft must be quark or a gluon
        if not cls.is_gluon(hard_leg, model) and not cls.is_quark(hard_leg, model):
            return None
        # The soft leg must be a gluon
        if not cls.is_gluon(soft_leg, model):
            return None
        # This current is implemented here.
        # Return the specific color charge to instantiate this kernel with,
        # in the form of a the name of the group factor to retrieve upon initialization.
        if cls.is_gluon(hard_leg, model):
            # This is a 'g > g g' soft-collinear splitting
            color_charge = 'CA'
        else:
            # This is a 'q > g q' soft-collinear splitting
            color_charge = 'CF'
        return {'color_charge': color_charge}

    def evaluate_subtraction_current(self, current,
                                     PS_point,
                                     reduced_process=None,
                                     leg_numbers_map=None,
                                     hel_config=None,
                                     mapping_variables={}
                                     ):
        """Evaluates this current and return the corresponding instance of
        SubtractionCurrentResult. See documentation of the mother function for more details."""

        if not hel_config is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s" % self.__class__.__name__ +
                " does not support helicity assignment.")

        if leg_numbers_map is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s" % self.__class__.__name__ +
                " requires the leg_number_map.")

        ss = current.get('singular_structure')

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        children = (ss.legs[0].n, ss.substructures[0].legs[0].n)
        parent = leg_numbers_map.inv[frozenset(children)]

        # Read mapping variables
        try:
            y        = mapping_variables['y']
            Q        = mapping_variables['Q']
            jacobian = mapping_variables['jacobian']
        except:
            raise CurrentImplementationError(
                "Subtraction current %s" % self.__class__.__name__ +
                "requires mapping variables y, Q and jacobian" )

        # Include the counterterm only up to y_0
        if y > self.y_0:
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Compute kinematic variables
        na = PS_point[parent]
        nb = Q
        kin_variables = dict()
        mappings.FinalCollinearVariables.get(PS_point, children, na, nb, kin_variables)
        pC = PS_point[children[0]] + PS_point[children[1]]

        # The hard child was placed first in children,
        # so 'z' is the energy fraction of the non-soft particle
        z = kin_variables['z%d' % ss.legs[0].n]
        s12 = pC.square()

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations': [None],
            'color_correlations': [None],
            'values': {(0, 0): {'finite': None}}
        })

        evaluation['values'][(0, 0)]['finite'] = self.color_charge * (2. * z) / (1. - z)

        # Now add the normalization factors
        norm = 4.0 * math.pi * alpha_s * (2.0 / s12)
        if self.divide_by_jacobian: norm /= jacobian
        norm *= (1-y)**(self.d_0_prime-2)

        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result
