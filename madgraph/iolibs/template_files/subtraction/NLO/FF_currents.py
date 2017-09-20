################################################################################
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
################################################################################
"""Implementation of NLO type of currents."""

import os
import sys
import math

import madgraph.core.subtraction as subtraction
import madgraph.integrator.phase_space_generators as PS_utils
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
CS_utils = utils.CS_utils

# The class below isjust a template useful for copy-pasting when starting
# a new implementation
class Template_current(utils.VirtualCurrentImplementation):
    """ Implements the XXXX current."""

    def __init__(self, *args, **opts):
        super(Template_current, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Returns None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning 
        a dictionary, it specifies potential options that must passed upon instantiating
        this implementation for the current given in argument. """
        return None
    
    def evaluate_subtraction_current(self,  current, 
                                            PS_point, 
                                            reduced_process = None,
                                            leg_numbers_map = None,
                                            hel_config = None,
                                            mapping_variables = {}
                                     ):
        """ Evaluates this current and return the corresponding instance of
        SubtractionCurrentResult. See documentation of the mother function for more details."""
        raise NotImplementedError
        
class NLO_FF_QCD_collinear_qqx(utils.VirtualCurrentImplementation):
    """ Implements the GluonToQQbar collinear NLO current."""

    def __init__(self, *args, **opts):
        super(NLO_FF_QCD_collinear_qqx, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Returns None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning 
        a dictionary, it specifies potential options that must passed upon instantiating
        this implementation for the current given in argument. """

        squared_orders = current.get('squared_orders')

        # First check that we indeed have a pure NLO QCD current
        if squared_orders['QCD'] != 2 or \
           any(squared_orders[order] != 0 for order in squared_orders if order!='QCD'):
            return None

        # Now check that it is tree-level
        if current.get('n_loops')!=0:
            return None

        # Make sure we don't need to sum over the quantum number of the mother leg
        if not current.get('resolve_mother_spin_and_color'):
            return None
        
        # Finally check that the singular structure and PDG matches
        singular_structure = current.get('singular_structure')

        # It should be a collinear type of structure
        if singular_structure.name()!='C':
            return None

        # It should not have nested structures
        if singular_structure.substructures:
            return None

        # It should consist in exactly two legs going collinear
        if len(singular_structure.legs)!=2:
            return None

        for leg in singular_structure.legs:
            if leg.state != subtraction.SubtractionLeg.FINAL:
                return None
            if abs(leg.pdg) not in range(1,7):
                return None
            if model.get_particle(leg.pdg).get('mass').upper()!='ZERO':
                return None
        
        # We now know that this current is implemented here, so we return
        # an empty dictionary which could potentially have contained specific
        # options to specify upon instantiating this class.
        return {}

    def get_cache_and_result_key(self,  current, 
                                        PS_point,
                                        reduced_process=None,
                                        leg_numbers_map=None,
                                        hel_config = None,
                                        mapping_variables = {},
                                        **opts
                      ):
        """ If this subtraction current depends on more than just the PS point and
        alpha_s, then complement the cache key here with the additional necessary information. """
        

        cache_key, result_key = super(NLO_FF_QCD_collinear_qqx, self).get_cache_and_result_key(
            current, PS_point,
            reduced_process=reduced_process, leg_numbers_map=leg_numbers_map, 
            hel_config=hel_config, mapping_variables=mapping_variables)
        
        # cache_key['another_call_specifier'] = 'set_it_here'

        return cache_key, result_key

    def evaluate_subtraction_current(self,  current, 
                                            PS_point, 
                                            reduced_process = None,
                                            leg_numbers_map = None,
                                            hel_config = None,
                                            mapping_variables = {}
                                     ):
        """ Now evalaute the current and return the corresponding instance of
        SubtractionCurrentResult. See documentation of the mother function for more details.
        Implementation according to Eq.12 of https://arxiv.org/pdf/hep-ph/9810389.pdf."""

        if not hel_config is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                            "%s does not support helicity assignment."%self.__class__.__name__) 

        if leg_numbers_map is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                      "%s requires the leg_number_map."%self.__class__.__name__)

        result = utils.SubtractionCurrentResult()

        ss = current.get('singular_structure')
        
        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        children_numbers = tuple(leg.n for leg in ss.legs)
        parent_number    = leg_numbers_map.inv[frozenset(children_numbers)]
        
        kin_variables = CS_utils.get_massless_collinear_CS_variables(
                PS_point, parent_number, children_numbers, mapping_variables=mapping_variables)
        z       = kin_variables['za%d'%ss.legs[0].n]
        kT_vec  = kin_variables['nt%d'%ss.legs[0].n]
        s12     = kin_variables['s%d'%parent_number]

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ,
                                      [(parent_number,( tuple(kT_vec), )),],
                                    ],
            'color_correlations'  : [ None ],
            'values'              : { (0,0): { 'finite' : None },
                                      (1,0): { 'finite' : None },
                                  }
          }
        )

        # The two lines below implement the g_\mu\nu part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_i={1,2} \epsilon_i^\mu \epsilon_i^{\star\nu} = g^{\mu\nu} + longitudinal terms
        # are irrelevant because ward identities evaluates them to zero anyway.
        evaluation['values'][(0,0)]['finite'] = 1.
        # The extra overall minus sign comes from the fact that kT_vec has been normalized
        # with kT_vec.square() which is negative.
        evaluation['values'][(1,0)]['finite'] = -4.*z*(1.-z)
                
        # Now add the normalization factors
        norm = 4.*math.pi*alpha_s*(2./s12)*self.TR
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm
        
        result.add_result(evaluation, 
                          hel_config=hel_config, 
                          squared_orders=tuple(sorted(current.get('squared_orders').items()))
                         )
 
        return result

class NLO_FF_QCD_collinear_gq(utils.VirtualCurrentImplementation):
    """ Implements the NLO_FF_QCD_collinear_gq current."""

    def __init__(self, *args, **opts):
        super(NLO_FF_QCD_collinear_gq, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Returns None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning 
        a dictionary, it specifies potential options that must passed upon instantiating
        this implementation for the current given in argument. """
        
        squared_orders = current.get('squared_orders')

        # First check that we indeed have a pure NLO QCD current
        if squared_orders['QCD'] != 2 or \
           any(squared_orders[order] != 0 for order in squared_orders if order!='QCD'):
            return None

        # Now check that it is tree-level
        if current.get('n_loops')!=0:
            return None

        # For this current, it doesn't matter if one needs to sum over the 
        # quantum numbers of the mother leg
        # if not current.get('resolve_mother_spin_and_color'):
        #    return None
        
        # Finally check that the singular structure and PDG matches
        singular_structure = current.get('singular_structure')

        # It should be a collinear type of structure
        if singular_structure.name()!='C':
            return None

        # It should not have nested structures
        if singular_structure.substructures:
            return None

        # It should consist in exactly two legs going collinear
        if len(singular_structure.legs)!=2:
            return None

        for leg in singular_structure.legs:
            if leg.state != subtraction.SubtractionLeg.FINAL:
                return None
            if model.get_particle(leg.pdg).get('mass').upper()!='ZERO':
                return None
        
        # Check that it is indeed a quark and a gluon going collinear
        if set([abs(leg.pdg) for leg in singular_structure.legs]) not in \
                                [set([21,q_pdg]) for q_pdg in range(1,7)]:
            return None            
        
        # We now know that this current is implemented here, so we return
        # an empty dictionary which could potentially have contained specific
        # options to specify upon instantiating this class.
        return {}

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
            raise CurrentImplementationError("Subtraction current implementation "+
                            "%s does not support helicity assignment."%self.__class__.__name__) 

        if leg_numbers_map is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                      "%s requires the leg_number_map."%self.__class__.__name__)

        result = utils.SubtractionCurrentResult()

        ss = current.get('singular_structure')
        
        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        children_numbers = tuple(leg.n for leg in ss.legs)
        parent_number    = leg_numbers_map.inv[frozenset(children_numbers)]
        
        kin_variables = CS_utils.get_massless_collinear_CS_variables(
                PS_point, parent_number, children_numbers, mapping_variables=mapping_variables)
        z       = kin_variables['za%d'%ss.legs[0].n]
        s12     = kin_variables['s%d'%parent_number]

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [ None ],
            'values'              : { (0,0): { 'finite' : None } }
          }
        )
        
        if ss.legs[0].pdg == 21:
            # z is the energy fraction of the gluon, therefore use P_{gq}
            kernel = (1.+(1.-z)**2)/z
        else:
            # z is the energy fraction of the quark, therefore use P_{qg}
            kernel = (1.+z**2)/(1.-z)

        evaluation['values'][(0,0)]['finite'] = kernel
            
        # Now add the normalization factors
        norm = 4.0*math.pi*alpha_s*(2.0/s12) * self.CF
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm
        
        result.add_result(evaluation, 
                          hel_config=hel_config, 
                          squared_orders=tuple(sorted(current.get('squared_orders').items()))
                         )
 
        return result

class NLO_FF_QCD_collinear_gg(utils.VirtualCurrentImplementation):
    """ Implements the NLO_FF_QCD_collinear_gg current."""

    def __init__(self, *args, **opts):
        super(NLO_FF_QCD_collinear_gg, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Returns None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning 
        a dictionary, it specifies potential options that must passed upon instantiating
        this implementation for the current given in argument. """
        
        squared_orders = current.get('squared_orders')

        # First check that we indeed have a pure NLO QCD current
        if squared_orders['QCD'] != 2 or \
           any(squared_orders[order] != 0 for order in squared_orders if order!='QCD'):
            return None

        # Now check that it is tree-level
        if current.get('n_loops')!=0:
            return None

        # Make sure we don't need to sum over the quantum number of the mother leg
        if not current.get('resolve_mother_spin_and_color'):
            return None
        
        # Finally check that the singular structure and PDG matches
        singular_structure = current.get('singular_structure')

        # It should be a collinear type of structure
        if singular_structure.name()!='C':
            return None

        # It should not have nested structures
        if singular_structure.substructures:
            return None

        # It should consist in exactly two legs going collinear
        if len(singular_structure.legs)!=2:
            return None

        # They should be massless final states
        for leg in singular_structure.legs:
            if leg.state != subtraction.SubtractionLeg.FINAL:
                return None
            if model.get_particle(leg.pdg).get('mass').upper()!='ZERO':
                return None
        
        # Check that it is indeed two gluons going collinear
        if [abs(leg.pdg) for leg in singular_structure.legs] != [21,21]:
            return None
        
        # We now know that this current is implemented here, so we return
        # an empty dictionary which could potentially have contained specific
        # options to specify upon instantiating this class.
        return {}
    
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
            raise CurrentImplementationError("Subtraction current implementation "+
                            "%s does not support helicity assignment."%self.__class__.__name__) 

        if leg_numbers_map is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                      "%s requires the leg_number_map."%self.__class__.__name__)

        result = utils.SubtractionCurrentResult()

        ss = current.get('singular_structure')
        
        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        children_numbers = tuple(leg.n for leg in ss.legs)
        parent_number    = leg_numbers_map.inv[frozenset(children_numbers)]
        
        kin_variables = CS_utils.get_massless_collinear_CS_variables(
                PS_point, parent_number, children_numbers, mapping_variables=mapping_variables)
        z       = kin_variables['za%d'%ss.legs[0].n]
        kT_vec  = kin_variables['nt%d'%ss.legs[0].n]
        s12     = kin_variables['s%d'%parent_number]

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ,
                                      [(parent_number,( tuple(kT_vec), )),],
                                    ],
            'color_correlations'  : [ None ],
            'values'              : { (0,0): { 'finite' : None },
                                      (1,0): { 'finite' : None },
                                  }
          }
        )

        # The two lines below implement the g_\mu\nu part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_i={1,2} \epsilon_i^\mu \epsilon_i^{\star\nu} = g^{\mu\nu} + longitudinal terms
        # are irrelevant because ward identities evaluates them to zero anyway.
        evaluation['values'][(0,0)]['finite'] = ( (z/(1.-z)) + ((1.-z)/z) )
        # The extra overall minus sign comes from the fact that kT_vec has been normalized
        # with kT_vec.square() which is negative.
        evaluation['values'][(1,0)]['finite'] = 2.*z*(1.-z)
                
        # Now add the normalization factors
        norm = 4.*math.pi*alpha_s*(2./s12)*2.*self.CA
        
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm
        
        result.add_result(evaluation, 
                          hel_config=hel_config, 
                          squared_orders=tuple(sorted(current.get('squared_orders').items()))
                         )
 
        return result

class NLO_QCD_soft_gluon(utils.VirtualCurrentImplementation):
    """ Implements the soft gluon Eikonel current.
    See Eq.4.12-4.13 of ref. https://arxiv.org/pdf/0903.1218.pdf"""

    def __init__(self, *args, **opts):
        super(NLO_QCD_soft_gluon, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False
        
    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Returns None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning 
        a dictionary, it specifies potential options that must passed upon instantiating
        this implementation for the current given in argument. """

        squared_orders = current.get('squared_orders')

        # First check that we indeed have a pure NLO QCD current
        if squared_orders['QCD'] != 2 or \
           any(squared_orders[order] != 0 for order in squared_orders if order!='QCD'):
            return None

        # Now check that it is tree-level
        if current.get('n_loops')!=0:
            return None

        # Make sure we don't need to sum over the quantum number of the mother leg
        if not current.get('resolve_mother_spin_and_color'):
            return None
        
        # Finally check that the singular structure and PDG matches
        singular_structure = current.get('singular_structure')

        # It should be a soft type of structure
        if singular_structure.name()!='S':
            return None

        # It should not have nested structures
        if singular_structure.substructures:
            return None

        # It should consist in exactly one legs going soft
        if len(singular_structure.legs)!=1:
            return None

        for leg in singular_structure.legs:
            if leg.state != subtraction.SubtractionLeg.FINAL:
                return None
            if abs(leg.pdg) not in [21]:
                return None
            if model.get_particle(leg.pdg).get('mass').upper()!='ZERO':
                return None
        
        # We now know that this current is implemented here, so we return
        # an empty dictionary which could potentially have contained specific
        # options to specify upon instantiating this class.
        return {}
   
    @classmethod
    def eikonal(cls, PS_point, i, j, r):
        """ Eikonal factor for soft particle with number 'r' emitted from 'i' and reconnecting
        to 'j'."""
        
        return ( PS_point[i].dot(PS_point[j]) ) / ( 
                 PS_point[i].dot(PS_point[r])*PS_point[j].dot(PS_point[r])
                )
    
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
            raise CurrentImplementationError("Subtraction current implementation "+
                            "%s does not support helicity assignment."%self.__class__.__name__) 

        if leg_numbers_map is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                      "%s requires the leg_number_map."%self.__class__.__name__)
        
        if reduced_process is None:
            raise CurrentImplementationError("Subtraction current implementation "+
                                      "%s requires a reduced_process."%self.__class__.__name__)
        
        result = utils.SubtractionCurrentResult()

        ss = current.get('singular_structure')
        
        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Retrieve kinematic variables from the specified PS point
        soft_leg_number = ss.legs[0].n
        # Use the momenta map, in case it has been remapped.
        # Although for the soft current it's typically not the case
        soft_leg_number   = leg_numbers_map.inv[frozenset([soft_leg_number,])]
        
        # Now find all colored leg numbers in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color')==1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [ ],
            'values'              : { }
          })
        
        # Normalization factors
        norm = -8.*math.pi*alpha_s

        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b) and add the corrsponding 
        # contributions to this current
        for i, a in enumerate(all_colored_parton_numbers):
            # Use the symmetricity of the color-correlation and soft current (a,b) <-> (b,a)
            for b in all_colored_parton_numbers[i+1:]:
                evaluation['color_correlations'].append( (a, b) )
                # Write the eikonal for that pair
                evaluation['values'][(0,color_correlation_index)] = {
                    'finite': norm * self.eikonal(PS_point, a, b, soft_leg_number)
                }
                color_correlation_index += 1
        
        result.add_result(evaluation, 
                          hel_config=hel_config, 
                          squared_orders=tuple(sorted(current.get('squared_orders').items()))
                         )
 
        return result