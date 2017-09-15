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
"""All virtual mother classes to IR-subtraction related hard-coded 
implementations"""

import madgraph.various.misc as misc
import madgraph.integrator.phase_space_generators as PS_utils

#===============================================================================
# Subtraction current evaluation and result
#===============================================================================
class SubtractionCurrentEvaluation(dict):
    """ A container class to specify the output of the evaluation of a given 
    current."""

    # All residues of the poles for the epsilon expansion can be specified with a
    # key names 'eps^n'.
    
    # values will list the weight from this current for each particular pairing
    # of the spin_correlations and color_correlations
    main_layer_result_order  = ['spin_correlations','color_correlations','values']
    sub_layer_result_order   = ['finite','return_code','accuracy','eps'] 
    result_order = main_layer_result_order+sub_layer_result_order

    main_layer_result_format = {'spin_correlations':'%s','color_correlations':'%s','values':'%s'}
    sub_layer_result_format  = {'tree':'%.15e','finite':'%.15e','eps':'%.15e',
                                'return_code':'%d','accuracy': '%.2g'}

    result_format = dict(it for it in main_layer_result_format.items()+sub_layer_result_format.items())

    def __init__(self, *args, **opts):
        super(SubtractionCurrentEvaluation, self).__init__(*args, **opts)      

    @classmethod
    def get_max_length_attribute(cls):
        return max(len(k) for k in cls.result_order)

    def get_result_order(self, result_name):
        """ Returns an index that specifies the ordering of the results to be displayed, for the 
        purpose of sorting."""
        
        index = 0
        if result_name.startswith('eps'):
            index += int(result_name.split('^')[1])
            result_name = 'eps'
        
        return index + ((100*self.result_order.index(result_name)) if result_name in self.result_order else 100000)

    def nice_string(self, max_len_attribute=-1):
        """ Formats nicely the output of a particular subtraction current evaluation."""
        
        res = []
        length = (max_len_attribute if max_len_attribute>0 else max(len(k) for k in self.keys()))
        template = '%%-%ds = %%s'%length
        subtemplate = '%%-%ds = %%s'%max((length-5),10)

        sorted_result_keys = sorted(self.keys(), key=lambda el: self.get_result_order(el))

        def format_result(key, res):

            # Special handling for the subresult for each spin/lorentz connection pairing
            if key=='values':
                if len(res)==0:
                    return 'N/A'
                formatted_res = ['']
                for spin_lorentz_pair in res:
                    formatted_res.append(misc.bcolors.BLUE+'     %s:'%str(spin_lorentz_pair)+misc.bcolors.ENDC)
                    sorted_subresult_keys = sorted(res[spin_lorentz_pair].keys(), 
                                                    key=lambda el: self.get_result_order(el))
                    for subkey in sorted_subresult_keys:
                        formatted_res.append(('       -> %s'%subtemplate)%(subkey, 
                                                format_result(subkey, res[spin_lorentz_pair][subkey])))
                return '\n'.join(formatted_res)

            if res is None or key=='accuracy' and res < 0.:
                return 'N/A'
            if key.startswith('eps'):
                key = 'eps'
            formatted_res = self.result_format[key]%res if result_key in self.result_format else str(res)
            if isinstance(res,float) and res > 0.:
                formatted_res = ' %s'%formatted_res
            return formatted_res

        res.append(misc.bcolors.BLUE+'Result:'+misc.bcolors.ENDC)
        for result_key in sorted_result_keys:
            res.append(('  -> %s'%template)%(result_key, format_result(result_key, self[result_key])))
        
        return '\n'.join(res)
    
    def __str__(self):
        return self.nice_string()

class SubtractionCurrentResult(dict):
    """ A class to store th e different results of current evaluation call for one specific PS point / scale."""

    def __init(self, *args, **opts):
        super(SubtractionCurrentResult, self).__init__(*args, **opts)        
        
    def nice_string(self):
        """ Print out all the results available in a nice form."""
        # First lists the result with the least amount of attributes specified
        sorted_keys = sorted(self.keys(), key=lambda k:[el[1] for el in k].count(None))

        res = []
        max_len_attribute = max(max(len(el[0]) for el in k) for k in sorted_keys)
        max_len_attribute = max(max_len_attribute, SubtractionCurrentEvaluation.get_max_length_attribute())
        template = '%%-%ds : %%s'%max_len_attribute
        for i,k in enumerate(sorted_keys):
            res.append(misc.bcolors.GREEN+'Entry #%d:'%(i+1)+misc.bcolors.ENDC)
            for name, value in k:
                if value is None:
                    continue
                res.append(('  -> %s'%template)%(name, str(value)))
            if any(el[1] is None for el in k):
                res.append('  Unspecified or summed over properties:')
                res.append('  -> %s'%(', '.join(el[0] for el in k if el[1] is None)))
            res.extend('  %s'%line for line in self[k].nice_string(
                                            max_len_attribute=max_len_attribute).split('\n'))
        
        return '\n'.join(res)

    def __str__(self):
        return self.nice_string()

    def get_result(self, **opts):
        """ Attemp ts to recycle the result from previous computations.
        Returns None otherwise. Opts are typically:
           helicities: A tuple of integer, None means summed over.
           squared_orders: a tuple repr. of the usual dict., None means all contributions.
        """
        key_opts = {'hel_config'           : None,
                    'squared_orders'       : None}
        key_opts.update(opts)
        
        try:
            result = self[tuple(sorted(opts.items()))]
#            misc.sprint('Recycled a result.')
            return result
        except KeyError:
            return None
    
    def add_result(self, value, **opts):
        """ Ad d a result to the current record."""
        
        key_opts = {'hel_config'           : None,
                    'squared_orders'       : None}
        key_opts.update(opts)
        
        self[tuple(sorted(opts.items()))] = value

#===============================================================================
# CurrentImplementationError
#===============================================================================
class CurrentImplementationError(Exception):
    """Exception raised if an exception is triggered in implementation of the currents.""" 

#===============================================================================
# VirtualCurrentImplementation
#===============================================================================
class VirtualCurrentImplementation(object):
    """A virtual class defining what a current implementation must specify"""
     
    def __init__(self, model, **opts):
        """ Saves some general properities or quantities useful for the current evaluation."""
        
        # A property variable specifying if helicity assignment is possible with
        # this implementation
        self.supports_helicity_assignment = True

        self.model = model
        # Extract some constants from the UFO model if present, otherwise take default values
        try:
            model_param_dict = self.model.get('parameter_dict')
        except:
            model_param_dict = {}
        try:
            self.TR = model_param_dict['TR']
        except:
            self.TR = 0.5
        try:
            self.CF = model_param_dict['CF']
        except:
            self.CF = 4.0/3.0
        try:
            self.NC = model_param_dict['NC']
        except:
            self.NC = 3.0

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ Returns None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements. When returning 
        a dictionary, it specifies potential options that must passed upon instanciating
        the current implementation for the current given in argument. """
        
        # This virtual class of course does not implement any current.
        return None 

    def get_cache_and_result_key(self,  current, 
                                        PS_point,
                                        reduced_process=None,
                                        leg_numbers_map=None,
                                        hel_config = None,
                                        mapping_variables = {},
                                        **opts
                                ):
        """ Generates a key for the cache dictionary. Make sure that everything that can
        lead to a different evaluation of the current ends up in the key so that two evaluations
        that lead to different results cannot be wrongly recycled.
        The minimal key is simply the PS point and alpha_s.
        Do not bother with helicity configurations and squared orders since those are handled in 
        the subtraction current result structure directly.
        The result key is trivial and is direclty related to the structure of the SubtractionCurrentResult."""

        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']

        cache_key = {'PS_point': tuple(sorted([(k, tuple(value)) for (k, value) in PS_point.items()])),
                     'alpha_s' : model_param_dict['aS']
                     }
        
        result_key = {'hel_config':hel_config,
                      'squared_orders': tuple(sorted(current.get('squared_orders').items()))}

        return cache_key, result_key
        # Alternatively, use: 
        #return None, result_key 
        # to disable the caching system

    def evaluate_subtraction_current(self,  current, 
                                            PS_point,
                                            reduced_process = None,
                                            leg_numbers_map = None,
                                            hel_config = None,
                                            mapping_variables = {}
                                     ):
        """ Returns an instance of SubtractionCurrentResult, with SubtractionCurrentEvaluation as keys which store 
        various output of this evaluation of the subtraction current, most importantly its 
        finite part and other coefficients of its epsilon expansion.
        The model provides quantum numbers of legs in the reduced process (useful for soft currents)
        and the bidirectional leg_numbers_map provides the number of the mother leg which is useful
        to provide complete spin-correlation specification.
        The value of alpha_s and mu_r can be retrieved from the model.
        Seee the documentation of the class method 'apply_permutations' of VirtualMEAccessor in the module
        madgraph.core.contributions for the specification of the syntax of spin_correlations and
        color_correlations.
        Mapping variables can be passed and possibly reused in the subtraction current if they are the same.
        """
        ## Example:
        #
        # subtraction_current_result = SubtractionCurrentResult()
        #
        # one_eval_for_QCD_splitting = SubtractionCurrentEvaluation({
        #   'spin_correlations'   : [ 
        #                             [ (3,[(1.1,1.2,1.3,1.4),(2.1,2.2,2.3,2.4)]),
        #                               (4,[(1.1,1.2,1.3,1.4),(2.1,2.2,2.3,2.4)]) ],
        #                             [ (3,[(5.1,5.2,5.3,5.4),(6.1,6.2,6.3,6.4)]),
        #                               (4,[(8.1,8.2,8.3,8.4),(7.1,7.2,7.3,7.4)]) ],
        #                           ]
        #   'color_correlations'  : [ (2,4),
        #                             (2,5)
        #                           ],
        #   'values'              : { (0,0): { 'finite' : 32.33,
        #                                      'eps^-1' : 21.2,
        #                                      'eps^-2' : 42.6,
        #                                    },
        #                             (1,0): { 'finite' : 52.33,
        #                                      'eps^-1' : 71.2,
        #                                      'eps^-2' : 82.6,
        #                                    },
        #                             (1,1): { 'finite' : 32.33,
        #                                      'eps^-1' : 21.2,
        #                                      'eps^-2' : 42.6,
        #                                    }
        #                           }
        #   })
        #
        ##  In the example above there was no contribution from the first spin correlation paired with the second
        ##  color correlations, i.e. the key (0,1). We stress that most of the time the structure will be significantly
        ##  simple from the example above (i.e. in principle elementary current only have either spin or color correlations,
        ##  not both.)
        #
        # subtraction_current_result.add_result(one_eval_for_QCD_squared_order, hel_config=hel_config, squared_orders={'QCD':2,'QED':0})
        #
        # another_eval_for_EW_splitting = SubtractionCurrentEvaluation({
        #   'spin_correlations'   = [etc...] 
        #   })
        # subtraction_current_result.add_result(another_eval_for_EW_splitting, hel_config=hel_config, squared_orders={'QCD':0,'QED':2})
        #
        # return subtraction_current_result

        raise NotImplemented

class DefaultCurrentImplementation(VirtualCurrentImplementation):
    """ This default implementation class will be used with a warning and *only* if none of the other implementation matches and 
    the function 'does_implement_this_current' of this class evalauted to true.
    This is typically useful for debugging and one has not completed the implementation of all currents but already wants to test
    a subset of them."""
    
    def __init__(self, *args, **opts):
        super(DefaultCurrentImplementation, self).__init__(*args, **opts)
        self.supports_helicity_assignment = True

    @classmethod
    def does_implement_this_current(cls, current, model):
        """ For production, it is preferable to turn off this default implementation by simply returnin None below instead of {}."""
        
        # Simply uncomment the line below to de-activate this default implementation.
        # return None
        return {}
    
    def get_cache_and_result_key(self, *args, **opts):
        return super(DefaultCurrentImplementation, self).__init__(*args, **opts)

    def evaluate_subtraction_current(self, current, PS_point,
                                            reduced_process = None,
                                            leg_numbers_map = None,
                                            hel_config = None,
                                            mapping_variables = {}
                                     ):
        """ Simply return 1.0 for this current default implementation."""
        
        subtraction_current_result = SubtractionCurrentResult()
        
        subtraction_current_eval = SubtractionCurrentEvaluation({
                'spin_correlations'   : [1],
                'color_correlations'  : [1],
                'values'              : {(0,0): {'finite' : '0.0',
                                                 'eps^-1' : '0.0',
                                                 'eps^-2' : '0.0'}}
            })
        
        subtraction_current_result.add_result(subtraction_current_eval, 
                                              hel_config=hel_config, 
                                              squared_orders=current.get('squared_orders'))
        
        return subtraction_current_result  

#===============================================================================
# CS util
#===============================================================================
class CS_utils(object):
    """ A container function for useful class methods for NLO currents ala CS."""

    @classmethod
    def get_massless_collinear_CS_variables(cls, 
                    PS_point, parent_number, children_numbers, mapping_variables={}):
        """ Returns the collinear splitting variables following Catani-Grazzini
        conventions (Eq.6 of https://arxiv.org/pdf/hep-ph/9810389.pdf).
        This is very similar to the function get_collinear_variables() of
        phase_space_generators.ElementaryMappingCollinearFinal.get_collinear_variables"""
        
        kin_variables = dict( ('z%d'%n,0.0) for n in children_numbers[:-1] )
        kin_variables.update( dict( ('kt%d'%n,0.0) for n in children_numbers[:-1] ) )
        kin_variables['s%d'%parent_number] = 0.0
        
        # Attempt recycling the variables from the mapping
        missing_variable = False
        for var in kin_variables:
            try:
                kin_variables[var] = mapping_variables[var]
            except KeyError:
                missing_variable = True
                break
        if not missing_variable:
            return kin_variables

        # Retrieve the parent's momentum
        p = PS_utils.LorentzVector(PS_point[parent_number])
        # Compute the sum of momenta
        q = PS_utils.LorentzVector(4)
        for i in children_numbers:
            q += PS_utils.LorentzVector(PS_point[i])
        # Pre-compute scalar products
        q2 = q.square()
        kin_variables['s%d'%parent_number] = q2
        pq = p.dot(q)
        # Compute all kinematic variables
        for i in children_numbers[:-1]:
            pi = PS_utils.LorentzVector(PS_point[i])
            ppi = p.dot(pi)
            qpi = q.dot(pi)
            zi = 2*qpi/q2 - ppi/pq
            kti = pi + (q2*ppi - pq*qpi)/(pq**2) * p - ppi/pq * q
            kin_variables['z%d'%i] = zi
            kin_variables['kt%d'%i] = kti
        return kin_variables
