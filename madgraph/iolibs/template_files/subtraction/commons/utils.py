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
from __builtin__ import classmethod
"""All virtual mother classes to IR-subtraction related hard-coded implementations"""

import math
import madgraph.various.misc as misc
import madgraph.integrator.phase_space_generators as PS_utils
from madgraph.core.base_objects import EpsilonExpansion

#=========================================================================================
# Useful constants
#=========================================================================================
class Constants(object):
    """Constants used throughout the implementation of the counterterms."""
    
    # Epsilon expansion constants
    EulerGamma = 0.57721566490153286061
    SEpsilon =  EpsilonExpansion({ 
         0 : 1., 
         1 : -EulerGamma + math.log(4.*math.pi),
         2 : 0.5*(EulerGamma**2-2.*EulerGamma*math.log(4.*math.pi)+math.log(4.*math.pi)**2)
    })
    
    # SU(3) group constants
    TR = 0.5
    NC = 3.0
    CF = (NC ** 2 - 1) / (2 * NC)
    CA = NC

#=========================================================================================
# Subtraction current evaluation and result
#=========================================================================================

class SubtractionCurrentEvaluation(dict):
    """Container class for the output of the evaluation of a current."""

    # All epsilon expansion poles can be specified with key names 'eps^n'.
    
    # values will list the weight from this current for each particular pairing
    # of the spin_correlations and color_correlations
    main_layer_result_order  = ['spin_correlations','color_correlations',
                                'Bjorken_rescalings','reduced_kinematics','values']
    sub_layer_result_order   = ['finite','return_code','accuracy','eps'] 
    result_order = main_layer_result_order+sub_layer_result_order

    main_layer_result_format = {'spin_correlations':'%s','color_correlations':'%s',
                                'Bjorken_rescalings':'%s','reduced_kinematics':'%s','values':'%s'}
    sub_layer_result_format  = {'finite':'%.15e','eps':'%.15e',
                                'return_code':'%d','accuracy': '%.2g'}

    result_format = dict(it for it in main_layer_result_format.items()+sub_layer_result_format.items())

    def __init__(self, *args, **opts):
        super(SubtractionCurrentEvaluation, self).__init__(*args, **opts)      

    @classmethod
    def get_max_length_attribute(cls):
        return max(len(k) for k in cls.result_order)

    def get_result_order(self, result_name):
        """Return an index that specifies the ordering of the results to be displayed,
        for the purpose of sorting."""
        
        index = 0
        if isinstance(result_name, int):
            index+=result_name
            result_name='eps'
        elif result_name.startswith('eps'):
            index += int(result_name.split('^')[1])
            result_name = 'eps'
        return index + ((100*self.result_order.index(result_name))
                        if result_name in self.result_order else 100000)

    def subresult_lines(self, subresult, subtemplate):
        """ Returns the string line describing the subresult in argument corresponding to 
        a paricular combination of color and Lorentz structure. The subtemplate specified
        will be used for rendering each line."""
        
        lines = []        
        sorted_subresult_keys = sorted(subresult.keys(), 
                                                 key=lambda el: self.get_result_order(el))
        for subkey in sorted_subresult_keys:
            if isinstance(subkey, int):
                key_name='eps^%d'%subkey
            else:
                key_name=subkey
            lines.append(subtemplate%(key_name, 
                                    self.format_result(subkey, subresult[subkey])))
        return lines

    def format_result(self, key, res):
        """ Format a string for the pair key-result specified."""
        
        # Special handling for the subresult for each spin/lorentz connection pairing
        if key=='values':
            if len(res)==0:
                return 'N/A'
            formatted_res = ['']
            for spin_lorentz_pair in res:
                formatted_res.append(misc.bcolors.BLUE+'     %s:'%str(spin_lorentz_pair)+misc.bcolors.ENDC)
                formatted_res.extend(['       -> %s'%line for 
                        line in self.subresult_lines(res[spin_lorentz_pair],'%-10s=%s')])
            return '\n'.join(formatted_res)

        if res is None or key=='accuracy' and res < 0.:
            return 'N/A'
        if isinstance(key, int):
            key = 'eps'
        elif key.startswith('eps'):
            key = 'eps'
        formatted_res = self.result_format[key]%res if key in self.result_format else str(res)
        if isinstance(res,float) and res > 0.:
            formatted_res = ' %s'%formatted_res
        return formatted_res

    def nice_string(self, max_len_attribute=-1):
        """Formats nicely the output of a particular subtraction current evaluation."""
        
        res = []
        length = (max_len_attribute if max_len_attribute>0 else max(len(k) for k in self.keys()))
        template = '%%-%ds = %%s'%length
        subtemplate = '%%-%ds = %%s'%max((length-5), 10)

        sorted_result_keys = sorted(self.keys(), key=lambda el: self.get_result_order(el))
        
        res.append(misc.bcolors.BLUE+'Result:'+misc.bcolors.ENDC)
        for result_key in sorted_result_keys:
            res.append(('  -> %s'%template)%(result_key, self.format_result(result_key, self[result_key])))
        
        return '\n'.join(res)
    
    def __str__(self):

        return self.nice_string()

    @classmethod
    def zero(cls, spin_correlations=None, color_correlations=None, Bjorken_rescalings=(None,None), reduced_kinematics=None):

        return SubtractionCurrentEvaluation({
            'spin_correlations'   : [ spin_correlations, ],
            'color_correlations'  : [ color_correlations, ],
            'Bjorken_rescalings'  : [ Bjorken_rescalings, ],
            'reduced_kinematics'  : [ reduced_kinematics, ],
            'values'              : {(0,0,0,0): { 'finite' : 0.0 }}
        })

class BeamFactorizationCurrentEvaluation(SubtractionCurrentEvaluation):
    """Container class for the output of the evaluation of beam factorization current.
    The main difference w.r.t the mother class is that values of the dictionary in the
    attribute 'values' of this class are themselves dictionaries representing the flavor
    matrix."""

    main_layer_result_order  = ['spin_correlations','color_correlations','Bjorken_rescalings','reduced_kinematics','values']
    main_layer_result_format = {'spin_correlations':'%s','color_correlations':'%s',
                                'Bjorken_rescalings':'%s', 'reduced_kinematics':'%s', 'values':'%s'}

    def subresult_lines(self, subresult, subtemplate):
        """ Returns the string line describing the subresult in argument corresponding to 
        a paricular combination of color and Lorentz structure. The subtemplate specified
        will be used for rendering each line. Here we must process a matrix output."""
        
        lines = []
        reduced_to_resolved_flavors = []
        for reduced_IS_flavor_PDG, values1 in subresult.items():
            for resolved_IS_flavor_PDGs, values2 in values1.items():
                reduced_to_resolved_flavors.append((reduced_IS_flavor_PDG, resolved_IS_flavor_PDGs))
        
        reduced_to_resolved_flavors.sort()
        for (reduced_IS_flavor_PDG, resolved_IS_flavor_PDGs) in reduced_to_resolved_flavors:
            if reduced_IS_flavor_PDG is None and resolved_IS_flavor_PDGs is None:
                lines.append('Diagonal flavor configuration')
            else:
                lines.append('Flavor configuration: %d -> (%s)'%( reduced_IS_flavor_PDG,
                    ','.join('%d'%pdg for pdg in resolved_IS_flavor_PDGs) +
                    (',' if len(resolved_IS_flavor_PDGs)==1 else '') ))

            values = subresult[reduced_IS_flavor_PDG][resolved_IS_flavor_PDGs]
            sorted_subresult_keys = sorted(values.keys(), 
                                                     key=lambda el: self.get_result_order(el))
            for subkey in sorted_subresult_keys:
                if isinstance(subkey, int):
                    key_name='eps^%d'%subkey
                else:
                    key_name=subkey
                lines.append('   %s'%(subtemplate%(key_name, 
                                             self.format_result(subkey, values[subkey]))))
        return lines

    def __init__(self, *args, **opts):
        super(SubtractionCurrentEvaluation, self).__init__(*args, **opts)      

    @classmethod
    def zero(cls, spin_correlations=None, color_correlations=None, Bjorken_rescalings=(None,None), reduced_kinematics=None):
        
        # None for reduced_IS_flavor_PDG and resolved_IS_flavor_PDGs means purely diagonal
        return SubtractionCurrentEvaluation({
            'spin_correlations'      : [ spin_correlations, ],
            'color_correlations'     : [ color_correlations, ],
            'Bjorken_rescalings'     : [ Bjorken_rescalings, ],
            'reduced_kinematics'     : [ reduced_kinematics, ],
            'values'                 : { (0,0,0,0): {
                    None : { # reduced_IS_flavor_PDG
                        None : { 'finite' : 0.0 } #resolved_IS_flavor_PDGs
                    }
                }
            }
        })

class SubtractionCurrentResult(dict):
    """A class to store the different results of current evaluation call
    for one specific PS point / scale."""

    def __init(self, *args, **opts):
        super(SubtractionCurrentResult, self).__init__(*args, **opts)        
        
    def nice_string(self):
        """Print out all the results available in a nice form."""
        # First lists the result with the least amount of attributes specified
        sorted_keys = sorted(self.keys(), key=lambda k:[el[1] for el in k].count(None))

        res = []
        max_len_attribute = max(max(len(el[0]) for el in k) for k in sorted_keys)
        max_len_attribute = max(
            max_len_attribute, SubtractionCurrentEvaluation.get_max_length_attribute())
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
        """Attempt to recycle the result from previous computations,
        and return None otherwise.
        Opts are typically:
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
        """Add a result to the current record."""
        
        key_opts = {'hel_config'           : None,
                    'squared_orders'       : None}
        key_opts.update(opts)
        
        self[tuple(sorted(opts.items()))] = value

    @staticmethod
    def zero(squared_orders=None, current=None, hel_config=None, **opts):
        """Return a 'zero' result."""

        if squared_orders is None and current is not None:
            sqo = tuple(sorted(current.get_squared_orders().items()))
        else:
            sqo = squared_orders
        subtraction_current_result = SubtractionCurrentResult()
        subtraction_current_result.add_result(
            SubtractionCurrentEvaluation.zero(**opts),
            hel_config=hel_config, squared_orders=sqo)
        return subtraction_current_result

#=========================================================================================
# CurrentImplementationError
#=========================================================================================

class CurrentImplementationError(Exception):
    """Exception raised if an exception is triggered in implementation of the currents.""" 
    pass

#=========================================================================================
# VirtualCurrentImplementation
#=========================================================================================

class VirtualCurrentImplementation(object):
    """A virtual class defining what a current implementation must specify"""

    # Class parameter to detect currents that vanish and do not need to be evaluated
    is_zero = False

    def __init__(self, model, **opts):
        """Save properties or quantities useful for the current evaluation."""
        
        # A property variable specifying if helicity assignment is possible with
        # this implementation
        self.supports_helicity_assignment = True

        self.model = model

    @classmethod
    def name(cls):
        """Extended name string to print the identity of current implementations."""

        return "Subtraction current implementation " + cls.__name__

    @classmethod
    def does_implement_these_currents(cls, currents, model):
        """Return None/a_dictionary depending on whether this particular current is
        part of what this particular current class implements.
        When returning a dictionary, it specifies potential options that must be passed
        upon instantiating the current implementation.
        """
        
        # This virtual class of course does not implement any current.
        raise NotImplementedError("The function 'does_implement_these_currents' must be implemented in class %s"%(cls.name()))

    def get_cache_and_result_key(
        self, currents_block, higher_PS_point=None, reduced_process=None, hel_config=None, **opts ):
        """Generate a key for the cache dictionary.
        Make sure that everything that can lead to a different evaluation of the current
        ends up in the key so that two evaluations that lead to different results
        cannot be wrongly recycled.
        The minimal key is simply the PS point and alpha_s.
        Do not bother with helicity configurations and squared orders
        since those are handled in the subtraction current result structure directly.
        The result key is trivial and is directly related to the structure
        of the SubtractionCurrentResult.
        """

        # The caching system is computationally heavy and is not so useful since we now
        # anticipate having a rust backend. We therefore disable it in
        # Also, for now do not differentiate squared orders.
        squared_orders = tuple(sorted(currents_block.get_squared_orders().items()))
        result_key = {'hel_config': hel_config, 'squared_orders': squared_orders}


        # Returning None for the cache key effectively disables the cache system.
        # A particular implementation of the subtraction scheme is free to modify this
        # if it wants to recover a proper cache (for example as done below).
        return None, result_key

        singular_structure_str = tuple(crt['singular_structure'].__str__(
            print_n=True, print_pdg=True, print_state=True) for crt in currents_block)
        higher_PS_point_tuple = tuple(sorted([
            (k, tuple(value)) for (k, value) in higher_PS_point.items() ]))
        if reduced_process:
            reduced_process_hash = tuple(
                (leg['number'], leg['id']) for leg in reduced_process['legs'] )
        else:
            reduced_process_hash = None
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']

        cache_key = {
            'higher_PS_point': higher_PS_point_tuple,
            'alpha_s': alpha_s,
            'singular_structure': singular_structure_str,
            'reduced_process': reduced_process_hash
        }

        return cache_key, result_key


#=========================================================================================
# DefaultCurrentImplementation
#=========================================================================================

class DefaultCurrentImplementation(VirtualCurrentImplementation):
    """This default implementation class will be used with a warning
    *only* if none of the other implementation matches
    and the function 'does_implement_this_current' of this class evaluated to true.
    This is typically useful for development when one has not completed
    the implementation of all currents but already wants to test a subset of them.
    """

    # Include these counterterms even if they are zero,
    # uncomment the following line to discard them instead
    # They must typically be discarded otherwise their lack of default mapping
    # will trigger a crash in test_IR_limits.
    # Coming up with a default mapping could be possible though so left as a todo
    # TODO: implement a default generic mapping
    is_zero = False
    
    def __init__(self, *args, **opts):

        super(DefaultCurrentImplementation, self).__init__(*args, **opts)
        self.supports_helicity_assignment = True

    @classmethod
    def does_implement_these_currents(cls, currents, model):
        """For production, it is preferable to turn off this default implementation
        by simply returning None below instead of {}."""
        
        # Simply uncomment the line below to de-activate this default implementation.
        # return None
        return {}
    
    def evaluate_subtraction_current(self, currents_block,
        higher_PS_point=None, momenta_dict=None, reduced_process=None,
        hel_config=None, xis=None, Q=None, **opts):
        """Simply return 0 for this current default implementation."""

        misc.sprint("WARNING: The default subtraction current implementation is used for the currents block:\n%s"%currents_block+
            "This will typically result in a crash of the evaluation of the integrand at run-time since no reduced kinematics are defined.\n"
            "Please specify an actual implementation for this current in the subtraction scheme.")

        return SubtractionCurrentResult.zero(current=currents_block, hel_config=hel_config)
