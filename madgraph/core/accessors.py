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
"""Classes for the implemention of accessors of Matrix elements and currents."""

import copy
import time
import logging
import sys
import importlib
import os
import shutil
import collections
import numpy as np
pjoin = os.path.join

import madgraph.core.base_objects as base_objects
import madgraph.integrator.phase_space_generators as PS_utils
from madgraph.interface.madevent_interface import MadLoopInitializer
import madgraph.various.misc as misc
import madgraph.core.subtraction as subtraction
import madgraph.interface.ME7_interface as ME7_interface
import madgraph.integrator.ME7_integrands as ME7_integrands
from madgraph import InvalidCmd, MadGraph5Error, DTYPE
from madgraph.iolibs.files import cp, ln, mv

logger = logging.getLogger('madevent7')

##########################################################################################
# ProcessKey, to be used as keys in the MEAccessorDict and the processes_map of contributions
##########################################################################################
class ProcessKey(object):
    """ We store here the relevant information for accessing a particular ME."""
    
    cache_active = False
    
    def get_key_for_cache(self,
                 process,  PDGs, sort_PDGs, vetoed_attributes, allowed_attributes, opts):
        """ Generates the look-up key for the dynamically created dictionary 'process_key_cache'
        of the objects for which a ProcessKey is generated."""
        
        key = [ tuple(PDGs), 
                sort_PDGs, 
                tuple(vetoed_attributes),
                tuple(allowed_attributes) if allowed_attributes else None ]
        if opts:
            key.extend((k,opts[k]) for k in sorted(opts.keys()))
        return tuple(key)

    def __init__(self, 
                # The user can initialize an ProcessKey directly from a process, or leave it empty
                # in which case the default argument of a base_objects.Process() instance will be considered.      
                process=None,
                # Instead, or on top of the above, one can decide to overwrite what set of PDGs this key 
                # will it correspond to
                PDGs = [],
                # Decide whether to store sorted PDGs or the original order provided by the process or other
                # intputs from this __init__
                sort_PDGs = True,
                # Exclude certain process attributes when creating the key for a certain process.
                vetoed_attributes = ['model','legs','uid','has_mirror_process','legs_with_decays'],
                # Specify only selected attributes to end up in the ProcessKey.
                # If None, this filter is deactivated.
                allowed_attributes = None,
                # Finally the user can overwrite the relevant attributes of the process passed in argument
                # if he so wants. Remember however that only attributes with *hashable* values are allowed
                # to be specified here.
                **opts):
        

        # Try to use the dynamically created caching system for the ProcessKey
        # which is a bit time consuming to generated
        if self.cache_active and (not process is None) and hasattr(process, 'process_key_cache'):

            key = self.get_key_for_cache(process,  PDGs, sort_PDGs, 
                                               vetoed_attributes, allowed_attributes, opts)
            if key in process.process_key_cache:
                cached_process_key = process.process_key_cache[key]  
                self.key_dict = cached_process_key.key_dict
                self.canonical_key = cached_process_key.canonical_key
                return
            
        # Initialize a dictionary which will be used to form the final tuple encoding all the information for 
        # this particular entry
        self.key_dict = {}
        
        # Cache of the key_dict converted into a canonical key. It will be set upon
        # calling 'get_canonical_key' for the first time.
        self.canonical_key = None
        
        # PDGs
        if (allowed_attributes is None or 'PDGs' in allowed_attributes) and (not 'PDGs' in vetoed_attributes):
            if PDGs:
                self.key_dict['PDGs'] = PDGs
            elif 'legs' in opts:
                self.key_dict['PDGs'] = ( tuple([l.get('id') for l in opts['legs'] if not l['state']]), 
                                          tuple([l.get('id') for l in opts['legs'] if l['state']]) )
            elif process:
                self.key_dict['PDGs'] = ( tuple(process.get_initial_ids()), 
                                          tuple(process.get_final_ids_after_decay()) )
            else:
                raise MadGraph5Error("When creating an ProcessKey, it is mandatory to specify the PDGs, either"+
                                     " via the options 'PDGs', 'legs' or 'process' (with precedence in this order).")
            
            if sort_PDGs:
                self.key_dict['PDGs'] = ( tuple(sorted(self.key_dict['PDGs'][0])), tuple(sorted(self.key_dict['PDGs'][1])) )
        
        # Now make sure to instantiate a default process if the user didn't select one
        if not process:
            process = base_objects.Process()

        # And dynamically check if one can indeed create a naive hash for the dict and list in attribute of Process.
        def hash_dict(a_dict, name):
            if not all(isinstance(el, collections.Hashable) for el in a_dict.keys()) or \
               not all(isinstance(el, collections.Hashable) for el in a_dict.values()):
                raise MadGraph5Error("Found unhashable dictionary '%s' as attribute of Process in ProcessKey."%name+
                  "Either define its hash in an ad-hoc way or veto it from taking part in the key of MEAccessorDict.")
            return tuple(sorted(a_dict.items()))

        def hash_list(a_list, name):
            if not all(isinstance(el, collections.Hashable) for el in a_list):
                raise MadGraph5Error("Found unhashable list '%s' as attribute of Process in ProcessKey."%name+
                  " Either define its hash in an ad-hoc way or veto it from taking part in the key of MEAccessorDict.")
            return tuple(a_list)
            
        # Now loop over all attributes of the Process and also additional options that the user might have specified
        # and do not correspond to Process attribute values to be overwritten
        for proc_attr in set(process.keys()+opts.keys()):
            # Check the veto rules
            if proc_attr in vetoed_attributes:
                continue
            # Check the allowance rule
            if not allowed_attributes is None and not proc_attr in allowed_attributes:
                continue
            # Take the value from the user_provided options if given
            # This is because we want the option specifically specified by the user to have precedence over the process
            # attributes and also because the user might specify options which are *not* process attributes.
            if proc_attr in opts:
                value = opts[proc_attr]
            else:
                value = process.get(proc_attr)

            # Special treatment for some attributes:
            if proc_attr == 'legs_with_decays':
                self.key_dict['legs_with_decays'] = tuple((l.get('id'), l.get('number'), l.get('state')) for l in value)
                
            elif proc_attr == 'decay_chains':
                # Group all the hashes of the processes in decay_chains and store them here.
                # BUT BEWARE THAT THE PDGs in self.key_dict only refer to the core production process then.
                self.key_dict['decay_chains'] = tuple( ProcessKey(proc).get_canonical_key() for proc in value)
            
            elif proc_attr == 'required_s_channels':
                # Group all the hashes of the processes in decay_chains and store them here.
                # BUT BEWARE THAT THE PDGs in self.key_dict only refer to the core production process then.
                self.key_dict['required_s_channels'] = tuple( tuple(pdg for pdg in pdg_list) for pdg_list in value)
            
            elif proc_attr == 'singular_structure':
                self.key_dict['singular_structure'] = value.get_canonical_representation(track_leg_numbers=False)

            # Let us not worry about WEIGHTED orders that are added automatically added when doing process matching
            # Also ignore squared order constraints that are not == as those are added automatically to improve
            # diagram efficiency and should not be used anyway to discriminate different processes since they may
            # overlap and cannot be constructed unambiguously in reduced processes for counterterms
            elif proc_attr in ['orders']:
                self.key_dict[proc_attr] = hash_dict(dict((k, v) for k, v in value.items() 
                                                             if k!='WEIGHTED'), proc_attr)
            elif proc_attr in ['sqorders_types', 'squared_orders']:
                self.key_dict[proc_attr] = hash_dict(dict((k, v) for k, v in value.items() if \
                        ( k!='WEIGHTED' and ( isinstance(process, subtraction.Current) or
                                              (k not in process['sqorders_types']) or 
                                              process['sqorders_types'][k]=='==') )
                    ), proc_attr)
    
            elif proc_attr == 'parent_subtraction_leg':
                parent_subtraction_leg = process[proc_attr]
                self.key_dict['singular_structure'] = (parent_subtraction_leg.pdg, parent_subtraction_leg.state)

            elif proc_attr == 'split_orders':
                # The ordering does not matter for split orders
                self.key_dict[proc_attr] = hash_list(sorted(value), proc_attr)
    
            # Now generically treat other attributes
            elif isinstance(value, (int, str, tuple, bool)):
                self.key_dict[proc_attr] = value
            elif isinstance(value, list):
                self.key_dict[proc_attr] = hash_list(value, proc_attr)
            elif isinstance(value, dict):
                self.key_dict[proc_attr] = hash_dict(value, proc_attr)
            else:
                raise MadGraph5Error("The attribute '%s' to be considered as part of a key for the MEAccessorDict"%proc_attr
                 +" is not of a basic type or list / dictionary, but '%s'. Therefore consider vetoing this"%type(value)+
                 " attribute or adding a dedicated ad-hoc rule for it in ProcessKey.")

        # Assign a caching system automatically
        if self.cache_active and (not process is None):
            key = self.get_key_for_cache(process,  PDGs, sort_PDGs, 
                                               vetoed_attributes, allowed_attributes, opts)
            if hasattr(process, 'process_key_cache'):
                process.process_key_cache[key] = self
            else:
                process.process_key_cache = {key : self}

    def set(self, key, value):
        """ Modify an entry in the key_dict created."""
        if self.cache_active:
            raise MadGraph5Error("ProcessKeys instances cannot be modified while their dynamic caching is active.")

        if key not in self.key_dict:
            raise MadGraph5Error("Key '%s' was not found in the key_dict created in ProcessKey."%key)
        if not isinstance(value, (int, str, bool, tuple)):
            raise MadGraph5Error("Values for the key_dict created in ProcessKey should be of type (int, str, bool, tuple).")
        
        self.key_dict[key] = value
        # Force the canonical key to be recomputed
        if hasattr(self, 'canonical_key'):
            self.canonical_key = None

    def get_canonical_key(self, force=False):
        """ Simply uses self.key_dict to return a hashable canonical representation of this object."""
        if not force and self.canonical_key:
            return self.canonical_key
        
        self.canonical_key = tuple(sorted(self.key_dict.items()))        
        return self.canonical_key

################################################################################
# MEAccessor mother class
################################################################################
class VirtualMEAccessor(object):
    """ A class wrapping the access to one particular MatrixEleemnt."""
    
    # A class variables that decides whether caching is active or not
    cache_active = True
    
    def __new__(cls, *args, **opts):
        """ Factory for the various Accessors, depending on the process a different 
        class will be selected.
        The typical arguments passed are:
           process, f2py_module_path, slha_card_path
        """
        if cls is VirtualMEAccessor:
            # Use possibly customized CurrentAccessors for the access to subtraction currents
            if isinstance(args[0], subtraction.Current):
                target_type = 'CurrentAccessor'
            # Check if the input match what is needed for a PythonAccessor namely that the first
            # two are ProcessInstance, (f2py_module_path, f2py_module_name)
            elif len(args)<2 or not isinstance(args[0],base_objects.Process) or \
               not isinstance(args[1],tuple) or not len(args[1])==2 or not \
               all(isinstance(path, str) for path in args[1]):
                # Check if the first second argument is key declared in MEAccessor_classes_map.
                # This way, custom contributions can declare the type of MEAccessor they want to
                # instantiate.
                if len(args)>=2 and isinstance(args[1],base_objects.Process) and \
                   args[2] in MEAccessor_classes_map:
                    target_type = args[2]
                else:
                    target_type = 'Unknown'
            else:
                process = args[0]
                f2py_module_path = args[1]
                target_type = 'PythonAccessor'
            
            target_class = MEAccessor_classes_map[target_type]
            if not target_class:
                raise MadGraph5Error(
                    "Cannot find a class implementation for target MEAccessor type '%s'."%target_type)
            else:
                return target_class.__new__(target_class, *args, **opts)
        else:
            return super(VirtualMEAccessor, cls).__new__(cls, *args, **opts)

    def __init__(self, process, mapped_pdgs = [], root_path='', **opts):
        
        """ Initialize an MEAccessor from a process, allowing to possibly overwrites n_loops.
        and mapped_pdgs to indicate other list of pdgs corresponding to processes that can be 
        obtained by this same one.
        root_path is the process output root path from which all other paths are relative."""
        
        # Save the initialization arguments and options to facilitate the generation and use of
        # a the dump. root_path should not be saved since this will be updated whenever the
        # instance is reconstructed from a dump.
        self.initialization_inputs = {'args':[], 'opts':{}}
        self.initialization_inputs['args'].append(process)
        self.initialization_inputs['opts'].update(
            {'mapped_pdgs': mapped_pdgs,
             'root_path': root_path})
        self.initialization_inputs['opts'].update(opts)
                
        # Define the attributes for this particular ME.
        self.pdgs_combs  = [ ( tuple(process.get_initial_ids()), tuple(process.get_final_ids()) ) ] + mapped_pdgs

        self.ME_dict_key = ProcessKey(process)
        
        if not os.path.isabs(root_path):
            raise MadGraph5Error("The initialization of VirtualMEAccessor necessitates an "+
                                                        "absolute path for root_path.")
        self.root_path = root_path
        
        # Finally instantiate a cache for this MEaccessor.
        self.cache = MEAccessorCache()

    @classmethod
    def activate_cache(cls):
        """ Activates the caching of results during __call___ """
        cls.cache_active = True
        
    @classmethod
    def deactivate_cache(cls):
        """ Activates the caching of results during __call___ """
        cls.cache_active = False

    def generate_dump(self, **opts):
        """ Generate a serializable dump of self, which can later be used, along with some more 
        information, in initialize_from_dump in order to regenerate the object. 
        This can be overloaded by the daughter classes."""
        
        dump = {}
        # We opt here for a dictionary of relevant attribute to be used in __init__ during
        # as well as the class of the accessor
        dump['class'] = self.__class__
        dump['args'] = self.initialization_inputs['args']
        dump['opts'] = self.initialization_inputs['opts']
        
        return dump

    @classmethod
    def initialize_from_dump(cls, dump, root_path, *args, **opts):
        """ Regenerate this instance from its dump and possibly other information."""
        
        MEAccessorClass = dump['class']
        assert (cls==MEAccessorClass)
        MEAccessor_instance = MEAccessorClass(*dump['args'], **dump['opts'])
        # Make sure to override the root_path with the newly provided one
        MEAccessor_instance.root_path = root_path        
        return MEAccessor_instance

    def get_canonical_key_value_pairs(self):
        """ Return all the canonical keys that should point to this ME. 
        This is intended to be used for the MEAccessorDict."""
        
        key_value_pairs = []
        for pdgs_comb in self.pdgs_combs:
            sorted_pdgs = (tuple(sorted(pdgs_comb[0])), tuple(sorted(pdgs_comb[1])))
            self.ME_dict_key.set('PDGs',sorted_pdgs)
            # Store both this accessor and the original ordering of the pdgs_comb
            # which is lost in the key since it has been sorted there.
            value = (self, pdgs_comb)
            key_value_pairs.append( (self.ME_dict_key.get_canonical_key(), value) )

        return key_value_pairs
    
    def __call__(self, PS_point, permutation=None, defining_pdgs_order=None, squared_orders = None, **opts):
        """ The evaluation will be implemented by the daughter classes, but we can already here apply 
        the specified permutation and record the defining_pdgs_order (i.e. which of the mapped process
        was the user really interested in, it is one of those present in self.pdgs_combs) 
        in a dictionary that we can return to the daughter class."""

        permuted_PS_point, permuted_opts = self.apply_permutations(permutation, PS_point=PS_point, **opts)
        if isinstance(squared_orders, dict):
            permuted_opts['squared_orders'] = tuple(sorted(squared_orders.items()))
        else:
            permuted_opts['squared_orders'] = squared_orders
        permuted_opts['defining_pdgs_order'] = defining_pdgs_order

        # Return the inputs properly treated
        return permuted_PS_point, permuted_opts

    def synchronize(self, **opts):
        """ Synchronizes this accessor with the possibly updated value of parameter cards and ME source code.
        Must be defined by daughter classes."""
        
        raise NotImplemented

    @classmethod
    def apply_permutations(cls, permutation, PS_point = [],
                                 spin_correlation  = [], 
                                 color_correlation = [], 
                                 hel_config = tuple([]),
                                 **opts):
        
        """ Processes the options when calling the evaluation of the matrix_element and apply the permutation
        to all options that need to be permuted. Additional options provided via opts will simply be propagated,
        but not permuted. The format for the three options above is:
        
        spin_correlation :
        -> Spin correlation defined as a list of 2-tuple of the following form
                    
           [ (leg_IDA, (vec_A1, vec_A2, vec_A3,...)), 
             (leg_IDB, (vec_B1, vec_B2, vec_B3,..)), 
             (leg_IDC, (vec_C1, vec_C2, vec_C3,...)), ... ]
            
            So basically this means replace the list of possible helcity polarization vectors for each 
            leg_IDA, leg_IDB, leg_IDC, etc.. ... (these are integer) with the corresponding list of 4-vectors,
            (vec_A1, vec_A2, vec_A3,...), (vec_B1, vec_B2, vec_B3,..), etc... where vec_i are simply 4-tuples of floats.

#        -> Color connections are specified by providing what SU(3) generator chain to apply in 
#          between the amplitude product. The format is a list of 2-tuple of the form ('generator_name', (indices)).
#          For instance
#             [ ('T',(-100,3,-1)), ('T',(-200,-1,-2)), ('T',(-100,-2,-3)), ('T',(-200,-3,4)), ('f',(7,-300,-400)), ('f',(-400,-300,8)), ...]
#          Indicates:
#            T^a_(3 i) T^b_(i j) T^a_(j k) T^b_(k 4) f^(7, c, d) f^(d, c, 8) 

        -> Color connections are specified at NLO by a simple list of 2-tuple, containing the leg number of each of the
           two legs color-connected.
           For instance
             [ (2, 4), (5, 1), etc...]
           Indicates the user wants
             <M| T2 T4 |M>, <M| T5 T1 |M>, etc...
           With T == t^a_{ij} for quarks and T = i f^{abc} for gluons.
           /!\ : This notation will require an extension for NNLO and beyond.
        
        -> Helicity configuration to be considered, this can be a tuple like
            (-1, -1, +1 ,-1, +1)
          For the helicity configuration --+-+
          If a spin-connection is specified for a specific particle, the helicity specification will be dropped.
        """
        
        # No permutation, just return local copies of the arguments
        if permutation is None:
            all_opts = {'spin_correlation': tuple(spin_correlation) ,
                        'color_correlation': tuple(color_correlation),
                        'hel_config': tuple(hel_config)}
            all_opts.update(opts)
            return PS_point, all_opts
        
        
        # Apply the permutations while retaining a canonical representation for each attribute

        # The permutation is how to go *from* the user-provided order to the order assumed by the underlying ME.
        # So When performing the permutation like below, we must use the reversed_permutation.
        reversed_permutation = dict((v,k) for (k,v) in permutation.items())
        permuted_PS_point = [PS_point[reversed_permutation[i]] for i in range(len(PS_point))] if len(PS_point)>0 else None

        permuted_spin_correlation = tuple(sorted([ (permutation[leg_ID-1]+1, vectors) for leg_ID, vectors 
                                in spin_correlation ], key=lambda el:el[0] )) if spin_correlation else None
        
        permuted_color_correlation = tuple([ tuple(sorted([permutation[ind1-1]+1, permutation[ind2-1]+1])) 
                                    for (ind1, ind2) in color_correlation ]) if color_correlation else None
        
        permuted_hel_config = tuple( hel_config[permutation[i]] for i in range(len(hel_config)) ) \
                                                                          if hel_config else None

        all_opts = {'spin_correlation' : permuted_spin_correlation,
                    'color_correlation' : permuted_color_correlation,
                    'hel_config' : permuted_hel_config}
        all_opts.update(opts)

        return permuted_PS_point, all_opts

class MEEvaluation(dict):
    """ A very basic class to store the various output of an ME evaluation output."""
    
    # All residues of the poles for the epsilon expansion can be specified with a
    # key names 'eps^n'.
    
    result_order  = ['tree','finite','return_code','accuracy','eps']    
    result_format = {'tree':'%.15e','finite':'%.15e','eps':'%.15e',
                     'return_code':'%d','accuracy': '%.2g'}
    
    def __init__(self, *args, **opts):
        super(MEEvaluation, self).__init__(*args, **opts)      

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
        """ Formats nicely the output of a particular ME evaluation."""
        
        res = []
        template = '%%-%ds = %%s'%(max_len_attribute if max_len_attribute>0 else
                                                    max(len(k) for k in self.keys()))
        sorted_result_keys = sorted(self.keys(), key=lambda el: self.get_result_order(el))

        def format_result(key, res):
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

class MEResult(dict):
    """ A class to store the different results of MatrixElement call for one specific PS point / scale."""

    def __init(self, *args, **opts):
        super(MEResult, self).__init__(*args, **opts)        
        
    def nice_string(self):
        """ Print out all the results available in a nice form."""
        # First lists the result with the least amount of attributes specified
        sorted_keys = sorted(self.keys(), key=lambda k:[el[1] for el in k].count(None))

        res = []
        max_len_attribute = max(max(len(el[0]) for el in k) for k in sorted_keys)
        max_len_attribute = max(max_len_attribute, MEEvaluation.get_max_length_attribute())
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
        """ Attempts to recycle the result from previous computations.
        Returns None otherwise. Opts are typically:
           helicities: A tuple of integer, None means summed over.
           squared_orders: a tuple repr. of the usual dict., None means all contributions.
           color_correlation: same convention as described in function apply_permutation.
           spin_correlation: same convention as described in function apply_permutation.
        """
        key_opts = {'hel_config'           : None,
                    'squared_orders'       : None,
                    'color_correlation'    : None,
                    'spin_correlation'     : None}
        key_opts.update(opts)
        
        try:
            result = self[tuple(sorted(opts.items()))]
            return result
        except KeyError:
            return None
    
    def add_result(self, value, **opts):
        """ Add a result to the current record.
        value is typically a dictionary with at leas the
        key 'finite'."""
        assert(isinstance(value, dict) and 'finite' in value), "Dictionaries with at least the key 'finite' "+\
                                                 "should be added as MatrixElement results, not %s."%str(value)
        key_opts = {'hel_config'           : None,
                    'squared_orders'       : None,
                    'color_correlation'    : None,
                    'spin_correlation'     : None}
        key_opts.update(opts)
        
        self[tuple(sorted(opts.items()))] = value

    def get_inverse_permuted_copy(self, permutation=None):
        """ Apply the inverse permutation of the one specified to a copy of self and return."""
        
        # Nothing to do in this case, then just add a copy of self
        if permutation is None:
            return MEResult(self)
        
        inverse_perm = dict((v,k) for k,v in permutation.items())
        returned_copy = MEResult()
        for attributes, value in self.items():
            attr = dict(attributes)
            # Make sure to temporarily redefine the color_correlators as list of color_correlator
            # since this is what apply_permutations expects
            if attr['color_correlation']:
                attr['color_correlation'] = [attr['color_correlation'],]
            _, inverse_permuted_attributes = VirtualMEAccessor.apply_permutations(inverse_perm, **attr)
            # Now recast the color_correlation as a single element:
            if inverse_permuted_attributes['color_correlation']:
                inverse_permuted_attributes['color_correlation']=inverse_permuted_attributes['color_correlation'][0]
            returned_copy[tuple(sorted(inverse_permuted_attributes.items()))] = value
        
        return returned_copy
            
class MEAccessorCache(dict):
    """ A simple class that implements a caching mechanism for Matrix Element call and results."""

    def __init__(self, min_n_entries=1, max_n_entries=1, max_cache_size=100000):
        """ One can limit the cache size with a combination of number of elements it contains
        and its size in bytes. The min_n_entries option allows to set a minimum numebr of elements
        to keep, irrespectively of the max_cache_size value."""
        self.NO_RESULT = MEResult()
        # List of entries added, keeping their order.
        self.added_entries = []
        if min_n_entries < 0:
            raise MadGraph5Error('The option min_n_entries of a MEAccessorCache should be positive.')
        self.min_n_entries = min_n_entries
        self.max_n_entries = max_n_entries
        self.max_cache_size = max_cache_size        
        if max_n_entries < 0 and max_cache_size < 0:
            raise MadGraph5Error('At least max_n_entries or max_cache_size should be definite positive'+
                                              ' so as to keep the MEAccessorCache cache size in check.')

    def get_result(self, **opts):
        """ Attempts to recycle a computation performed with opts as input. Returns an empty MEResult if
        not found. Opts are typically:
           PS_point : A list of Lorentz5DVectors
           alpha_s  : alpha_s value (float)
           mu_r     : renormalization_scale (float)
        """
        key_opts = {'PS_point' : None,
                    'alpha_s'  : None,
                    'mu_r'     : None}
        key_opts.update(opts)

        try:
            ME_result = self[tuple(sorted(key_opts.items()))]
            return ME_result
        except KeyError:
            return self.NO_RESULT

    def is_cache_too_large(self):
        """ Checks whether the current cache size is too large."""
        
        n_items = len(self)
        if n_items<=self.min_n_entries:
            return False
        
        if self.max_n_entries > 0 and n_items > self.max_n_entries:
            return True
        
        if self.max_cache_size > 0 and sys.getsizeof(self) > self.max_cache_size:
            return True
        
        return False
        
    def check_cache_size(self):
        """ Controls the size of the cache and possibly remove entries if too large."""
        
        while self.is_cache_too_large():
            del self[self.added_entries.pop(0)]
        
    def add_result(self, result=None, **opts):
        """ Add a result to the cache, making sure it does not extend its target maximum size.
        """
        key_opts = {'PS_point' : None,
                    'alpha_s'  : None,
                    'mu_r'     : None}
        key_opts.update(opts)
        key = tuple(sorted(key_opts.items()))
        
        if key in self:
            if result:
                self[key].update(result)
            return self[key]

        if result is None:
            new_result = MEResult()
        else:
            new_result = result
            
        self[key] = new_result
        self.added_entries.append(key)
        
        # Now check Cache size:
        self.check_cache_size()
        
        return self[key]

class SubtractionCurrentAccessorCache(MEAccessorCache):
    # We can reuse completely the MEAccessorCache for the purpose of the SubtractionCurrentAccessorCache.
    pass

class SubtractionCurrentAccessor(VirtualMEAccessor):
    """ A class wrapping the access to a particular set of mapped subtraction currents."""

    cache_active = False

    def __init__(self, defining_current,
                       relative_module_path,
                       current_implementation_class_name,
                       relative_generic_module_path, 
                       instantiation_options, 
                       mapped_process_keys=[], 
                       root_path='', 
                       model=None,
                       **opts):
        """ Initialize a SubtractionCurrentAccessor from a current and mapped process keys to indicate 
        other currents that can be obtained by this same accessor.
        root_path is the parent directory of the SubtractionCurrent directory, from which all
        module paths are relative."""
        
        # Save the initialization arguments and options to facilitate the generation and use of
        # a the dump. root_path should not be saved since this will be updated whenever the
        # instance is reconstructed from a dump.
        self.initialization_inputs = {'args':[], 'opts':{}}
        self.initialization_inputs['args'].append(defining_current)
        self.initialization_inputs['args'].append(relative_module_path)
        self.initialization_inputs['args'].append(current_implementation_class_name)
        self.initialization_inputs['args'].append(relative_generic_module_path)
        self.initialization_inputs['args'].append(instantiation_options)
        self.initialization_inputs['opts'].update(
            {'mapped_process_keys': mapped_process_keys,
             'root_path': root_path})
        self.initialization_inputs['opts'].update(opts)
                
        # Now define the attributes for this particular subtraction current accessor.
        self.defining_current = defining_current
        self.mapped_current_keys  = mapped_process_keys
        self.relative_module_path = relative_module_path
        self.relative_generic_module_path = relative_generic_module_path
        self.current_implementation_class_name = current_implementation_class_name
        self.instantiation_options = instantiation_options
        
        if not os.path.isabs(root_path):
            raise MadGraph5Error("The initialization of SubtractionCurrentAccessor necessitates an "+
                                                        "absolute path for root_path.")
        self.root_path = root_path
        
        # Finally instantiate a cache for this SubtractionCurrentAccessor. 
        # (Here this class basically reuses the MEAccessorCache one)
        self.cache = SubtractionCurrentAccessorCache()

        # Now load the modules and specify the result classes        
        self.module = self.load_module(self.relative_module_path)
        self.generic_module = self.load_module(self.relative_generic_module_path)
        self.evaluation_class = self.generic_module.SubtractionCurrentEvaluation
        self.result_class = self.generic_module.SubtractionCurrentResult
        
        # Instantiate the subtraction_current_instance with the specified model
        self.model = model
        self.subtraction_current_instance = None
        self.synchronize(model=model) 

    def synchronize(self, model=None, **opts):
        """ Synchronizes this accessor with the possibly updated model and value."""
        if not model is None:
            self.model = model
            self.subtraction_current_instance = getattr(self.module, 
                self.current_implementation_class_name)(model, **self.instantiation_options)

    def generate_dump(self, **opts):
        """ Generate a serializable dump of self, which can later be used, along with some more 
        information, in initialize_from_dump in order to regenerate the object. 
        This can be overloaded by the daughter classes."""
        
        dump = {}
        # We opt here for a dictionary of relevant attribute to be used in __init__ during
        # as well as the class of the accessor
        dump['class'] = self.__class__
        dump['args'] = self.initialization_inputs['args']
        dump['opts'] = self.initialization_inputs['opts']
        
        return dump

    def compile(self, *args, **opts):
        """ For basic python subtraction current there is not compilation necessary for now."""
        pass
    
    @classmethod
    def initialize_from_dump(cls, dump, root_path, model, *args, **opts):
        """ Regenerate this instance from its dump and possibly other information."""
        
        SubtractionCurrentAccessorClass = dump['class']
        assert (cls==SubtractionCurrentAccessorClass)
        # Make sure to override the root_path with the newly provided one
        if 'root_path' in dump['opts']:
            dump['opts']['root_path'] = root_path
        # And also add the current active model to the constructor options
        SubtractionCurrentAccessor_instance = SubtractionCurrentAccessorClass(*dump['args'], model=model, **dump['opts'])
        return SubtractionCurrentAccessor_instance

    def nice_string(self):
        """ Summary of the details of this Subtraction current accessor."""
        res = []
        res.append("%s: %s @ '%s'"%(self.__class__.__name__,self.relative_module_path, self.current_implementation_class_name))
        res.append('Defining subtraction current: %s'%str(self.defining_current))
        return '\n'.join(res)

    def get_canonical_key_value_pairs(self):
        """ Return all the canonical keys that should point to this ME. 
        This is intended to be used for the MEAccessorDict."""
        
        key_value_pairs = []
        for mapped_current_key in self.mapped_current_keys:
            # There is no remapping of the inputs to be done for currents, so we just store
            # None as a second entry to the value in the accessor dictionary.
            value = (self, None)
            key_value_pairs.append( (mapped_current_key.get_canonical_key(), value) )

        return key_value_pairs
    
    def check_inputs_validity(self, opts, current):
        """ Check the validity of the inputs of the call to this current accessor."""
        
        # Work on a copy of th eoption dictionary, filtered from irrelevant keys
        new_opts = dict((k,v) for k,v in opts.items() if k not in 
                                                            ['permutation','process_pdgs'])
        
        if 'hel_config' in opts and opts['hel_config']:
            new_opts['hel_config'] = tuple(opts['hel_config'])
            # In this case it is not optimal to check the validity of the helicity configuration, so 
            # we limit ourselves to checking if it supports helicity assignment
            if not self.subtraction_current_instance.supports_helicity_assignment:
                raise MadGraph5Error("The following subtraction current accessor:\n%s"%(
                                    self.nice_string() ) + "\ndoes not support helicity assignment.")

        squared_orders = None
        if 'squared_orders' in opts:
            new_opts.pop('squared_orders')
            if isinstance(opts['squared_orders'], dict):
                squared_orders = tuple(sorted(opts['squared_orders'].items()))
            elif isinstance(opts['squared_orders'], list):
                squared_orders = tuple(sorted(opts['squared_orders']))
            if squared_orders:
                # This information must be passed via the 'squared_orders' attribute of the current, so we
                # simply make sure that it is identical if specified and then remove it
                if squared_orders != tuple(sorted(current.get('squared_orders').items())):
                    raise MadGraph5Error("The following subtraction current accessor:"+\
                        "\n%s\ncannot provide squared orders %s."%(
                            self.nice_string(), str(squared_orders)))

        return new_opts

    def call_subtraction_current(self, *args, **opts):
        """ Wrapper around the actual call of the subtraction currents, so as to be
        able to easily time it with a profiler."""
        
        return self.subtraction_current_instance.evaluate_subtraction_current(*args, **opts)

    def __call__(self, current, PS_point, **opts):
        """ Evaluation of the subtraction current. """
        
        if self.subtraction_current_instance is None:
            raise MadGraph5Error("This subtraction current accessor\n'%s'\nhas not been properly initialized."%(
                                                                                            self.nice_string()))
        
        # Parse options and check their validity
        call_opts = self.check_inputs_validity(opts, current)
    
        # Set the arguments of the call
        call_args = [current, PS_point]
        
        is_cache_active = opts.get('cache_active', self.cache_active)
        if is_cache_active:
            # Now obtain the cache key directly from the current implementation
            cache_key, result_key = self.subtraction_current_instance.get_cache_and_result_key(
                                                                   *call_args, **call_opts)
        else:
            cache_key, result_key = None, None
    
        # Attempt to recycle the result
        if not cache_key is None:
            recycled_call = self.cache.get_result(**cache_key)
            recycled_result = recycled_call.get_result(**result_key)
            if recycled_result:
                return self.evaluation_class(recycled_result), self.result_class(recycled_call)

        all_evaluations = self.call_subtraction_current(*call_args, **call_opts)
        
        if len(all_evaluations)==1:
            return all_evaluations.values()[0], all_evaluations
        
        # If there are several evaluations we need the result_key even in the absence
        # of caching, in which case it must be recomputed here.
        if result_key is None:
            cache_key, result_key = self.subtraction_current_instance.get_cache_and_result_key(
                                                                   *call_args, **call_opts)
        
        evaluation_asked_for = all_evaluations.get_result(**result_key)
        if not evaluation_asked_for:
            raise MadGraph5Error("Could not obtain result '%s' from evaluation:\n%s"%(str(result_key), str(all_evaluations)))
        
        # Update the cache with the new results produced
        if not cache_key is None:
            all_evaluations = self.cache.add_result(all_evaluations, **cache_key)
        
        # Return both the specific evaluation asked for and all results available for this cache_key
        return self.evaluation_class(evaluation_asked_for), self.result_class(all_evaluations)

    def load_module(self, module_path):
        """ Loads a particular subtraction module given its path"""
        
        # Make sure to temporarily adjust the environment
        added_path = False
        if self.root_path not in sys.path:
            added_path = True
            sys.path.insert(0, self.root_path)
        try:
            subtraction_current_module = importlib.import_module(module_path)
        except ImportError as e:
            raise MadGraph5Error("Could not load subtraction current module '%s' in '%s'."%(module_path,self.root_path)+
                                 "\nThe import error is: %s"%e.message)

        if added_path:
            sys.path.pop(sys.path.index(self.root_path))

        return subtraction_current_module

    @classmethod
    def apply_permutations(cls, *args, **opts):
        """ apply_permutations should never be called in the context of currents because there is not remapping
        of the inputs to be done. This is performed upstream."""
        raise MadGraph5Error("Function 'apply_permutations' should never have been called in the context of "+
                             "subtraction currents.")

class F2PYMEAccessor(VirtualMEAccessor):
    """ A class wrapping the access to one particular MatrixEleemnt wrapped with F2PY """
    
    cache_active = False
    
    def __new__(cls, process, f2py_module_path, slha_card_path, **opts):
        """ Factory for the various F2PY accessors, depending on the process a different class
        will be selected."""

        if cls is F2PYMEAccessor:
            target_class = None
            n_loops = process.get('n_loops')
            if process.get('perturbation_couplings') and n_loops==1:
                target_class = F2PYMEAccessorMadLoop
            elif n_loops==0:
                target_class = F2PYMEAccessor
            else:
                raise MadGraph5Error("Could not determine the type of F2PYMEAccessor suited for "+
                                     " %s with %d loops."%(process.nice_string(), n_loops))

            return super(F2PYMEAccessor, cls).__new__(target_class, 
                                                process, f2py_module_path, slha_card_path, **opts)
        else:
            return super(F2PYMEAccessor, cls).__new__(cls, 
                                                process, f2py_module_path, slha_card_path, **opts)
    
    def __init__(self, process, f2py_module_path, slha_card_path, **opts):
        """ Initialize a F2PYMEAccessor with the path to the f2py_module.
        The path f2py_module_path is a 2-tuple of the form 
           (absolute_path_of_module_location, module_python_path) 
        for instance:
           ('/usr/john/Documents/MyProcOutput', 'LO_udx_epve_1.matrix_1_udx_epve_py') 
        """
        
        compile_if_necessary = True
        if 'compile_if_necessary' in opts:
            compile_if_necessary = opts.pop('compile_if_necessary')

        super(F2PYMEAccessor, self).__init__(process, **opts)
        
        # Add the additional inputs to back up to be used later for the dump
        self.initialization_inputs['args'].extend([f2py_module_path, slha_card_path])
        
        # Save the location of the SLHA card path
        if not os.path.isabs(slha_card_path):
            raise MadGraph5Error("The initialization of F2PYMEAccessor necessitates an "+
                                                        "absolute path for slha_card_path.")
        self.slha_card_path = os.path.relpath(slha_card_path, self.root_path)

        self.f2py_module_path = (os.path.relpath(f2py_module_path[0],self.root_path),f2py_module_path[1])
        
        # Process directory and name
        self.proc_dir = pjoin(*(self.f2py_module_path[1].split('.')[:-1]))
        name = self.f2py_module_path[1].split('.')[-1]
        assert(name.startswith('matrix_') and name.endswith('_py')), "Non-standard f2py'ed so library name: %s"%name
        self.proc_name = name[7:-3]
        
        shared_library_path = pjoin(f2py_module_path[0],*(f2py_module_path[1].split('.')))+'.so'
        if not os.path.isfile(shared_library_path):
            if compile_if_necessary:
                self.compile()
            else:
                # We are at the generation stage and the f2py module has never been compiled yet.
                # This accessor is therefore only a placeholder for now.
                self.f2py_module = None
                return

        self.f2py_module = self.load_f2py_module(f2py_module_path)

        # Try to guess the process prefix if not defined
        possible_proc_prefix_path = pjoin(self.root_path,self.f2py_module_path[0],'proc_prefix.txt')
        if 'proc_prefix' in opts:
            self.proc_prefix = opts['proc_prefix']
        elif os.path.isfile(possible_proc_prefix_path):
            self.proc_prefix = open(possible_proc_prefix_path).read()
        elif hasattr(self.f2py_module, 'smatrix'):
            self.proc_prefix = ''
        else:
            candidates = list(set(
                [attribute[:-7] for attribute in dir(self.f2py_module) if attribute.endswith('smatrix')]+\
                [attribute[:-14] for attribute in dir(self.f2py_module) if attribute.endswith('sloopmatrixhel')]))
            if len(candidates)>1:
                raise MadGraph5Error("Cannot automatically detect process prefix in f2py module %s @ '%s'."%
                                     (self.f2py_module_path[1], self.f2py_module_path[0])+
                                     "\n. Possible options are: '%s'."%str(candidates))
            self.proc_prefix = candidates[0]
        
        # Sanity check
        if not self.has_function('loopme_accessor_hook') and not self.has_function('me_accessor_hook'):
            raise MadGraph5Error("The specified f2pymodule %s @ '%s' , with proc_prefix = '%s'"%
                (self.f2py_module_path[1], self.f2py_module_path[0], self.proc_prefix)+
                " does not seem to define the subroutine 'me_accessor_hook' or 'loopme_accessor_hook'.\n"+
                "Check the sanity of the proc_prefix value.")
            
        self.synchronize(from_init=True)
    
    def compile(self,mode='auto'):
        """ Compiles the source code associated with this MatrixElement accessor."""
        
        root_dir = pjoin(self.root_path, self.proc_dir)
        source_dir = pjoin(root_dir, 'Source')
        Pdir = pjoin(root_dir, 'SubProcesses','P%s'%self.proc_name)
        if not os.path.isdir(Pdir):
            raise InvalidCmd("The expected subprocess directory %s could not be found."%Pdir)
        
        if os.path.isfile(pjoin(Pdir, 'matrix_%s_py.so'%self.proc_name)) and mode=='never':
            return

        if mode=='always':
             misc.compile(arg=['clean'], cwd=Pdir)
             misc.compile(arg=['clean'], cwd=pjoin(source_dir,'DHELAS'))
             
        if not os.path.isfile(pjoin(Pdir, 'matrix_%s_py.so'%self.proc_name)):
            logger.info("Compiling directory %s ..."%(pjoin(self.proc_dir, 
                                                     'SubProcesses','P%s'%self.proc_name)))                    
            misc.compile(arg=['../lib/libdhelas.a',], cwd=source_dir)
            misc.compile(arg=['../lib/libmodel.a',], cwd=source_dir)
            misc.compile(arg=['matrix_%s_py.so'%self.proc_name, 'MENUM=_%s_'%self.proc_name], cwd=Pdir)
        
        if not os.path.isfile(pjoin(Pdir, 'matrix_%s_py.so'%self.proc_name)):
            raise InvalidCmd("The f2py compilation of SubProcess '%s' failed.\n"%Pdir+
                "Try running 'MENUM=_%s_ make matrix_%s_py.so' by hand in this directory."\
                                                          %(self.proc_name,self.proc_name))

        if not os.path.isfile(pjoin(root_dir, 'matrix_%s_py.so'%self.proc_name)):
            # Refresh the soft link
            ln( pjoin(Pdir, 'matrix_%s_py.so'%self.proc_name), starting_dir = root_dir )

    def synchronize(self, ME7_options = None, from_init=False, compile='auto', **opts):
        """ Synchronizes this accessor with the possibly updated value of parameter cards and ME source code.
        Must be defined by daughter classes."""

        # Reset the initialization to False
        self.module_initialized = False

        # Recompile and reload the module if synchronize was not call from the __init__ function
        if not from_init:
            self.compile(mode=compile)
            self.f2py_module = reload(self.f2py_module)

        # Now gather various properties about the Matrix Elements
        # Get all helicities orderin
        self.helicity_configurations = dict((tuple(hels),i+1) for i, hels in enumerate(
                                              self.get_function('get_helicity_definitions')() ))
        # Reversed map
        self.id_to_helicity = dict((value,key) for (key,value) in self.helicity_configurations.items())
        
        # Available spin correlations
        if not self.has_function('get_max_n_spin_corr_legs'):
            self.spin_correlations = None
            self.max_spin_corr_vectors = None
        else:
            self.spin_correlations = 'N'*self.get_function('get_max_n_spin_corr_legs')()+'LO'
            self.max_spin_corr_vectors = self.get_function('get_max_n_spin_corr_vectors')()
        
        # Available color correlations
        if not self.has_function('get_n_color_correlators'):
            self.color_correlations = None
            self.id_to_color_correlation = None
        else:
            # Initialize a map of the color connectors
            self.color_correlations = dict( 
                (tuple(sorted(self.get_function('get_color_correlator_for_id')(i_correlator))),i_correlator)
                              for i_correlator in range(1,self.get_function('get_n_color_correlators')()+1))
            # Reversed map
            self.id_to_color_correlation = dict((value,key) for (key,value) in self.color_correlations.items())

        # Available squared orders
        # 'None' means summed over.
        self.squared_orders = {None: 0}
        self.split_order_names = []
        if self.has_function('get_nsqso_born'):
            n_squared_orders = self.get_function('get_nsqso_born')()
            self.split_order_names = [''.join(name).strip() for name in self.get_function('get_split_order_names')()]
            for i_sq_order in range(1,n_squared_orders+1):
                key = tuple(sorted([(self.split_order_names[i], value) for i, value in enumerate(
                                            self.get_function('get_squared_orders_for_soindex')(i_sq_order) )]))
                self.squared_orders[key] = i_sq_order
        
        # Reversed map
        self.id_to_squared_order = dict((value,key) for (key,value) in self.squared_orders.items())

    def nice_string(self):
        """ Summary of the details of this ME accessor."""
        res = []
        res.append("%s: %s @ '%s'"%(self.__class__.__name__,self.f2py_module_path[1], self.f2py_module_path[0]))
        res.append('%-40s:   %s'%('Number of hel. configs',len(self.helicity_configurations)))
        res.append('%-40s:   %s'%('Spin correlators available','None' if self.spin_correlations is None \
                                                                                        else self.spin_correlations))
        res.append('%-40s:   %s'%('Number of color correlators available','None' if self.color_correlations is None \
                                                                                        else len(self.color_correlations)))
        res.append('%-40s:   %s'%('Squared orders available', 'Only summed' if self.squared_orders.keys()==[None] else \
                                                                str([k for k in self.squared_orders.keys() if not k is None])))
        return '\n'.join(res)

    def has_function(self, function_name):
        """ Simply tests if the current f2py module has the desired function."""
        
        return hasattr(self.f2py_module,'%s%s'%(self.proc_prefix, function_name))

    def get_function(self, function_name):
        """ Returns the specified f2py function."""
        
        if not self.has_function(function_name):
            raise MadGraph5Error("The loaded f2py module '%s' @ '%s' does not have function %s."%(
                                        self.f2py_module_path[1], self.f2py_module_path[0], function_name))
            
        return getattr(self.f2py_module, '%s%s'%(self.proc_prefix, function_name))

    def load_f2py_module(self, module_path):
        """ Loads a f2py module given its path and returns it."""

        # Make sure to temporarily adjust the environment
        added_path = False
        if module_path[0] not in sys.path:
            added_path = True
            sys.path.insert(0, pjoin(self.root_path,module_path[0]))
        try:
            f2py_module = importlib.import_module(module_path[1])
        except:
            # Try again after having recompiled this module
            self.compile(mode='always')
            f2py_module = importlib.import_module(module_path[1])
        if added_path:
            sys.path.pop(sys.path.index(pjoin(self.root_path,module_path[0])))
        return f2py_module
    
    def format_momenta_for_f2py(self, p):
        """ fortran/C-python do not order table in the same order.
        Also, remove the mass component of the momenta. """
        
        return np.array(p, dtype=DTYPE).transpose()
        
#        new_p = []        
#        for i in range(4):
#            new_p.append([0]*len(p))
#        for i, onep in enumerate(p):
#            for j, x in enumerate(onep):
#                if j==4: continue
#                new_p[j][i] = x
#        return new_p

    def get_squared_order_entry(self, squared_order):
        """ Returns the squared order index corresponding to the squared order in argument.
        Returns None if not available."""
        try:
            return self.squared_orders[squared_order]
        except KeyError:
            return None

    def check_inputs_validity(self, opts):
        """ Check the validity of the inputs."""
        
        new_opts = { 'squared_orders'   : tuple(sorted(opts['squared_orders'])) if opts['squared_orders'] else None,
                     'spin_correlation' : tuple( (sc[0], tuple(sc[1]) ) for sc in opts['spin_correlation']) if \
                                                                                      opts['spin_correlation'] else None,
                     'color_correlation': tuple(tuple(tuple(cc) for cc in opts['color_correlation'])) if \
                                                                                     opts['color_correlation'] else None,
                     'hel_config'       : tuple(opts['hel_config']) if opts['hel_config'] else None
                   }
    
        if new_opts['squared_orders']:
            if self.get_squared_order_entry(new_opts['squared_orders']) is None:
                raise MadGraph5Error("The following accessor:\n%s\ncannot provide squared orders %s."%(
                        self.nice_string(), str(new_opts['squared_orders'])))
        
        if new_opts['spin_correlation']:
            if not self.spin_correlations or len(new_opts['spin_correlation']) > self.spin_correlations.count('N'):
                raise MadGraph5Error("%sLO spin-correlations not supported in accessor:\n%s."%(
                                                    'N'*len(new_opts['spin_correlation']) , self.nice_string() ) )
                
        if new_opts['color_correlation']:
            if not self.color_correlations or any(cc not in self.color_correlations for cc in new_opts['color_correlation']):
                raise MadGraph5Error("Color correlation %s not supported in accessor:\n%s."%(
                                                       str(new_opts['color_correlation']) , self.nice_string() ) )
        if new_opts['hel_config']:
            if new_opts['hel_config'] not in self.helicity_configurations:
                raise MadGraph5Error("Helicity configuration %s not present in accessor:\n%s."%(
                                                         str(new_opts['hel_config']) , self.nice_string() ) )
        
        return new_opts
        
    def clean_ME_settings(self, opts):
        """ Clean up possible Matrix Elements setting for a particular call."""
        
        if opts['spin_correlation'] and self.spin_correlations:
            # By default, do no compute any spin correlators
            self.get_function('reset_spin_correlation_vectors')()
        
        if opts['color_correlation'] and self.color_correlations:
            # By default, do no compute any color correlators
            self.get_function('set_color_correlators_to_consider')(0,0)
        
    def setup_ME(self, opts):
        """ Setup some MatrixElement steering variables according to user-defined options."""
        
        if opts['spin_correlation']:
            for (legID, spin_correlation_vectors) in opts['spin_correlation']:
                sc_vectors = [[0.0,0.0,0.0,0.0]]*self.max_spin_corr_vectors
                for i, sc_vector in enumerate(spin_correlation_vectors):
                    sc_vectors[i] = list(sc_vector)
                self.get_function('set_spin_correlation_vectors')(legID, len(spin_correlation_vectors), sc_vectors)

        if opts['color_correlation']:
            # If one is using shortcuts with -1 entries to indicate all sums over a particular leg,
            # then use 'set_color_correlators_to_consider' otherwise use 'add_color_correlators_to_consider'
            if len(opts['color_correlation'])==1:
                self.get_function('set_color_correlators_to_consider')\
                                    (opts['color_correlation'][0][0],opts['color_correlation'][0][1])
            else:
                for color_correlator in opts['color_correlation']:
                    self.get_function('add_color_correlators_to_consider')(color_correlator[0],color_correlator[1])
        
    
    def call_tree_ME(self, func, *args, **opts):
        """ Wrapper around the actual call of the tree-level matrix element, so as to be
        able to easily time it with a profiler."""

        return func(*args, **opts)
    
    def __call__(self, PS_point, alpha_s, mu_r=91.188, return_all_res=False, **opts):
        """ Actually performs the f2py call. """

        permutation = opts['permutation']

        # The mother class takes care of applying the permutations for the generic options
        PS_point, opts = VirtualMEAccessor.__call__(self, PS_point, **opts)

        is_cache_active = opts.get('cache_active', self.cache_active)

        PS_point = self.format_momenta_for_f2py(PS_point)
    
        new_opts = self.check_inputs_validity(opts)
        
        # tuple(tuple(p) for p in PS_point)
        # Make the numpy array hashable
        PS_point.flags.writeable = False
        this_call_key = { 'PS_point' : hash(PS_point.data),
                          'alpha_s'  : alpha_s }

        # We can only recycle results where color correlations are either not specified or only one is specified.    
        no_multiple_color_connections = \
            new_opts['color_correlation'] is None or len(new_opts['color_correlation'])==1 and \
                                       all(cc>0 for cc in new_opts['color_correlation'][0]) 
        # always return all results if more than one color connection asked for
        return_all_res = return_all_res or not no_multiple_color_connections
        
        if is_cache_active and no_multiple_color_connections:
            result_key = dict(new_opts)
            result_key['color_correlation'] = None if not new_opts['color_correlation'] else new_opts['color_correlation'][0]
            
            recycled_call = self.cache.get_result(**this_call_key)
            recycled_result = recycled_call.get_result(**result_key)
            if recycled_result:
                if return_all_res:
                    return MEEvaluation(recycled_result), recycled_call.get_inverse_permuted_copy(permutation)
                else:
                    return MEEvaluation(recycled_result), None

        # If/When grouping several processes in the same f2py module (so as to reuse the model for example),
        # we will be able to use the information of self.process_pdgs to determine which one to call.
        # misc.sprint(" I was called from :",self.process_pdgs)
        if not self.module_initialized:
            self.get_function('initialise')(pjoin(self.root_path, self.slha_card_path))
            self.module_initialized = True
        
        # Setup Matrix Element code variables for the user-defined options
        self.setup_ME(new_opts)

        # Actual call to the matrix element
        main_output = self.call_tree_ME(
            self.get_function('me_accessor_hook'),
            PS_point,
            (-1 if not new_opts['hel_config'] else self.helicity_configurations[new_opts['hel_config']]),
            alpha_s)

        # Gather additional newly generated output_data to be returned and placed in the cache.
        output_datas = self.gather_output_datas(main_output, new_opts, return_all_res)

       # Make sure to clean up specification of various properties for that particular call
        self.clean_ME_settings(new_opts)

        if is_cache_active:
            ME_result = self.cache.add_result(**this_call_key)
        else:
            ME_result = MEResult()
        
        if is_cache_active or return_all_res:
            for output_data in output_datas:
                ME_result.add_result(output_data[1], **output_data[0])
            
        if not return_all_res:
            return output_datas[0][1], None

        # Now recover the main result the user expects. If he did not specify a specific
        # single color_correlators, we will chose to return the result without color_correlation.
        main_result_key = dict(new_opts)
        if new_opts['color_correlation'] and len(new_opts['color_correlation'])==1 and not \
                                        any(sc<=0 for sc in new_opts['color_correlation'][0]):
            main_result_key['color_correlation'] = new_opts['color_correlation'][0]
        else:
            main_result_key['color_correlation'] = None

        main_result = MEEvaluation(ME_result.get_result(**main_result_key))
        
        # Now return a dictionary containing the expected result anticipated by the user given the specified options,
        # along with a copy of the ME_result dictionary storing all information available at this point for this call_key
        return main_result, ME_result.get_inverse_permuted_copy(permutation)

    def is_color_correlation_selected(self, color_correlator, color_correlation_specified):
        """ Check if a particular spin_correlator is among those specified by the user."""

        if (color_correlation_specified[0][0]< 0 or color_correlation_specified[0][0]==color_correlator[0]) and \
           (color_correlation_specified[0][1]< 0 or color_correlation_specified[0][1]==color_correlator[1]):
            return True
        
        return False

    def gather_output_datas(self, main_output, user_opts, return_all_res):
        """ Gather additional newly generated output_data to be returned and placed in the cache.
            This functions returns output_datas  which is a list of 2-tuples of the form:
                  ( dictionary_describing_data, dictionary_of_data ) 
            where 'dictionary_describing_data' typically contains a particular spin_correlation describer, 
            color_correlation, etc.. and 'dictionary_of_data' typically contains the keys ['finite', 'eps**x', 'ML_return_code', etc...]  
            or just 'finite' for tree-level.
        """

        output_datas = []
        
        # Basic output template dictionary. We can already assign hel_config and spin_correlation since nothing
        # can be computed on top.
        output_key_template = {'hel_config'           : user_opts['hel_config'],
                               'squared_orders'       : None,
                               'color_correlation'    : None,
                               'spin_correlation'     : user_opts['spin_correlation']}
        
        output_result_template = MEEvaluation({'finite' : 0.})
        
        if not return_all_res:
            if user_opts['squared_orders'] is None:
                i_sqso = 0
            else:
                output_key_template['squared_orders'] = user_opts['squared_orders']
                i_sqso = self.squared_orders[user_opts['squared_orders']]
            if user_opts['color_correlation'] is None:
                output_result_template['finite'] = main_output[i_sqso]
            else:
                output_key_template['color_correlation'] = user_opts['color_correlation'][0]
                color_correlated_mes = self.get_function('get_color_correlated_me')()
                output_result_template['finite'] = color_correlated_mes[\
                    self.color_correlations[user_opts['color_correlation'][0]]-1][i_sqso]

            return [ (output_key_template, output_result_template) ]

        # First add the general squared order results
        for i_sqso, squared_order in self.id_to_squared_order.items():
            output_result = MEEvaluation(output_result_template)
            output_result['finite'] = main_output[i_sqso]
            # Now add this piece of data to the list to be added to the MEResult
            output_key = dict(output_key_template)
            output_key['squared_orders'] = squared_order
            output_datas.append( (output_key, output_result) )
        
        # Now add the color-correlated results.
        # This is a rank-2 array, with first index labeling the squared orders and the second labelling
        # the color correlators
        if user_opts['color_correlation']:
            color_correlated_mes = self.get_function('get_color_correlated_me')()
            for i_cc, color_correlator in self.id_to_color_correlation.items():
                if not self.is_color_correlation_selected(color_correlator, user_opts['color_correlation']):
                    continue
                output_key_template['color_correlation'] = color_correlator
                for i_sqso, squared_order in self.id_to_squared_order.items():
                    output_key = dict(output_key_template)
                    output_key['squared_orders'] = squared_order
                    output_result = MEEvaluation(output_result_template)
                    output_result['finite'] = color_correlated_mes[i_cc-1][i_sqso]
                    output_datas.append( (output_key, output_result) )

        return output_datas
    
class F2PYMEAccessorMadLoop(F2PYMEAccessor):
    """ Specialization of F2PYMEAccessor for MadLoop."""

    cache_active = False

    def __init__(self, process, f2py_module_path, slha_card_path, madloop_resources_path=None, **opts):
        """ Use the MadLoop resources path for MadLoop outputs """
  
        super(F2PYMEAccessorMadLoop, self).__init__(process, 
                                                 f2py_module_path, slha_card_path, **opts)
        
        # Add the additional inputs to the back-up to be used later for the dump 
        #  (if not already added by the mother)
        if not 'madloop_resources_path' in self.initialization_inputs['opts']:
            self.initialization_inputs['opts']['madloop_resources_path'] = madloop_resources_path

        if madloop_resources_path:
            # Save the location of MadLoop resources directory path
            if not os.path.isabs(madloop_resources_path):
                raise MadGraph5Error("The initialization of F2PYMEAccessorMadLoop necessitates an "+
                                                            "absolute path for madloop_resources_path.")
            self.madloop_resources_path = os.path.relpath(madloop_resources_path, self.root_path)
        else:
            self.madloop_resources_path = None

        if self.f2py_module is None:
            # We are at the generation stage and the f2py module has never been compiled yet.
            # This accessor is therefore only a placeholder for now.
            return
    
        # Obtain loop_squared_orders. Split names are the same as those already determined at tree-level
        # 'None' means summed over. 
        self.loop_squared_orders = {None: 0}
        if self.has_function('get_nsqso_loop'):
            n_squared_orders = self.get_function('get_nsqso_loop')()
            for i_sq_order in range(1,n_squared_orders+1):
                key = tuple(sorted([(self.split_order_names[i], value) for i, value in enumerate(
                                            self.get_function('ml5get_squared_orders_for_soindex')(i_sq_order) )]))
                self.loop_squared_orders[key] = i_sq_order
        
        # Reversed map
        self.id_to_loop_squared_order = dict((value,key) for (key,value) in self.loop_squared_orders.items())
        
        # List all available squared order combinations to consider and which entries to consider for each
        all_entries = {'tree':True, 'finite':True, 'eps^-1':True, 'eps^-2':True, 'accuracy':True}
        tree_entries = {'tree':True, 'finite':False, 'eps^-1':False, 'eps^-2':False, 'accuracy':False}
        loop_entries = {'tree':False, 'finite':True, 'eps^-1':True, 'eps^-2':True, 'accuracy':True}
        self.all_squared_orders = [  (all_entries,0,None) ] + \
            [ (tree_entries, sqso[0], sqso[1]) for sqso in self.id_to_squared_order.items() if not sqso[1] is None] + \
            [ (loop_entries, sqso[0], sqso[1]) for sqso in self.id_to_loop_squared_order.items() if not sqso[1] is None]
        
        # Returned array dimension
        self.res_dim = self.get_function('get_answer_dimension')()


    def synchronize(self, ME7_options = None, from_init=False, refresh_filters='auto', compile='auto',
                                                                                        proc_dirs_initialized=[], **opts):
        """ Synchronizes this accessor with the possibly updated value of parameter cards and ME source code.
        Must be defined by daughter classes."""
        
        ret_value = super(F2PYMEAccessorMadLoop, self).synchronize(ME7_options, from_init=from_init,**opts)
        
        # Never do anything when called from the __init__ function
        if from_init or not ME7_options or refresh_filters=='never':
            return ret_value

        # (Re-)generate the Helicity (and loop) filters if needed.        
        # Check if the specified madloop_resources_path is standard
        Pdir  = pjoin(self.root_path, self.proc_dir)
        # Let us not re-reinitialize process dirs that have already been re-initialized from other ME accesors.
        if Pdir in proc_dirs_initialized:
            return ret_value

        default_ML5_resources_dir = os.path.normpath(pjoin(Pdir,'SubProcesses','MadLoop5_resources'))
        if os.path.normpath(pjoin(self.root_path,self.madloop_resources_path)) != default_ML5_resources_dir:
            logger.info("Non-standard MadLoop resources directory provided. The automatic"+
                        " setup of the loop and helicity filter is therefore skipped.")
            return ret_value
        
        # Check if it needs initialization
        if MadLoopInitializer.need_MadLoopInit(Pdir, subproc_prefix='P', force_initialization=(refresh_filters=='always')):
            # This will run in parallel the initialization of *all* the processes in this 'SubProcesses' directory.
            # Not only the filters relevant to *this* accessor will be affected then.
            MadLoopInitializer.init_MadLoop(Pdir, subproc_prefix='P', MG_options=ME7_options)
            # Flag this directory as initialized so that other ME accessor do not immediately re-initialize it,
            # as it could be the case if refresh_filters is set to 'always'.
            proc_dirs_initialized.append(Pdir)
        return ret_value

    def nice_string(self):
        """ Additional details for this loop MEaccessor."""
        
        res = [super(F2PYMEAccessorMadLoop, self).nice_string()]
        res.append('%-40s:   %s'%('Loop squared orders available', 'Only summed' if self.loop_squared_orders.keys() is None else \
                                                                str([k for k in self.loop_squared_orders.keys() if not k is None])))
        return '\n'.join(res)

    def get_squared_order_entry(self, squared_order):
        """ Returns (i,j) where j is the squared order index corresponding to the squared order in argument.
        'i' is zero if this is a Born (tree) squared-order specification and one if it is a loop one.
        Returns None if not available."""
        try:
            return (1,self.loop_squared_orders[squared_order])
        except KeyError:
            try:
                return (0,self.squared_orders[squared_order])
            except KeyError:
                return None

    def setup_ME(self, opts):
        """ Setup some extra MatrixElement steering variables according to user-defined options."""
        
        super(F2PYMEAccessorMadLoop, self).setup_ME(opts)

        if opts['squared_orders']:
            sq_orders = self.get_squared_order_entry(opts['squared_orders'])
            if sq_orders and sq_orders[0]==1:
                self.get_function('set_couplingorders_target')(sq_orders[1])

    def clean_ME_settings(self, opts):
        """ Additional clean-up of steering variables for Matrix Elements after it was called."""
        
        super(F2PYMEAccessorMadLoop, self).clean_ME_settings(opts)

        # Make sure to compute all squared orders available by default
        self.get_function('set_couplingorders_target')(-1)      
        
    def call_loop_ME(self, func, *args, **opts):
        """ Wrapper around the actual call of the loop-level matrix element, so as to be
        able to easily time it with a profiler."""
        
        return func(*args, **opts)

    # Disable stdout for this function so as not to be flooded by MadLoop
    @ME7_interface.wrap_with_ME7RunEnvironment(silence=True, accessor_optimization=False, loggers = [])     
    def __call__(self, PS_point, alpha_s, mu_r, **opts):
        """ Actually performs the f2py call.
        """

        permutation = opts['permutation']
                
        # The mother class takes care of applying the permutations for the generic options
        PS_point, opts = VirtualMEAccessor.__call__(self, PS_point, **opts)
        new_opts = self.check_inputs_validity(opts)

        required_accuracy = -1.0
        if 'required_accuracy' in opts:
            required_accuracy = opts['required_accuracy']  
        
        this_call_key = { 'PS_point' : tuple(tuple(p) for p in PS_point), 
                          'alpha_s'  : alpha_s,
                          'mu_r'     : mu_r,
                          'required_accuracy' : required_accuracy}

        # We can only recycle results where color correlations are either not specified or only one is specified.
        if new_opts['color_correlation'] is None or \
                len(new_opts['color_correlation'])==1 and all(cc>0 for cc in new_opts['color_correlation'][0]):
            result_key = dict(new_opts)
            result_key['color_correlation'] = None if not new_opts['color_correlation'] else new_opts['color_correlation'][0]
            
            recycled_call = self.cache.get_result(**this_call_key)
            recycled_result = recycled_call.get_result(**result_key)
            if recycled_result:
                return MEEvaluation(recycled_result), recycled_call.get_inverse_permuted_copy(permutation)

        # If/When grouping several processes in the same f2py module (so as to reuse the model for example),
        # we will be able to use the information of self.process_pdgs to determine which one to call.
        # misc.sprint(" I was called from :",self.process_pdgs)
        if not self.module_initialized:
            # These functions are not prefixed, so we should not ask to get them via the accessor self.get_function
            self.f2py_module.initialise(pjoin(self.root_path, self.slha_card_path))
            if self.madloop_resources_path:
                    self.f2py_module.initialise_madloop_path(pjoin(self.root_path, self.madloop_resources_path))
            self.module_initialized = True
        
        
        # Setup Matrix Element code variables for the user-defined options
        self.setup_ME(new_opts)

        # Actual call to the matrix element
        evals, estimated_accuracies, return_code = self.call_loop_ME(
            self.get_function('loopme_accessor_hook'),
            self.format_momenta_for_f2py(PS_point), 
            (-1 if not new_opts['hel_config'] else self.helicity_configurations[new_opts['hel_config']]),                              
            alpha_s, mu_r, required_accuracy)

        main_output = {'evals': evals,
                       'estimated_accuracies': estimated_accuracies,
                       'return_code': return_code}   

        # Gather additional newly generated output_data to be returned and placed in the cache.
        output_datas = self.gather_output_datas(main_output, new_opts)
        
        ME_result = self.cache.add_result(**this_call_key)
        for output_data in output_datas:
            ME_result.add_result(output_data[1], **output_data[0])
        
        # Now recover the main result the user expects. If he did not specify a specific single color_correlators, we will
        # chose to return the result without color_correlation.
        main_result_key = dict(new_opts)
        if new_opts['color_correlation'] and len(new_opts['color_correlation'])==1 and not \
                                                                    any(sc<=0 for sc in new_opts['color_correlation'][0]):
            main_result_key['color_correlation'] = new_opts['color_correlation'][0]
        else:
            main_result_key['color_correlation'] = None
        main_result = MEEvaluation(ME_result.get_result(**main_result_key))
        
        # Make sure to clean up specification of various properties for that particular call
        self.clean_ME_settings(new_opts)
        
        # Now return a dictionary containing the expected result anticipated by the user given the specified options,
        # along with a copy of the ME_result dictionary storing all information available at this point for this call_key
        return main_result, ME_result.get_inverse_permuted_copy(permutation)

    def gather_output_datas(self, main_output, user_opts):
        """ Gather additional newly generated output_data to be returned and placed in the cache.
            This functions returns a 'output_datas' which is a list of 2-tuples of the form:
                  ( dictionary_describing_data, dictionary_of_data ) 
            where 'dictionary_describing_data' typically contains a particular spin_correlation describer, 
            color_correlation, etc.. and 'dictionary_of_data' typically contains the keys ['finite', 'eps^-1', 'ML_return_code', etc...]  
            or just 'finite' for tree-level.
        """

        output_datas = []

        # Basic output template dictionary. We can already assign hel_config and spin_correlation since nothing
        # can be computed on top.
        output_key_template = {'hel_config'           : user_opts['hel_config'],
                               'squared_orders'        : None,
                               'color_correlation'    : None,
                               'spin_correlation'     : user_opts['spin_correlation']}
        
        output_result_template = MEEvaluation(
                                 {'tree'        : None,
                                  'finite'      : None,
                                  'eps^-1'      : None,
                                  'eps^-2'      : None,
                                  'accuracy'    : None,
                                  'return_code' : main_output['return_code']})
    
        i_sqso_user = self.get_squared_order_entry(user_opts['squared_orders'])
        if not i_sqso_user is None:
            i_sqso_user = i_sqso_user[1]
  
        # First add the general tree squared order results
        for entries_to_consider, i_sqso, squared_order in self.all_squared_orders:
            # Do not include results for squared order that were selected out by having invoked
            # set_couplingorders_target. Also do not include the sum if a squared order selection was performed.
            if i_sqso_user and (squared_order is None or i_sqso > i_sqso_user):
                continue
            output_result = MEEvaluation(output_result_template)
            if entries_to_consider['tree']:
                output_result['tree']     = main_output['evals'][0][i_sqso]
            if entries_to_consider['finite']:
                output_result['finite']   = main_output['evals'][1][i_sqso]
            if entries_to_consider['eps^-1']:    
                output_result['eps^-1']    = main_output['evals'][2][i_sqso]
            if entries_to_consider['eps^-2']:
                output_result['eps^-2'] = main_output['evals'][3][i_sqso]
            if entries_to_consider['accuracy']:
                output_result['accuracy'] = main_output['estimated_accuracies'][i_sqso]
            # Now add this piece of data to the list to be added to the MEResult
            output_key = dict(output_key_template)
            output_key['squared_orders'] = squared_order
            output_datas.append( (output_key, output_result) )
        
        # Now add the color-correlated results.
        if user_opts['color_correlation']:
            color_correlated_mes = self.get_function('get_color_correlated_me')()            
            for i_cc, color_correlator in self.id_to_color_correlation.items():
                if not self.is_color_correlation_selected(color_correlator, user_opts['color_correlation']):
                    continue
                output_key_template['color_correlation'] = color_correlator
                for entries_to_consider, i_sqso, squared_order in self.all_squared_orders:
                    # Do not include results for squared order that were selected out by having invoked
                    # set_couplingorders_target. Also do not include the sum if a squared order selection was performed.
                    if i_sqso_user and (squared_order is None or i_sqso > i_sqso_user):
                        continue
                    output_result = MEEvaluation(output_result_template)
                    if entries_to_consider['tree']:
                        output_result['tree']     = color_correlated_mes[i_cc-1][0][i_sqso]
                    if entries_to_consider['finite']:
                        output_result['finite']   = color_correlated_mes[i_cc-1][1][i_sqso]
                    if entries_to_consider['eps^-1']:    
                        output_result['eps^-1']    = color_correlated_mes[i_cc-1][2][i_sqso]
                    if entries_to_consider['eps^-2']:
                        output_result['eps^-2'] = color_correlated_mes[i_cc-1][3][i_sqso]
                    if entries_to_consider['accuracy']:
                        output_result['accuracy'] = main_output['estimated_accuracies'][i_sqso]    
                    # Now add this piece of data to the list to be added to the MEResult
                    output_key = dict(output_key_template)
                    output_key['squared_orders'] = squared_order
                    output_datas.append( (output_key, output_result) )            

        return output_datas

class MEAccessorDict(dict):
    """ A class for nicely wrapping the access to the (possibly spin- and color- correlated) matrix elements of 
    various processes (mostly via f2py)."""
    
    cache_active = False
    
    def __init__(self, *args, **opts):
        """ Initialize an allMEAccessor. """
        super(MEAccessorDict, self).__init__(*args, **opts)

    def generate_dump(self):
        """ Generate a serializable dump of self, which can later be used, along with some more 
        information, in initialize_from_dump in order to regenerate the object."""
        
        all_MEAccessor_dumps = []
        
        # Remove the module attribute of the MEAccessor as this can't be pickled.
        for (ME_accessor, defining_pdgs_order) in self.values():
            dump = ME_accessor.generate_dump()
            if dump not in all_MEAccessor_dumps:
                all_MEAccessor_dumps.append(dump)
        
        final_dump = {}
        final_dump['class'] = self.__class__
        final_dump['all_MEAccessor_dumps']  = all_MEAccessor_dumps
        
        return final_dump
    
    @classmethod
    def initialize_from_dump(cls, dump, root_path, model):
        """ Initialize self from a dump and possibly other information necessary for reconstructing this
        contribution."""
        
        new_MEAccessorDict_instance = dump['class']()
        
        all_MEAccessors = [accessor_dump['class'].initialize_from_dump(accessor_dump, root_path, model) for 
                           accessor_dump in  dump['all_MEAccessor_dumps']]
        
        new_MEAccessorDict_instance.add_MEAccessors(all_MEAccessors)
        # For now everything is dumped, so nothing needs to be done.
        return new_MEAccessorDict_instance

    def __getitem__(self, key):
        return self.get_MEAccessor(key)

    def get_MEAccessor(self, key, pdgs=None):
        """ Provides access to a given ME, provided its ProcessKey. See implementation of this ProcessKey to
        understand how the properties of the process being accessed are stored.
        The user can specify here a particular order in which the particles/flavors are provided. When doing so
        the corresponding permutation will be computed and provided in the call_key returned 
        (which corresponds to the options to provide when calling the returned MEaccessor).
        
        Notice that pdgs specifies not only the desired order of pdgs, but also the flavors, i.e if
        
        {'u u~ > g g' : (MEAccessorInstanceA, 'u u~ > g g'), 'c c~ > g g' : (MEAccessorInstanceA, 'c c~ > g g')}
        
        and you ask for pdgs = 'c~ c > g g', then this function will return 
           MEAccessorInstanceA
        with the "call_key":
           {'permutation': [1,0,2,3], 'process_pdgs': 'c c~ > g g'}
        """

        if isinstance(key, subtraction.Current):
            # Automatically convert the process to a ProcessKey
            accessor_key = key.get_key()
        elif isinstance(key, base_objects.Process):
            # Automatically convert the process to a ProcessKey
            accessor_key = ProcessKey(process=key, PDGs=pdgs if pdgs else [])
        elif isinstance(key, ProcessKey):
            accessor_key = key
        else:
            raise MadGraph5Error("Key passed to get_MEAccessor should always be of type ProcessKey or base_objects.Process")
            
        try:
            (ME_accessor, defining_pdgs_order) = super(MEAccessorDict, self).__getitem__(accessor_key.get_canonical_key())
        except KeyError:
            raise MadGraph5Error("This collection of matrix elements does not contain process %s."%str(accessor_key.key_dict))
        
        if pdgs is None:
            # None indicates that no permutation will be necessary to apply
            return ME_accessor, {'permutation': None, 'process_pdgs': defining_pdgs_order}
        
        # First deduce the inverse mapping, that is the mapping to apply *TO the order assumed in the ME* to the user inputs.
        permutation = {}
        # Create a look_up list from the user-provided list of PDGs, whose elements will progressively be set to zero
        # as they are being mapped
        look_up_list = [list(pdgs[0]), list(pdgs[1])]
        
        # Map the initial states
        for i, pdg in enumerate(defining_pdgs_order[0]):
            try:
                permutation[i] = look_up_list[0].index(pdg)
                # Set the element that just got mapped to 0 in the look_up_list so that it will not be reused
                look_up_list[0][permutation[i]] = 0
            except ValueError:
                raise MadGraph5Error("Cannot map two PDGs list: %s and %s"%(str(defining_pdgs_order[0]), str(pdgs[0])))
        
        # Map final states now
        n_initial = len(defining_pdgs_order[0])
        for i, pdg in enumerate(defining_pdgs_order[1]):
            try:
                permutation[i+n_initial] = look_up_list[1].index(pdg)+n_initial
                # Set the element that just got mapped to 0 in the look_up_list so that it will not be reused
                look_up_list[1][permutation[i+n_initial]-n_initial] = 0
            except ValueError:
                raise MadGraph5Error("Cannot map two PDGs list: %s and %s"%(str(defining_pdgs_order[1]), str(pdgs[1])))        

        # Now inverse the mapping so as to obtain the permutations to apply *TO the user inputs* in order to get 
        # to the *order assumed in the ME*
        permutation = dict((v,k) for (k,v) in permutation.items())

        # Now return the corresponding ME_accessor, along with the options to provide when calling it, aka call_key
        # Remember that defining_pdgs_order is the list of PDGs of the original ordering of the ME_accessor, with the correct flavor
        # information , as specified  by the user via the option 'pdgs'.
        return ME_accessor, {'permutation': permutation, 'process_pdgs': defining_pdgs_order}
                
    def format_color_correlation(self, process, color_correlations):
        """ Synchronize the numbers in the color correlation specifier to the leg_number in the process."""

        leg_number_to_pos_dict = {}
        for leg_pos, leg in enumerate(process.get_initial_legs()+process.get_final_legs()):
            leg_number_to_pos_dict[leg.get('number')] = leg_pos+1
        
        new_color_correlations = []
        for color_correlation in color_correlations:
            new_color_correlations.append(tuple(leg_number_to_pos_dict[n] for n in color_correlation))
        
        return new_color_correlations
            
    def format_spin_correlation(self, process, spin_correlations):
        """ Synchronize the numbers in the spin correlation specifier to the leg_number in the process."""        
        
        leg_number_to_pos_dict = {}
        for leg_pos, leg in enumerate(process.get_initial_legs()+process.get_final_legs()):
            leg_number_to_pos_dict[leg.get('number')] = leg_pos+1
        
        new_spin_correlations = []
        for spin_correlation in spin_correlations:
            new_spin_correlations.append( ( leg_number_to_pos_dict[spin_correlation[0]], spin_correlation[1] ) )

        return new_spin_correlations
        
    def __call__(self, *args, **opts):
        """ Quick access to directly calling a Matrix Element. The first argument should always be the MEdictionary key
        and in the options the user can specify the desired_pdgs_order which will be used by the MEDictionary to figure
        out which permutation to apply and which flavours to specify.
        """

        assert (len(args)>0 and isinstance(args[0], (ProcessKey, base_objects.Process, subtraction.Current))), "When using the shortcut "+\
            "__call__ method of MEAccessorDict, the first argument should be an instance of a ProcessKey or base_objects.Process or"+\
            " subtraction.Current."

        if self.cache_active and hasattr(args[0], 'accessor'):
            if isinstance(args[0], subtraction.Current):
                return args[0].accessor(*args, **opts)
            elif isinstance(args[0], base_objects.Process):
                try:
                    accessor, call_options = args[0].accessor[tuple(sorted(opts.items()))]
                    # Remove the process isntance
                    call_args = list(args[1:])
                    # Format PS point
                    call_args[0] = args[0].format_PS_point_for_ME_call(call_args[0])
                    return accessor(*call_args, **call_options)
                except KeyError:
                    pass

        # Now store the me_accessor_key and remove it from the arguments to be passed to the MEAccessor call.
        me_accessor_key = args[0]
        
        desired_pdgs_order = None
        call_options = dict(opts)
        pdgs_specified = False
        if 'pdgs' in call_options:
            pdgs_specified = True
            desired_pdgs_order = call_options.pop('pdgs')
        
        specified_process_instance = None
        if not isinstance(args[0], subtraction.Current) and isinstance(args[0], base_objects.Process):
            # The user called this MEAccessorDictionary with a specific instance of a Process (not current), therefore
            # we must enforce the pdgs ordering specified in it (if not overwritten by the user).
            specified_process_instance = args[0]
            if not pdgs_specified:
                desired_pdgs_order = specified_process_instance.get_cached_initial_final_pdgs()
            
        ME_accessor, call_key = self.get_MEAccessor(me_accessor_key, pdgs=desired_pdgs_order)
        call_options.update(call_key)

        # Now for subtraction current accessors, we must pass the current as first argument
        call_args = list(args)
        if isinstance(ME_accessor, SubtractionCurrentAccessor):
            if not isinstance(call_args[0], subtraction.Current):
                raise MadGraph5Error("SubtractionCurrentAccessors must be called from the accessor dictionary with "+
                                     "an instance of a current as first argument.")
        else:
            # For Matrix element, we don't need to provide the specific process since this is done via the PDG list
            call_args = call_args[1:]
            # If the user specified a process instance, he might also have passed the PS point as a dictionary,
            # so we must transform it here into a flatlist:
            PS_point = call_args[0]
            if specified_process_instance and isinstance(PS_point, dict):
                call_args[0] = specified_process_instance.format_PS_point_for_ME_call(PS_point)
            # Also, if spin and color correlation are specified, we must change their ordering
            # according to the leg numbers
            if specified_process_instance and 'color_correlation' in call_options and call_options['color_correlation']:
                call_options['color_correlation'] = self.format_color_correlation(specified_process_instance, 
                                                                                  call_options['color_correlation'])
            if specified_process_instance and 'spin_correlation' in call_options and call_options['spin_correlation']:    
                call_options['spin_correlation'] = self.format_spin_correlation(specified_process_instance, 
                                                                                  call_options['spin_correlation'])
        
        if self.cache_active:
            if isinstance(args[0], subtraction.Current):
                args[0].accessor = ME_accessor
            elif isinstance(args[0], base_objects.Process):
                key = tuple(sorted(opts.items()))
                if hasattr(args[0],'accessor'):
                    args[0].accessor[key] = (ME_accessor, call_options)
                else:
                    args[0].accessor = {key:(ME_accessor, call_options)}
        return ME_accessor(*call_args, **call_options)
    
    def add_MEAccessor(self, ME_accessor, allow_overwrite=False):
        """ Add a particular ME accessor to the collection of available ME's."""
        if not isinstance(ME_accessor, VirtualMEAccessor):
            raise MadGraph5Error("MEAccessorDict can only be assigned values inheriting from VirtualMEAccessor.")
        
        for key, value in ME_accessor.get_canonical_key_value_pairs():
            if key in self and not allow_overwrite:
                raise MadGraph5Error("Attempting to assign two MEAccessors to the same key in MEAccessorDict.")
            self[key] = value

    def add_MEAccessors(self, ME_accessor_list, allow_overwrite=False):
        """ Add a list of ME_accessors."""
        for ME_accessor in ME_accessor_list:
            self.add_MEAccessor(ME_accessor, allow_overwrite=allow_overwrite)

    def synchronize(self, *args, **opts):
        """ Synchronizes all the accessors with the possibly updated value of parameter cards and ME source code."""
        
        # We want to make sure that we don't refresh the filter of all processes in a given directory output
        # more than one time. So we always pass here the addional option proc_dir_initialized which gets filled
        # in as the different accessor initialize the MadLoop filters.
        opts['proc_dirs_initialized'] = []
        for (ME_accessor, defining_pdgs_order) in self.values():
            ME_accessor.synchronize(*args, **opts)

def activate_cache():
    """ Activate caches of various accessor classes. This is safe to do during an actual MC run."""
    ProcessKey.cache_active = True
    MEAccessorDict.cache_active = True

def deactivate_cache():
    """ Deactivate caches of various accessor classes. This should be deactivated during process generation."""
    ProcessKey.cache_active = False
    MEAccessorDict.cache_active = False

# MEAccessor classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Accessor class.
# For instance 'Unknown' can be mapped to a user-defined class.
# Notice that this map must be placed after the MEAccessor daughter classes have been declared.
MEAccessor_classes_map = {'PythonAccessor': F2PYMEAccessor,
                          'CurrentAccessor': SubtractionCurrentAccessor,
                          'Unknown': None}
