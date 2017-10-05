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
"""Classes for implementions Contributions.
Contributions are abstract layer, between the madgraph_interface and amplitudes.
They typically correspond to the Born, R, V, RV, RR, VV, etc.. pieces of higher
order correction computations.
"""
import array
import copy
import itertools
import time
import logging
import sys
import importlib
import os
import shutil
import collections
pjoin = os.path.join

import madgraph.core.base_objects as base_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.loop.loop_diagram_generation as loop_diagram_generation
import madgraph.integrator.phase_space_generators as PS_utils
import madgraph.interface.madevent_interface as madevent_interface
import madgraph.core.helas_objects as helas_objects
import madgraph.loop.loop_helas_objects as loop_helas_objects
import madgraph.iolibs.group_subprocs as group_subprocs
import madgraph.iolibs.helas_call_writers as helas_call_writers
import madgraph.iolibs.export_v4 as export_v4
import madgraph.various.misc as misc
import madgraph.core.subtraction as subtraction
import madgraph.interface.ME7_interface as ME7_interface
from madgraph import InvalidCmd, MadGraph5Error
from madgraph.iolibs.files import cp, ln, mv

logger = logging.getLogger('contributions')

##########################################################################################
# ProcessKey, to be used as keys in the MEAccessorDict and the processes_map of contributions
##########################################################################################
class ProcessKey(object):
    """ We store here the relevant information for accessing a particular ME."""
    
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
                # if he so wants.
                **opts):
              
        # Initialize a dictionary which will be used to form the final tuple encoding all the information for 
        # this particular entry
        self.key_dict = {}

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

    def set(self, key, value):
        """ Modify an entry in the key_dict created."""
        
        if key not in self.key_dict:
            raise MadGraph5Error("Key '%s' was not found in the key_dict created in ProcessKey."%key)
        if not isinstance(value, (int, str, bool, tuple)):
            raise MadGraph5Error("Values for the key_dict created in ProcessKey should be of type (int, str, bool, tuple).")
        
        self.key_dict[key] = value
        # Force the canonical key to be recomputed
        if hasattr(self, 'canonical_key'):
            delattr(self, 'canonical_key')

    def get_canonical_key(self, force=False):
        """ Simply uses self.key_dict to return a hashable canonical representation of this object."""
        if not force and hasattr(self,'canonical_key'):
            return self.canonical_key
        
        self.canonical_key = tuple(sorted(self.key_dict.items()))        
        return self.canonical_key

################################################################################
# MEAccessor mother class
################################################################################
class VirtualMEAccessor(object):
    """ A class wrapping the access to one particular MatrixEleemnt."""
    
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
            return list(PS_point), all_opts
        
        
        # Apply the permutations while retaining a canonical representation for each attribute

        # The permutation is how to go *from* the user-provided order to the order assumed by the underlying ME.
        # So When performing the permutation like below, we must use the reversed_permutation.
        reversed_permutation = dict((v,k) for (k,v) in permutation.items())
        permuted_PS_point = [PS_point[reversed_permutation[i]] for i in range(len(PS_point))] if PS_point else None

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
#            misc.sprint('Recycled a result.')
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
        
        new_opts = dict(opts)
        new_opts['hel_config'] = tuple(opts['hel_config']) if ('hel_config' in opts and opts['hel_config']) else None
        
        # Make sure to drop useless specifier only relevant to MEAccessors
        for irrelevant_opt in ['permutation','process_pdgs']:
            new_opts.pop(irrelevant_opt)

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

        if new_opts['hel_config']:
            # In this case it is not optimal to check the validity of the helicity configuration, so 
            # we limit ourselves to checking if it supports helicity assignment
            if not self.subtraction_current_instance.supports_helicity_assignment:
                raise MadGraph5Error("The following subtraction current accessor:\n%s"%(
                                    self.nice_string() ) + "\ndoes not support helicity assignment.")
        

        return new_opts

    def __call__(self, current, PS_point, **opts):
        """ Evaluation of the subtraction current. """
        
        if self.subtraction_current_instance is None:
            raise MadGraph5Error("This subtraction current accessor\n'%s'\nhas not been properly initialized."%(
                                                                                            self.nice_string()))
        
        # Parse options and check their validity
        call_opts = self.check_inputs_validity(opts, current)
        
        # Set the arguments of the call
        call_args = [current, PS_point]
        
        # Now obtain the cache key directly from the current implementation
        cache_key, result_key = self.subtraction_current_instance.get_cache_and_result_key(*call_args, **call_opts)
        
        # Attempt to recycle the result
        if not cache_key is None:
            recycled_call = self.cache.get_result(**cache_key)
            recycled_result = recycled_call.get_result(**result_key)
            if recycled_result:
                return self.evaluation_class(recycled_result), self.result_class(recycled_call)

        all_evaluations = self.subtraction_current_instance.evaluate_subtraction_current(*call_args, **call_opts)
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
        Pdir = pjoin(self.root_path,self.proc_dir, 'SubProcesses','P%s'%self.proc_name)
        if not os.path.isdir(Pdir):
            raise InvalidCmd("The expected subprocess directory %s could not be found."%Pdir)
        
        if os.path.isfile(pjoin(Pdir, 'matrix_%s_py.so'%self.proc_name)) and mode=='never':
            return
        
        if mode=='always':
             misc.compile(arg=['clean'], cwd=Pdir)
        if not os.path.isfile(pjoin(Pdir, 'matrix_%s_py.so'%self.proc_name)):
            logger.debug("Compiling directory %s ..."%(pjoin(self.proc_dir, 'SubProcesses','P%s'%self.proc_name)))
            misc.compile(arg=['matrix_%s_py.so'%self.proc_name, 'MENUM=_%s_'%self.proc_name], cwd=Pdir)
        
        if not os.path.isfile(pjoin(Pdir, 'matrix_%s_py.so'%self.proc_name)):
            raise InvalidCmd("The f2py compilation of SubProcess '%s' failed.\n"%Pdir+
                "Try running 'MENUM=_%s_ make matrix_%s_py.so' by hand in this directory."%(self.proc_name,self.proc_name))

    def synchronize(self, ME7_options = None, from_init=False, compile='auto', **opts):
        """ Synchronizes this accessor with the possibly updated value of parameter cards and ME source code.
        Must be defined by daughter classes."""

        # Reset the initialization to False
        self.module_initialized = False

        # Recompile and reload the module if synchronize was not call from the __init__ function
        if not from_init:
            self.compile(mode=compile)
            with misc.Silence(active=False and (logger.level>logging.DEBUG)):
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
        # Use logger.level >= loggign.DEBUG and not > so as to forward
        # the loading text only if the user *really* wants it by setting logger.level = 0 :)
        with misc.Silence(active=False and (logger.level>=logging.DEBUG)):
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
        new_p = []
        for i in range(4):
            new_p.append([0]*len(p))
        for i, onep in enumerate(p):
            for j, x in enumerate(onep):
                if j==4: continue
                new_p[j][i] = x
        return new_p

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
        
    def clean_ME_settings(self):
        """ Clean up possible Matrix Elements setting for a particular call."""
        
        if self.spin_correlations:
            # By default, do no compute any spin correlators
            self.get_function('reset_spin_correlation_vectors')()
        
        if self.color_correlations:
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
        
        
    def __call__(self, PS_point, alpha_s, mu_r=91.188, **opts):
        """ Actually performs the f2py call. """

        permutation = opts['permutation']

        # The mother class takes care of applying the permutations for the generic options
        PS_point, opts = VirtualMEAccessor.__call__(self, PS_point, **opts)
    
        new_opts = self.check_inputs_validity(opts)
        
        this_call_key = { 'PS_point' : tuple(tuple(p) for p in PS_point), 
                          'alpha_s'  : alpha_s }
        
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
            with misc.Silence(active=(logger.level>logging.DEBUG)):
                self.get_function('initialise')(pjoin(self.root_path, self.slha_card_path))
            self.module_initialized = True
        
        # Setup Matrix Element code variables for the user-defined options
        self.setup_ME(new_opts)

        # Actual call to the matrix element
#        start = time.time() #TOBECOMMENTED
        with misc.Silence(active=(logger.level>logging.DEBUG)):
            main_output  = self.get_function('me_accessor_hook')(
                    self.format_momenta_for_f2py(PS_point), 
                    (-1 if not new_opts['hel_config'] else self.helicity_configurations[new_opts['hel_config']]),
                    alpha_s)
#        misc.sprint("The tree call took %10g"%(time.time()-start)) #TOBECOMMENTED

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
        self.clean_ME_settings()
        
        # Now return a dictionary containing the expected result anticipated by the user given the specified options,
        # along with a copy of the ME_result dictionary storing all information available at this point for this call_key
        return main_result, ME_result.get_inverse_permuted_copy(permutation)

    def is_color_correlation_selected(self, color_correlator, color_correlation_specified):
        """ Check if a particular spin_correlator is among those specified by the user."""

        if (color_correlation_specified[0][0]< 0 or color_correlation_specified[0][0]==color_correlator[0]) and \
           (color_correlation_specified[0][1]< 0 or color_correlation_specified[0][1]==color_correlator[1]):
            return True
        
        return False

    def gather_output_datas(self, main_output, user_opts):
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

    def __init__(self, process, f2py_module_path, slha_card_path, madloop_resources_path=None, **opts):
        """ Use the MadLoop resources path for MadLoop outputs """
  
        super(F2PYMEAccessorMadLoop, self).__init__(process, f2py_module_path, slha_card_path, **opts)

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
        ML_initializer = madevent_interface.MadLoopInitializer
        if ML_initializer.need_MadLoopInit(Pdir, subproc_prefix='P', force_initialization=(refresh_filters=='always')):
            # This will run in parallel the initialization of *all* the processes in this 'SubProcesses' directory.
            # Not only the filters relevant to *this* accessor will be affected then.
            ML_initializer.init_MadLoop(Pdir, subproc_prefix='P', MG_options=ME7_options)
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

    def clean_ME_settings(self):
        """ Additional clean-up of steering variables for Matrix Elements after it was called."""
        
        super(F2PYMEAccessorMadLoop, self).clean_ME_settings()

        # Make sure to compute all squared orders available by default
        self.get_function('set_couplingorders_target')(-1)      
        
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
            with misc.Silence(active=(logger.level>logging.DEBUG)):
                self.f2py_module.initialise(pjoin(self.root_path, self.slha_card_path))
                if self.madloop_resources_path:
                        self.f2py_module.initialise_madloop_path(pjoin(self.root_path, self.madloop_resources_path))
            self.module_initialized = True
        
        
        # Setup Matrix Element code variables for the user-defined options
        self.setup_ME(new_opts)

        # Actual call to the matrix element
#        start = time.time() #TOBECOMMENTED
        with misc.Silence(active=(logger.level>logging.DEBUG)):
            evals, estimated_accuracies, return_code = self.get_function('loopme_accessor_hook')(
            self.format_momenta_for_f2py(PS_point), 
            (-1 if not new_opts['hel_config'] else self.helicity_configurations[new_opts['hel_config']]),                              
            alpha_s, mu_r, required_accuracy)
#        misc.sprint("The loop call took %10g"%(time.time()-start)) #TOBECOMMENTED
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
        self.clean_ME_settings()
        
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
    
    def format_PS_point_for_ME_call(self, PS_point, process):
        """ From a dictionary formatted PS point and a process, returns the PS point as a flat list, ordered as
        the legs in the process."""
        
        formatted_PS_point = PS_utils.LorentzVectorList()
        for leg in process.get_initial_legs()+process.get_final_legs():
            try:
                formatted_PS_point.append(PS_point[leg.get('number')])
            except KeyError:
                raise MadGraph5Error("Cannot find leg #%d in the following PS point specifications:\n%s"%(
                                                                                leg.get('number'),str(PS_point)))
        return formatted_PS_point

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
            "subtraction.Current."

        # Now store the me_accessor_key and remove it from the arguments to be passed to the MEAccessor call.
        me_accessor_key = args[0]

        desired_pdgs_order = None
        call_options = dict(opts)
        if 'pdgs' in call_options:
            desired_pdgs_order = call_options.pop('pdgs')
        
        specified_process_instance = None
        if not isinstance(args[0], subtraction.Current) and isinstance(args[0], base_objects.Process):
            # The user called this MEAccessorDictionary with a specific instance of a Process (not current), therefore
            # we must enforce the pdgs ordering specified in it (if not overwritten by the user).
            specified_process_instance = args[0]
            if desired_pdgs_order is None:
                desired_pdgs_order = (tuple(specified_process_instance.get_initial_ids()),
                                      tuple(specified_process_instance.get_final_ids_after_decay()))
            
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
                call_args[0] = self.format_PS_point_for_ME_call(PS_point, specified_process_instance)
            # Also, if spin and color correlation are specified, we must change their ordering
            # according to the leg numbers
            if specified_process_instance and 'color_correlation' in call_options and call_options['color_correlation']:
                call_options['color_correlation'] = self.format_color_correlation(specified_process_instance, 
                                                                                  call_options['color_correlation'])
            if specified_process_instance and 'spin_correlation' in call_options and call_options['spin_correlation']:    
                call_options['spin_correlation'] = self.format_spin_correlation(specified_process_instance, 
                                                                                  call_options['spin_correlation'])
        
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

#===============================================================================
# Contribution mother class
#===============================================================================
class Contribution(object):
    """ An high-level object that stores all the relevant information for a 
    particular contribution to a cross-section calculation, including higher
    order corrections.
    We do not inherit from PhysicsObject as we want to depart from the dict
    structure imposed by this mother."""
    
    def __new__(cls, contribution_definition, cmd_interface, **opt):

        if cls is Contribution:
            target_type = 'Unknown'
            if contribution_definition.correction_order == 'LO':
                if contribution_definition.n_loops == 0 and \
                  contribution_definition.n_unresolved_particles == 0:
                    target_type = 'Born'
                elif contribution_definition.n_loops == 1 and \
                  contribution_definition.n_unresolved_particles == 0:
                    target_type = 'LoopInduced_Born'                    
            elif contribution_definition.correction_order == 'NLO':
                if contribution_definition.n_loops == 1 and \
                   contribution_definition.n_unresolved_particles == 0:
                    target_type = 'Virtual'
                elif contribution_definition.n_loops == 0 and \
                     contribution_definition.n_unresolved_particles == 1:
                    target_type = 'SingleReals'
            elif contribution_definition.correction_order == 'NNLO':
                if contribution_definition.n_loops == 0 and \
                   contribution_definition.n_unresolved_particles == 2:
                    target_type = 'DoubleReals'
                else:
                    raise MadGraph5Error("Some NNLO type of contributions are not implemented yet.")
            elif contribution_definition.correction_order == 'NNNLO':
                if contribution_definition.n_loops == 0 and \
                   contribution_definition.n_unresolved_particles == 3:
                    target_type = 'TripleReals'
                else:
                    raise MadGraph5Error("Some NNNLO type of contributions are not implemented yet.")

            else:
                target_type = 'Unknown'
            target_class = Contribution_classes_map[target_type]
            if not target_class:
                raise MadGraph5Error("Could not determine the class for contribution of type '%s' to be added for"%target_type+
                                     " the contribution definiton:\n%s"%str(contribution_definition.nice_string()))
            return super(Contribution, cls).__new__(target_class, contribution_definition, cmd_interface, **opt)
        else:
            return super(Contribution, cls).__new__(cls, contribution_definition, cmd_interface, **opt)
    
    def __init__(self, contribution_definition, cmd_interface, **opts):
        """ Instantiates a particular contribution."""
        
        self.contribution_definition = contribution_definition
        
        # Now set what additional export options need to be set
        self.additional_exporter_options = {'color_correlators' : None,
                                            'spin_correlators'  : None}
        # Below correlators_needed will be 1 if we are in an NLO-type of contribution (i.e. correction_order='NLO')
        # within an NNLO general computation (i.e. overall_correction_order='NNLO').
        # In this case we indeed expect to need NLO-type of correlators.
        correlators_needed = self.contribution_definition.overall_correction_order.count('N') - \
                             self.contribution_definition.correction_order.count('N')

        ##############################################################################################################
        ##############################################################################################################
        ###                                                 TEMPORARY HACK                             
        ### For testing purposes, one can force to always include NLO types of correlators in all matrix elements
        ### outputs simply with the line below
        correlators_needed = max(correlators_needed,1)
        ###
        ##############################################################################################################
        ##############################################################################################################
        
        if correlators_needed > 0:
            self.additional_exporter_options['color_correlators'] ='N'*correlators_needed+'LO'
            # Also force MG5_aMC option 'loop_color_flows' to be true as it is necessary to compute color_correlators.
            # We can do so in a neat way here by simply adding this option to self.additional_exporter_options since they
            # will overwrite the interface option when the exporter will be instantiated.
            self.additional_exporter_options['loop_color_flows'] = True
            ##############################################################################################################
            ##############################################################################################################
            ###                                                 TEMPORARY HACK
            ### Since NNLO color correlators are not available yet and we already want to be able to tinker with NNLO  outputs
            ### we force here the color correlators to be at most NLO type here. This should be removed eventually of course.
            ###
            self.additional_exporter_options['color_correlators'] ='N'*min(correlators_needed,1)+'LO'
            ####
            ##############################################################################################################
            ##############################################################################################################
            self.additional_exporter_options['spin_correlators']  ='N'*correlators_needed+'LO'
                                        
        self.amplitudes              = diagram_generation.AmplitudeList()
        self.all_matrix_elements     = helas_objects.HelasMultiProcess()
        self.exporter                = None
        
        # Things borrowed from the cmd_interface
        self.options                 = dict(cmd_interface.options)
        self.options['_model_v4_path'] = cmd_interface._model_v4_path
        self.options['_mgme_dir']    = cmd_interface._mgme_dir
        self.model                   = cmd_interface._curr_model
        
        # Initialize an IR subtraction module if necessary
        self.IR_subtraction = None
        if self.contribution_definition.n_unresolved_particles > 0:
            self.IR_subtraction = subtraction.IRSubtraction(
                self.model,
                coupling_types=self.contribution_definition.correction_couplings,
                n_unresolved=self.contribution_definition.n_unresolved_particles )
        
        # The following two attributes dicate the type of Exporter which will be assigned to this contribution
        self.output_type             = 'default'
        self.format                  = 'standalone'
        self.export_dir              = 'None'
        if self.options['group_subprocesses'] == 'Auto':
            self.collect_mirror_procs = True
            self.group_subprocesses   = True
        else:
            self.collect_mirror_procs   = self.options['group_subprocesses']
            self.group_subprocesses     = self.options['group_subprocesses']

        # Options relevant only for LO diagram generation
        self.ignore_six_quark_processes = False
        self.diagram_filter             = False 
        self.optimize                   = False
        # Flag to indicate whether this contribution supports decay chains
        # Typically only LO contributions support decay chain amplitudes
        self.supports_decay_chain       = False
        
        # Specifies if this contribution output requires its own MODEL in Source
        self.requires_its_own_model     = False
        
        # Specify the MultiProcessClass to use to generate amplitudes
        self.MultiProcessClass          = diagram_generation.MultiProcess
        
        # Add default phase-space topology specifiers
        self.processes_to_topologies = None
        self.topologies_to_processes = None
        
    def add_integrated_counterterm(self, integrated_CT_properties):
        """ By default, do not support adding integrated counterterms."""
        
        raise MadGraph5Error("The contribution of type %s cannot receive"%type(self)+
                                            " contributions from integrated counterterms.")

    def set_export_dir(self, prefix):
        """ Assigns an export directory name."""
        dir_name = self.contribution_definition.get_shell_name()
        # Use the name of the first process since we can't access the name of the ProcessDefinition
        dir_name += '_%s'%self.amplitudes[0].get('process').shell_string(
            schannel=False, forbid=False, main=False, pdg_order=False, print_id = False)
        dir_name += '_%d'%self.contribution_definition.process_definition.get('id')
        export_dir = pjoin(prefix, dir_name)
        if os.path.exists(export_dir):
            raise MadGraph5Error("The following contribution:\n"+self.nice_string()+
                "\ncannot be exported at location:\n"+export_dir+
                "\nsince this directory already exists.")
        self.export_dir = export_dir

    def initialize_exporter(self, cmd_interface, noclean, group_subprocesses=True):
        """ Initialize the exporter that will be associated to that particular contribution.
        noclean specifies what to do in case the output directory already exists and group_subprocesses
        whether the exporter should attempt to group identical subprocesses.
        """
        self.set_export_dir(cmd_interface._export_dir)
        self.exporter = export_v4.ExportV4Factory(
                cmd_interface, noclean, output_type=self.output_type, group_subprocesses=group_subprocesses,
                curr_amps = self.amplitudes,
                export_dir = self.export_dir,
                format = self.format,
                additional_options = self.additional_exporter_options)

    def copy_template(self, model):
        """ Copy the template structure for that contribution. Quite often, this limits itself to aksing its
        exporter to do this."""
        ret_value =  self.exporter.copy_template(model)
         # Make sure that this contribution's output as an __init__ file
        if not os.path.isfile(pjoin(self.export_dir, '__init__.py')):
            open(pjoin(self.export_dir, '__init__.py'),'w').write('')
        return ret_value

    def pass_information_from_cmd(self, cmd_interface):
        """ Pass information from the command_interface to this contribution. Most of the time, this only amounts
        to passing information to the active exporter."""
        return self.exporter.pass_information_from_cmd(cmd_interface)

    def generate_matrix_elements(self, group_processes=True):
        """Generate the Helas matrix elements before exporting. Uses the main function argument 
        'group_processes' to decide whether to use group_subprocess or not."""

        cpu_time1 = time.time()
        ndiags = 0
        
        if group_processes:
            dc_amps = diagram_generation.DecayChainAmplitudeList(\
                [amp for amp in self.amplitudes if isinstance(amp, diagram_generation.DecayChainAmplitude)])
            non_dc_amps = diagram_generation.AmplitudeList(\
                [amp for amp in self.amplitudes if not isinstance(amp, diagram_generation.DecayChainAmplitude)])
            subproc_groups = group_subprocs.SubProcessGroupList()
            matrix_elements_opts = {'optimized_output': self.options['loop_optimized_output']}
            
            if non_dc_amps:
                subproc_groups.extend(group_subprocs.SubProcessGroup.group_amplitudes(\
                  non_dc_amps, self.exporter.grouped_mode, matrix_elements_opts=matrix_elements_opts))

            if dc_amps:
                dc_subproc_group = group_subprocs.DecayChainSubProcessGroup.\
                        group_amplitudes(dc_amps, self.exporter.grouped_mode, matrix_elements_opts=matrix_elements_opts)
                subproc_groups.extend(dc_subproc_group.generate_helas_decay_chain_subproc_groups())

            ndiags = sum([len(m.get('diagrams')) for m in subproc_groups.get_matrix_elements()])
            self.all_matrix_elements = subproc_groups
            # assign a unique id number to all groups
            uid = 0
            for group in subproc_groups:
                uid += 1 # update the identification number
                for me in group.get('matrix_elements'):
                    me.get('processes')[0].set('uid', uid)

        else: # Not grouped subprocesses
            generation_mode = {}
            
            # The conditional statement tests whether we are dealing with a loop induced process.
            if isinstance(self.amplitudes[0], loop_diagram_generation.LoopAmplitude):
                generation_mode['optimized_output']=self.options['loop_optimized_output']
                HelasMultiProcessClass = loop_helas_objects.LoopHelasProcess
                compute_loop_nc = True
            else:
                HelasMultiProcessClass = helas_objects.HelasMultiProcess
                compute_loop_nc = False
            
            # Explicitly set combine_matrix_elements to False, since the grouping of
            # subprocesses has been turned off
            self.all_matrix_elements = HelasMultiProcessClass(
                self.amplitudes, compute_loop_nc=compute_loop_nc, 
                matrix_element_opts=generation_mode,
                combine_matrix_elements=False)
            
            ndiags = sum([len(me.get('diagrams')) for me in 
                                            self.all_matrix_elements.get_matrix_elements()])
            # assign a unique id number to all process
            uid = 0
            for me in self.all_matrix_elements.get_matrix_elements()[:]:
                uid += 1 # update the identification number
                me.get('processes')[0].set('uid', uid)

        cpu_time2 = time.time()

        return ndiags, cpu_time2 - cpu_time1

    def set_helas_model(self):
        """ Instantiate properly the helas model """

        if self.exporter.exporter == 'cpp':       
            self.helas_model = helas_call_writers.CPPUFOHelasCallWriter(self.model)
        elif self.options['_model_v4_path']:
            assert self.exporter.exporter == 'v4'
            self.helas_model = helas_call_writers.FortranHelasCallWriter(self.model)
        else:
            assert self.exporter.exporter == 'v4'
            self.helas_model = helas_call_writers.FortranUFOHelasCallWriter(self.model)

    def generate_code(self):
        """ Assuming the Helas Matrix Elements are now generated, we can write out the corresponding code."""

        matrix_elements = self.all_matrix_elements.get_matrix_elements()
        
        calls=0
        cpu_time_start = time.time()
        if isinstance(self.all_matrix_elements, group_subprocs.SubProcessGroupList) and \
                                                                                self.exporter.grouped_mode:
            modified, self.all_matrix_elements = self.exporter.modify_grouping(self.all_matrix_elements)
            for me_number, me in enumerate(self.all_matrix_elements):
                calls = calls + self.exporter.generate_subprocess_directory(me, self.helas_model, me_number)

        else: # Non-grouped mode
            for nb, me in enumerate(matrix_elements[:]):
                new_calls = self.exporter.generate_subprocess_directory(me, self.helas_model, nb)
                if isinstance(new_calls, int):
                    if new_calls ==0:
                        matrix_elements.remove(me)
                    else:
                        calls = calls + new_calls
        
        return calls, time.time() - cpu_time_start

    def get_maximum_order(self, order):
        """ Returns the maximum value of the coupling order specified among all the diagrams in this
        contribution."""
        
        max_order = 0
        for amplitude in self.amplitudes:
            max_order = max(max_order, amplitude.get('diagrams').get_max_order(order))
        return max_order

    def export(self, nojpeg=False, group_processes=True, args=[]):
        """ Perform the export duties, that include generation of the HelasMatrixElements and 
        the actual output of the matrix element code by the exporter."""

        # Update self.group_subprocesses with the specified value
        self.group_subprocesses = group_processes
        
        # Check if the current contribution supports Decay chain
        if not self.supports_decay_chain and any(isinstance(amp, 
                                        diagram_generation.DecayChainAmplitude) for amp in self.amplitudes):
            raise MadGraph5Error("This contribution:\n%s\nshould "%self.nice_string()+
                                 "not have been used with decay chain amplitudes.")
        
        # Sort amplitudes according to number of diagrams,
        # to get most efficient multichannel output
        self.amplitudes.sort(lambda a1, a2: a2.get_number_of_diagrams() - a1.get_number_of_diagrams())

        # Check if matrix elements are already generated
        if self.all_matrix_elements.get_matrix_elements():
            ndiags, delta_time1 = 0, 0.
        else:
            ndiags, delta_time1 = self.generate_matrix_elements(group_processes)
        
        self.set_helas_model()
        
        ncalls, delta_time2 = self.generate_code()

        matrix_elements = self.all_matrix_elements.get_matrix_elements()

        # Replace the amplitudes with the actual amplitudes from the
        # matrix elements, which allows proper diagram drawing also of
        # decay chain processes
        self.amplitudes = diagram_generation.AmplitudeList([me.get('base_amplitude') for me in matrix_elements])
        if ncalls:
            logger.info("Generated output code with %d helas calls in %0.3f s" % (ncalls, delta_time1+delta_time2))
            
        # Now investigate phase-space topologies and identify the list of kinematic configurations present
        # and, for each process key in the process map, what configuration it includes and with which defining diagrams.
        self.set_phase_space_topologies()
        
        # The printout below summarizes the topologies identified and the information gathered about the
        # matrix elements and their diagrams
#        misc.sprint('='*50)
#        processes_map = self.get_processes_map()
#        for process_key, value in self.processes_to_topologies.items():
#            misc.sprint(processes_map[process_key][0].nice_string())
#            misc.sprint('This process has the following topologies:',value['topologies'])
#            misc.sprint('And its diagrams correspond to the following topologies:',value['diagrams_topology'])
#        misc.sprint('='*50)
#        for topology, value in self.topologies_to_processes.items():
#            misc.sprint('The topology with key:',topology)
#            misc.sprint('Applies to the following processes:',
#                '\n'.join(processes_map[proc_key][0].nice_string() for proc_key in value['process_keys']))
#            misc.sprint('And is characterized by the following s-channels:\n%s\nand t-channels:\n%s'%
#                    ( ', '.join('%s > %d'%(
#                        ' '.join('%d'%leg['number'] for leg in vertex['legs'][:-1]),
#                        vertex['legs'][-1]['number']) for vertex in value['s_and_t_channels'][0]), 
#                      ', '.join('%s > %d'%(
#                        ' '.join('%d'%leg['number'] for leg in vertex['legs'][:-1]),
#                        vertex['legs'][-1]['number']) for vertex in value['s_and_t_channels'][1])
#                    ) )
#        misc.sprint('='*50)
#        stop
        
        # In principle information can be passed back to ME7_exporter with the return value
        return {}

    def set_phase_space_topologies(self):
        """ Investigate phase-space topologies and identify the list of kinematic configurations present
        and, for each process key in the process map, what configuration it includes and with which defining diagrams."""
        
        processes_map = self.get_processes_map()
        
        # Store the kinematic configurations identified in a dictionary, with format
        # key: (subprocess_group_number, configuration_number)
        # values: {'s_and_t_channels':s_and_t_channels, 'process_keys': ProcessKeysIncludingIt)
        topologies_to_processes = {}
        
        # This is the inverse dictionary, with keys being processes keys and values listing the 
        # topologies keys (subprocess_group, config_number) relevant for it.
        # key: ProcessKey
        # values: {'topologies' : [(topology_key1, representative_diagram_number_for_it), ...], 
        #          'diagrams_topology': [list of the topology_key for each diagram],
        processes_to_topologies = {}

        if isinstance(self.all_matrix_elements, group_subprocs.SubProcessGroupList):
            for subprocess_group_number, subproc_group in enumerate(self.all_matrix_elements):
                subproc_diagrams_for_config = subproc_group.get('diagrams_for_configs')
                # This 'subproc_diagrams_for_config' a list of length equal to the number of configurations:
                #
                # subproc_diagrams_for_config = [diags_for_config1, diags_for_config2, etc...]
                #
                # And for each configuration, the diags_for_configI is a simple list returned by
                # the function 'get_subproc_diagrams_for_config' which is of length equal to the 
                # number of matrix elements in that Subprocess group and contains, for each matrix element, 
                # the diagram number that is representative of configI.
                #
                # Example from the first SubProcessGroup of 'p p > e+ e- j'
                # [[1, 1, 3, 3], [2, 2, 4, 4], [3, 3, 1, 1], [4, 4, 2, 2]]
                # Means that there is 4 configurations, and the subrocess group as 4 matrix elements.
                # The representative diagrams for the first configuration are, for each of the 4 matrix elements: 
                # 1, 1, 3 and 3.
                # On the other hand, in the 3rd matrix elements, the diagrams that are representative of the 
                # configurations are:
                # diag #3 is representative of Config 1
                # diag #4 is representative of Config 2
                # diag #1 is representative of Config 3
                # diag #2 is representative of Config 4
                # Notice that if no diagrams in a ME corresponds to a particular config, then the number 
                # 0 will be given.
                matrix_elements = subproc_group.get('matrix_elements')
                # Get initial and final number of legs
                (nexternal, ninitial) = subproc_group.get_nexternal_ninitial()
                for config_number, representative_diag_numbers in enumerate(subproc_diagrams_for_config):
                    # Skip a topology that has no representative diagram at all.
                    if set(representative_diag_numbers)==set([0]):
                        continue
                    topology_key = (subprocess_group_number, config_number)
                    topologies_to_processes[topology_key] = {
                            's_and_t_channels' : None,
                            'process_keys'     : [] }
                    representative_diagrams = []
                    for me_number, diag_number in enumerate(representative_diag_numbers):
                        matrix_element = matrix_elements[me_number]
                        process_key = ProcessKey(matrix_element.get('processes')[0],sort_PDGs=False).get_canonical_key()
                        # Check if that matrix element has a diagram matching this
                        if diag_number!=0:
                            topologies_to_processes[topology_key]['process_keys'].append(process_key)
                            representative_diagrams.append(matrix_element.get('diagrams')[diag_number-1])
                            # Add the identified process <-> topology relation to the processes_to_topologies
                            # dictionary
                            if not process_key in processes_to_topologies:
                                processes_to_topologies[process_key] = {
                                    'topologies' : [],
                                    'diagrams_topology': []}
                            proc_to_topo = processes_to_topologies[process_key]
                            proc_to_topo['topologies'].append( (topology_key, diag_number) )
                    # Set the s and t channels pertaining to this topology
                    topologies_to_processes[topology_key]['s_and_t_channels'] = \
                            self.extract_s_and_t_channels(nexternal, ninitial, representative_diagrams)
                            
                # Now fill in the entry 'diagrams_topology' of each process_key in processes_to_topologies.
                diagrams_maps = subproc_group.get('diagram_maps')
                # The 'diagrams_maps' is a dictionary with a number of values equal to the number of matrix elements
                # in that particular subprocess group. So:
                # diagram_map = {me_number: diagram_map_for_that_ME, etc...]
                # Each diagram_map_for_that_ME is then a list of length of the number of diagrams in the 
                # I^th matrix element, which contains, for each diagram, the configuration ID; that is 
                # the position that this diagram corresponds to in 'mapping diagrams' (with index starting at 1)
                for me_number, diagram_numbers in diagrams_maps.items():
                    matrix_element = matrix_elements[me_number]
                    process_key = ProcessKey(matrix_element.get('processes')[0],sort_PDGs=False).get_canonical_key()
                    processes_to_topologies[process_key]['diagrams_topology'] = [ 
                        (subprocess_group_number, topo_ID-1) for topo_ID in diagram_numbers ]
                
        else: # Non-grouped mode, here each diagram of each process is a topology
            for me_number, matrix_element in enumerate(self.all_matrix_elements.get_matrix_elements()):
                (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()
                process_key = ProcessKey(matrix_element.get('processes')[0],sort_PDGs=False).get_canonical_key()
                processes_to_topologies[process_key] = {
                                    'topologies' : [],
                                    'diagrams_topology': []}
                for diag_number, diagram in enumerate(matrix_element.get('diagrams')): 
                    # Notice that here the convention is to assign a key to the topology that is not
                    # (subproc_group_number, proc_number) but (me_number, diag_number) instead
                    topology_key = (me_number, diag_number)
                    topologies_to_processes[topology_key] = {
                            's_and_t_channels' : self.extract_s_and_t_channels(nexternal, ninitial, [diagram]),
                            'process_keys'     : [process_key,] }
                    processes_to_topologies[process_key]['topologies'].append(topology_key)
                    processes_to_topologies[process_key]['diagrams_topology'].append( (topology_key, (diag_number+1)) )
                
                                  
        # Now finally save the generate topology dictionaries as attributes of this class
        self.topologies_to_processes = topologies_to_processes
        self.processes_to_topologies = processes_to_topologies

    def extract_s_and_t_channels(self, nexternal, ninitial, helas_diagrams):
        """ Given the list of helas diagrams representative of a given topology, return the list of s- and t-
        channels they contain."""
        
        # Obtain an unused PDG code
        new_pdg = self.model.get_first_non_pdg()
        
        stchannels = []
        empty_verts = []
        for h in helas_diagrams:
            if h:
                # get_s_and_t_channels gives vertices starting from
                # final state external particles and working inwards
                stchannels.append(h.get('amplitudes')[0].\
                                  get_s_and_t_channels(ninitial, self.model, new_pdg))
            else:
                stchannels.append((empty_verts, None))

        # For t-channels, just need the first non-empty one
        tchannels = [t for s,t in stchannels if t != None][0]
        
        # For s_and_t_channels (to be used later) use only first config
        return ([ s for s,t in stchannels if t != None ][0], tchannels)

    def add_ME_accessors(self, all_MEAccessors, root_path):
        """ Adds all MEAccessors for the matrix elements generated as part of this contribution."""
        
        MEAccessors = []
        for process_key, (defining_process, mapped_processes) in self.get_processes_map().items():
            # The PDGs of the hashed representations correspond to entry [0][0]
            mapped_process_pdgs = [ (proc.get_initial_ids(), proc.get_final_ids()) 
                                                            for proc in mapped_processes ]
            f2py_module_path = pjoin(self.export_dir,'matrix_%s_py.so'%
                                                   self.process_dir_name(defining_process))
            if not os.path.exists(f2py_module_path):
                raise MadGraph5Error("Cannot find the compiled f2py module for %s at %s"%
                                                            (defining_process.nice_string(),f2py_module_path))
            f2py_load_path = (os.path.dirname(self.export_dir), 
                        '%s.matrix_%s_py'%( os.path.basename(self.export_dir), 
                                                self.process_dir_name(defining_process) ) )
            slha_card_path = pjoin(self.export_dir,'Cards','param_card.dat')
            if os.path.isdir(pjoin(self.export_dir,'SubProcesses','MadLoop5_resources')):
                madloop_resources_path = pjoin(self.export_dir,'SubProcesses','MadLoop5_resources')
            else:
                madloop_resources_path = ''
                
            MEAccessors.append(VirtualMEAccessor(
                defining_process, 
                f2py_load_path, 
                slha_card_path,
                madloop_resources_path=madloop_resources_path,
                mapped_pdgs = mapped_process_pdgs, 
                root_path=root_path
                ) )

        # Only allow overwriting accessors if processes were not grouped
        all_MEAccessors.add_MEAccessors(MEAccessors, 
                                          allow_overwrite = (not self.group_subprocesses) )

    def can_processes_be_integrated_together(self, processA, processB):
        """ Investigates whether processA and processB can be integrated together for this contribution."""

        # For now only make sure external masses are identical
        if processA.get_external_masses(self.model) != processB.get_external_masses(self.model):
            return False
        
        return True
    
    def get_integrands_for_process_map(self, process_map, model, run_card, all_MEAccessors, ME7_configuration):
        """ Returns all the integrands implementing this contribution for the specified process_map.
        The instance of MEAccessorDict is necessary so as to be passed to the integrand instances.
        """
        return [ ME7_interface.ME7Integrand(model, run_card,
                                       self.contribution_definition,
                                       process_map,
                                       self.topologies_to_processes,
                                       self.processes_to_topologies,
                                       all_MEAccessors,
                                       ME7_configuration)
               ]
        
    def get_integrands(self, *args):
        """ Returns all the integrands implementing this contribution.
        The *args are passed to the integrand instances.
        """
        
        # A list of processes maps we will distribute to each integrand
        integrand_process_maps = []

        # Regroup the general process map into smaller sub-process maps where one is guaranteed
        # that all processes in these submaps can be integrated together. 
        for process_key, (defining_process, mapped_processes) in self.get_processes_map().items():
            found_partner = False
            for integrand_process_map in integrand_process_maps:
                if all(self.can_processes_be_integrated_together(defining_process,
                                         proc[0]) for proc in integrand_process_map.values()):
                    integrand_process_map[process_key] = (defining_process, mapped_processes)
                    found_partner = True
                    break
            if not found_partner:
                integrand_process_maps.append({process_key:(defining_process, mapped_processes)})
            
        all_integrands = []
        for integrand_process_map in integrand_process_maps:
            all_integrands.extend(self.get_integrands_for_process_map(integrand_process_map, *args))
        
        return all_integrands
        
    def add_content_to_global_ME7_resources(self, global_ME7_dir, **opts):
        """ Possibly add content to the global ME7 resources pool."""
        pass

    def remove_superfluous_content(self):
        """ At the end of the export of this contributions, remove extra files not desired/necessary. """
        
        dirs_to_remove = [ pjoin(self.export_dir, 'Cards') ]
        
        files_to_remove = [ pjoin(self.export_dir, 'Source','make_opts') ]
        
        for dir in dirs_to_remove:
            if os.path.isdir(dir):
                shutil.rmtree(dir)
            elif os.path.islink(dir):
                os.remove(dir)                
        
        if file in files_to_remove:
            if os.path.isfile(file) or os.path.islink(file):
                os.remove(file)                

    def finalize(self, flaglist=[], interface_history=[]):
        """ Finalize the output of the code necessary for this contribution."""
        
        # WARNING: It would be ideal to use the same DHELAS output for all contributions.
        # This should be done eventually, but there are issues with coef_specs, complex vs real momenta,
        # etc.. So for this first test we stick with a local DHELAS for each contributions
        
        # This stores the wanted couplings which should be exported by the overall ME7 exporter.
        global_wanted_couplings = []
        
        # Handling of the model.
        if self.options['_model_v4_path']:
            logger.info('Copy %s model files to directory %s' % \
                            (os.path.basename(self.options['_model_v4_path']), self.export_dir))
            self.exporter.export_model_files(self.options['_model_v4_path'])
            self.exporter.export_helas(pjoin(self._mgme_dir,'HELAS'))        
        else:
            # wanted_lorentz are the lorentz structures which are actually used in the 
            # wavefunctions and amplitudes in these processes
            wanted_lorentz = self.all_matrix_elements.get_used_lorentz()
            if self.requires_its_own_model:
                wanted_couplings = self.all_matrix_elements.get_used_couplings()
            else:
                # This prevents 'convert_model' to write out the Source/MODEL directory
                wanted_couplings = []
                # This tells the ME7 exporter to link it to the global model instead
                global_wanted_couplings = self.all_matrix_elements.get_used_couplings()
            self.exporter.convert_model(self.model, wanted_lorentz, wanted_couplings)

        # Dedicated finalize function.
        self.exporter.finalize(self.all_matrix_elements, interface_history, self.options, flaglist)
        
        return global_wanted_couplings

    def link_global_ME7_resources(self, global_ME7_dir):
        """ Remove some of the local resources and instead link the global ones from the ME7 exporters."""
        
        # Link to the global model if necessary
        if not self.requires_its_own_model:
            ln(pjoin(global_ME7_dir, 'Source','MODEL'), starting_dir=pjoin(self.export_dir, 'Source'))
            ln(pjoin(global_ME7_dir, 'lib','libmodel.a'), starting_dir=pjoin(self.export_dir, 'lib'))

        # Link the card directory
        ln(pjoin(global_ME7_dir, 'Cards'), starting_dir=self.export_dir)

        # Use a global make_opt file
        ln(pjoin(global_ME7_dir, 'Source','make_opts'), starting_dir=pjoin(self.export_dir, 'Source'))

    def make_model_symbolic_link(self):
        """ Links some of the MODEL resources to the other folders exported within this contribution."""
        self.exporter.make_model_symbolic_link()
    
    def get_inverse_processes_map(self, force=False):
        """ Returns a dictionary with the keys (with ordered PDGs) of all processes present
        in this contribution and the single value being their corresponding "defining" process key.
        If force is True, then the map will be reconstructed entirely and the cache updated."""
        
        # Sort PDGs only when processes have been grouped
        sort_PDGs = self.group_subprocesses
        
        if not self.all_matrix_elements.get_matrix_elements():
            raise MadGraph5Error("The processes reverse map can only be built *after*"+
                                                 " the matrix element have been exported.")
        
        if hasattr(self, 'inverse_processes_map') and not force:
            return self.inverse_processes_map
            
        process_map = self.get_processes_map(force=force)
        inverse_map = {}
        for process_key, (defining_process, mapped_processes) in process_map.items():
            defining_process_key = ProcessKey(defining_process,sort_PDGs=sort_PDGs).get_canonical_key()
            if defining_process_key in inverse_map:
                raise MadGraph5Error("The following process key appears twice in "+
                                                  "the processes map: %s"%str(process_key))
            inverse_map[defining_process_key] = process_key
            for mapped_process in mapped_processes:
                mapped_process_key = ProcessKey(mapped_process,sort_PDGs=sort_PDGs).get_canonical_key()
                if mapped_process_key in inverse_map:
                    raise MadGraph5Error("The following mapped process key appears twice "+
                                             "in the processes map: %s"%mapped_process_key)                    
                inverse_map[mapped_process_key] = process_key
        
        self.inverse_processes_map = inverse_map
        
        return self.inverse_processes_map
                    

    def get_processes_map(self, force=False):
        """ Returns a dictionary with keys obtained from ProcessKey and values being a tuple with format
                 (defining_process, [mapped processes])
            where defining_process is an actual instance of Process and [mapped processes] is a list 
            of processes instances mapped to that defining process.
        """
        if not self.amplitudes:
            return {}
        
        # Attempt reusing a cached version
        if hasattr(self, 'processes_map') and not force:
            if (self.processes_map[0]['had_amplitudes'] == bool(self.amplitudes)) and \
               (self.processes_map[0]['had_matrix_elements'] == bool(
                                          self.all_matrix_elements.get_matrix_elements())):
                return self.processes_map[1]
        
        all_defining_procs = [amp.get('process') for amp in self.amplitudes]
        
        # The values are (defining_process_instance, list_of_mapped_process_instances)
        all_defining_procs = dict( ( ProcessKey(proc,sort_PDGs=False).get_canonical_key()
                                            , (proc, []) ) for proc in all_defining_procs)
        
        if not self.all_matrix_elements.get_matrix_elements():
            # Cache the process map
            self.processes_map = ({
                'had_amplitudes':bool(self.amplitudes),
                'had_matrix_elements':bool(self.all_matrix_elements.get_matrix_elements())},
                all_defining_procs)
            return all_defining_procs
        
        for me in self.all_matrix_elements.get_matrix_elements():
            # List all the processes in this Matrix eElement
            all_procs_pdgs = dict((ProcessKey(proc,sort_PDGs=False).get_canonical_key(), proc) 
                                                                for proc in me.get('processes'))
            # Add the attribute 'has_mirror_process' to these processes instance
            for proc in all_procs_pdgs.values():
                proc.set('has_mirror_process', me.get('has_mirror_process'))
            # Look through all the processes identified in the amplitudes
            for proc in all_defining_procs:
                # Make sure it is found in the processes mapped in this Matrix element
                if not proc in all_procs_pdgs:
                    continue
                # Also add 'has_mirror_process' to the defining process instance
                all_defining_procs[proc][0].set('has_mirror_process', me.get('has_mirror_process'))
                # Now add them to list of mapped processes for this Matrix element
                all_defining_procs[proc][1].extend([all_procs_pdgs[p] for p in all_procs_pdgs
                                            if p!=proc and p not in all_defining_procs[proc][1]])



        # Cache the process map
        self.processes_map = ({
                'had_amplitudes':bool(self.amplitudes),
                'had_matrix_elements':bool(self.all_matrix_elements.get_matrix_elements())},
                all_defining_procs)
        
        return all_defining_procs        

    def process_dir_name(self, process):
        """ Given a specified process, return the directory name in which it is output."""
        # This function is essentially useful in order to be shadowed by daughter classes     
        return process.shell_string()       
    
    def compile(self):
        """ Compiles the f2py shared library to provide easy access to this contribution."""
        
        processes_map = self.get_processes_map()
        i_process=0
        for process_key, (defining_process, mapped_processes) in processes_map.items():
            i_process += 1
            name = self.process_dir_name(defining_process)
            Pdir = pjoin(self.export_dir, 'SubProcesses', 'P%s'%name)
            logger.debug("[%d/%d] Compiling %s in '%s'..."%(
                        i_process, len(processes_map.keys()),
                        defining_process.nice_string().replace('Process','process'),
                        os.path.basename(Pdir) ))
            if not os.path.isdir(Pdir):
                raise MadGraph5Error("The expected subprocess directory %s could not be found."%Pdir)
            misc.compile(arg=['matrix_%s_py.so'%name, 'MENUM=_%s_'%name], cwd=Pdir)
            if not os.path.isfile(pjoin(Pdir, 'matrix_%s_py.so'%name)):
                raise InvalidCmd("The f2py compilation of subprocess '%s' failed.\n"%Pdir+
                                 "Try running 'make matrix2py.so' by hand in this directory.")
            ln(pjoin(Pdir, 'matrix_%s_py.so'%name), starting_dir=self.export_dir)
    
    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        contribution."""
        GREEN = '\033[92m'
        ENDC = '\033[0m'        
        return GREEN+'  %s'%defining_process.nice_string(print_weighted=False).\
                                                               replace('Process: ','')+ENDC
    
    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this contribution."""
        return []
    
    def nice_string(self, format=0):
        """ Nice string representation of self."""
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        res = ['%-30s:   %s'%('contribution_type',type(self))]
        res.extend([self.contribution_definition.nice_string()])
        if not self.topologies_to_processes is None:
            res.append('%-30s:   %d'%('Number of topologies', 
                                                    len(self.topologies_to_processes.keys())))
        res.extend(self.get_additional_nice_string_printout_lines())
        if self.amplitudes and not self.all_matrix_elements.get_matrix_elements():
            if format < 1:
                res.append('Amplitudes generated for the following processes: %d'%
                                                                    len(self.amplitudes) )
            else:
                res.append('Amplitudes generated for the following processes:')
                for amp in self.amplitudes:
                    res.append(GREEN+'  %s'%amp.get('process').nice_string(print_weighted=False).\
                                                                            replace('Process: ','')+ENDC)
        elif self.amplitudes and self.all_matrix_elements.get_matrix_elements():
            processes_map = self.get_processes_map()
            if format < 1:
                res.append('Generated and mapped processes for this contribution: %d (+%d mapped)'%
                           ( len(processes_map.keys()),
                             len(sum([v[1] for v in processes_map.values()],[])) ) )
            else:
                res.append('Generated and mapped processes for this contribution:')
                for process_key, (defining_process, mapped_processes) in processes_map.items():
                    res.append(self.get_nice_string_process_line(process_key, defining_process, format=format))                
                    for mapped_process in mapped_processes:
                        res.append(BLUE+u'   \u21b3  '+mapped_process.nice_string(print_weighted=False)\
                                                                            .replace('Process: ','')+ENDC)
        else:
            res.append(BLUE+'No amplitudes generated yet.'+ENDC)
            
        return '\n'.join(res).encode('utf-8')
        
    def generate_amplitudes(self, force=False):
        """ Generates the relevant amplitudes for this contribution."""
        
        # First check if the amplitude was not already generated
        if self.amplitudes and not force:
            return

        myproc = self.MultiProcessClass(self.contribution_definition.process_definition,
                    collect_mirror_procs = self.collect_mirror_procs,
                    ignore_six_quark_processes = self.ignore_six_quark_processes,
                    optimize=self.optimize, diagram_filter=self.diagram_filter)
    
        for amp in myproc.get('amplitudes'):
            if amp not in self.amplitudes:
                self.amplitudes.append(amp)
            else:
                logger.warning('Duplicate process found in contribution '+
                               '%s. Sanity check needed.'%str(type(self)))
 
class Contribution_B(Contribution):
    """ Implements the handling of a Born-type of contribution."""
    
    def __init__(self, contribution_definition, cmd_interface, **opts):
        super(Contribution_B,self).__init__(contribution_definition, cmd_interface, **opts)
        self.ignore_six_quark_processes = self.options['ignore_six_quark_processes'] if \
            "ignore_six_quark_processes" in self.options else []
        if 'diagram_filter' in opts:
            self.diagram_filter = opts['diagram_filter']
        if 'optimize' in opts:
            self.optimize = opts['optimize']
        self.supports_decay_chain       = True

class Contribution_LIB(Contribution_B):
    """ Implements the handling of loop-induced Born-type of contribution."""
    pass

class Contribution_R(Contribution):
    """ Implements the handling of a single real-emission type of contribution."""

    def __init__(self, contribution_definition, cmd_interface, **opts):
        """ Instantiates a real-emission contribution with additional attributes."""
        
        super(Contribution_R, self).__init__(contribution_definition, cmd_interface, **opts)
        # Add a default empty list of counterterms
        self.counterterms = None

    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this contribution."""
        res = []
        if self.counterterms:
            res.append('%-30s:   %d'%('Number of local counterterms', 
               len([1 for CT in sum(self.counterterms.values(),[]) if CT.is_singular()]) ))
        return res
        
    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        contribution."""
        
        GREEN = '\033[92m'
        ENDC = '\033[0m'        
        res =  GREEN+'  %s'%defining_process.nice_string(print_weighted=False).\
                                                               replace('Process: ','')+ENDC

        if not self.counterterms:
            return res                                                               

        if format < 2:
            if process_key in self.counterterms:
                res += ' | %d local counterterms'%len([ 1 for CT in 
                                      self.counterterms[process_key] if CT.is_singular() ])
            else:
                res += ' | 0 local counterterm'
                
        else:
            long_res = [' | with the following local counterterms:']
            for CT in self.counterterms[process_key]:
                if CT.is_singular():
                    if format==2:
                        long_res.append( '   | %s'%CT.get_singular_structure_string(
                                            print_n=True, print_pdg=False, print_state=False )  )
                    elif format==3:
                        long_res.append( '   | %s'%CT.get_singular_structure_string(
                                            print_n=True, print_pdg=True, print_state=True )  )
                    elif format>3:
                        long_res.append( '   | %s'%str(CT))
            res += '\n'.join(long_res)

        return res

    def generate_all_counterterms(self, group_processes=True):
        """ Generate all counterterms associated to the processes in this contribution."""
        
        self.counterterms = {}
        all_integrated_counterterms = []
        for process_key, (defining_process, mapped_processes) in self.get_processes_map().items():
            local_counterterms, integrated_counterterms = \
                                 self.IR_subtraction.get_all_counterterms(defining_process)
            self.counterterms[process_key] = local_counterterms
            # Add a additional information about the integrated_counterterms that will
            # be necessary for their evaluations in the corresponding lower-multiplicity
            # contribution
            
            # First, a dictionary with keys being a one- or two- tuple of the initial
            # states part of this mapping and the values being the number of mapped
            # processes that possess these initial states
            flavors_combinations = []
            # Notice that the defining flavor combination will always be placed first
            for process in [defining_process,]+mapped_processes:
                initial_pdgs = process.get_initial_ids()
                final_pdgs = tuple(process.get_final_ids_after_decay())
                flavors_combinations.append( (
                    tuple(initial_pdgs),
                    final_pdgs
                ) )
                if process.get('has_mirror_process'):
                    flavors_combinations.append( (
                        tuple(reversed(initial_pdgs)),
                        final_pdgs
                    ) )

            for integrated_counterterm in integrated_counterterms:
                all_integrated_counterterms.append({
                    'integrated_counterterm'    :   integrated_counterterm,
                    'flavors_combinations'      :   flavors_combinations
                })

        return all_integrated_counterterms
    
    @classmethod
    def remove_counterterms_with_no_reduced_process(cls, all_MEAccessors, counterterms):
        """ Given the list of available reduced processes encoded in the MEAccessorDict 'all_MEAccessors'
        given in argument, remove all the counterterms whose underlying reduced process does not exist."""


        for counterterm in list(counterterms):
            # Of course don't remove any counterterm that would be the real-emission itself.
            if not counterterm.is_singular():
                continue
            try:
                all_MEAccessors.get_MEAccessor(counterterm.process)
            except MadGraph5Error:
                # This means that the reduced process could not be found and
                # consequently, the corresponding counterterm must be removed.
                # Example: C(5,6) in e+ e- > g g d d~
                counterterms.remove(counterterm)

    def export(self, *args, **opts):
        """ Overloads export so as to export subtraction currents as well."""
        ret_value = super(Contribution_R, self).export(*args, **opts)

        # Fish out the group_processes option as it could be used when attempting to
        # generate all currents.
        
        if 'group_processes' in opts:
            group_processes = opts['group_processes']
        else:
            group_processes = True

        integrated_counterterms = self.generate_all_counterterms(
                                                          group_processes=group_processes)
      
        # Add the integrated counterterms to be passed to the exporter
        ret_value.update({'integrated_counterterms': integrated_counterterms})
        
        return ret_value

    def get_all_necessary_subtraction_currents(self, all_MEAccessors):
        """ Given the counterterms in place and the currents already accessible in the 
        all_MEAccessors, return what subtraction currents are needed."""
        
        all_currents = []
        for process_key, counterterms in self.counterterms.items():
            for current in subtraction.IRSubtraction.get_all_currents(counterterms):
                # Retain only a single copy of each needed current.
                # We must remove the leg information since this is information is irrelevant
                # for the selection of the hard-coded current implementation to consider.
                copied_current = current.get_copy(('squared_orders','singular_structure'))
                copied_current.discard_leg_numbers()
                if copied_current not in all_currents:
                    all_currents.append(copied_current)

        # Now further remove currents that are already in all_MEAccessors
        all_currents = [current
                        for current in all_currents
                        if current.get_key().get_canonical_key() not in all_MEAccessors]
        
        return all_currents

    @classmethod
    def add_current_accessors(cls, model, all_MEAccessors, root_path, currents_to_consider):
        """  Generates and add all subtraction current accessors to the MEAccessorDict."""

        # Now generate the computer code and exports it on disk for the remaining new currents
        current_exporter = subtraction.SubtractionCurrentExporter(model, root_path)
        mapped_currents = current_exporter.export(currents_to_consider)
        
        logger.debug("The following subtraction current implementation are exported:\n%s"%\
                    ( '\n'.join((" > %-35s for representative current '%s'"%("'%s'"%class_name, str(current_properties['defining_current']))
                                if class_name!='DefaultCurrentImplementation' else " > %-35s for a total of %d currents."%
                                ("'DefaultCurrentImplementation'",len(current_properties['mapped_process_keys'])))
                            for  (module_path, class_name, _), current_properties in mapped_currents.items() ) ))

        all_current_accessors = []  
        # Finally instantiate the CurrentAccessors corresponding to all current implementations identified and needed
        for (module_path, class_name, _), current_properties in mapped_currents.items():
            all_current_accessors.append(VirtualMEAccessor(
                current_properties['defining_current'],
                module_path,
                class_name,
                '%s.subtraction_current_implementations_utils'%current_exporter.main_module_name, 
                current_properties['instantiation_options'], 
                mapped_process_keys=current_properties['mapped_process_keys'], 
                root_path=root_path,
                model=model
            ))
        
        all_MEAccessors.add_MEAccessors(all_current_accessors)

    def add_ME_accessors(self, all_MEAccessors, root_path):
        """ Adds all MEAccessors for the matrix elements and currents generated as part of this contribution."""
        
        # Get the basic accessors for the matrix elements
        super(Contribution_R, self).add_ME_accessors(all_MEAccessors, root_path)

        for process_key, counterterms in self.counterterms.items():
            # Remove counterterms with non-existing underlying Born processes
            self.remove_counterterms_with_no_reduced_process(all_MEAccessors, counterterms)

        # Obtain all necessary currents
        currents_to_consider = self.get_all_necessary_subtraction_currents(all_MEAccessors)

        self.add_current_accessors(self.model, all_MEAccessors, root_path, currents_to_consider)
     
    def get_integrands_for_process_map(self, process_map, model, run_card, all_MEAccessors, ME7_configuration):
        """ Returns all the integrands implementing this contribution for the specified process_map.
        The instance of MEAccessorDict is necessary so as to be passed to the integrand instances.
        """
        
        relevant_counterterms = {}
        for process_key in process_map:
            relevant_counterterms[process_key] = self.counterterms[process_key]

        return [ ME7_interface.ME7Integrand(model, run_card,
                                       self.contribution_definition,
                                       process_map,
                                       self.topologies_to_processes,
                                       self.processes_to_topologies,
                                       all_MEAccessors,
                                       ME7_configuration,
                                       counterterms=relevant_counterterms)
               ]
        
class Contribution_RR(Contribution_R):
    """ Implements the handling of a double real-emission type of contribution."""
    pass

class Contribution_RRR(Contribution_R):
    """ Implements the handling of a triple real-emission type of contribution."""
    pass

class Contribution_V(Contribution):
    """ Implements the handling of a virtual type of contribution."""

    def __init__(self, contribution_definition, cmd_interface, **opts):
        """ Bring in the couple of modifications necessary for this type of contributions."""
        super(Contribution_V,self).__init__(contribution_definition, cmd_interface, **opts)
        # Make sure to adjust the MultiProcessClass to be used
        self.MultiProcessClass          = loop_diagram_generation.LoopMultiProcess
        self.output_type                = 'madloop'
        
        # Store values being integration counterterms and their properties (as a dictionary), 
        # with the ProcessKey they are attached to as keys of this attribute's dictionary.
        self.integrated_counterterms    = {}
    
    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this contribution."""
        res = []
        if self.integrated_counterterms:
            res.append('%-30s:   %d'%('Nb. of integrated counterterms', 
                                      len(sum(self.integrated_counterterms.values(),[]))))
        return res
    
    @classmethod
    def remove_counterterms_with_no_reduced_process(cls, all_MEAccessors, CT_properties_list):
        """ Given the list of available reduced processes encoded in the MEAccessorDict 'all_MEAccessors'
        given in argument, remove all the counterterms whose underlying reduced process does not exist."""

        for CT_properties in list(CT_properties_list):
            counterterm = CT_properties['integrated_counterterm']
            # Of course don't remove any counterterm that would be the real-emission itself.
            if not counterterm.is_singular():
                continue
            try:
                all_MEAccessors.get_MEAccessor(counterterm.process)
            except MadGraph5Error:
                # This means that the reduced process could not be found and
                # consequently, the corresponding counterterm must be removed.
                # Example: C(5,6) in e+ e- > g g d d~
                CT_properties_list.remove(CT_properties)
    
    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        contribution."""
        GREEN = '\033[92m'
        ENDC = '\033[0m'        
        res = GREEN+'  %s'%defining_process.nice_string(print_weighted=False).\
                                                               replace('Process: ','')+ENDC

        if not self.integrated_counterterms:
            return res

        if format<2:
            if process_key in self.integrated_counterterms:
                res += ' | %d integrated counterterms'%len(self.integrated_counterterms[process_key])
            else:
                res += ' | 0 integrated counterterm'
                
        else:
            long_res = [' | with the following integrated counterterms:']
            for CT_properties in self.integrated_counterterms[process_key]:
                CT = CT_properties['integrated_counterterm']
                if format==2:
                    long_res.append( '   | %s'%CT.get_singular_structure_string(
                                        print_n=True, print_pdg=False, print_state=False )  )
                elif format==3:
                    long_res.append( '   | %s'%CT.get_singular_structure_string(
                                        print_n=True, print_pdg=True, print_state=True )  )
                elif format==4:
                    long_res.append( '   | %s'%str(CT))
                elif format>4:
                    long_res.append( '   | %s'%str(CT))
                    for key, value in CT_properties.items():
                        if not key in ['integrated_counterterm', 'matching_process_key']:
                            long_res.append( '     + %s : %s'%(key, str(value)))

            res += '\n'.join(long_res)

        return res

    def generate_matrix_elements(self, group_processes=True):
        """Generate the Helas matrix elements before exporting. Uses the main function argument 
        'group_processes' to decide whether to use group_subprocess or not."""

        cpu_time1 = time.time()
        ndiags = 0

        generation_mode = {'optimized_output': self.options['loop_optimized_output']}

        self.all_matrix_elements = loop_helas_objects.LoopHelasProcess(
            self.amplitudes,
            compute_loop_nc = True,
            matrix_element_opts = generation_mode,
            optimized_output = self.options['loop_optimized_output'],
            combine_matrix_elements=group_processes
        )
        ndiags = sum([len(me.get('diagrams')) for me in self.all_matrix_elements.get_matrix_elements()])
        
        # assign a unique id number to all process
        uid = 0 
        for me in self.all_matrix_elements.get_matrix_elements():
            uid += 1 # update the identification number
            me.get('processes')[0].set('uid', uid)

        cpu_time2 = time.time()
        return ndiags, cpu_time2 - cpu_time1        

    def add_content_to_global_ME7_resources(self, global_ME7_dir, **opts):
        """ Pass the MadLoopParams.dat card to the global ME7 resources."""
        super(Contribution_V, self).add_content_to_global_ME7_resources(global_ME7_dir, **opts)
        
        destination = pjoin(global_ME7_dir,'Cards','MadLoopParams.dat')
        if not os.path.exists(destination):
            shutil.copyfile(pjoin(self.export_dir,'Cards','MadLoopParams.dat'), destination)

    def set_helas_model(self):
        """ Instantiate properly the helas model """

        assert self.exporter.exporter == 'v4'
        assert (not self.options['_model_v4_path'])
        self.helas_model = helas_call_writers.FortranUFOHelasCallWriter(self.model)

    def get_all_necessary_integrated_currents(self, all_MEAccessors):
        """ Given the list of currents already available encoded in the all_MEAccessors,
        generate all the integrated currents that must be exported."""
        
        all_currents = []
        for CT_properties in sum(self.integrated_counterterms.values(),[]):
            counterterm = CT_properties['integrated_counterterm']
            # For now we only support a basic integrated counterterm which is not broken
            # down in subcurrents but contains a single CountertermNode with a single
            # current in it that contains the whole singular structure describing this
            # integrated counterterm.
            if len(counterterm.nodes) != 1 or len(counterterm.nodes[0].nodes) != 0:
                raise MadGraph5Error(
                    """For now, MadEvent7 only support simple integrated
                    counterterms that consists of single current encompassing the full
                    singular structure that must be analytically integrated over.""")
            integrated_current = counterterm.nodes[0].current
            assert(isinstance(integrated_current, subtraction.IntegratedCurrent))
            
            # Retain only a single copy of each needed current.
            # We must remove the leg information since this is information is irrelevant
            # for the selection of the hard-coded current implementation to consider.
            copied_current = integrated_current.get_copy(('squared_orders','singular_structure'))
            copied_current.discard_leg_numbers()
            if copied_current not in all_currents:
                all_currents.append(copied_current)

        # Now further remove currents that are already in all_MEAccessors
        all_currents = [current for current in all_currents if 
                        current.get_key().get_canonical_key() not in all_MEAccessors]
        
        return all_currents

    @classmethod
    def add_current_accessors(cls, model, all_MEAccessors, 
                                           root_path, currents_to_consider, *args, **opts):
        """Generate and add all integrated current accessors to the MEAccessorDict.
        For now we can recycle the implementation of the Contribution_R class. """

        return Contribution_R.add_current_accessors(model, 
                          all_MEAccessors, root_path, currents_to_consider, *args, **opts)

    def get_integrands_for_process_map(self, process_map, model, run_card, all_MEAccessors, ME7_configuration):
        """ Returns all the integrands implementing this contribution for the specified process_map.
        The instance of MEAccessorDict is necessary so as to be passed to the integrand instances.
        """
        
        relevant_counterterms = {}
        if self.integrated_counterterms:
            for process_key in process_map:
                relevant_counterterms[process_key] = self.integrated_counterterms[process_key]

        return [ ME7_interface.ME7Integrand(model, run_card,
                                       self.contribution_definition,
                                       process_map,
                                       self.topologies_to_processes,
                                       self.processes_to_topologies,
                                       all_MEAccessors,
                                       ME7_configuration,
                                       integrated_counterterms=relevant_counterterms)
               ]

    def add_ME_accessors(self, all_MEAccessors, root_path):
        """ Adds all MEAccessors for the matrix elements and currents generated as part 
        of this contribution."""
        
        # Get the basic accessors for the matrix elements
        super(Contribution_V, self).add_ME_accessors(all_MEAccessors, root_path)

        for process_key, CT_properties in self.integrated_counterterms.items():
            # Remove integrated counterterms with non-existing underlying Born processes
            self.remove_counterterms_with_no_reduced_process(all_MEAccessors, CT_properties)

        # Obtain all necessary currents
        currents_to_consider = self.get_all_necessary_integrated_currents(all_MEAccessors)

        self.add_current_accessors(self.model, all_MEAccessors, root_path, currents_to_consider)

    def add_integrated_counterterm(self, integrated_CT_properties):
        """ Virtual contributions can receive integrated counterterms and they will
        be stored in the attribute list self.integrated_counterterms."""
        
        # Extract quantities from integrated_counterterm_properties
        integrated_counterterm = integrated_CT_properties['integrated_counterterm']
        # flavors_combinations = integrated_CT_properties['flavors_combinations']
        
        # Sort PDGs only when processes have been grouped
        sort_PDGs = self.group_subprocesses
        
        # First access the inverse list of processes key, i.e. a dictionary with process keys
        # for all processes in this contribution with the value being the corresponding
        # defining process key
        inverse_processes_map = self.get_inverse_processes_map()
        
        # Make sure the dictionary self.integrated_counterterms is initialized
        # with one key for each mapped process
        if len(self.integrated_counterterms)==0:
            processes_map = self.get_processes_map()
            for key in processes_map.keys():
                self.integrated_counterterms[key] = []

        # Obtain the ProcessKey of the reduced process of the integrated counterterm
        # Use ordered PDGs since this is what is used in the inverse_processes_map
        # Also verwrite some of the process properties to make it match the 
        # loop process definitions in this contribution
        counterterm_reduced_process_key = ProcessKey(integrated_counterterm.process,
            sort_PDGs=sort_PDGs, n_loops=1, NLO_mode='virt', 
            perturbation_couplings=self.contribution_definition.correction_couplings)\
                                                                       .get_canonical_key()
        
        # misc.sprint(inverse_processes_map.keys(), len(inverse_processes_map.keys()))
        # misc.sprint(counterterm_reduced_process_key)
        if counterterm_reduced_process_key not in inverse_processes_map:
            # The reduced process of the integrated counterterm is not included in this
            # contribution so we return here False, letting the ME7 exporter know that
            # we cannot host this integrated counterterm.
            return False
        
        # Now we can simply add this integrated counterterm in the group of the
        # defining process to which the reduced process is mapped
        defining_key = inverse_processes_map[counterterm_reduced_process_key]
        
        integrated_counterterm_properties = dict(integrated_CT_properties)
        
        # also store what was the matching virtual process key, not the defining one
        # Fow not, it is not useful, but it can be added later if necessary.
        # integrated_counterterm_properties['matching_process_key'] = counterter_reduced_process_key
        
        self.integrated_counterterms[defining_key].append(integrated_counterterm_properties)
        
        # Let the ME7 exporter that we could successfully host this integrated counterterm
        return True

    def process_dir_name(self, process):
        """ Given a specified process, return the directory name in which it is output."""
        
        # For MadLoop ME's there is extra prefixes to avoid symbol and resource file clashes.
        return "%d_%s_%s"%(process.get('id'), process.get('uid'),
                                                      process.shell_string(print_id=False))


    def generate_code(self):
        """ Assuming the Helas Matrix Elements are now generated, we can write out the corresponding code."""

        matrix_elements = self.all_matrix_elements.get_matrix_elements()
            
        cpu_time_start = time.time()

        # Pick out the matrix elements in a list
        matrix_elements = self.all_matrix_elements.get_matrix_elements()
        # MadLoop standalone fortran output
        calls = 0
        for me in matrix_elements:
            # Choose the group number to be the unique id so that the output prefix
            # is nicely P<proc_id>_<unique_id>_
            calls = calls + self.exporter.generate_loop_subprocess(me, 
                self.helas_model,
                group_number = me.get('processes')[0].get('uid'),
                proc_id = None,
                config_map=None,
                unique_id=me.get('processes')[0].get('uid'))

        # If all ME's do not share the same maximum loop vertex rank and the
        # same loop maximum wavefunction size, we need to set the maximum
        # in coef_specs.inc of the HELAS Source. The SubProcesses/P* directory
        # all link this file, so it should be properly propagated
        if self.options['loop_optimized_output'] and len(matrix_elements)>1:
            max_lwfspins = [m.get_max_loop_particle_spin() for m in \
                                                            matrix_elements]
            max_loop_vert_ranks = [me.get_max_loop_vertex_rank() for me in \
                                                            matrix_elements]
            if len(set(max_lwfspins))>1 or len(set(max_loop_vert_ranks))>1:
                self.exporter.fix_coef_specs(max(max_lwfspins),\
                                                   max(max_loop_vert_ranks))
        
        return calls, time.time() - cpu_time_start


class Contribution_RV(Contribution_R, Contribution_V):
    """ Implements the handling of a real-virtual type of contribution."""
    pass

class ContributionList(base_objects.PhysicsObjectList):
    """ A container for storing a list of contributions."""
    
    contributions_natural_order = [
        ('LO',    (Contribution_B, Contribution_LIB) ),
        ('NLO',   (Contribution_R, Contribution_V) ),
        ('NNLO',  (Contribution_RR, ) ),
        ('NNNLO', (Contribution_RRR, ) )
    ]
    
    def is_valid_element(self, obj):
        """Test if object obj is a valid instance of Contribution."""
        return isinstance(obj, Contribution)
    
    def get_contributions_of_order(self, correction_order):
        """ Returns a list of all contributions of a certain correction_order in argument."""
        return ContributionList([contrib for contrib in self if
                contrib.contribution_definition.correction_order==correction_order])

    def get_contributions_of_type(self, correction_classes):
        """ Returns a list of all contributions that are direct instances of certain classes."""
        if not isinstance(correction_classes, tuple):
            if isinstance(correction_classes, list):
                correction_classes = tuple(correction_classes)
            else:
                correction_classes = (correction_classes,)                
        return ContributionList([contrib for contrib in self if isinstance(contrib, correction_classes)])

    def nice_string(self, format=0):
        """ A nice representation of a list of contributions. 
        We can reuse the function from ContributionDefinitions."""
        return base_objects.ContributionDefinitionList.contrib_list_string(self, format=format)

    def sort_contributions(self):
        """ Sort contributions according to the order dictated by the class_attribute contributions_natural_order"""
        
        new_order = []
        for correction_order, contribution_types in self.contributions_natural_order:
            for contrib_type in contribution_types:
                selected_contribs = self.get_contributions_of_order(correction_order).\
                                        get_contributions_of_type(contrib_type)
                new_order.extend(selected_contribs)
                for contrib in selected_contribs:
                    self.pop(self.index(contrib))
    
        # Finally all remaining contributions of unknown types
        new_order.extend(self)
        self[:] = new_order

    def apply_method_to_all_contribs(self, method, log=None, method_args = [], method_opts = {}):
        """ Apply a given method to all contributions nicely sorted."""
        
        # Keep track of the return values
        return_values = []
        
        remaining_contribs = list(self)
        for correction_order, contribution_types in self.contributions_natural_order:
            for contrib_type in contribution_types:
                selected_contribs = self.get_contributions_of_order(correction_order).\
                                        get_contributions_of_type(contrib_type)
                if selected_contribs:
                    if log: logger.info('%s the %d contribution%s of type %s (%s)...'%(log,
                        len(selected_contribs), 's' if len(selected_contribs)>1 else '',
                        contrib_type.__name__, correction_order))
                for i, contrib in enumerate(selected_contribs):
                    if log: logger.info('%s (%s) %d/%d'%
                        (contrib_type.__name__, correction_order, i+1, len(selected_contribs)))
                    try:
                        contrib_function = getattr(contrib, method)
                    except AttributeError:
                        raise MadGraph5Error("The contribution\n%s\n does not have function '%s' defined."%(
                                                                                    contrib.nice_string(), method))
                    return_values.append((contrib, contrib_function(*method_args, **method_opts)))
                    remaining_contribs.pop(remaining_contribs.index(contrib))

        if remaining_contribs:
            if log: logger.info('%s the %d contribution%s of unknown type...'%(log,
                        len(remaining_contribs), 's' if len(remaining_contribs)>1 else ''))                    
        for i, contrib in enumerate(remaining_contribs):
            if log: logger.info('%s %d/%d'%(contrib_type.__name__, i+1, len(remaining_contribs)))
            try:
                contrib_function = getattr(contrib, method)
            except AttributeError:
                raise MadGraph5Error("The contribution\n%s\n does not have function '%s' defined."%(
                                                                            contrib.nice_string(), method))
            return_values.append((contrib, contrib_function(*method_args, **method_opts)))
        
        return return_values


# MEAccessor classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Accessor class.
# For instance 'Unknown' can be mapped to a user-defined class.
# Notice that this map must be placed after the MEAccessor daughter classes have been declared.
MEAccessor_classes_map = {'PythonAccessor': F2PYMEAccessor,
                          'CurrentAccessor': SubtractionCurrentAccessor,
                          'Unknown': None}

# Contribution classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own class of Contribution.
# Notice that this Contribution must be placed after all the Contribution daughter classes have been declared.
Contribution_classes_map = {'Born': Contribution_B,
                            'LoopInduced_Born': Contribution_LIB,
                            'Virtual': Contribution_V,
                            'SingleReals': Contribution_R,
                            'DoubleReals': Contribution_RR,
                            'TripleReals': Contribution_RRR,
                            'Unknown': None}