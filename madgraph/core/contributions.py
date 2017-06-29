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
pjoin = os.path.join

import madgraph.core.base_objects as base_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.loop.loop_diagram_generation as loop_diagram_generation
import madgraph.core.helas_objects as helas_objects
import madgraph.loop.loop_helas_objects as loop_helas_objects
import madgraph.iolibs.group_subprocs as group_subprocs
import madgraph.iolibs.helas_call_writers as helas_call_writers
import madgraph.iolibs.export_v4 as export_v4
import madgraph.various.misc as misc
import madgraph.interface.ME7_interface as ME7_interface
from madgraph import InvalidCmd, MadGraph5Error
from madgraph.iolibs.files import cp, ln, mv

logger = logging.getLogger('madgraph.contributions')

################################################################################
# MEAccessor mother class
################################################################################
class VirtualMEAccessor(object):
    """ A class wrapping the access to one particular MatrixEleemnt."""
    
    def __init__(self, process, n_loops=None, mapped_pdgs = [], root_path='', **opts):
        
        """ Initialize an MEAccessor from a process, allowing to possibly overwrites n_loops.
        and mapped_pdgs to indicate other list of pdgs corresponding to processes that can be 
        obtained by this same one.
        root_path is the process output root path from which all other paths are relative."""
        
        # Define the attributes for this particular ME.
        self.pdgs_combs  = [ ( tuple(process.get_initial_ids()), tuple(process.get_final_ids()) ) ] + mapped_pdgs
        if n_loops:
            self.n_loops = n_loops
        else:
            self.n_loops = 1 if process.get('perturbation_couplings') else 0
        self.perturbations = process.get('perturbation_couplings')
        
        if not os.path.isabs(root_path):
            raise MadGraph5Error("The initialization of VirtualMEAccessor necessitates an "+
                                                        "absolute path for root_path.")
        self.root_path = root_path
        
        # This permutation will determine how to interpret user-inputs during the call
        self.permutation = []
        # And in case the matrix-element evaluation depends on the process employed, 
        # (this can be the case when grouping several process under the same f2py module)
        # The MEAccessorDict will also specify this one.
        self.process_pdgs = tuple([])
        # User inputs to be used during the evaluation:
        self.PS_point         = []
        self.squared_orders   = {}
        self.spin_correlation = tuple([])
        self.color_connection = tuple([])
        self.hel_config       = tuple([])
        
    def get_canonical_key_value_pairs(self):
        """ Return all the canonical keys that should point to this ME. 
        This is intended to be used for the MEAccessorDict."""
        
        key_value_pairs = []
        for pdgs_comb in self.pdgs_combs:
            sorted_pdgs = (tuple(sorted(pdgs_comb[0])), tuple(sorted(pdgs_comb[1])))
            key = ( sorted_pdgs, self.n_loops, tuple(sorted(self.perturbations)) )
            value = (self, pdgs_comb)
            key_value_pairs.append( (key, value) )

        return key_value_pairs
    
    def __call__(self, *args, **opts):
        """ Make it explicitely virtual since this must be implemented by the daughter classes."""
        raise NotImplementedError("The __call__ function must be implemented in the daughter "+
                                  " classes of VirtualMEAccessor.")

    def apply_permutations(self, PS_point,
                                 squared_orders = {}, 
                                 spin_correlation = tuple([]), 
                                 color_connection = tuple([]), 
                                 hel_config = tuple([])):
        """ Return the matrix element evaluation, provided its arguments, whost format is:
        
        -> Spin correlation defined as a 2-tuple of lists of (4-vectors, index), for instance
             ( [(k1, 2), (k3, 3)], [(k1, 2), (k4, 1)] )
          Indicates that the vector k1 must be used on both side of the amplitude to saturate 
          the lorentz indices of leg #2 whereas k3 should be used to saturated leg #3 on the left
          and k4 leg#1 on the right.
          Notes that k_i here are simply 4-tuples of floats.

        -> Color connections are specified by providing what SU(3) generator chain to apply in 
          between the amplitude product. The format is a 2-tuple of lists providing first the outermost
          fundamental indices and then the summed adoint onces. For instance
             ( [3, 6], [-1, -2, -1, -2] )
          Indicates:
            T^a_(3 i) T^b_(i j) T^a_(j k) T^b_(k 6)  
        
        -> Helicity configuration to be considered, this can be a tuple like
            (-1, -1, +1 ,-1, +1)
          For the helicity configuration --+-+
        
        Empty arguments imply an explicit summation.
        """

        # Only take care of permuting the inputs if necessary. The daughter class will take care of 
        # actually evaluating the ME.
        self.PS_point = [PS_point[self.permutation[i]] for i in range(len(PS_point))]
        self.squared_orders = dict(squared_orders)
        
        if not spin_correlation:
            self.spin_correlation = tuple([])
        else:
            self.spin_correlation = ( [ (k, self.permutation[i])  for k, i in spin_correlation[0] ],
                                      [ (k, self.permutation[i])  for k, i in spin_correlation[1] ] )
        if not color_connection:
            self.color_connection = tuple([])
        else:
            self.color_connection = ( [self.permutation[i] for i in color_connection[0] if i > 0],
                                      [self.permutation[i] for i in color_connection[1] if i > 0]   )
            
        if not hel_config:
            hel_config = tuple([])
        else:
            self.hel_config = tuple( hel_config[self.permutation[i]] for i in range(len(hel_config)) )

class F2PYMEAccessor(VirtualMEAccessor):
    """ A class wrapping the access to one particular MatrixEleemnt wrapped with F2PY"""
    
    def __new__(cls, process, f2py_module_path, slha_card_path, **opts):
        """ Factory for the various F2PY accessors, depending on the process a different class
        will be selected."""
        if cls is F2PYMEAccessor:
            target_class = None
            if 'n_loops' in opts:
                n_loops = opts['n_loops']
            else:
                n_loops = 1 if process.get('perturbation_couplings') else 0
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

        # Keep track whether initialization is necessary or not
        self.module_initialized = False

        # Save the location of the SLHA card path
        if not os.path.isabs(slha_card_path):
            raise MadGraph5Error("The initialization of F2PYMEAccessor necessitates an "+
                                                        "absolute path for slha_card_path.")
        self.slha_card_path = os.path.relpath(slha_card_path, self.root_path)

        self.f2py_module_path = f2py_module_path
        self.f2py_module = self.load_f2py_module(f2py_module_path)

    def load_f2py_module(self, module_path):
        """ Loads a f2py module given its path and returns it."""
        
        # Make sure to temporarily adjust the environment
        if module_path[0] not in sys.path:
            sys.path.insert(0, module_path[0])
        f2py_module = importlib.import_module(module_path[1])
        sys.path.pop(sys.path.index(module_path[0]))
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
    
    def __call__(self, PS_point, alpha_s, mu_r=91.188, **opts):
        """ Actually performs the f2py call. """
        
        # The mother class takes care of applying the permuation
        self.apply_permutations(PS_point, **opts)
        
        # For now, only support basic information
        if self.squared_orders:
            raise MadGraph5Error("F2PYMEAccessor does not yet support the specification of squared orders.")
        if self.spin_correlation:
            raise MadGraph5Error("F2PYMEAccessor does not yet support the specification of spin_correlations.")
        if self.color_connection:
            raise MadGraph5Error("F2PYMEAccessor does not yet support the specification of color_connections.")
        if self.hel_config:
            raise MadGraph5Error("F2PYMEAccessor does not yet support the specification of helicity configurations.")

        # If/When grouping several processes in the same f2py module (so as to reuse the model for example),
        # we will be able to use the information of self.process_pdgs to determine which one to call.
        # misc.sprint(" I was called from :",self.process_pdgs)
        
        if not self.module_initialized:
            self.f2py_module.initialise(pjoin(self.root_path, self.slha_card_path))
            self.module_initialized = True
        
        # Most basic access for now
        output_data ={}
        output_data['finite'] = self.f2py_module.get_me( self.format_momenta_for_f2py(self.PS_point), alpha_s, 0)
        
        return output_data
    
class F2PYMEAccessorMadLoop(F2PYMEAccessor):
    """ Specialization of F2PYMEAccessor for MadLoop."""

    def __init__(self, process, f2py_module_path, slha_card_path, madloop_resources_path=None, **opts):
        """ Use the MadLoop resources path for MadLoop outputs """
  
        super(F2PYMEAccessorMadLoop, self).__init__(process, f2py_module_path, slha_card_path, **opts)

        if madloop_resources_path:
            # Save the location of MadLoop resources directory path
            if not os.path.isabs(madloop_resources_path):
                raise MadGraph5Error("The initialization of F2PYMEAccessorMadLoop necessitates an "+
                                                            "absolute path for madloop_resources_path.")
            self.madloop_resources_path = os.path.relpath(madloop_resources_path, self.root_path)
        else:
            self.madloop_resources_path = None

    def __call__(self, PS_point, alpha_s, mu_r, **opts):
        """ Actually performs the f2py call.
        """
        
        # The mother class takes care of applying the permuation
        self.apply_permutations(PS_point, **opts)
        
        # For now, only support basic information
        if self.squared_orders:
            raise MadGraph5Error("F2PYMEAccessorMadLoop does not yet support the specification of squared orders.")
        if self.spin_correlation:
            raise MadGraph5Error("F2PYMEAccessorMadLoop does not yet support the specification of spin_correlations.")
        if self.color_connection:
            raise MadGraph5Error("F2PYMEAccessorMadLoop does not yet support the specification of color_connections.")
        if self.hel_config:
            raise MadGraph5Error("F2PYMEAccessorMadLoop does not yet support the specification of helicity configurations.")

        # If/When grouping several processes in the same f2py module (so as to reuse the model for example),
        # we will be able to use the information of self.process_pdgs to determine which one to call.
        # misc.sprint(" I was called from :",self.process_pdgs)
        
        if not self.module_initialized:
            self.f2py_module.initialise(pjoin(self.root_path, self.slha_card_path))
            if self.madloop_resources_path:
                self.f2py_module.initialise_madloop_path(pjoin(self.root_path, self.madloop_resources_path))
            self.module_initialized = True
        
        # Most basic access for now
        output_data ={}
        finite_loop_me, return_code = self.f2py_module.get_me( 
                                                    self.format_momenta_for_f2py(self.PS_point), alpha_s, mu_r, 0)
        output_data['finite'] = finite_loop_me
        output_data['ML_return_code'] = return_code        
        
        return output_data
            
class MEAccessorDict(dict):
    """ A class for nicely wrapping the access to the (possibly spin- and color- correlated) matrix elements of 
    various processes (mostly via f2py)."""
    
    def __init__(self):
        """ Initialize an allMEAccessor. """
        
        # Define the general characteristic of the current query
        self.curr_n_loops           = 0
        self.curr_perturbations     = []
        
    def set_queries_details(self, n_loops=None, perturbations=None):
        """ Specifies global details about the queries for ME performed in this MEAccessorDict."""
        if not n_loops is None:
            self.curr_n_loops = n_loops
        if not perturbations is None:
            self.curr_perturbations = tuple(sorted(perturbations))

    def __getitem__(self, key):
        return self.get_MEAccessor(key)

    def get_MEAccessor(self, pdgs, **opts):
        """ Provides access to a given ME, provided its PDGs and possibly other specifications such as
        n_loops and perturbations, given in options. The PDGs should be provided in the following form:
           (initial_states_pdgs_in_a_tuple, final_states_pdgs_in_a_tuple)
        """

        # Propagate additional characteristics that the user might pass
        self.set_queries_details(**opts)
        
        # Build the canonical key corresponding to the request for the look-up
        canonical_key = ( ( tuple(sorted(pdgs[0])), tuple(sorted(pdgs[1])) ), 
                                                            self.curr_n_loops, self.curr_perturbations )
        
        try:
            (ME_accessor, defining_pdgs_order) = super(MEAccessorDict, self).__getitem__(canonical_key)
        except KeyError:
            raise MadGraph5Error("This collection of matrix elements does not contain process %s."%str(canonical_key))
        
        # Figure out the permutations to apply *TO the user inputs* in order to get to the *order assumed in the ME*
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

        # Now assign this permutation to the ME_accessor found and return it
        ME_accessor.permutation = permutation
        ME_accessor.process_pdgs = defining_pdgs_order
        return ME_accessor

    def add_MEAccessor(self, ME_accessor):
        """ Add a particular ME accessor to the collection of available ME's."""
        if not isinstance(ME_accessor, VirtualMEAccessor):
            raise MadGraph5Error("MEAccessorDict can only be assigned values inheriting from VirtualMEAccessor.")
        
        for key, value in ME_accessor.get_canonical_key_value_pairs():
            if key in self:
                raise MadGraph5Error("Attempting to assign two MEAccessors to the same key in MEAccessorDict.")
            self[key] = value

    def add_MEAccessors(self, ME_accessor_list):
        """ Add a list of ME_accessors."""
        for ME_accessor in ME_accessor_list:
            self.add_MEAccessor(ME_accessor)

    def modify_root_path_in_accessors(self, new_root_path):
        """ Modifies the root_path of all accessors stored."""
        
        if not os.path.isabs(new_root_path):
            raise MadGraph5Error("Accessor root paths must be absolute.")

        for accessor in self.values():
            accessor.root_path = new_root_path

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
            target_class = None
            if contribution_definition.correction_order == 'LO':
                if contribution_definition.n_loops == 0 and \
                  contribution_definition.n_unresolved_particles == 0:
                   target_class = Contribution_B
                elif contribution_definition.n_loops == 0 and \
                  contribution_definition.n_unresolved_particles == 0:
                   target_class = Contribution_LIB                    
            elif contribution_definition.correction_order == 'NLO':
                if contribution_definition.n_loops == 1 and \
                   contribution_definition.n_unresolved_particles == 0:
                    target_class = Contribution_V
                elif contribution_definition.n_loops == 0 and \
                     contribution_definition.n_unresolved_particles == 1:
                    target_class = Contribution_R
            elif contribution_definition.correction_order == 'NNLO':
                if contribution_definition.n_loops == 0 and \
                   contribution_definition.n_unresolved_particles == 2:
                    target_class = Contribution_RR
                else:
                    raise MadGraph5Error("Some NNLO type of contributions are not implemented yet.")                
            if not target_class:
                raise MadGraph5Error("Could not determine the type of contribution to be added for"+
                                     " the contribution definiton:\n%s"%str(contribution_definition.nice_string()))
            return super(Contribution, cls).__new__(target_class, contribution_definition, cmd_interface, **opt)
        else:
            return super(Contribution, cls).__new__(cls, contribution_definition, cmd_interface, **opt)
    
    def __init__(self, contribution_definition, cmd_interface, **opts):
        """ Instantiates a particular contribution."""
        
        self.contribution_definition = contribution_definition
        self.amplitudes              = diagram_generation.AmplitudeList()
        self.all_matrix_elements     = helas_objects.HelasMultiProcess()
        self.exporter                = None
        
        # Things burrowed from the cmd_interface
        self.options                 = dict(cmd_interface.options)
        self.options['_model_v4_path'] = cmd_interface._model_v4_path
        self.options['_mgme_dir']    = cmd_interface._mgme_dir
        self.model                   = cmd_interface._curr_model
        
        # The following two attributes dicate the type of Exporter which will be assigned to this contribution
        self.output_type             = 'default'
        self.format                  = 'standalone'
        self.export_dir              = 'None'
        if self.options['group_subprocesses'] == 'Auto':
            self.collect_mirror_procs = True
        else:
            self.collect_mirror_procs   = self.options['group_subprocesses']

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
                format = self.format)

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
            
            self.all_matrix_elements = HelasMultiProcessClass(
                self.amplitudes, compute_loop_nc=compute_loop_nc, matrix_element_opts=generation_mode)
            
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

    def export(self, nojpeg=False, group_processes=True, args=[]):
        """ Perform the export duties, that include generation of the HelasMatrixElements and 
        the actual output of the matrix element code by the exporter."""
        
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

    def get_MEAccessors(self, root_path):
        """ Returns all MEAccessors for the matrix elemements generated as part of this contribution."""
        
        MEAccessors = []
        for proc_hash, (defining_process, mapped_processes) in self.get_processes_map().items():
            # The PDGs of the hashed representations correspond to entry [0][0]
            mapped_process_pdgs = [ (proc.get_initial_ids(), proc.get_final_ids()) for proc in mapped_processes ]
            f2py_module_path = pjoin(self.export_dir,'matrix_%s_py.so'%defining_process.shell_string())
            if not os.path.exists(f2py_module_path):
                raise MadGraph5Error("Cannot find the compiled f2py module for %s at %s"%
                                                            (defining_process.nice_string(),f2py_module_path))
            f2py_load_path = (os.path.dirname(self.export_dir), 
                    '%s.matrix_%s_py'%( os.path.basename(self.export_dir), defining_process.shell_string() ) )
            slha_card_path = pjoin(self.export_dir,'Cards','param_card.dat')
            if os.path.isdir(pjoin(self.export_dir,'SubProcesses','MadLoop5_resources')):
                madloop_resources_path = pjoin(self.export_dir,'SubProcesses','MadLoop5_resources')
            else:
                madloop_resources_path = ''
                
            MEAccessors.append(F2PYMEAccessor(
                defining_process, 
                f2py_load_path, 
                slha_card_path,
                madloop_resources_path=madloop_resources_path,
                n_loops=self.contribution_definition.n_loops,
                mapped_pdgs = mapped_process_pdgs, 
                root_path=root_path
                ) )
            
        return MEAccessors

    def get_integrands(self, model, run_card, all_MEAccessors, ME7_configuration):
        """ Returns all (though typically only one) integrand implementing this contribution.
        The instance of MEAccessorDict is necessary so as to be passed to the integrand instances.
        """

        return [ ME7_interface.ME7Integrand(model, run_card,
                                           self.contribution_definition,
                                           self.get_processes_map(),
                                           all_MEAccessors,
                                           ME7_configuration) ]
        
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
        
    def get_processes_map(self, force=False):
        """ Returns a dictionary with a tuple of PDGS and other information, using the format provided by
                 base_objects.Process.get_hashable_representation(self)
            as keys and values being a tuple with format
                 (defining_process, [mapped processes])
            where defining_process is an actual instance of Process and [mapped processes] is a list of elements
            of the same format as the keys.
        """
        if not self.amplitudes:
            return {}
        
        # Attempt reusing a cached version
        if hasattr(self, 'processes_map') and not force:
            if (self.processes_map[0]['had_amplitudes'] == bool(self.amplitudes)) and \
               (self.processes_map[0]['had_matrix_elements'] == bool(self.all_matrix_elements.get_matrix_elements())):
                return self.processes_map[1]
        
        all_defining_procs = [amp.get('process') for amp in self.amplitudes]
        all_defining_procs = dict( ( proc.get_hashable_representation()
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
            all_procs_pdgs = dict((proc.get_hashable_representation(), proc) for proc in me.get('processes'))
            # Look through all the processes identified in the amplitudes
            for proc in all_defining_procs:
                # Make sure it is found in the processes mapped in this Matrix element
                if not proc in all_procs_pdgs:
                    continue
                # Now add them to list of mapped processes for this Matrix element
                all_defining_procs[proc][1].extend([all_procs_pdgs[p] for p in all_procs_pdgs
                                            if p!=proc and p not in all_defining_procs[proc][1]])
        
        # Cache the process map
        self.processes_map = ({
                'had_amplitudes':bool(self.amplitudes),
                'had_matrix_elements':bool(self.all_matrix_elements.get_matrix_elements())},
                all_defining_procs)
        
        return all_defining_procs        
    
    def compile(self):
        """ Compiles the f2py shared library to provide easy access to this contribution."""
        
        for proc_hash, (defining_process, mapped_processes) in self.get_processes_map().items():
            name = defining_process.shell_string()
            Pdir = pjoin(self.export_dir, 'SubProcesses', 'P%s'%name)
            if not os.path.isdir(Pdir):
                raise MadGraph5Error("The expected subprocess directory %s could not be found."%Pdir)
            misc.compile(arg=['matrix_%s_py.so'%name, 'MENUM=_%s_'%name], cwd=Pdir)
            if not os.path.isfile(pjoin(Pdir, 'matrix_%s_py.so'%name)):
                raise InvalidCmd("The f2py compilation of subprocess '%s' failed.\n"%Pdir+
                                 "Try running 'make matrix2py.so' by hand in this directory.")
            ln(pjoin(Pdir, 'matrix_%s_py.so'%name), starting_dir=self.export_dir)
    
    def nice_string(self):
        """ Nice string representation of self."""
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        res = ['%-30s:   %s'%('contribution_type',type(self))]
        res.extend([self.contribution_definition.nice_string()])
        if self.amplitudes and not self.all_matrix_elements.get_matrix_elements():
            res.append('Amplitudes generated for the following processes:')
            for amp in self.amplitudes:
                res.append(GREEN+'  %s'%amp.get('process').nice_string().replace('Process: ','')+ENDC)
        elif self.amplitudes and self.all_matrix_elements.get_matrix_elements():
            res.append('Generated and mapped processes for this contribution:')
            for process_pdgs, (defining_process, mapped_processes) in self.get_processes_map().items():
                res.append(GREEN+'  %s'%defining_process.nice_string().replace('Process: ','')+ENDC)
                for mapped_process in mapped_processes:
                    res.append(BLUE+u'   \u21b3  '+mapped_process.nice_string().replace('Process: ','')+ENDC)
        else:
            res.append(BLUE+'No amplitudes generated yet.'+ENDC)                
        return '\n'.join(res)
        
    def generate_amplitudes(self):
        """ Generates the relevant amplitudes for this contribution."""

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
    pass

class Contribution_RR(Contribution):
    """ Implements the handling of a double real-emission type of contribution."""
    pass

class Contribution_V(Contribution):
    """ Implements the handling of a virtual type of contribution."""

    def __init__(self, contribution_definition, cmd_interface, **opts):
        """ Bring in the couple of modifications necessary for this type of contributions."""
        super(Contribution_V,self).__init__(contribution_definition, cmd_interface, **opts)
        # Make sure to adjust the MultiProcessClass to be used
        self.MultiProcessClass   = loop_diagram_generation.LoopMultiProcess
        self.output_type         = 'madloop'

    def generate_matrix_elements(self, group_processes=True):
        """Generate the Helas matrix elements before exporting. Uses the main function argument 
        'group_processes' to decide whether to use group_subprocess or not."""

        cpu_time1 = time.time()
        ndiags = 0
        
        self.all_matrix_elements = loop_helas_objects.LoopHelasProcess(self.amplitudes,
                                        optimized_output = self.options['loop_optimized_output'])
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

    def generate_code(self):
        """ Assuming the Helas Matrix Elements are now generated, we can write out the corresponding code."""

        matrix_elements = self.all_matrix_elements.get_matrix_elements()
            
        cpu_time_start = time.time()

        # Pick out the matrix elements in a list
        matrix_elements = self.all_matrix_elements.get_matrix_elements()
        # MadLoop standalone fortran output
        calls = 0
        for me in matrix_elements:
            calls = calls + self.exporter.generate_subprocess_directory(me, self.helas_model)
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

class ContributionList(base_objects.PhysicsObjectList):
    """ A container for storing a list of contributions."""
    
    contributions_natural_order = [
        ('LO',   (Contribution_B, Contribution_LIB) ),
        ('NLO',  (Contribution_R, Contribution_V) ),
        ('NNLO', (Contribution_RR, ) )
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

    def nice_string(self):
        """ A nice representation of a list of contributions. 
        We can reuse the function from ContributionDefinitions."""
        return base_objects.ContributionDefinitionList.contrib_list_string(self)

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
                    contrib_function(*method_args, **method_opts)
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
            contrib_function(*method_args, **method_opts)
