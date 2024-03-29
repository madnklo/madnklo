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
import madgraph.core.accessors as accessors
from accessors import ProcessKey
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
import madgraph.integrator.ME7_integrands as ME7_integrands
from madgraph import InvalidCmd, MadGraph5Error
from madgraph.iolibs.files import cp, ln, mv

logger = logging.getLogger('contributions')

# The initial state collinear counterterms may be combined because the
# BeamCurrents that implements them return the full backward evolution flavor
# matrix, so that several initial state collinear counterterms that differ only
# by the resulting backward evolved flavor (e.g. 'd g <- d' and 'g d~ <- d') could
# be combined. This option is not implemented as of now, so the only valid value
# for the global parameter below is False.
_COMBINE_INITIAL_STATE_INTEGRATED_COLLINEAR_COUNTERTERMS = False

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
            correction_order = contribution_definition.correction_order.count('N')
            beam_factorization_count = 0
            if contribution_definition.is_beam_active('beam_one'):
                beam_factorization_count += 1
            if contribution_definition.is_beam_active('beam_two'):
                beam_factorization_count += 1
            n_loops = contribution_definition.n_loops
            n_unresolved_particles = contribution_definition.n_unresolved_particles
            # Beam factorization contributions are automatically of type RV because
            # they must both generate local counterterms (for the form factors) and
            # accept integrated ISR ones.
            if beam_factorization_count > 0:
                if contribution_definition.correlated_beam_convolution:
                    target_type = 'BeamSoft'
                else:
                    target_type = 'RealVirtual'                
            elif n_loops == 0 and n_unresolved_particles == 0:
                target_type = 'Born'
            elif n_loops == 1 and n_unresolved_particles == 0:
                if correction_order < 1:
                    target_type = 'LoopInduced_Born'
                else:
                    target_type = 'Virtual'
            elif n_loops == 0 and n_unresolved_particles == 1:
                target_type = 'SingleReals'
            elif n_loops == 0 and n_unresolved_particles == 2:
                target_type = 'DoubleReals'
            elif n_loops == 0 and n_unresolved_particles == 3:
                target_type = 'TripleReals'                
            else:
                raise MadGraph5Error("Some %s type of contributions are not implemented yet."%
                                                  contribution_definition.correction_order)

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
        ### For testing purposes, one can force to always include NNLO types of correlators in all matrix elements
        ### outputs simply with the line below
        #correlators_needed = max(correlators_needed,2)
        ###
        ##############################################################################################################
        ##############################################################################################################
        
        if correlators_needed > 0:
            self.additional_exporter_options['color_correlators'] ='N'*correlators_needed+'LO'
            # Also force MG5_aMC option 'loop_color_flows' to be true as it is necessary to compute color_correlators.
            # We can do so in a neat way here by simply adding this option to self.additional_exporter_options since they
            # will overwrite the interface option when the exporter will be instantiated.
            self.additional_exporter_options['loop_color_flows'] = True

            self.additional_exporter_options['color_correlators'] ='N'*correlators_needed+'LO'
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
        self.track_leg_numbers = False
        if self.contribution_definition.overall_correction_order.count('N') > 0:
            self.IR_subtraction = subtraction.IRSubtraction(
                self.model,
                coupling_types     = self.contribution_definition.correction_couplings,
                beam_types         = self.contribution_definition.get_beam_types(),
                subtraction_scheme = self.options['subtraction_scheme'] )
            if self.IR_subtraction.subtraction_scheme_module is not None:
                self.track_leg_numbers = self.IR_subtraction.subtraction_scheme_module.are_current_instances_for_specific_leg_numbers
        
        # The following two attributes dictate the type of Exporter which will be assigned to this contribution
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
        if 'diagram_filter' in opts:
            self.diagram_filter         = opts['diagram_filter']
        else:
            self.diagram_filter         = False
        self.optimize                   = False
        # Flag to indicate whether this contribution supports decay chains
        # Typically only LO contributions support decay chain amplitudes
        self.supports_decay_chain       = False
        
        if 'loop_filter' in opts:
            self.loop_filter            = opts['loop_filter']
        else:
            self.loop_filter            = None

        # Specifies if this contribution output requires its own MODEL in Source
        self.requires_its_own_model     = False
        
        # The following two attributes dictate the type of Exporter which will be assigned to this contribution
        self.output_type             = self.get_output_type()
        self.format                  = self.get_output_format()
        
        # Specify the MultiProcessClass to use to generate amplitudes
        self.MultiProcessClass       = self.get_multi_process_class()
        
        # Add default phase-space topology specifiers
        self.processes_to_topologies = {}
        self.topologies_to_processes = {}
        
        # Add an empty processes map, specifying that it is obtained before anything
        # was generated yet
        self.processes_map = ({
                'had_amplitudes'        :   False,
                'had_matrix_elements'   :   False
            }, {})
        
    def get_multi_process_class(self):
        """ Function returning the appropriate multiprocess class to instantiate during
        the amplitude generation for this particular contribution."""
        
        if self.output_type == 'madloop':
            # Make sure to adjust the MultiProcessClass to be used in the case of a 'madloop' output
            return loop_diagram_generation.LoopMultiProcess
        elif self.output_type == 'default':
            return diagram_generation.MultiProcess
        else:
            raise MadGraph5Error(
                "Output type '%s' not reckognized by the base Contribution class."%self.output_type+
                " It must be implemented explicitely in the daughter class.")

    def get_output_type(self):
        """ Function returning the appropriate output type to be used when outputting
        the amplitudes from this particular contribution."""
        
        # Decide whether a MadLoop output is required based on the process definition
        # underlying this contribution.
        if (self.contribution_definition.process_definition.get('NLO_mode') in ['virt','sqrvirt']):
            return 'madloop'
        else:
            return 'default'

    def get_output_format(self):
        """ Function returning the appropriate output format to be used when outputting
        the amplitudes from this particular contribution."""
        return 'standalone'
        
    def add_content_to_global_ME7_resources(self, global_ME7_dir, **opts):
        """ Possibly add content to the global ME7 resources pool, depending on the output 
        format requested."""
        
        if self.format == 'standalone' and self.output_type == 'madloop':        
            destination = pjoin(global_ME7_dir,'Cards','MadLoopParams.dat')
            if not os.path.exists(destination):
                shutil.copyfile(pjoin(self.export_dir,'Cards','MadLoopParams.dat'), destination)


    def set_helas_model(self):
        """ Instantiate properly the helas model according to the exporter and the
        output type."""

        if self.exporter.exporter == 'cpp':       
            self.helas_model = helas_call_writers.CPPUFOHelasCallWriter(self.model)
        else:
            assert self.exporter.exporter == 'v4'
            if self.output_type == 'madloop':
                assert (not self.options['_model_v4_path'])
                self.helas_model = helas_call_writers.FortranUFOHelasCallWriter(self.model)
            else:
                if self.options['_model_v4_path']:
                    self.helas_model = helas_call_writers.FortranHelasCallWriter(self.model)
                else:
                    self.helas_model = helas_call_writers.FortranUFOHelasCallWriter(self.model)

    def add_integrated_counterterm(self, integrated_CT_properties):
        """By default, do not support adding integrated counterterms."""
        
        raise MadGraph5Error(
            "The contribution of type %s cannot receive" % type(self) +
            " contributions from integrated counterterms." )

    def combine_initial_state_counterterms(self, *args, **opts):
        """ By default, does nothing. For inheritance purposes only."""
        pass

    def set_export_dir(self, prefix):
        """Assigns an export directory name."""

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
            cmd_interface, noclean, output_type=self.output_type,
            group_subprocesses=group_subprocesses,
            curr_amps=self.amplitudes,
            export_dir=self.export_dir,
            format=self.format,
            additional_options=self.additional_exporter_options)

    def copy_template(self, model):
        """Copy the template structure for that contribution.
        Quite often, this limits itself to asking its exporter to do this.
        """

        ret_value =  self.exporter.copy_template(model)
         # Make sure that this contribution's output as an __init__ file
        if not os.path.isfile(pjoin(self.export_dir, '__init__.py')):
            open(pjoin(self.export_dir, '__init__.py'),'w').write('')
        return ret_value

    def pass_information_from_cmd(self, cmd_interface):
        """Pass information from the command_interface to this contribution. Most of the time, this only amounts
        to passing information to the active exporter."""

        return self.exporter.pass_information_from_cmd(cmd_interface)

    def generate_matrix_elements(self, group_processes=True):
        """ Generate the matrix element in the manner appropriate to the selected output type."""
        
        if self.output_type == 'madloop':
            return self.generate_loop_matrix_elements(group_processes=group_processes)
        else:
            return self.generate_tree_matrix_elements(group_processes=group_processes)
    
    def generate_loop_matrix_elements(self, group_processes=True):
        """Generate the *loop* Helas matrix elements before exporting. Uses the main function argument 
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

    def generate_tree_matrix_elements(self, group_processes=True):
        """Generate the *tree* Helas matrix elements before exporting. Uses the main function argument 
        'group_processes' to decide whether to use group_subprocess or not."""

        cpu_time1 = time.time()
        ndiags = 0
        
        if group_processes:
            dc_amps = diagram_generation.DecayChainAmplitudeList(
                [amp for amp in self.amplitudes if isinstance(amp, diagram_generation.DecayChainAmplitude)])
            non_dc_amps = diagram_generation.AmplitudeList(
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

    def generate_code(self):
        """ Assuming the Helas Matrix Elements are now generated, we can write out the corresponding code."""

        matrix_elements = self.all_matrix_elements.get_matrix_elements()
        calls=0
        cpu_time_start = time.time()
        
        if self.output_type == 'madloop':
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
        else:
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
    
    def subtract(self, ignore_integrated_counterterms=False, **opts):
        """ Generate and export all necessary subtraction counterterterms and currents
        (including beam factorization counterterms as well [bulk ones are already generated
        at this stage] ones.). Nothing to do in this base class."""
        
        # Assign an empty list of countrerterms for now to this contribution
        self.counterterms = {}

        # No information need to be returned at this stage since nothing was generated
        return {}
        
    def export(self, nojpeg=False, group_processes=True, args=[], **opts):
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
        
        # In principle information can be passed back to ME7Exporter with the return value
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
        
        # This is not needed for factorization contributions as they always use processes already
        # in other contributions.
        if self.has_beam_factorization():
            return MEAccessors

        for process_key, (defining_process, mapped_processes) in self.get_processes_map().items():
            # The PDGs of the hashed representations correspond to entry [0][0]
            mapped_process_pdgs = [ (proc.get_initial_ids(), proc.get_final_ids()) for proc in mapped_processes ]
            proc_dir = pjoin(self.export_dir, 'SubProcesses', 'P%s' % self.process_dir_name(defining_process) )
            if not os.path.isdir(proc_dir):
                raise MadGraph5Error(
                    "Cannot find the process directory '%s' for process %s." %
                    (proc_dir, defining_process.nice_string()) )

            f2py_load_path = (
                os.path.dirname(self.export_dir),
                '%s.matrix_%s_py' % (
                    os.path.basename(self.export_dir),
                    self.process_dir_name(defining_process) ) )
            slha_card_path = pjoin(self.export_dir,'Cards','param_card.dat')
            if os.path.isdir(pjoin(self.export_dir,'SubProcesses','MadLoop5_resources')):
                madloop_resources_path = pjoin(self.export_dir,'SubProcesses','MadLoop5_resources')
            else:
                madloop_resources_path = ''
                
            MEAccessors.append(accessors.VirtualMEAccessor(
                defining_process, 
                f2py_load_path,
                slha_card_path,
                madloop_resources_path=madloop_resources_path,
                mapped_pdgs=mapped_process_pdgs,
                root_path=root_path,
                compile_if_necessary=False,
            ) )

        # Only allow overwriting accessors if processes were not grouped
        all_MEAccessors.add_MEAccessors(
            MEAccessors, allow_overwrite=(not self.group_subprocesses) )

    def add_current_accessors(
        self, model, all_MEAccessors, root_path, current_set, currents_to_consider, CT_type='local'):
        """Generate and add all subtraction current accessors to the MEAccessorDict.
        The optioin CT_type is only present in order to give a more precise printout
        """

        # Generate the computer code and export it on disk for the remaining new currents
        current_exporter = subtraction.SubtractionCurrentExporter(model, root_path, current_set)
        mapped_currents = current_exporter.export(currents_to_consider)
        # Print to the debug log which currents were exported
        log_string = "The following subtraction %s current implementations are exported "%CT_type+\
                     "in contribution '%s':\n"%self.short_name()
        default_implementation_string = ''
        already_listed = []
        for (module_path, class_name, _), current_properties in mapped_currents.items():
            if class_name != 'DefaultCurrentImplementation':
                if class_name not in already_listed:
                    already_listed.append(class_name)
                    quote_class_name = "'%s'" % class_name
                    defining_current_str = str(current_properties['defining_current'])
                    line_pars = (quote_class_name, defining_current_str)
                    log_string += " > %-35s for representative current '%s'\n" % line_pars
            else:
                quote_default_name = "'DefaultCurrentImplementation'"
                number_of_currents = len(current_properties['mapped_process_keys'])
                line_pars = (quote_default_name, number_of_currents)
                default_implementation_string = " > %-35s for a total of %d currents.\n" % line_pars

        log_string +=default_implementation_string
        if len(mapped_currents)>0:
            logger.debug(log_string)
        # Instantiate the CurrentAccessors corresponding
        # to all current implementations identified and needed
        all_current_accessors = []

        for (subtraction_scheme, class_identifier, _), current_properties in mapped_currents.items():
            all_current_accessors.append(accessors.VirtualMEAccessor(
                current_properties['defining_current'],
                subtraction_scheme,
                class_identifier,
                '%s.commons.utils'%current_exporter.main_module_name,
                current_properties['instantiation_options'],
                mapped_process_keys=current_properties['mapped_process_keys'],
                root_path=root_path,
                model=model
            ))
        # Add MEAccessors
        all_MEAccessors.add_MEAccessors(all_current_accessors,)
        return mapped_currents

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
        return [ ME7_integrands.ME7Integrand(model, run_card,
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

    def remove_superfluous_content(self):
        """At the end of the export of this contributions,
        remove extra files not desired/necessary.
        """
        
        dirs_to_remove  = [ pjoin(self.export_dir, 'Cards') ]
        files_to_remove = [ pjoin(self.export_dir, 'Source','make_opts') ]
        
        for dir in dirs_to_remove:
            if os.path.isdir(dir):
                shutil.rmtree(dir)
            elif os.path.islink(dir):
                os.remove(dir)
        for f in files_to_remove:
            if os.path.isfile(f) or os.path.islink(f):
                os.remove(f)

    def finalize(self, flaglist=[], interface_history=[]):
        """ Finalize the output of the code necessary for this contribution."""
        
        # WARNING: It would be ideal to use the same DHELAS output for all contributions.
        # This should be done eventually, but there are issues with coef_specs,
        # complex vs real momenta, etc..
        # So for this first test we stick with a local DHELAS for each contributions
        
        # This stores the wanted couplings which should be exported by the overall ME7 exporter.
        global_wanted_couplings = []
        
        # For some contributions, like the beam factorization ones, there is no additional
        # export to be performed.
        if self.has_beam_factorization():
            return global_wanted_couplings

        # Handling of the model
        if self.options['_model_v4_path']:
            logger.info(
                'Copy %s model files to directory %s' %
                (os.path.basename(self.options['_model_v4_path']), self.export_dir) )
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
        finalize_options = dict(self.options)
        finalize_options['no_compilation'] = True
        finalize_options['ME7_output'] = True
        self.exporter.finalize(
            self.all_matrix_elements, interface_history, finalize_options, flaglist)
        
        return global_wanted_couplings

    def link_global_ME7_resources(self, global_ME7_dir):
        """ Remove some of the local resources and instead link the global ones from the ME7 exporters."""

        # For some contributions, like the beam factorization ones, there is no additional
        # export to be performed.
        if self.has_beam_factorization(): 
            return
        
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
        
        # For some contributions, like the beam factorization ones, there is no additional
        # export to be performed.
        if self.has_beam_factorization():
            return        
        
        self.exporter.make_model_symbolic_link()
    
    def get_inverse_processes_map(self, force=False):
        """ Returns a dictionary with the keys (with ordered PDGs) of all processes present
        in this contribution and the single value being a 2-tuple with their corresponding 
        "defining" process key as first entry and the corresponding mapped process instance as 
        a second entry.
        If force is True, then the map will be reconstructed entirely and the cache updated."""
        
        # Sort PDGs only when processes have been grouped
        sort_PDGs = self.group_subprocesses

        if hasattr(self, 'inverse_processes_map') and not force:
            return self.inverse_processes_map
            
        process_map = self.get_processes_map(force=force)
        inverse_map = {}
        for process_key, (defining_process, mapped_processes) in process_map.items():
            defining_process_key = ProcessKey(defining_process,sort_PDGs=sort_PDGs).get_canonical_key()
            if defining_process_key in inverse_map:
                raise MadGraph5Error("The following process key appears twice in "+
                                                  "the processes map: %s"%str(process_key))
            inverse_map[defining_process_key] = (process_key, defining_process)
            for mapped_process in mapped_processes:
                mapped_process_key = ProcessKey(mapped_process,sort_PDGs=sort_PDGs).get_canonical_key()
                if mapped_process_key in inverse_map:
                    raise MadGraph5Error(
                        "The following mapped process key:\n%s\n"%str(mapped_process_key)+
                        "of %s appears twice in the processes map."%
                               (mapped_process.nice_string().replace('Process','process')))                    
                inverse_map[mapped_process_key] = (process_key, mapped_process)
        
        self.inverse_processes_map = inverse_map
        
        return self.inverse_processes_map
                    

    def get_processes_map(self, force=False):
        """ Returns a dictionary with keys obtained from ProcessKey and values being a tuple with format
                 (defining_process, [mapped processes])
            where defining_process is an actual instance of Process and [mapped processes] is a list 
            of processes instances mapped to that defining process.
        """

        # Attempt reusing a cached version
        if hasattr(self, 'processes_map') and not force:
            if (self.processes_map[0]['had_amplitudes'] == bool(self.amplitudes)) and \
               (self.processes_map[0]['had_matrix_elements'] == bool(
                                          self.all_matrix_elements.get_matrix_elements())):
                return self.processes_map[1]
        
        if not self.amplitudes:
            return {}
        
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

        # Insure that all processes have leg numbers that are consecutive
        for process_key, (defining_process, mapped_processes) in all_defining_procs.items():
            for p in [defining_process,]+mapped_processes:
                for i, leg in enumerate(p.get('legs')):
                    leg.set('number', i+1)

        # Cache the process map
        self.processes_map = ({
                'had_amplitudes':bool(self.amplitudes),
                'had_matrix_elements':bool(self.all_matrix_elements.get_matrix_elements())},
                all_defining_procs)
        
        return all_defining_procs        

    def process_dir_name(self, process):
        """ Given a specified process, return the directory name in which it is output."""
        # This function is essentially useful in order to be shadowed by daughter classes
        if self.output_type == 'madloop':
            return "%d_%s_%s"%(process.get('id'), process.get('uid'),
                                                      process.shell_string(print_id=False))
        else:
            return process.shell_string()

    def compile(self):
        """ Compiles the f2py shared library to provide easy access to this contribution."""
        
        # This is not needed for factorization contributions as they always use processes already
        # in other contributions.
        if self.has_beam_factorization():
            return
        
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

        return GREEN+'  %s'%(
            defining_process.nice_string(print_weighted=False).replace('Process: ',''))+ENDC
    
    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this contribution."""
        return []
    
    def short_name(self):
        """ Returns a short-hand notation for that contribution."""
        return self.contribution_definition.short_name()
    
    def nice_string(self, format=0):
        """ Nice string representation of self."""
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        res = ['< %s%s%s >'%(BLUE,self.short_name(),ENDC)]
        res.append('%-30s:   %s'%('contribution_type',type(self)))
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
        elif self.amplitudes and self.all_matrix_elements.get_matrix_elements() or \
                                                             self.has_beam_factorization():
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
    
    def has_beam_factorization(self):
        """ Checks whether this contribution has as any beam factorization active."""

        return self.contribution_definition.has_beam_factorization()
    
    def import_topologies_from_contrib(self, process_key, contrib):
        """ Add to the current topologies of this contribution, the ones defined for the
        process key of the contribution passed in argument."""
        
        # TODO a proper combination needs to be performed here. For now, I simply add them
        # to the attribute processes_to_topologies and topologies_to_processes, but that
        # doesn't work because the topology_key used is based on the 
        #    (subprocess_group_number, config_number)
        # which will clash between the various contributions added. So matching configurations
        # must be combined and the different ones must be added with a topology key upgraded
        # so as to make them unique (for example using a combination of:
        #    contrib.contribution_definiton.process_definition.get('id')
        # and
        #    contrib.contribution_definiton.n_loops
        #.which is best done when generating the topologies in the first place anyway.
        # But they will need to be combined anyway.
        #
        # For now we implement a minimal behaviour below which is supposed to work for
        # NLO with a simple 'add process' command, where we can basically copy the topologies
        # directly from the Born
        self.processes_to_topologies[process_key] = contrib.processes_to_topologies[process_key]
        
        for (topology_key, processes) in contrib.topologies_to_processes.items():
                self.topologies_to_processes[topology_key] = processes

    def add_beam_factorization_processes_from_contribution(self, contrib, beam_factorization_order):
        """ 'Import' the processes_map and all generated attributes of the contribution
        passed in argument to assign them to this beam factorization contributions, and
        generate all necessary "bulk" BeamCurrents as we go along.
        """
        
        beam_factorization_configuration = (
           (not self.contribution_definition.beam_factorization['beam_one'] is None) and \
            self.contribution_definition.beam_factorization['beam_one']['active'],
           (not self.contribution_definition.beam_factorization['beam_two'] is None) and \
            self.contribution_definition.beam_factorization['beam_two']['active'] )

        # If both beams must be assigned a factorization term, then make sure that the
        # beam_factorization_order is at least two
        if beam_factorization_configuration == (True,True) and beam_factorization_order < 2:
            # We have nothing to add from this contribution.
            return 
        
        # Make sure the source contribution does have a factorizing beam for the active
        # ones of this contribution.
        if ( (beam_factorization_configuration[0] and 
              (contrib.contribution_definition.beam_factorization['beam_one'] is None)) or
             (beam_factorization_configuration[1] and 
               (contrib.contribution_definition.beam_factorization['beam_two'] is None)) ):
            return
        
        # Now fetch the process map to load from
        source_processes_map = contrib.get_processes_map()
        processes_map = self.processes_map[1]
        # To the new contribution instances, also assign references to the generated attributes
        # of the original contributions (no deep copy necessary here)
        for process_key, (defining_process, mapped_processes) in source_processes_map.items():
            
            # First add this process to the current processes map
            if process_key in processes_map:
                raise MadGraph5Error('MadNkLO attempts to add beam factorization terms for'+
                    ' the %s a second time. This should never happen.'%(
                             defining_process.nice_string().replace('Process','process')))
            processes_map[process_key] = (
                  defining_process.get_copy(), [mp.get_copy() for mp in mapped_processes] )
            
            # Create the corresponding BeamCurrents
            all_beam_currents = []
                        
            beam_one_type = contrib.contribution_definition.beam_factorization['beam_one']['beam_type']
            beam_two_type = contrib.contribution_definition.beam_factorization['beam_two']['beam_type']
            beam_one_PDGs = contrib.contribution_definition.beam_factorization['beam_one']['beam_PDGs']
            beam_two_PDGs = contrib.contribution_definition.beam_factorization['beam_two']['beam_PDGs']
            # The double for-loop below creates all compatible assignments of n_loop to the beam
            # factorization currents applied to the first and second beam. So (2,1),(1,2) at N^3LO for instance
            # In principle it would also be possible here to further differentiate P^(0) x P^(0) from P^(1), but
            # we do not do this for now.
            for factorization_order_beam_one in range(beam_factorization_order+1):
                if factorization_order_beam_one > 0 and not beam_factorization_configuration[0]:
                    break
                for factorization_order_beam_two in range(beam_factorization_order-factorization_order_beam_one+1):
                    if factorization_order_beam_one == factorization_order_beam_two == 0:
                        continue
                    if factorization_order_beam_two > 0 and not beam_factorization_configuration[1]:
                        break
                    initial_legs = { l.get('number'): l for l in defining_process.get_initial_legs() }
                    all_beam_currents.append({
                        'beam_one' : None if factorization_order_beam_one == 0 else
                          subtraction.BeamCurrent({
                            'beam_type'         : beam_one_type,
                            'beam_PDGs'         : beam_one_PDGs,
                            'distribution_type' : 'bulk',
                            'n_loops'           : factorization_order_beam_one-1,
                            'perturbation_couplings' : self.contribution_definition.correction_couplings,
                            'squared_orders'    : { order : factorization_order_beam_one*2 for
                                   order in self.contribution_definition.correction_couplings },
                            'sqorders_types'    : { order : '<=' for order in 
                                            self.contribution_definition.correction_couplings },
                            'singular_structure': subtraction.SingularStructure(
                                 substructures = [ subtraction.BeamStructure(
                                                                 legs = [initial_legs[1],]), ]
                            ),
                          }),
                        'beam_two' : None if factorization_order_beam_two == 0 else
                          subtraction.BeamCurrent({
                            'beam_type'         : beam_two_type,
                            'beam_PDGs'         : beam_two_PDGs,
                            'distribution_type' : 'bulk',
                            'n_loops'           : factorization_order_beam_two-1,
                            'perturbation_couplings' : self.contribution_definition.correction_couplings,
                            'squared_orders'    : { order : factorization_order_beam_two*2 for
                                   order in self.contribution_definition.correction_couplings },
                            'sqorders_types'    : { order : '<=' for order in 
                                            self.contribution_definition.correction_couplings },
                            'singular_structure': subtraction.SingularStructure(
                                 substructures = [ subtraction.BeamStructure(
                                                                 legs = [initial_legs[2],]), ]
                            ),
                          })
                    })

            # Add these beam currents to the processes just added to this contribution
            processes_map[process_key][0].set('beam_factorization', all_beam_currents)
            for mp in processes_map[process_key][1]:
                mp.set('beam_factorization', all_beam_currents)

            # Also add the corresponding topologies
            self.import_topologies_from_contrib(process_key, contrib)

    def generate_all_counterterms(self, 
                        ignore_integrated_counterterms=False, group_processes=True, **opts):
        """ Generate all counterterms associated to the processes in this contribution.
        Since we are in the mother class Contribution here, we only consider generating the
        beam factorization counterterm of type 'local_CT'."""
        
        # Initialise local and integrated counterterms to an empty list
        self.counterterms = {}
        all_integrated_counterterms = []
        
        # Nothing done by default
        return all_integrated_counterterms
        
    def generate_amplitudes(self, force=False):
        """ Generates the relevant amplitudes for this contribution and the construction
        of the instances of currents for the beam factorization terms."""
        
        # First check if the amplitude was not already generated
        if self.amplitudes and not force:
            return

        myproc = self.MultiProcessClass(self.contribution_definition.process_definition,
                    collect_mirror_procs = self.collect_mirror_procs,
                    ignore_six_quark_processes = self.ignore_six_quark_processes,
                    optimize=self.optimize, diagram_filter=self.diagram_filter,
                    loop_filter=self.loop_filter)
    
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
        if 'optimize' in opts:
            self.optimize = opts['optimize']
        self.supports_decay_chain       = True

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
        res =  GREEN+'  %s'%(defining_process.nice_string(print_weighted=False).\
                                                              replace('Process: ',''))+ENDC

        if not self.counterterms:
            return res                                                               

        if format < 2:
            if process_key in self.counterterms:
                res += ' | %d local counterterms'%len([ 1
                    for CT in self.counterterms[process_key] if CT.is_singular() ])
            else:
                res += ' | 0 local counterterm'
                
        else:
            long_res = [' | with the following local counterterms:']
            for CT in self.counterterms[process_key]:
                if CT.is_singular():
                    if format==2:
                        long_res.append( '   | %s' % str(CT))
                    elif format==3:
                        long_res.append( '   | %s' % CT.__str__(
                            print_n=True, print_pdg=True, print_state=True ) )
                    elif format>3:
                        long_res.append(CT.nice_string("   | "))
            res += '\n'.join(long_res)

        return res

    def generate_all_counterterms(self, ignore_integrated_counterterms=False, group_processes=True):
        """ Generate all counterterms associated to the processes in this contribution."""

        all_integrated_counterterms = []

        for process_key, (defining_process, mapped_processes) in self.get_processes_map().items():

            local_counterterms, integrated_counterterms =  self.IR_subtraction.get_all_counterterms(
                defining_process, self.contribution_definition.n_unresolved_particles,
                ignore_integrated_counterterms=ignore_integrated_counterterms)
            
            if ignore_integrated_counterterms:
                logger.warning(
"""The integrated counterpart of the local counterterms of contribution '%s' for %s are ignored as per user's request.
The resulting output must therefore be used for debugging only as it will not yield physical results."""
                %(self.short_name(),defining_process.nice_string(print_weighted = False).replace('Process','process')))


            # Make sure that the non-singular counterterm which correspond to the real-emission
            # matrix elements themselves are not added here
            local_counterterms = [ct for ct in local_counterterms if ct.is_singular()]

            # Set the reduced flavor map of all the counterterms
            for CT in integrated_counterterms:
                CT.set_reduced_flavors_map(defining_process, mapped_processes, self.IR_subtraction)
            for CT in local_counterterms:
                CT.set_reduced_flavors_map(defining_process, mapped_processes, self.IR_subtraction)

#            misc.sprint('Local CTs for %s'%(defining_process.nice_string()))
#            for lc in local_counterterms:
#                misc.sprint(str(lc))
#            misc.sprint('Integrated CTs for %s'%(defining_process.nice_string()))
#            for ic in integrated_counterterms:
#                misc.sprint(str(ic))

            self.counterterms[process_key] = local_counterterms

            # Now Add a additional information about the integrated_counterterms that will
            # be necessary for their evaluations in the corresponding lower-multiplicity
            # contribution
            
            # First, a list of two tuples ( (in_pdgs, out_pdgs), (in_leg_numbers, out_leg_numbers) )
            # listing all combinations of pdgs and leg numbers for this contribution
            flavors_combinations = []
            
            # Notice that the defining flavor combination will always be placed first
            for process in [defining_process,]+mapped_processes:
                initial_numbers = tuple( leg.get('number') for leg in process.get_initial_legs() )
                final_numbers   = tuple( leg.get('number') for leg in process.get_final_legs() )
                leg_numbers = ( initial_numbers, final_numbers  )
                initial_pdgs = process.get_initial_ids()
                final_pdgs = tuple(process.get_final_ids_after_decay())
                flavors_combinations.append( ( ( tuple(initial_pdgs),final_pdgs ), leg_numbers ) )
#               Mirrored processes are not accounted for at this stage, as these are processed at the
#               integrand evaluation level where mirrored copy of events are generated.
#               Notice that this mirroring shouldn't be done simply as in the line below.
##
##                if process.get('has_mirror_process'):
##                    flavors_combinations.append( 
##                           ( ( tuple(reversed(initial_pdgs)), final_pdgs ), leg_numbers ) )
##
#               because if one takes the following process for instance:
#                     u d~ > u d~ a
#               then the counterterm C(1,3) is only valid if leg #1 is u and leg #2 is d~. 
#               So swapping u(1) d~(2) for u(2) d~(1) would require propagating the indices
#               substitution in the entire counterterm (most importantly, its singular structure)
#               and create a copy of it. This pretty heavy (and the idea of the 'use_mirror_process'
#               option is precisely to avoid this work), so this is entirely bypassed
#               by implementing this symmetrisation during the integrand evaluation.
            
#            misc.sprint('doing process:',defining_process.nice_string())
#            misc.sprint('flavors_combinations',flavors_combinations)
            # Keep track of all integrated current topologies considered so as to never
            # repeat the same one twice, because it will be distributed over all possible
            # assignment of kinematics in the reduced virtual process.
            # For instance, if the real emission process is:
            #   e+(1) e-(2) > d(3) d~(4) d(5) d~(6) g(7)
            # Then the counterterms C(3,7) and C(5,7) have the same topology and C(5,7) should
            # be discarded because the counterterm C(3,7) will be distributed on both
            # d-quark #3 and #5 of the following reduced virtual process:
            #    e+(1) e-(2) > d(3) d~(4) d(5) d~(6)
            # We must therefore first distribute all counterterms over identical topologies
            integrated_current_topologies = {}
            for integrated_counterterm in integrated_counterterms:
                # The integrated current topology can directly be encoded as its ProcessKey
                # where the leg numbers are removed
                copied_current = integrated_counterterm.get_integrated_current().\
                                          get_copy(('squared_orders','singular_structure'))
                # The two initial states are always distinguishable, so we should not 
                # group them into a single topology
                integrated_current_topology = copied_current.get_key().get_canonical_key()
                try:
                    integrated_current_topologies[integrated_current_topology].append(
                                                                    integrated_counterterm)             
#                    misc.sprint("Skipping integrated counterterm:\n%s\n"%str(integrated_counterterm)+
#                                "because it has the same topology has:\n%s"%
#                                    str(integrated_currents_topologies[integrated_current_topology][0]))
                except KeyError:
                    # Register the new integrated current topology
                    integrated_current_topologies[integrated_current_topology] = [integrated_counterterm,]

            # For each topology, count the number of counterterms that are associated to
            # various list of external leg numbers. Within a given topology, one expects
            # that there is the same number of counterterm repetition for each combination of external legs
            # Example: resolved process is
            #   e+(1) e-(2) > g(3) g(4) g(5) d(6) d~(7)
            # All the following counterterms belong to the same topology, which we group here according
            # to which list of external leg numbers they include:
            #   (3,4,6) : C(C(S(3),4),6), C(C(S(4),3),6) 
            #   (3,5,6) : C(C(S(3),5),6), C(C(S(5),3),6) 
            #   (4,5,6) : C(C(S(4),5),6), C(C(S(5),4),6)
            # And we see that the repetition is now 2, and it is the same for each group.
            # On the other hand the repetition would have been only 1 for the topology
            # of the form  C(C(3,4),6) (i.e. no soft) since C(C(4,3),6) is not generated.
            repetitions_for_topology = {}
            for integrated_current_topology, integrated_counterterms in integrated_current_topologies.items():
                repetitions={}
                for integrated_counterterm in integrated_counterterms:
                    # Gather what are the singular leg numbers for each singular structure involved
                    singular_external_leg_numbers_list = []
                    for node in integrated_counterterm.nodes:
                        singular_legs = node.current['singular_structure'].get_all_legs()
                        # Sorted tuple of the numbers of all legs that are external
                        singular_external_leg_numbers_list.append( 
                            tuple(sorted([l.n for l in singular_legs if 
                            integrated_counterterm.momenta_dict[l.n] == frozenset([l.n,])])) )
                    singular_external_leg_numbers_tuple = tuple(sorted(singular_external_leg_numbers_list))
                    try:
                        repetitions[singular_external_leg_numbers_tuple] += 1
                    except KeyError:
                        repetitions[singular_external_leg_numbers_tuple] = 1
                n_repetitions = set(repetitions.values())
                if len(n_repetitions) != 1:
                    raise MadGraph5Error(
"All the following integrated counterterms share the same topology:\n"+
'\n'.join(str(CT) for CT in integrated_counterterms)+'\n'+
'But the repetition numbers derived from them is not identical in each group of singular legs:\n'+
'\n'.join('%s -> %d'%(str(k), v) for k,v in repetitions.items())+'\n'+
'This is not expected to happen.')
                repetitions_for_topology[integrated_current_topology] = repetitions.values()[0]

            for integrated_current_topology, integrated_counterterms in integrated_current_topologies.items():
                # We only include the first integrated counterterm matching this topology
                # The others are included effectively by the fact that this countertm will
                # be distributed in all possible permutations of the virtual process flavors
                integrated_counterterm = integrated_counterterms[0]
    
#                misc.sprint('Doing counterterm:',str(integrated_counterterm))
                # List the reduced flavors combinations. The keys are the tuple of the reduced
                # flavors combination and the value is the multiplication factor.
                reduced_flavors_combinations    = {}
                resolved_flavors_combinations   = {}
                # One subtlety occurs because of the backward flavor evolution in the initial
                # states which will be registered in the event since it must be convoluted with
                # the correct PDFs.
                # Consider the following process:
                #             g(1) s(2) > s(3) c(4) c~(5) u(6) u~(7)
                #   + mapped  g(1) q(2) > q(3) qprime(4) qprime~(5) u(6) u~(7)
                # and the counterterm C(2,3)C(4,5) which yields the following reduced mapped flavor:
                #             g g > g u u~             
                # with a multiplicity of naively n_f**2. However, the actual events that will be generated
                # will be the following:
                #             g s > g u u~
                #             g d > g u u~
                #             ...
                # It is then clear that each of this flavor configuration in the ME7 event should not be
                # multiplied by n_f**2 but instead just n_f. For this reason, the function 
                # generate_all_counterterms of the contribution class also produces the dictionary 
                # 'reduced_flavors_with_resolved_initial_states_combinations' which, in the example above, 
                # would be (nf=5):
                #     {  (21,3),(21,2,-2) : 5,
                #        (21,1),(21,2,-2) : 5,
                #        ...
                #     }
                # which will allow to divide the weight of each of these flavor configurations at the very
                # by the corresponding multiplicity factor n_f.
                reduced_flavors_with_resolved_initial_states_combinations = {}
                
                # Each integrated counterterm contributes proportionally to the ratio
                # of the symmetry factors of the resolved process and reduced one,
                # further divided by the symmetry factor of each group external leg
                # part of a singular structure.
                # For instance if the resolved process is:
                #    e+(1) e-(2) > d(3) d~(4) d(5) d~(6) g(7) g(8) g(9)
                # And the counterterm is C(3,4), C(5,6) and C(7,8), then the reduced process is
                #    e+(1) e-(2) > g(10) g(11) g(12) g(9)
                # (ignore the fact that the reduced ME is zero) and it will contribute with
                # the weight:
                #   S_t = (S_r / S_s) / S_u
                # Where S_r is the symmetry factor of the resolved process:
                #   S_r = 2!*2!*3! = 24
                # And S_s the product of the symmetry factor each group of singular legs:
                #   S_s = (1!*1!)*(1!*1!)*(2!) = 2
                # And S_u the symmetry factor of the unresolved process
                #   S_u = 4! = 24
                # So that the final resulting S_t is 1/2. The "symmetry_factors" attribute
                # below stores the symmetry factor S_t for each reduced flavor.
                # There should always be a unique value of S_t for each reduced_flavor
                symmetry_factors                = {}
                for flavor_combination, leg_numbers in flavors_combinations:
                    # Access the reduced flavors
                    reduced_flavors = integrated_counterterm.get_reduced_flavors(
                                                    defining_flavors = flavor_combination )
                    # Create a map from leg number to flavor
                    number_to_flavor_map = dict( 
                        [ ( number, flavor_combination[0][i] ) for i, number in enumerate(leg_numbers[0]) ] +
                        [ ( number, flavor_combination[1][i] ) for i, number in enumerate(leg_numbers[1]) ]
                    )
#                    misc.sprint('From resolved flavor:',number_to_flavor_map)
#                    misc.sprint('I get:',reduced_flavors)

                    resolved_symmetry_factor = misc.symmetry_factor(list(flavor_combination[1]))
                    reduced_symmetry_factor = misc.symmetry_factor(list(reduced_flavors[1]))
                    # Now identify the groups of external legs relevant for each singular structure
                    # of this integrated counterterm
                    singular_legs_symmetry_factor = 1
                    for node in integrated_counterterm.nodes:
                        singular_legs = node.current['singular_structure'].get_all_legs()
                        # List the numbers of all legs that are external
                        singular_external_leg_pdgs = [number_to_flavor_map[l.n] for 
                            l in singular_legs if l.state==subtraction.SubtractionLeg.FINAL]
                        singular_legs_symmetry_factor *= misc.symmetry_factor(singular_external_leg_pdgs)

#                    misc.sprint('Counterterm: %s'%str(integrated_counterterm))
#                    misc.sprint('    -> multiplicity                   = %d'%len(integrated_counterterms))
#                    misc.sprint('    -> CT with same topologies        = %s'%(' | '.join(str(CT) for CT in integrated_counterterms)) )                    
#                    misc.sprint('    -> resolved_flavors               = %s'%str(flavor_combination))
#                    misc.sprint('    -> reduced                        = %s'%str(reduced_flavors))
#                    misc.sprint('    -> resolved_symmetry_factor       = %d'%resolved_symmetry_factor)
#                    misc.sprint('    -> reduced_symmetry_factor        = %d'%reduced_symmetry_factor)
#                    misc.sprint('    -> singular_legs_symmetry_factor  = %f'%singular_legs_symmetry_factor)
#                    misc.sprint('    -> counterterm_repetition_factor  = %f'%repetitions_for_topology[integrated_current_topology])
                        
                    # Now correct for the repetition number in that topology
                    # Notice it should always be an integer, the float here is just so that
                    # one can check this explicitly.
                    singular_legs_symmetry_factor /= float(repetitions_for_topology[integrated_current_topology])

                    overall_symmetry_factor = resolved_symmetry_factor/\
                               float(reduced_symmetry_factor*singular_legs_symmetry_factor)
                    try:
                        reduced_flavors_combinations[reduced_flavors] += 1
                    except KeyError:
                        reduced_flavors_combinations[reduced_flavors] = 1
                    resolved_flavors_combinations[flavor_combination] = reduced_flavors
                    
                    # Now also build a multiplicity factor for combination of resolved initial-state
                    # flavors and reduced final-state ones.
                    reduced_flavors_with_resolved_initial_states = (flavor_combination[0], reduced_flavors[1])
                    try:
                        reduced_flavors_with_resolved_initial_states_combinations[
                                            reduced_flavors_with_resolved_initial_states] += 1
                    except KeyError:
                        reduced_flavors_with_resolved_initial_states_combinations[
                                            reduced_flavors_with_resolved_initial_states] = 1
                    
                    if reduced_flavors in symmetry_factors:
                        # Sanity check, all reduced flavors should share the same S_t
                        if overall_symmetry_factor != symmetry_factors[reduced_flavors]:
                            raise MadGraph5Error(
"The overall symmetry factors derived for counterterm %s and resolved "%str(integrated_counterterm)+
"flavors %s is different than one previously derived from a different "%str(flavor_combination)+
"flavor combination.  %f =!= %f"%(overall_symmetry_factor, symmetry_factors[reduced_flavors]))
                    else:
                        symmetry_factors[reduced_flavors] = overall_symmetry_factor
                
#                misc.sprint('Finally having reduced_flavors_combinations=',reduced_flavors_combinations)
#                misc.sprint('For reduced process = ',integrated_counterterm.process.nice_string())
                # Final sanity check, which should hold as far as the developer anticipated
                if len(set(symmetry_factors.values()))!=1:
                    raise MadGraph5Error(
"For the counterterm %s with the following resolved and reduced flavors:\n"%str(integrated_counterterm)+
'\n'.join('%s => %s'%(str(k),str(v)) for k, v in resolved_flavors_combinations.items())+'\n'+
"The list of overall symmetry factors computed contains the following different values:\n"+
'\n'.join('%s => %f'%(str(k),v) for k, v in symmetry_factors.items())+'\n'+
"This is thought never to happen.")

                all_integrated_counterterms.append({
                    'integrated_counterterm'             :   integrated_counterterm,
                    'resolved_flavors_combinations'      :   resolved_flavors_combinations,
                    'reduced_flavors_combinations'       :   reduced_flavors_combinations,
                    'reduced_flavors_with_resolved_initial_states_combinations' :
                                        reduced_flavors_with_resolved_initial_states_combinations,
                    'symmetry_factor'                    :   symmetry_factors.values()[0],
                    'multiplicity'                       :   len(integrated_counterterms)
                })

        return all_integrated_counterterms
    
    @classmethod
    def remove_counterterms_with_no_reduced_process(cls, all_MEAccessors, counterterms):
        """Given the list of available reduced processes encoded in the MEAccessorDict 'all_MEAccessors'
        given in argument, remove all the counterterms whose underlying reduced process does not exist."""

        all_removed = True
        for counterterm_info in list(counterterms):
            # Integrated counterterm come as dictionaries with additional metadata
            if isinstance(counterterm_info, dict):
                counterterm = counterterm_info['integrated_counterterm']
            else:
                counterterm = counterterm_info
            # Of course don't remove any counterterm that would be the real-emission itself.
            if not counterterm.is_singular():
                continue
            try:
                all_MEAccessors.get_MEAccessor(counterterm.process)
                all_removed = False
            except MadGraph5Error:
                # This means that the reduced process could not be found and
                # consequently, the corresponding counterterm must be removed.
                # Example: C(5,6) in e+ e- > g g d d~
                # Useful debug lines below
#                misc.sprint('-'*50)
#                import pprint
#                misc.sprint('Process key:\n%s'%(
#                    pprint.pformat(dict(ProcessKey(counterterm.process).get_canonical_key()))))
#                misc.sprint('-'*50)
#                misc.sprint('All process keys in the accessor:\n')
#                for key in all_MEAccessors:
#                    misc.sprint(pprint.pformat(key))
#                misc.sprint('-'*50)
                counterterms.remove(counterterm_info)
        
        return all_removed

    def subtract(self, ignore_integrated_counterterms=False, **opts):
        """ Generate and export all necessary subtraction counterterterms and currents
        (including beam factorization integrated counterterms [bulk ones are already generated
        at this stage])."""

        ret_value = super(Contribution_R,self).subtract(ignore_integrated_counterterms=False, **opts)

        # Fish out the group_processes option as it could be used when attempting to
        # generate all currents.
        if 'group_processes' in opts:
            group_processes = opts['group_processes']
        else:
            group_processes = True

        integrated_counterterms = self.generate_all_counterterms(
            group_processes=group_processes,
            ignore_integrated_counterterms=ignore_integrated_counterterms)
      
        # Add the integrated counterterms to be passed to the exporter
        ret_value.update({'integrated_counterterms': integrated_counterterms})
        
        return ret_value

    def get_all_necessary_local_currents(self, all_MEAccessors):
        """ Given the counterterms in place and the currents already accessible in the 
        all_MEAccessors, return what local currents are needed."""

        all_currents = []
        for process_key, counterterms in self.counterterms.items():
            for current in subtraction.IRSubtraction.get_all_currents(counterterms,track_leg_numbers=self.track_leg_numbers):
                # Retain only a single copy of each needed current.
                copied_current = current.get_copy(('squared_orders','singular_structure'))
                if copied_current not in all_currents or self.track_leg_numbers:
                    all_currents.append(copied_current)

        # Also add the beam_factorization currents from process itself if it has any
        for process_key, (defining_process, mapped_processes) in self.get_processes_map().items():
            for beam_factorization_currents in defining_process['beam_factorization']:
                for bfc in beam_factorization_currents.values():
                    if bfc is not None:
                        all_currents.append(bfc)

        # Now further remove currents that are already in all_MEAccessors


        # We can remove the leg information only if for the subtraction scheme considered it is irrelevant
        # for the selection of the hard-coded current implementation to consider.
        all_necessary_currents = {}
        for current in all_currents:
            current_key = current.get_key(track_leg_numbers=self.track_leg_numbers).get_canonical_key()
            if current_key not in all_MEAccessors and current_key not in all_necessary_currents:
                all_necessary_currents[current_key] = current

        return all_necessary_currents.values()

    def add_ME_accessors(self, all_MEAccessors, root_path):
        """ Adds all MEAccessors for the matrix elements and currents generated as part of this contribution."""

        # Get the basic accessors for the matrix elements
        super(Contribution_R, self).add_ME_accessors(all_MEAccessors, root_path)

        # Only initialize all_counterterms_removed to True if there is CT to removed to begin with
        all_counterterms_removed = any(len(CTs)>0 for CTs in self.counterterms.values())    
        for process_key, counterterms in self.counterterms.items():
            # Remove counterterms with non-existing underlying Born processes
            if not self.remove_counterterms_with_no_reduced_process(all_MEAccessors, counterterms):
                all_counterterms_removed = False
        # Sanity check
        if all_counterterms_removed:
            logger.warning("All local counterterms from contribution '%s' were removed because"%self.short_name()+
                " their corresponding reduced processes was not found in the list of processes exported.\n"+
                " This is likely wrong and originates from a mismatch in the process definitions of the various"+
                " contributions building this higher order computation.")

        # Obtain all necessary currents
        current_set = self.options['subtraction_scheme']

        currents_to_consider = self.get_all_necessary_local_currents(all_MEAccessors)

        self.add_current_accessors(
            self.model, all_MEAccessors, root_path, current_set, currents_to_consider, CT_type='local' )
        
    def remove_local_counterterms_set_to_zero(self, all_ME_accessors):
        """ Remove all local counterterms involving currents whose implementation has set
        them to zero."""
        for process_key, counterterms in self.counterterms.items():
            for counterterm in list(counterterms):
                # misc.sprint("Considering CT %s" % str(counterterm))
                if counterterm.is_singular():
                    for current in counterterm.get_all_currents():
                        accessor, _ = all_ME_accessors[current.get_key(track_leg_numbers=self.track_leg_numbers)]
                        if accessor.subtraction_current_instance.is_zero:
                            # misc.sprint("Zero current found in CT %s" % str(counterterm))
                            counterterms.remove(counterterm)
                            break

    def get_integrands_for_process_map(self, process_map, model, run_card, all_MEAccessors, ME7_configuration):
        """ Returns all the integrands implementing this contribution for the specified process_map.
        The instance of MEAccessorDict is necessary so as to be passed to the integrand instances.
        """
        
        relevant_counterterms = {}
        self.remove_local_counterterms_set_to_zero(all_MEAccessors)

        if self.counterterms:
            for process_key in process_map:
                relevant_counterterms[process_key] = self.counterterms[process_key]

        identified_sectors = self.IR_subtraction.get_sectors(
            self.contribution_definition,
            process_map,
            counterterms=relevant_counterterms,
        )

        return [
            ME7_integrands.ME7Integrand(
                model, run_card, self.contribution_definition,
                process_map, self.topologies_to_processes, self.processes_to_topologies,
                all_MEAccessors, ME7_configuration,
                sectors = identified_sectors,
                counterterms = relevant_counterterms
            )
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
    
    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        contribution."""
        GREEN = '\033[92m'
        ENDC = '\033[0m'        
        res = GREEN+'  %s'%(defining_process.nice_string(print_weighted=False)
                                                             .replace('Process: ',''))+ENDC

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
                    long_res.append( '   | %s' % str(CT))
                elif format==3:
                    long_res.append( '   | %s' % CT.__str__(
                        print_n=True, print_pdg=True, print_state=True ))
                elif format==4:
                    long_res.append(CT.nice_string("   | "))
                elif format>4:
                    long_res.append(CT.nice_string("   | "))
                    for key, value in CT_properties.items():
                        if not key in ['integrated_counterterm', 'matching_process_key']:
                            long_res.append( '     + %s : %s'%(key, str(value)))

            res += '\n'.join(long_res)

        return res

    def get_all_necessary_integrated_currents(self, all_MEAccessors):
        """ Given the list of currents already available encoded in the all_MEAccessors,
        generate all the integrated currents that must be exported."""
        
        all_currents = []
        for CT_properties in sum(self.integrated_counterterms.values(),[]):
            counterterm = CT_properties['integrated_counterterm']
            integrated_current = counterterm.get_integrated_current()
            
            # Retain only a single copy of each needed current.
            # We must remove the leg information since this is information is irrelevant
            # for the selection of the hard-coded current implementation to consider.
            copied_current = integrated_current.get_copy(('squared_orders','singular_structure'))
            if copied_current not in all_currents or self.track_leg_numbers:
                all_currents.append(copied_current)

        # We can remove the leg information only if for the subtraction scheme considered it is irrelevant
        # for the selection of the hard-coded current implementation to consider.
        all_necessary_currents = {}
        for current in all_currents:
            current_key = current.get_key(track_leg_numbers=self.track_leg_numbers).get_canonical_key()
            if current_key not in all_MEAccessors and current_key not in all_necessary_currents:
                all_necessary_currents[current_key] = current

        return all_necessary_currents.values()

    @classmethod
    def remove_counterterms_with_no_reduced_process(cls, all_MEAccessors, counterterms):
        """Given the list of available reduced processes encoded in the MEAccessorDict 'all_MEAccessors'
        given in argument, remove all the counterterms whose underlying reduced process does not exist.
        For now we can recycle the implementation of the Contribution_R class.
        """
        return Contribution_R.remove_counterterms_with_no_reduced_process(all_MEAccessors, counterterms)

    def remove_integrated_counterterms_set_to_zero(self, all_ME_accessors):
        """ Remove all integrated counterterms involving currents whose implementation has set
        them to zero."""

        for process_key, counterterms in self.integrated_counterterms.items():
            for integrated_counterterm_properties in list(counterterms):
                integrated_counterterm = integrated_counterterm_properties['integrated_counterterm']
                # misc.sprint("Considering integrated CT %s" % str(integrated_counterterm))
                if integrated_counterterm.is_singular():
                    for current in integrated_counterterm.get_all_currents():
                        accessor, _ = all_ME_accessors[current.get_key(track_leg_numbers=self.track_leg_numbers)]
                        if accessor.subtraction_current_instance.is_zero:
                            # misc.sprint("Zero current found in integrated CT %s" % str(integrated_counterterm))
                            logger.debug("Zero current found in integrated CT %s, removing it now."%str(integrated_counterterm))
                            counterterms.remove(integrated_counterterm_properties)
                            break

    def get_integrands_for_process_map(self, process_map, model, run_card, all_MEAccessors, ME7_configuration):
        """ Returns all the integrands implementing this contribution for the specified process_map.
        The instance of MEAccessorDict is necessary so as to be passed to the integrand instances.
        """
        
        relevant_counterterms = {}
        self.remove_integrated_counterterms_set_to_zero(all_MEAccessors)
        if self.integrated_counterterms:
            for process_key in process_map:
                relevant_counterterms[process_key] = self.integrated_counterterms[process_key]

        if self.IR_subtraction is None:
            identified_sectors = None
        else:
            identified_sectors = self.IR_subtraction.get_sectors(
                self.contribution_definition,
                process_map,
                integrated_counterterms=relevant_counterterms,
            )

        return [ ME7_integrands.ME7Integrand(model, run_card,
                                       self.contribution_definition,
                                       process_map,
                                       self.topologies_to_processes,
                                       self.processes_to_topologies,
                                       all_MEAccessors,
                                       ME7_configuration,
                                       sectors = identified_sectors,
                                       integrated_counterterms=relevant_counterterms)
               ]

    def add_ME_accessors(self, all_MEAccessors, root_path):
        """ Adds all MEAccessors for the matrix elements and currents generated as part 
        of this contribution."""
        
        # Get the basic accessors for the matrix elements
        super(Contribution_V, self).add_ME_accessors(all_MEAccessors, root_path)

        # Only initialize all_counterterms_removed to True if there is CT to removed to begin with
        all_counterterms_removed = any(len(CTs)>0 for CTs in self.integrated_counterterms.values())
        for process_key, CT_properties in self.integrated_counterterms.items():
            # Remove integrated counterterms with non-existing underlying Born processes
            if not self.remove_counterterms_with_no_reduced_process(all_MEAccessors, CT_properties):
                all_counterterms_removed = False
 
        # Sanity check
        if all_counterterms_removed:
            logger.warning("All integrated counterterms from contribution '%s' were removed because"%self.short_name()+
                " their corresponding reduced processes was not found in the list of processes exported.\n"+
                " This is likely wrong and originates from a mismatch in the process definitions of the various"+
                " contributions building this higher order computation.")
        
        # Obtain all necessary currents
        current_set = self.options['subtraction_scheme']
        currents_to_consider = self.get_all_necessary_integrated_currents(all_MEAccessors)
        self.add_current_accessors(
            self.model, all_MEAccessors, root_path, current_set, currents_to_consider, CT_type='integrated' )

    @classmethod
    def get_basic_permutation(cls, origin_pdg_orders, destination_pdg_orders):
        """ Figures out the permutation to apply to map the origin pdg orders
        (a 2-tuple of initial and final state pdgs) to the destination pdg orders (same format)."""

        # Create a look_up list from the user-provided list of PDGs, whose elements will progressively be set to zero
        # as they are being mapped
        look_up_list = list(list(a_tuple) for a_tuple in origin_pdg_orders)
        
        basic_permutation = {}
        # Map the initial states
        for i, pdg in enumerate(destination_pdg_orders[0]):
            try:
                basic_permutation[i] = look_up_list[0].index(pdg)
                # Set the element that just got mapped to 0 in the look_up_list so that it will not be reused
                look_up_list[0][basic_permutation[i]] = 0
            except ValueError:
                raise MadGraph5Error("Cannot map two PDGs list: %s and %s"%(
                                str(destination_pdg_orders[0]), str(origin_pdg_orders[0])))
        
        # Map final states now
        n_initial = len(destination_pdg_orders[0])
        for i, pdg in enumerate(destination_pdg_orders[1]):
            try:
                basic_permutation[i+n_initial] = look_up_list[1].index(pdg)+n_initial
                # Set the element that just got mapped to 0 in the look_up_list so that it will not be reused
                look_up_list[1][basic_permutation[i+n_initial]-n_initial] = 0
            except ValueError:
                raise MadGraph5Error("Cannot map two PDGs list: %s and %s"%(
                                str(destination_pdg_orders[1]), str(origin_pdg_orders[1])))   
        
        return basic_permutation
        
    @classmethod
    def distribute_parent_flavors(cls, parent_flavors, reduced_process_pdgs):
        """ Returns a list of mapping (aka. permutations) that assign the parent_flavors,
        which is a 2-tuple of the following format:
             ( { initial_state_"singular"_flavor_pdg : [leg_positions] }, 
               { final_state_"singular"_flavor_pdg   : [leg_positions] } )
        to all possible places in the reduced_process_pdgs_list which is simply
             ( (initial_state_pdgs), (final_state_pdgs) )
        For instance if the parent flavors to distribute are only final and read:
            ( {}, {1:[2,7], 21:[6,11]} )
        and the reduced_process_pdgs are:
           (11, -11), (1,1,-1,-1,21,1,1,-1,-1,21,21)
        Then the list of permuations returned correspond to all possible ways of 
        assigning d-quarks at positions 2 and 7 and gluons at position 6 and 11 given the
        flavor content of the reduced process.
        """

        n_initial = len(reduced_process_pdgs[0])
        initial_state_parent_flavors = parent_flavors[0]
        final_state_parent_flavors = parent_flavors[1]

        # First build the list of positions of each flavor in the reduced process
        initial_flavor_positions = {}
        for i, flavor in enumerate(reduced_process_pdgs[0]):
            try:
                initial_flavor_positions[flavor].append(i)
            except KeyError:
                initial_flavor_positions[flavor] = [i,]
        final_flavor_positions = {}
        for i, flavor in enumerate(reduced_process_pdgs[1]):
            try:
                final_flavor_positions[flavor].append(i+n_initial)
            except KeyError:
                final_flavor_positions[flavor] = [i+n_initial,]

        # Now for each parent flavor, build all ordered N-list of matching flavors in the
        # counterterm reduced process, where N is the number of parents of that flavor.
        # Example; if the real-emission process is
        #    e+(1) e-(2) > d(3) d~(4) d(5) d~(6) g(7)  
        # And the integrated counterterm is C(3,4)C(5,6), then the reduced process is:
        #    e+(1) e-(2) > g(8) g(9) g(7) 
        # And we need to distribute this integragrated counterterm of *all* assignements 
        # of the two (N=2) gluon momenta of leg 8 and 9 in the virtual process:
        #    e+(1) e-(2) > g(3) g(4) g(5)
        # Which leads to the following permutations: (remember that 8,9,7 will be automatically
        # considered as consecutive leg numbers in 'get_reduced_kinematics')
        #  (3,4,5), (3,5,4), (4,3,5), (4,5,3), (5,3,4), (5,4,3)
        # We do this for initial and final states separately
        initial_mappings = [{}]
        for flavor, parent_positions in initial_state_parent_flavors.items():
            picks = list( itertools.permutations(initial_flavor_positions[flavor], len(parent_positions)) )
            initial_mappings_for_this_flavor = []
            for pick in picks:
                mapping = {}
                # List all flavor position not selected for being parents
                unused_flavor_positions = [pos for pos in initial_flavor_positions[flavor] if pos not in pick]
                # loop over all flavor positions to map
                for flavor_pos in initial_flavor_positions[flavor]:
                    if flavor_pos in parent_positions:
                        # if the position of this flavor is one of the 'singular' parent ones,
                        # then assign it to what is selected in that specific pick.
                        destination = pick[parent_positions.index(flavor_pos)]
                        mapping[flavor_pos] = destination
                    else:
                        # Otherwise just assign the next available flavor position not selected
                        mapping[flavor_pos] = unused_flavor_positions.pop(0)
                initial_mappings_for_this_flavor.append(mapping)
            # Now combine the mappings found for this flavor with all the ones already
            # found for other flavors
            new_initial_mappings = []
            for initial_mapping in initial_mappings:
                for mapping in initial_mappings_for_this_flavor:
                    copied_mapping = dict(initial_mapping)
                    copied_mapping.update(mapping)
                    new_initial_mappings.append(copied_mapping)
            initial_mappings = new_initial_mappings
        
        final_mappings = [{}]
        for flavor, parent_positions in final_state_parent_flavors.items():
            picks = list( itertools.permutations(final_flavor_positions[flavor], len(parent_positions)) )
            final_mappings_for_this_flavor = []
            for pick in picks:
                mapping = {}
                # List all flavor position not selected for being parents
                unused_flavor_positions = [pos for pos in final_flavor_positions[flavor] if pos not in pick]
                # loop over all flavor positions to map
                for flavor_pos in final_flavor_positions[flavor]:
                    if flavor_pos in parent_positions:
                        # if the position of this flavor is one of the 'singular' parent ones,
                        # then assign it to what is selected in that specific pick.
                        destination = pick[parent_positions.index(flavor_pos)]
                        mapping[flavor_pos] = destination 
                    else:
                        # Otherwise just assign the next available flavor position not selected 
                        mapping[flavor_pos] = unused_flavor_positions.pop(0)
                final_mappings_for_this_flavor.append(mapping)
            # Now combine the mappings found for this flavor with all the ones already
            # found for other flavors
            new_final_mappings = []
            for final_mapping in final_mappings:
                for mapping in final_mappings_for_this_flavor:
                    copied_mapping = dict(final_mapping)
                    copied_mapping.update(mapping)
                    new_final_mappings.append(copied_mapping)
            final_mappings = new_final_mappings
            
        # Now combine initial and final mappings
        all_mappings = []
        for initial_mapping in initial_mappings:
            for final_mapping in final_mappings:
                copied_mapping = dict(initial_mapping)
                copied_mapping.update(final_mapping)
                all_mappings.append(copied_mapping)
        
        return all_mappings

    def combine_initial_state_counterterms(self):
        """ The initial state beam factorization current, including the integrated collinear
        ones, return a flavor matrix that represents all possible backward evolution of the
        initial state. For instance, if the BF1 process is:
              d d~ > z
        then there will be two [C(1,3)] coming from the real processes d d~ > z g and g d~ > z
        which correspond to the two possible backward evolutions to a d-quark and a gluon
        respectively. Given that those are accounted for at once in the beam-factorization
        current via the flavor matrix, we must make sure to combine all such integrated
        initial-state counterterms so as to avoid double-counting.
        Alternatively, one can disable this combination and always backward evolve only to 
        the particular initial-state flavor specified in this particular integrated counterterm.
        
        WARNING: Eventually this entire construction of the 'allowed_backward_evolved_flavors'
        should be dropped in favour of further differentiation of the integrated collinear ISR
        current that should apply only to particular current (for instance the one 
        *only* backward-evolving to a gluon). When this will be done, we will be able to 
        drop this construction altogether.
        """

        for process_key, (process, mapped_processes) in self.get_processes_map().items():
            # It may be that a particular process does not have any counterterms (typically when
            # considering forcing a beam-type by hand at generation time for debugging purposes, for 
            # example:
            #    generate e+ e- > j j j / a QED=2 --NLO=QCD --beam_types=proton@(1,-1,21)
            # In this case, self.integrated_counterterms would be empty and this combination 
            # should be skipped.
            if process_key not in self.integrated_counterterms:
                continue
            for counterterm_characteristics in self.integrated_counterterms[process_key]:
                # One should perform here the combination of the values of each key in
                # self.integrated_counterterms.
            
                # By default include all possible bbackward evolutions.
                counterterm_characteristics['allowed_backward_evolved_flavors'] = \
                                                        {'beam_one': 'ALL', 'beam_two': 'ALL',}
                if _COMBINE_INITIAL_STATE_INTEGRATED_COLLINEAR_COUNTERTERMS:
                    raise MadGraph5Error('Combination of integrated initial-state collinear '+
                                                                    'counterterms not implemented yet')
                else:
                    # Since in this case we don't combine different initial state integrated
                    # collinear counterterms that differ only by the backward evolved flavor,
                    # we must restrict each one to the particular backward evolved flavor they
                    # correspond to.
                    integrated_counterterm = counterterm_characteristics['integrated_counterterm']
                    if integrated_counterterm.does_require_correlated_beam_convolution():
                        # In this case the beam currents has a diagonal flavor matrix anyway
                        continue
                    beam_number_map = {'beam_one': 0, 'beam_two': 1}
                    beam_currents = integrated_counterterm.get_beam_currents()
                    active_beams = set([])
                    # Check if the counterterm contains any current that could have a 
                    # flavor matrix attached to either beam_one of beam_two
                    for bc in beam_currents:
                        for beam_name in beam_number_map:
                            # If it is a pure integrated mass collinear counterterm, then
                            # no masking should be applied
                            if (bc[beam_name] is None) or all(ss.name()=='F' for ss in 
                                      bc[beam_name].get('singular_structure').decompose()):
                                continue
                            # Now we can add this beam as a one for which a masking of the
                            # allowed initial state flavours must be specified.
                            active_beams.add(beam_name)
                    # Now for each beam identified as potentially having a flavour matrix
                    # then identify which flavours must be allowed as allowed ones after
                    # the application of the flavor matrix
                    for beam_name in active_beams:
                        counterterm_characteristics['allowed_backward_evolved_flavors'][beam_name] = \
                            tuple(set(fl[0][beam_number_map[beam_name]] for fl in 
                                counterterm_characteristics['resolved_flavors_combinations']))

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
        # Also overwrite some of the process properties to make it match the 
        # loop process definitions in this contribution
        n_loops_of_reduced_process = self.contribution_definition.n_loops - \
            integrated_counterterm.get_n_loops_in_beam_factorization_of_host_contribution()
        
        counterterm_reduced_process_key = ProcessKey(integrated_counterterm.process,
            sort_PDGs   = sort_PDGs,
            n_loops     = n_loops_of_reduced_process).get_canonical_key()
        
        # misc.sprint(inverse_processes_map.keys(), len(inverse_processes_map.keys()))
        # misc.sprint(counterterm_reduced_process_key)

        # Now look for a matching process key using only a limited set of criteria
        potential_matches = []
        criteria_considered = ['PDGs','n_loops','id']
        for process_key in inverse_processes_map:
            if all(pk[1] == crpk[1] for (pk, crpk) in zip(process_key,
                           counterterm_reduced_process_key) if pk[0] in criteria_considered):
                # Now we can simply add this integrated counterterm in the group of the
                # defining process to which the reduced process is mapped
                potential_matches.append(inverse_processes_map[process_key])

        if len(potential_matches)==0:
            # The reduced process of the integrated counterterm is not included in this
            # contribution so we return here False, letting the ME7 exporter know that
            # we cannot host this integrated counterterm here.
            return False
        
        if len(potential_matches)>1:
            logger.warning(
                ('MadNkLO found the following processes:\n%s\nto be potential matches'+
                 'to this integrated counterterm:\n%s\nThe first matching process will be'+
                 ' selected, but this may indicate that the current list of criteria considered'+
                 ' (%s) should be extended.')%(
                     '\n'.join(p[1].nice_string() for p in potential_matches),
                     integrated_counterterm.nice_string(), ','.join(criteria_considered)) )

        defining_key, virtual_process_instance = potential_matches[0]
        
        integrated_counterterm_properties = dict(integrated_CT_properties)
       
        # Given the reduced flavors, find all ways of distributing them over the (mapped) flavors
        # of the virtual process. Each possible way gives a specific permutations of the input
        # kinematics that will need to be considered when evaluating the integrated counterterm.
        # These permutations are listed in the variable "permutations" below.
        # Example: 
        #   e+(1) e-(2) > d(3) d~(4) d(5) d~(6) g(7)
        # Then the counterterms C(3,7) and C(5,7) have the same topology and C(5,7) is
        # be discarded because the counterterm C(3,7) will be distributed on both
        # d-quark #3 and #5 of the following virtual process:
        #  e+(1) e-(2) > d(3) d~(4) d(5) d~(6)
        permutations = []
        
        # First compute the base_permutation which is one particular permutation that correctly
        # sends legs of a given flavor from the virtual process to those of the integrated
        # counterterm reduced process. For example, if the virtual process is:
        #     e+(1) e-(2) > d(3) d(4) d~(5) d~(6)
        # and the reduced process of the integrated counterterm is:
        #     e+(1) e-(2) > d(3) d~(4) d(5) d~(6)
        # Then the base_permuation will be
        #     {0:0, 1:1, 2:2, 3:4, 4:3, 5:5}
        origin_pdgs_order = (virtual_process_instance.get_initial_ids(), 
                             virtual_process_instance.get_final_ids_after_decay())
        destination_pdgs_order = (integrated_counterterm.process.get_initial_ids(),
                                  integrated_counterterm.process.get_final_ids_after_decay())
        # Obtain the corresponding permutation
        basic_permutation = self.get_basic_permutation(origin_pdgs_order,destination_pdgs_order)

        # If the integrated counterterm requires only a single active beam factorization then
        # this beam-factorization contribution should indeed host this counterterm only if:
        #  a) The initial state legs *does* require swapping and the factorized beam number of
        #     this integrated CT does *not* match the factorized beam number of thhis contribution
        #  b) The initial state legs does *not* require swapping and the factorized beam number 
        #     of this integrated CT *does* match the factorized beam number of this contribution
        active_beams = {k:self.contribution_definition.is_beam_active(k) for k in ['beam_one','beam_two'] }
        integrated_CT_active_beams = integrated_counterterm.get_necessary_beam_convolutions()
        if active_beams.values().count(True)==1 and len(integrated_CT_active_beams)==1:
            require_initial_state_swapping = ( basic_permutation[0] == 1 and 
                                               basic_permutation[1] == 0 )
            integrated_CT_active_beam = list(integrated_CT_active_beams)[0]
            if (require_initial_state_swapping and active_beams[integrated_CT_active_beam]==True) or \
                ((not require_initial_state_swapping) and active_beams[integrated_CT_active_beam]==False):
                # Then this integrated CT should not be hosted in this contribution
                return False
        
        # Now we need to distribute all flavors that are parents of the integrated current
        # in all possible ways. In the example below, the parent of C(3,7) is a d-quark and
        # we must assign it as both leg #3 and leg #5 of the following reduced integrated
        # counterterm process:
        #          e+(1) e-(2) > d(3) d~(4) d(5) d~(6)
        #
        # First, we must get what are the "singular" parent flavors, and their position
        # in the list of external legs of the reduced process.
        initial_state_parent_flavors = {}
        final_state_parent_flavors = {}
        n_initial = len(destination_pdgs_order[0])
        momenta_dict = integrated_counterterm.momenta_dict
        for i, leg in enumerate(integrated_counterterm.process.get_initial_legs()):
            if momenta_dict[leg.get('number')] != frozenset([leg.get('number'),]):
                try:
                    initial_state_parent_flavors[leg.get('id')].append(i)
                except KeyError:
                    initial_state_parent_flavors[leg.get('id')] = [i,]
        for i, leg in enumerate(integrated_counterterm.process.get_final_legs()):
            if momenta_dict[leg.get('number')] != frozenset([leg.get('number'),]):
                try:
                    final_state_parent_flavors[leg.get('id')].append(i+n_initial)
                except KeyError:
                    final_state_parent_flavors[leg.get('id')] = [i+n_initial,]

        # Now obtain all possible ways of distributing these parent flavors over the
        # flavors of the underlying reduced process.
        # Initial-state are always distinguishible (and they were treated as so when 
        # generating the list of integrated counterterms to dispatch) so that we must not
        # distribute over identical initial_state_parent_flavors:
        initial_state_parent_flavors = {}
        flavor_permutations = self.distribute_parent_flavors(
            (initial_state_parent_flavors, final_state_parent_flavors), destination_pdgs_order)

        # Now combine these permutations with the basic permutation to obtain the final
        # list of permutations to consider
        for flavor_permutation in flavor_permutations:
            combined_permuation = {}
            for index in basic_permutation.keys():
                combined_permuation[index] = basic_permutation[flavor_permutation.get(index, index)]
            permutations.append(combined_permuation)

#       Example of how to print out information of a specific contribution.
#        if isinstance(self, Contribution_RV) and False:
#            misc.sprint('-'*50)
#            misc.sprint(momenta_dict)
#            misc.sprint(initial_state_parent_flavors)
#            misc.sprint(final_state_parent_flavors)
#            misc.sprint(destination_pdgs_order)            
#            misc.sprint(integrated_counterterm.process.nice_string())
#            misc.sprint('For integrated counterterm %s, flavor permutations are: %s'%(str(integrated_counterterm),str(permutations)))
#            misc.sprint('-'*50)

        # Store the mapping to apply to the virtual ME inputs
        integrated_counterterm_properties['input_mappings'] = permutations
        
        # Sanity check that is anticpated to always hold. Powerful non-trivial check
        symmetry_factor = integrated_counterterm_properties['symmetry_factor']
        reproduced_symmetry_factor = integrated_counterterm_properties['multiplicity']/float(len(permutations))
        if symmetry_factor != reproduced_symmetry_factor:
            raise MadGraph5Error(
"The overall symmetry factor of integrated counterterm %s"%str(integrated_counterterm)+
" does not match when computed from symmetry factors (%f) and from CT multiplicity (%f=%d/%d)."
%(symmetry_factor,reproduced_symmetry_factor,integrated_counterterm_properties['multiplicity'],len(permutations))+
"\nFlavors assignments are:\n"+'\n'.join('%s => %s'%(str(k),str(v)) for k, v in 
               integrated_counterterm_properties['resolved_flavors_combinations'].items()))

        # For safety in case the statement above is commented, force the symmetry factor
        # to be the one recomputed here
        integrated_counterterm_properties['symmetry_factor'] = reproduced_symmetry_factor

        self.integrated_counterterms[defining_key].append(integrated_counterterm_properties)
        
        # Let the ME7 exporter that we could successfully host this integrated counterterm
        return True

class Contribution_RV(Contribution_R, Contribution_V):
    """ Implements the general handling of contribution with arbitrary number of reals, virtuals
    and beam factorization contributions."""
    
    def get_integrands_for_process_map(self, process_map, model, run_card, all_MEAccessors, ME7_configuration):
        """ Returns all the integrands implementing this contribution for the specified process_map.
        The instance of MEAccessorDict is necessary so as to be passed to the integrand instances.
        """
        
        relevant_counterterms = {}
        self.remove_local_counterterms_set_to_zero(all_MEAccessors)
        if self.counterterms:
            for process_key in process_map:
                relevant_counterterms[process_key] = self.counterterms[process_key]
        
        relevant_integrated_counterterms = {}
        self.remove_integrated_counterterms_set_to_zero(all_MEAccessors)
        if self.integrated_counterterms:
            for process_key in process_map:
                relevant_integrated_counterterms[process_key] = self.integrated_counterterms[process_key]

        if self.IR_subtraction is None:
            identified_sectors = None
        else:
            identified_sectors = self.IR_subtraction.get_sectors(
                self.contribution_definition,
                process_map,
                counterterms=relevant_counterterms,
                integrated_counterterms=relevant_integrated_counterterms
            )

        return [ ME7_integrands.ME7Integrand(model, run_card,
                   self.contribution_definition,
                   process_map,
                   self.topologies_to_processes,
                   self.processes_to_topologies,
                   all_MEAccessors,
                   ME7_configuration,
                   sectors=identified_sectors,
                   counterterms=relevant_counterterms,
                   integrated_counterterms=relevant_integrated_counterterms
                ) ]

    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this contribution."""
        res = Contribution_R.get_additional_nice_string_printout_lines(self, format=format)
        res.extend(Contribution_V.get_additional_nice_string_printout_lines(self, format=format))
        return res

    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        contribution."""
        GREEN = '\033[92m'
        ENDC = '\033[0m'        
        res = GREEN+'  %s'%(defining_process.nice_string(print_weighted=False)
                                                             .replace('Process: ',''))+ENDC

        if self.integrated_counterterms is None or self.counterterms is None:
            return res

        if format<2:
            if process_key in self.counterterms:
                res += ' | %d local counterterms'%len([ 1
                    for CT in self.counterterms[process_key] if CT.is_singular() ])
            else:
                res += ' | 0 local counterterm'                
            if process_key in self.integrated_counterterms:
                res += ' and %d integrated counterterms'%len(self.integrated_counterterms[process_key])
            else:
                res += ' and 0 integrated counterterm'
        else:
            long_res = [' | with the following local and integrated counterterms:']
            for CT in self.counterterms[process_key]:
                if CT.is_singular():
                    if format==2:
                        long_res.append( '   | %s' % str(CT))
                    elif format==3:
                        long_res.append( '   | %s' % CT.__str__(
                            print_n=True, print_pdg=True, print_state=True ) )
                    elif format>3:
                        long_res.append(CT.nice_string("   | "))

            for CT_properties in self.integrated_counterterms[process_key]:
                CT = CT_properties['integrated_counterterm']
                if format==2:
                    long_res.append( '   | %s' % str(CT))
                elif format==3:
                    long_res.append( '   | %s' % CT.__str__(
                        print_n=True, print_pdg=True, print_state=True ))
                elif format==4:
                    long_res.append(CT.nice_string("   | "))
                elif format>4:
                    long_res.append(CT.nice_string("   | "))
                    for key, value in CT_properties.items():
                        if not key in ['integrated_counterterm', 'matching_process_key']:
                            long_res.append( '     + %s : %s'%(key, str(value)))
            res += '\n'.join(long_res)

        return res


class Contribution_BS(Contribution_RV):
    """ A class implementing the Beam-Soft factorization contributions (BS) originating
    from soft counterterms in the colorful currents scheme which can recoil against the
    initial states, as when using the SoftBeamsRecoilNLOWalker mappings."""

    def __init__(self, contribution_definition, cmd_interface, **opts):
        """ Make sure that this contribution is instantiated by a contribution definition
        that enforces a correlated two-beams convolution."""
        
        if not contribution_definition.correlated_beam_convolution:
            raise MadGraph5Error("A soft beam convolution contribution must be instantiated"+
                " from a contribution definition that specifies a *correlated* convolution "+
                "of the two beams.")
        super(Contribution_RV,self).__init__(contribution_definition, cmd_interface, **opts)

    def add_beam_factorization_processes_from_contribution(self, contrib, beam_factorization_order):
        """ 'Import' the processes_map and all generated attributes of the contribution
        passed in argument to assign them to this beam factorization contributions. Contrary
        to the implementation of this function in the base class, here we do not instantiate
        'bulk' beam factorization currents because they are not needed. The BS contribution is
        only meant to be a receptacle for the integrated soft counterterm in the 'colorful'
        currents scheme when the 'ppToOneWalker' mappings scheme is used (where the soft
        recoils equally against both initial states).
        """

        # Make sure the source contribution does have a factorizing beam for the active
        # ones of this contribution.
        if ( contrib.contribution_definition.beam_factorization['beam_one'] is None ) or \
           ( contrib.contribution_definition.beam_factorization['beam_two'] is None ):
            return
    
        # Now fetch the process map to load from
        source_processes_map = contrib.get_processes_map()
        processes_map = self.processes_map[1]
        # To the new contribution instances, also assign references to the generated attributes
        # of the original contributions (no deep copy necessary here)
        for process_key, (defining_process, mapped_processes) in source_processes_map.items():
            
            # First add this process to the current processes map
            if process_key in processes_map:
                raise MadGraph5Error('MadNkLO attempts to add beam factorization terms for'+
                    ' the %s a second time. This should never happen.'%(
                             defining_process.nice_string().replace('Process','process')))
            processes_map[process_key] = (
                  defining_process.get_copy(), [mp.get_copy() for mp in mapped_processes] )
            
            # Contrary to what is done in the base class, no 'bulk' beam factorization 
            # current need to be added here
            
            # Also add the corresponding topologies
            self.import_topologies_from_contrib(process_key, contrib)

    def generate_all_counterterms(self, ignore_integrated_counterterms=False, group_processes=True, *args, **opts):
        """ Contrary to the other contributions the BS one has no physical contribution that require subtraction
        and is only a place-holder to receive the integrated soft counterterms. Therefore nothing needs to be
        done a this stage here. """
    
        for process_key, (defining_process, mapped_processes) in self.get_processes_map().items():
            self.counterterms[process_key] = []
        all_integrated_counterterms = []
        return all_integrated_counterterms

class ContributionList(base_objects.PhysicsObjectList):
    """ A container for storing a list of contributions."""
    
    contributions_natural_order = [
        ('LO',    (Contribution_B,) ),
        ('NLO',   (Contribution_R, Contribution_V, Contribution_RV) ),
        ('NNLO',  (Contribution_RR, Contribution_RV) ),
        ('NNNLO', (Contribution_RRR, Contribution_RV) )
    ]
    
    def is_valid_element(self, obj):
        """Test if object obj is a valid instance of Contribution."""
        return isinstance(obj, Contribution)
    
    def get_loop_induced_contributions(self):
        """ Return a list of loop-induced contributions."""
        return ContributionList([contrib for contrib in self if
                                        contrib.contribution_definition.is_loop_induced()])
    
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
        return ContributionList([contrib for contrib in self if type(contrib) in correction_classes])

    def nice_string(self, format=0):
        """ A nice representation of a list of contributions. 
        We can reuse the function from ContributionDefinitions."""
        return base_objects.ContributionDefinitionList.contrib_list_string(self, format=format)

    def get_max_correction_order(self):
        """ Returns the maximum correction order found  among these contributions."""
        
        return 'N'*(max(contrib.contribution_definition.correction_order.count('N')
                                                                 for contrib in self))+'LO'

    def sort_contributions(self):
        """ Sort contributions according to the order dictated by the class_attribute contributions_natural_order"""
        
        new_order = []
        for correction_order, contribution_types in self.contributions_natural_order:
            for contrib_type in contribution_types:
                selected_contribs = self.get_contributions_of_order(correction_order).\
                                        get_contributions_of_type(contrib_type)
                # Now sort contributions of this type using the ordering F1, F2, F1F2
                factorization_order = [(False,False),(True, False),(False, True), (True, True)]
                new_order.extend(sorted(selected_contribs,
                    key = lambda contrib: factorization_order.index((
                        contrib.contribution_definition.is_beam_active('beam_one'),
                        contrib.contribution_definition.is_beam_active('beam_two')
                    ))
                ))
                for contrib in selected_contribs:
                    self.pop(self.index(contrib))
    
        # Finally all remaining contributions of unknown types
        new_order.extend(self)
        self[:] = new_order

    def apply_method_to_all_contribs(self, method, log=None, method_args=[], method_opts={}):
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
                    # If method opts depend on the contribution, it can be passed as a callable
                    # that can be evaluated here
                    if not isinstance(method_opts, dict) and callable(method_opts):
                        method_options = method_opts(contrib)
                    else:
                        method_options = dict(method_opts)
                    return_values.append((contrib, contrib_function(*method_args, **method_options)))
                    remaining_contribs.pop(remaining_contribs.index(contrib))

        if remaining_contribs:
            if log: logger.info('%s the %d contribution%s of unknown type...'%(log,
                        len(remaining_contribs), 's' if len(remaining_contribs)>1 else ''))                    
        for i, contrib in enumerate(remaining_contribs):
            if log: logger.info('%s %d/%d'%(type(contrib).__name__, i+1, len(remaining_contribs)))
            try:
                contrib_function = getattr(contrib, method)
            except AttributeError:
                raise MadGraph5Error("The contribution\n%s\n does not have function '%s' defined."%(
                                                                            contrib.nice_string(), method))
            # If method opts depend on the contribution, it can be passed as a callable
            # that can be evaluated here
            if not isinstance(method_opts, dict) and callable(method_opts):
                method_options = method_opts(contrib)
            else:
                method_options = dict(method_opts)
            return_values.append((contrib, contrib_function(*method_args, **method_options)))
        
        return return_values

# Contribution classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own class of Contribution.
# Notice that this Contribution must be placed after all the Contribution daughter classes have been declared.
Contribution_classes_map = {'Born': Contribution_B,
                            'Virtual': Contribution_V,
                            'SingleReals': Contribution_R,
                            'RealVirtual': Contribution_RV,
                            'DoubleReals': Contribution_RR,
                            'TripleReals': Contribution_RRR,
                            'BeamSoft': Contribution_BS,
                            'Unknown': None}
