##########################################################################################
#
# Copyright (c) 2011 The MadGraph5_aMC@NLO Development team and Contributors
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
"""Collection of ME7 Integrand classes."""

import numpy as np
import collections
import glob
import logging
import math
import os
import random
import re
import math
import itertools

import stat
import subprocess
import sys
import time
import tarfile
import StringIO
import shutil
import copy

try:
    import readline
    GNU_SPLITTING = ('GNU' in readline.__doc__)
except:
    GNU_SPLITTING = True

try:
    import pyjet
    PYJET_AVAILABLE = True
except:
    PYJET_AVAILABLE = False    

# useful shortcut
pjoin = os.path.join

#root_path = pjoin(os.path.dirname(os.path.realpath( __file__ )), *([os.path.pardir]*3))
#sys.path.insert(0, root_path)

# Special logger for the Cmd Interface
logger = logging.getLogger('madevent7') # -> stdout
logger_stderr = logging.getLogger('madevent7.stderr') # ->stderr

import madgraph.core.base_objects as base_objects
import madgraph.core.subtraction as subtraction
import madgraph.interface.extended_cmd as cmd
import madgraph.interface.common_run_interface as common_run
import madgraph.interface.madevent_interface as madevent_interface
import madgraph.iolibs.files as files
import madgraph.iolibs.save_load_object as save_load_object
import madgraph.various.banner as banner_mod
import madgraph.various.cluster as cluster
import madgraph.various.misc as misc
import madgraph.various.lhe_parser as lhe_parser
import madgraph.integrator.integrands as integrands
import madgraph.integrator.integrators as integrators
import madgraph.integrator.walkers as walkers
import madgraph.integrator.phase_space_generators as phase_space_generators
import madgraph.integrator.pyCubaIntegrator as pyCubaIntegrator
import madgraph.integrator.vegas3_integrator as vegas3_integrator
#    import madgraph.various.histograms as histograms  # imported later to not slow down the loading of the code
import models.check_param_card as check_param_card
import models.model_reader as model_reader
import models.import_ufo as import_ufo

from madgraph.iolibs.files import ln    
from madgraph import InvalidCmd, MadGraph5Error, MG5DIR, ReadWrite



class MadEvent7Error(Exception):
    pass

class ZeroResult(MadEvent7Error):
    pass

#===============================================================================
# ME7Events
#===============================================================================
class ME7Event(object):
    """ A class to store a particular configuration produced by the integrands
    when evaluated and which has undefinite flavor configuration to begin with."""
    
    def __init__(self, PS_point, weights_per_flavor_configurations, 
                 is_mirrored                  = False,
                 host_contribution_definition = 'N/A',
                 counterterm_structure        = None,
                 Bjorken_x_rescalings         = (1.0,1.0) ):
        """ Initilize this particular event with a unique PS point and a dictionary of
        the form:
             ( (initial_state_flavors,), (final_state_flavors,) ) : weight
        to specify all possible flavor configurations applying to this event.
        For PDF counterterms and integrated collinear CT counterterms, we must
        provide the possibility of applying a rescaling to the Bjorken x that will be
        used to compute each of the two beams. By default this rescaling is of course one."""
        
        self.PS_point = PS_point.to_list()
        self.weights_per_flavor_configurations = weights_per_flavor_configurations
        self.is_mirrored = is_mirrored
        self.host_contribution_definition = host_contribution_definition
        self.counterterm_structure = counterterm_structure
        self.Bjorken_x_rescalings = Bjorken_x_rescalings
        
        # This entry specifies which flavor has been selected for this event.
        # None implies that it has not been selected yet.
        self.selected_flavors = None
        
    def get_copy(self):
        """ Returns a shallow copy of self, except for the 
        'weights_per_flavor_configurations' attribute. """
        
        return ME7Event(
            self.PS_point, dict(self.weights_per_flavor_configurations), 
            is_mirrored                  = self.is_mirrored,
            host_contribution_definition = self.host_contribution_definition,
            counterterm_structure        = self.counterterm_structure,
            Bjorken_x_rescalings         = self.Bjorken_x_rescalings
        )
        
    def set_Bjorken_rescalings(self, *args):
        """ Assigns the specified Bjorken x rescaling to this event."""

        self.Bjorken_x_rescalings = tuple(args)

    def apply_PDF_convolution(self, pdf_accessor, all_pdf, all_x, all_mu_f_squared):
        """ Convolute the weights_per_flavor_configurations."""
        
        # In order for the + distributions of the PDF counterterms and integrated
        # collinear ISR counterterms to hit the PDF only (and not the matrix elements or
        # observables functions), a change of variable is necessary: xb_1' = xb_1 * chsi1
        # In this case, the chsi1 to factor to apply is given by the Bjorken_x_rescalings
        # factor, which must be used both for rescaling the argument of the PDF as well
        # as bringing a factor 1/\chsi for the jacobian of the change of variable.
        for flavors in self.weights_per_flavor_configurations:
            PDFs = 1.
            for i, flavor in enumerate(flavors[0]):
                PDFs *= pdf_accessor(all_pdf[i], flavor, 
                                     all_x[i]*self.Bjorken_x_rescalings[i], all_mu_f_squared[i])
                PDFs *= self.Bjorken_x_rescalings[i]
            self.weights_per_flavor_configurations[flavors] *= PDFs

    def get_total_weight(self):
        """ Return the total weight of this event for all flavor considerations to consider."""

        summed_weight = None
        for wgt in self.weights_per_flavor_configurations.values():
            if summed_weight is None:
                summed_weight = copy.copy(wgt)
            else:
                summed_weight += wgt
        if summed_weight is None:
            return 0.
        if self.is_mirrored:
            return summed_weight*2.
        else:
            return summed_weight

    def filter_flavor_configurations(self, flavor_cut, **opts):
        """ Apply the flavor cut on this event for each flavor configuration, removing
        all those failings. If there is None left, this returns False."""
        
        new_weights_per_flavor_configurations = {}
        for flavors, weight in self.weights_per_flavor_configurations.items():
            if flavor_cut(self.PS_point, flavors, **opts):
                new_weights_per_flavor_configurations[flavors] = weight
        
        self.weights_per_flavor_configurations = new_weights_per_flavor_configurations
        
        return len(self.weights_per_flavor_configurations) != 0

    def apply_flavor_blind_cuts(self, flavor_blind_cut, *args, **opts):
        """ Apply the flavor-blind cut to this event, returning False if it failed."""
        
        return flavor_blind_cut(self.PS_point, *args, **opts)

    def select_a_flavor_configuration(self):
        """ Select a particular flavor configuration for this event which may be used
        by the flavor-sensitive observables. The selection is performed randomly with 
        a probability proportional to the weight of each flavor configuration."""

        weights_by_flavor_configurations = self.weights_per_flavor_configurations.items()
        index_selected=0
        abs_flavor_wgts_sum = sum(abs(value) for value in 
                                         self.weights_per_flavor_configurations.values())
        running_sum = abs(weights_by_flavor_configurations[index_selected][1])
        rv_flavor = random.random()*abs_flavor_wgts_sum
        while rv_flavor > running_sum:
            index_selected +=1
            running_sum += abs(weights_by_flavor_configurations[index_selected][1])
        self.selected_flavors = weights_by_flavor_configurations[index_selected][0]

    def select_epsilon_expansion_term(self, term_name):
        """ Select a particular EpsilonExpansion term in the weights of the flavor 
        matrix of this event. If the weights are float already, then they will be
        set to zero unless the term_name is 'finite'."""
        
        is_finite_specified = (term_name.lower()=='finite')
        for flavors, weight in self.weights_per_flavor_configurations.items():
            if isinstance(weight, float):
                if not is_finite_specified:
                    del self.weights_per_flavor_configurations[flavors]
            else:
                self.weights_per_flavor_configurations[flavors] = weight.get_term(term_name)

    def convolve_flavors(self, flavor_matrix, leg_index, initial_state=True):
        """ Convolves the flavor matrix in argument (of the form: 
            {   start_flavor_a :  {  end_flavor_1 : wgt_x,
                                     end_flavor_2 : wgt_y,
                                     [...]
                                 },
                start_flavor_b :  {  end_flavor_1 : wgt_z,
                                     end_flavor_2 : wgt_w,
                                     [...]
                                 },
                [...]
            }
            where 'wgt_x' can be EpsilonExpansion instances.
            and modifies self.weights_per_flavor_configurations with the convolved result.
            leg index specifies which leg must be convolved.
               (e.g. beam #1 correpsonds to leg_index=0 and initial_state = True."""
            
        def substitute_flavor(flavor_config, new_flavor):
            """ Substitute in the flavor config the convoluted flavor with the one in 
            argument and returns the corresponding new flavor configuration as a two tuples,
            for the initial and final states respectively."""
            return ( 
                tuple( (pdg if i!=leg_index else new_flavor) 
                       for i, pdg in enumerate(flavor_config[0]) ) if initial_state else flavor_config[0],
                tuple( (pdg if i!=leg_index else new_flavor) 
                       for i, pdg in enumerate(flavor_config[1]) ) if (not initial_state) else flavor_config[1]
            )

        new_flavor_configurations = {}
        for flavor_config, wgt in self.weights_per_flavor_configurations.items():
            start_flavor = flavor_config[0][leg_index] if initial_state else flavor_config[1][leg_index]
            if start_flavor not in flavor_matrix:
                new_flavor_configs = [(flavor_config, wgt),]
            else:
                new_flavor_configs = []
                for end_flavors, matrix_wgt in flavor_matrix[start_flavor].items():
                    new_flavor_configs.extend([
                    ( substitute_flavor(flavor_config, end_flavor), 
                                base_objects.EpsilonExpansion(matrix_wgt)*wgt )
                                                           for end_flavor in end_flavors ])

            for new_flavor_config, new_wgt in new_flavor_configs:
                try:
                    new_flavor_configurations[new_flavor_config] += new_wgt
                except KeyError:
                    new_flavor_configurations[new_flavor_config] = new_wgt
        
        # Now assign the newly create flavor configurations
        self.weights_per_flavor_configurations = new_flavor_configurations

    def __add__(self, other):
        """ overload the '+' operator."""
        new_event = self.get_copy()
        return new_event.__iadd__(other)
            
    def __iadd__(self, other):
        """ overload the '+=' operator."""
        if __debug__: assert(self.is_mirrored == other.is_mirrored)
        if __debug__: assert(self.PS_point == other.PS_point)
        if __debug__: assert(self.Bjorken_x_rescalings  == other.Bjorken_x_rescalings)
        for flavor_configuration, wgt in other.weights_per_flavor_configurations.items():
            try:
                self.weights_per_flavor_configurations[flavor_configuration] += wgt
            except KeyError:
                self.weights_per_flavor_configurations[flavor_configuration] = wgt
        return self

    def __mul__(self, multiplier):
        """ overload the '*' operator with a float or an EpsilonExpansion. """
        new_event = self.get_copy()
        return new_event.__imul__(multiplier)
    
    def __imul__(self, multiplier):
        """ overload the '*=' operator with a float or an EpsilonExpansion. """
        if __debug__: assert( isinstance(multiplier, (float, base_objects.EpsilonExpansion)) )
        for flavor_configuration, wgt in self.weights_per_flavor_configurations.items():
            try:
                self.weights_per_flavor_configurations[flavor_configuration] = multiplier*wgt
            except KeyError:
                self.weights_per_flavor_configurations[flavor_configuration] = multiplier*wgt
        return self

    def __str__(self):
        return self.nice_string()

    def nice_string(self):
        """ Returns a nice and short string representation of this event."""
        BLUE = misc.bcolors.BLUE
        GREEN = misc.bcolors.GREEN
        ENDC = misc.bcolors.ENDC
        res = []
        if self.counterterm_structure is None:
            res.append("%sEvent from %s contribution '%s'%s:"%(
                GREEN,
                self.host_contribution_definition.correction_order,
                self.host_contribution_definition.short_name(),
                ENDC))
        else:
            if isinstance(self.counterterm_structure, list):
                if len(self.counterterm_structure) > 1:
                    counterterm_structure_str = "( %s )"%(' | '.join(
                                             str(ct) for ct in self.counterterm_structure))
                else:
                    counterterm_structure_str = str(self.counterterm_structure)
            else:
                counterterm_structure_str = str(self.counterterm_structure)            
            res.append('%sCounterterm event from limit%s %s of %s contribution %s%s'%(
                GREEN,
                's' if isinstance(self.counterterm_structure, list) and len(self.counterterm_structure)>1 else '',
                counterterm_structure_str,
                self.host_contribution_definition.correction_order,
                self.host_contribution_definition.short_name(),
                ENDC))
        res.append('%s  Kinematic configuration:%s'%(BLUE, ENDC))
        res.extend('     %s'%line for line in str(self.PS_point).split('\n'))
        if not any(bs is None for bs in self.Bjorken_x_rescalings):
            res.append("     Bjorken x's scalings: %s"%(', '.join('%.16e'%bs for bs in self.Bjorken_x_rescalings)))
        res.append('%s  Flavors configurations:%s'%(BLUE, ENDC))
        for flavors in sorted(self.weights_per_flavor_configurations.keys()):
            wgt = self.weights_per_flavor_configurations[flavors]
            res.append('%s  |  %s'%('%-25s'%('     %s -> %s'%(
                ' '.join('%s'%pdg for pdg in flavors[0]),
                ' '.join('%s'%pdg for pdg in flavors[1]))), 
                wgt.__str__(format='.16e') if isinstance(wgt,base_objects.EpsilonExpansion) 
                                                                         else '%.16e'%wgt ))
        if not self.selected_flavors is None:
            res.append('%s  Selected flavors configuration : %s%s -> %s'%(BLUE,  ENDC,
             ' '.join('%s'%pdg for pdg in self.selected_flavors[0]), 
             ' '.join('%s'%pdg for pdg in self.selected_flavors[1]) ))
        
        return '\n'.join(res)

class ME7EventList(list):
    """ A class handling a collection of ME7Events, with helper function to 
    collectively act on them. This customized list is mainly intended to
    store a list of ME7Events that all originate from the same integrand 
    evaluation call."""
    
    def __init__(self, *args, **opts):
        super(ME7EventList,self).__init__(*args, **opts)

    def apply_PDF_convolution(self, *args,**opts):
        """ Convolute the weights_per_flavor_configurations of each event 
        with the PDF luminosity."""
        for event in self:
            event.apply_PDF_convolution(*args, **opts)

    def filter_flavor_configurations(self, flavor_cut, **opts):
        """ Apply the flavor cut on each event in this list, removing those with
        no valid flavor configuration left. This returns false if there is no events left."""
        
        new_event_list = []
        for event in self:
            if event.filter_flavor_configurations(flavor_cut, **opts):
                new_event_list.append(event)
        self[:] = new_event_list

        return len(self) != 0

    def filter_with_flavor_blind_cuts(self, cut_function, *args, **opts):
        """ Apply the kinematic flavor-blind cut on each event in this list, removing those with
        failing the cut. This returns false if there is no events left."""
        
        new_event_list = []
        for event in self:
            if event.apply_flavor_blind_cuts(cut_function, *args, **opts):
                new_event_list.append(event)
        self[:] = new_event_list

        return len(self) != 0
    
    def get_total_weight(self):
        """ Return the total weight over all events in this list."""
        if len(self) > 0:
            return sum(event.get_total_weight() for event in self)
        else:
            return 0.

    def select_a_flavor_configuration(self):
        """ Select a particular flavor configuration to be used in the flavor-sensitive
        observables for each event. """
        
        for event in self:
            event.select_a_flavor_configuration()
    
    def select_epsilon_expansion_term(self, term_name):
        """ Select a particular EpsilonExpansion term in the weights of the flavor 
        matrix of these events. If the weights are float already, then they will be
        set to zero unless the term_name is 'finite'."""
        
        for event in self:
            event.select_epsilon_expansion_term(term_name)
    
    def apply_observables(self, observable_list, xb_1=None, xb_2=None):
        """ Applies the observable list to the entire set of events of this list."""

        for event in self:
            observable_list.apply_observables(event, xb_1, xb_2)

    def __str__(self):
        return self.nice_string()

    def nice_string(self):
        """ Returns a nice string representation of this list of event."""
        
        BLUE = misc.bcolors.BLUE
        GREEN = misc.bcolors.GREEN
        ENDC = misc.bcolors.ENDC
        
        n_events = len(self)
        res=['%s>> Start of a list of %d event%s:%s'%(
                                        GREEN, n_events, 's' if n_events > 1 else '',ENDC)]
        res.append('')
        for i_event, event in enumerate(self):
             res.append('#%d '%(i_event+1) + event.nice_string())
             res.append('')
        res.append('%s>> End of a list of %d event%s.%s'%(
                                        GREEN, n_events, 's' if n_events > 1 else '',ENDC))
        return '\n'.join(res)

#===============================================================================
# ME7Integrand
#===============================================================================
class ME7Integrand(integrands.VirtualIntegrand):
    """ Specialization for multi-purpose integration with ME7."""
    
    # Maximum size of the cache for PDF calls
    PDF_cache_max_size = 1000
    
    def __new__(cls, model, 
                     run_card,
                     contribution_definition,
                     processes_map,
                     topologies_to_processes,
                     processes_to_topologies,
                     all_MEAccessors,
                     ME7_configuration, **opt):
        all_args = [model, run_card, contribution_definition, processes_map,
                    all_MEAccessors, ME7_configuration ]

        if cls is ME7Integrand:
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
                raise MadGraph5Error("Some %s type of ME7Integrands are not implemented yet."%
                                                  contribution_definition.correction_order)
            target_class = ME7Integrand_classes_map[target_type]

            if not target_class:
                raise MadGraph5Error("Could not determine the class of ME7Integrand of type '%s' to be added for"%target_type+
                                     " the contribution definiton:\n%s"%str(contribution_definition.nice_string()))

            return super(ME7Integrand, cls).__new__(target_class, *all_args, **opt)
        else:
            return super(ME7Integrand, cls).__new__(cls, *all_args, **opt)
    
    def __init__(self, model, 
                       run_card,
                       contribution_definition,
                       processes_map,
                       topologies_to_processes,
                       processes_to_topologies,
                       all_MEAccessors,
                       ME7_configuration):
        """ Initializes a generic ME7 integrand defined by the high-level abstract information
        given in arguments, whose explanations are provided in the comments of this initializer's body. """
        
        super(ME7Integrand, self).__init__()
        
        # Keep the initialization inputs to generate the dump
        # Only keep the relevant information, in and set to None the ones that will need to be updated by
        # ME7 interface directly.
        self.initialization_inputs = { 'model'                      : None,
                                       'run_card'                   : None,
                                       'contribution_definition'    : contribution_definition,
                                       'processes_map'              : processes_map,
                                       'topologies_to_processes'    : topologies_to_processes,
                                       'processes_to_topologies'    : processes_to_topologies,
                                       'all_MEAccessors'            : None,
                                       'ME7_configuration'          : None,
                                       'options'                    : {} }
        
        # The original ContributionDefinition instance at the origin this integrand 
        self.contribution_definition    = contribution_definition
        self.contribution_short_name    = contribution_definition.short_name()

        # The process map of the Contribution instance at the origin of this integrand.
        # The format is identical to the one generated from the function 'get_process_map' of a contribution.
        self.processes_map              = processes_map
        
        # Add information about the topology of the diagrams constituting the processes,
        # so as to be able to build efficient phase-space parametrizations. The format of these dictionaries
        # is specified in the function 'set_phase_space_topologies' of the class contributions.Contribution
        self.topologies_to_processes    = topologies_to_processes
        self.processes_to_topologies    = processes_to_topologies
        
        # An instance of accessors.MEAccessorDict providing access to all ME available as part of this
        # ME7 session.
        self.all_MEAccessors            = all_MEAccessors

        # Update and define many properties of self based on the provided run-card and model.
        self.synchronize(model, run_card, ME7_configuration)

    def has_integrated_counterterms(self):
        """ A property of the class, indicating whether it can contain integrated counterterms."""
        return False

    def has_local_counterterms(self):
        """ A property of the class, indicating whether it can contain local counterterms."""
        return False

    def get_additional_nice_string_printout_lines(self):
        """ Additional printout information for nice_string. 
        Meant to possibly be overloaded by daughter classes."""
        return []

    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        contribution. Can be overloaded by daughter classes."""
        GREEN = '\033[92m'
        ENDC = '\033[0m'        
        return GREEN+'  %s'%defining_process.nice_string(print_weighted=False).\
                                                               replace('Process: ','')+ENDC

    def nice_string(self, format=0):
        """ Nice string representation of self."""
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        res = ['< %s%s%s >'%(BLUE,self.contribution_definition.short_name(),ENDC)]
        res.append('%-30s:   %s'%('ME7Integrand_type',type(self)))
        res.extend([self.contribution_definition.nice_string()])
        if not self.topologies_to_processes is None:
            res.append('%-30s:   %d'%('Number of topologies', 
                                                    len(self.topologies_to_processes.keys())))
        res.extend(self.get_additional_nice_string_printout_lines())

        if format < 1:
            res.append('Generated and mapped processes for this contribution: %d (+%d mapped)'%
                       ( len(self.processes_map.keys()),
                         len(sum([v[1] for v in self.processes_map.values()],[])) ) )
        else:
            res.append('Generated and mapped processes for this contribution:')
            for process_key, (defining_process, mapped_processes) in self.processes_map.items():
                res.append(self.get_nice_string_process_line(process_key, defining_process, format=format))                
                for mapped_process in mapped_processes:
                    res.append(BLUE+u'   \u21b3  '+mapped_process.nice_string(print_weighted=False)\
                                                                        .replace('Process: ','')+ENDC)
            
        return '\n'.join(res).encode('utf-8')

    def __str__(self):
        return self.nice_string()

    def synchronize(self, model, run_card, ME7_configuration):
        """ Synchronize this integrand with the most recent run_card and model."""

        # The option dictionary of ME7
        self.ME7_configuration          = ME7_configuration
        
        # A ModelReader instance, initialized with the values of the param_card.dat of this run
        self.model                      = model
        if not isinstance(self.model, model_reader.ModelReader):
            raise MadGraph5Error("The ME7Integrand must be initialized with a ModelReader instance.")

        # A RunCardME7 instance, properly initialized with the values of the run_card.dat of this run
        self.run_card                   = run_card
        
        # Set external masses
        all_processes = [p[0] for p in self.processes_map.values()]
        self.masses = all_processes[0].get_external_masses(self.model)
        for proc in all_processes[1:]:
            this_proc_masses = proc.get_external_masses(self.model)
            if this_proc_masses != self.masses:
                raise MadGraph5Error("A contribution must entail processes with all the same external masses.\n"
                 "This is not the case; process\n%s\nhas masses '%s' while process\n%s\n has masses '%s'."%
                 (all_processes[0].nice_string(), self.masses, proc.nice_string(), this_proc_masses) )
        self.n_initial = len(self.masses[0])
        self.n_final = len(self.masses[1])
        
        if self.n_initial==1:
            raise InvalidCmd("MadEvent7 does not yet support decay processes.")
        
        if not (self.run_card['lpp1']==self.run_card['lpp2']==1) and \
           not (self.run_card['lpp1']==self.run_card['lpp2']==0):
            raise InvalidCmd("MadEvent7 does not support the following collider mode yet (%d,%d)."%\
                                                            (self.run_card['lpp1'], self.run_card['lpp2']))
        
        # Always initialize the basic flat PS generator. It can be overwritten later if necessary.
        simplified_beam_types = (
            0 if self.contribution_definition.beam_factorization['beam_one'] is None else 1,
            0 if self.contribution_definition.beam_factorization['beam_two'] is None else 1,
        )
        self.phase_space_generator = phase_space_generators.FlatInvertiblePhasespace(
            self.masses[0], self.masses[1],
            beam_Es             = (self.run_card['ebeam1'], self.run_card['ebeam2']),
            beam_types          = simplified_beam_types,
            is_beam_factorization_active = 
                        ( self.contribution_definition.is_beam_active('beam_one'),
                          self.contribution_definition.is_beam_active('beam_two') ) 
        )

        # Add a copy of the PS generator dimensions here.
        # Notice however that we could add more dimensions pertaining to this integrand only, and PS generation.
        # This is in particular true for discrete integration dimension like sectors, helicities, etc...
        integrand_dimensions = integrands.DimensionList(self.phase_space_generator.dimensions)
        self.set_dimensions(integrand_dimensions)
        self.dim_ordered_names = [d.name for d in self.get_dimensions()]
        self.dim_name_to_position = dict((name,i) for i, name in enumerate(self.dim_ordered_names))
        self.position_to_dim_name = dict((v,k) for (k,v) in self.dim_name_to_position.items())

        self.collider_energy = self.run_card['ebeam1'] + self.run_card['ebeam2']
        
        # Set the seed
        if self.run_card['iseed'] > 0:
            random.seed(self.run_card['iseed'])

        # Initialize the PDF, if necessary
        # Setup the PDF cache
        self.PDF_cache = {}
        self.PDF_cache_entries = []
        if self.run_card['lpp1']==0 and self.run_card['lpp2']==0:
            self.pdf = None 
            self.pdfsets = None
        else:
            if self.run_card['pdlabel'] != 'lhapdf':
                raise InvalidCmd("MadEvent7 does not support built-in PDFs.")
            # Load the PDFs
            if self.ME7_configuration['lhapdf']:
                lhapdf_config = self.ME7_configuration['lhapdf']
            else:
                lhapdf_config = misc.which('lhapdf-config')
            lhapdf = misc.import_python_lhapdf(lhapdf_config)
            if not lhapdf:
                raise MadGraph5Error("The python lhapdf API could not be loaded.")
            # Adjust LHAPDF verbosity to current logger's verbosity
            # Ask for logging.DEBUG-1 so as to only have lhapdf verbose if really desired.
            lhapdf.setVerbosity(1 if logger.level<=(logging.DEBUG-1) else 0)

            pdfsets_dir = subprocess.Popen([lhapdf_config,'--datadir'],\
                                           stdout=subprocess.PIPE).stdout.read().strip()
            lhapdf.pathsPrepend(pdfsets_dir)
            lhapdf_version = subprocess.Popen([lhapdf_config,'--version'],\
                                           stdout=subprocess.PIPE).stdout.read().strip()
            pdf_info = common_run.CommonRunCmd.get_lhapdf_pdfsets_list_static(pdfsets_dir, lhapdf_version)
            lhaid = self.run_card.get_lhapdf_id()
            if lhaid not in pdf_info:
                raise InvalidCmd("Could not find PDF set with lhaid #%d in %s."%(lhaid, pdfsets_dir))
            pdf_set_name = pdf_info[lhaid]['filename']
            if not os.path.isdir(pjoin(pdfsets_dir, pdf_set_name)):
                raise InvalidCmd("Could not find PDF set directory named "+
                    "'%s' in '%s'.\n"%(pdf_set_name, pdfsets_dir)+
                    "It can be downloaded from LHAPDF official online resources.")

            self.pdfsets = lhapdf.getPDFSet(pdf_info[lhaid]['filename'])
            # Pick the central PDF for now
            self.pdf = self.pdfsets.mkPDF(0)

        # Initialize access to a mean of running alpha_s
        if self.pdf:
            self.alpha_s_runner = self.pdf.alphasQ
        else:
            # We assume here that the model passed to this integrand for initialization started off
            # with its couplings (i.e. alpha_s) defined for Q**2= MZ**2.
            model_param_dict = model.get('parameter_dict')
            as_running_params = {'aS':0.118, 'mdl_MZ':91.188, 'mdl_MC':1.55, 'mdl_MB':4.7}
            for param in as_running_params:
                if param not in model_param_dict:
                    if param in ['mdl_MC', 'mdl_MB']:
                        # Leave these parameters to their default if not specified in the model
                        continue
                    else:
                        raise InvalidCmd("When not using PDFsets, MadEvent7 requires a model with the"+
                                    " parameter %s to be defined so as to be able to run alpha_S."%param)
                if model_param_dict[param] != 0.:
                    as_running_params[param] = model_param_dict[param]
            # For now always chose to run alpha_S at two loops.
            n_loop_for_as_running = 2
            self.alpha_s_runner = model_reader.Alphas_Runner(as_running_params['aS'], n_loop_for_as_running, 
                      as_running_params['mdl_MZ'], as_running_params['mdl_MC'], as_running_params['mdl_MB'])
            

    def generate_dump(self):
        """ Generate a serializable dump of self, which can later be used, along with some more 
        information, in initialize_from_dump in order to regenerate the object."""
        
        dump = {'class': self.__class__}
        dump.update(self.initialization_inputs)
        return dump
    
    @classmethod
    def initialize_from_dump(cls, dump, model, run_card, all_MEAccessors, ME7_configuration):
        """ Initialize self from a dump and possibly other information necessary for reconstructing this
        integrand."""
        
        return cls(  model, 
                     run_card,                             
                     dump['contribution_definition'],
                     dump['processes_map'],
                     dump['topologies_to_processes'],
                     dump['processes_to_topologies'],
                     all_MEAccessors,
                     ME7_configuration,
                     **dump['options'])
    
    def set_phase_space_generator(self, PS_generator):
        """ Overwrites current phase-space generator."""
        if not isinstance(PS_generator, phase_space_generators.VirtualPhaseSpaceGenerator):
            raise MadGraph5Error("Cannot assign to a MadEvent7 integrand a phase-space generator that "+
                                 " does not inherit from VirtualPhaseSpaceGenerator.")
        if PS_generator.nDimPhaseSpace() != self.phase_space_generator.nDimPhaseSpace():
            raise MadGraph5Error("A MadEvent7 integrand was assigned a phase-space generator with the"+
                                 " wrong number of integration dimensions: %d instead of %d"%
                (PS_generator.nDimPhaseSpace(),self.phase_space_generator.nDimPhaseSpace()))
        self.phase_space_generator = PS_generator

    def is_part_of_process_selection(self, process_list, selection=None):
        """ Checks whether any of the specified processes in the process_list provided matches the user's process
        selection. If not provided, returns True by default. 'selection' is a dictionary with the format:
           {'in_pdgs'  : ( (in_pgs1), (in_pdgs2), ...)
            'out_pdgs' : ( (out_pgs1), (out_pdgs2), ...) 
            'n_loops'  : n_loops }"""
        
        def pdg_list_match(target_list, selection_list):
            if len(target_list) != len(selection_list):
                return False
            targets = dict( (k, target_list.count(k)) for k in set(target_list) )
            found_it = False
            for sel in itertools.product(*selection_list):
                found_it = True
                for k, v in targets.items():
                    if sel.count(k) != v:
                        found_it = False
                        break
                if found_it:
                    break
            return found_it                
    
        for process in process_list:

            if (not selection['in_pdgs'] is None) and \
               (not pdg_list_match(process.get_initial_ids(), selection['in_pdgs'])):
                continue

            if (not selection['out_pdgs'] is None) and \
               (not pdg_list_match(process.get_final_ids_after_decay(), selection['out_pdgs'])):
                continue
            
            if (not selection['n_loops'] is None) and process.get('n_loops') != selection['n_loops']:
                continue
            return True

        return False

    def find_counterterms_matching_regexp(self,
                                    counterterms, limit_pattern=None, integrated_CT=False):
        """ Find all mappings that match a particular limits given in argument
        (takes a random one if left to None). This function is placed here given that
        it can be useful for both the ME7Integrnd_V and ME7_integrand_R."""

        # First select only the counterterms which are not pure matrix elements 
        # (i.e. they have singular structures) and also exclude here soft-collinear 
        # counterterms since they do not have an approach limit function.
        if not integrated_CT:
            selected_counterterms = [ ct for ct in counterterms if ct.is_singular() ]
        else:
            selected_counterterms = [ ct for ct in counterterms if ct['integrated_counterterm'].is_singular() ]

        if len(selected_counterterms)==0:
            return []
        
        limit_pattern_re = None
        if limit_pattern:
            if limit_pattern.lower() == 'soft':
                limit_pattern_re = re.compile(r'.*S.*')
            elif limit_pattern.lower() == 'collinear':
                limit_pattern_re = re.compile(r'.*C.*')
            elif limit_pattern.lower() == 'puresoft':
                limit_pattern_re = re.compile(r'^[S\d,\(\)]*$')
            elif limit_pattern.lower() == 'purecollinear':
                limit_pattern_re = re.compile(r'^[C\d,\(\)]*$')
            elif limit_pattern.lower() == 'all':
                limit_pattern_re = re.compile(r'.*')
            elif any(limit_pattern.startswith(start) for start in ['r"', "r'"]):
                limit_pattern_re = re.compile(eval(limit_pattern))
            else:
                # Check if a list of counterterms is specified
                try:
                    list_limit_pattern = eval(limit_pattern)
                    if not isinstance(list_limit_pattern, list):
                        raise BaseException
                except:
                    list_limit_pattern = [limit_pattern]
                new_list_limit_pattern = []
                for limit_pattern in list_limit_pattern:
                    if not integrated_CT and not limit_pattern.startswith('('):
                        # If not specified as a raw string, we take the liberty of adding 
                        # the enclosing parenthesis.
                        limit_pattern = '(%s,)' % limit_pattern
                    elif integrated_CT and not limit_pattern.startswith('['):
                        # If not specified as a raw string, we take the liberty of adding 
                        # the enclosing parenthesis.
                        limit_pattern = '[%s,]' % limit_pattern                    
                    # We also take the liberty of escaping the parenthesis
                    # since this is presumably what the user expects.
                    limit_pattern = limit_pattern.replace('(', '\(').replace(')', '\)')
                    limit_pattern = limit_pattern.replace('[', '\[').replace(']', '\]')
                    new_list_limit_pattern.append(limit_pattern)
                limit_pattern_re = re.compile(r'^(%s)$'%(
                    '|'.join(limit_pattern for limit_pattern in new_list_limit_pattern) ))    

        returned_counterterms = []
        if not limit_pattern:
            returned_counterterms.append(random.choice(selected_counterterms))
        else:
            for counterterm in selected_counterterms:
                ct = counterterm if not integrated_CT else counterterm['integrated_counterterm']
                if re.match(limit_pattern_re, str(ct)):
                    returned_counterterms.append(counterterm)

        return returned_counterterms

    def pass_all_cuts(self, PS_point, flavors, 
                          n_jets_allowed_to_be_clustered = None, xb_1 = None, xb_2 = None):
        """ Calls both the flavor-blind and flavor-sensitive cuts."""
        return self.pass_flavor_blind_cuts(PS_point, flavors, 
                    n_jets_allowed_to_be_clustered = n_jets_allowed_to_be_clustered, 
                    xb_1 = xb_1, xb_2 = xb_2 ) and \
                self.pass_flavor_sensitive_cuts(PS_point, flavors, 
                    xb_1 = xb_1, xb_2 = xb_2 )

    def pass_flavor_blind_cuts(self, PS_point, process_pdgs, 
                          n_jets_allowed_to_be_clustered = None, xb_1 = None, xb_2 = None):
        """ Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and 
        will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
        We consider here a two-level cuts system, this first one of which is flavour blind.
        The 'n_jets_allowed_to_be_clustered' is an option that allows to overwrite the
        maximum number of jets that can be clustered and which is by default taken to be:
            self.contribution_definition.n_unresolved_particles 
        This is useful when using this function to apply cuts to the reduced PS of the CTs.
        Finally, the options 'xb_1' and 'xb_2' allow to specify the boost bringing
        the PS point back from the c.o.m. to the lab frame, necessary if some cuts are not
        invariant under a boost along the beam axis."""

        debug_cuts = False
        # This is a temporary function anyway which should eventually be replaced by a full
        # fledged module for handling generation level cuts, which would also make use of fjcore.
        # return True

        from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList, LorentzVector
        
        # If the cuts depend on the boost to the lab frame in case of hadronic collision
        # then the quantity below can be used:
#        boost_vector_to_lab_frame = None
#        if (xb_1 is not None) and (xb_2 is not None) and (xb_1, xb_2)!=(1.,1.):
#            boost_vector_to_lab_frame = PS_point.get_boost_vector_to_lab_frame(xb_1, xb_2)
        
        # These cuts are not allowed to resolve flavour, but only whether a particle is a jet or not
        def is_a_jet(pdg):
            return abs(pdg) in range(1,self.run_card['maxjetflavor']+1)+[21]
        
        def is_a_lepton(pdg):
            return abs(pdg) in [11,13,15]
        
        def is_a_neutrino(pdg):
            return abs(pdg) in [12,14,16]
        
        def is_a_photon(pdg):
            return pdg==22
        
        if debug_cuts: logger.debug( "Processing flavor-blind cuts for process %s and PS point:\n%s"%(
            str(process_pdgs), LorentzVectorList(PS_point).__str__(n_initial=self.phase_space_generator.n_initial) ))

        if n_jets_allowed_to_be_clustered is None:
            n_jets_allowed_to_be_clustered = self.contribution_definition.n_unresolved_particles

        ###################################################################################
        # JET CLUSTERING AND CUTS
        ###################################################################################
        
        ptj_cut = self.run_card['ptj']
        drjj_cut = self.run_card['drjj']
        etaj_cut = self.run_card['etaj']
        
        if ptj_cut <= 0. and    \
           drjj_cut <= 0. and   \
           etaj_cut <= 0. :
            return True
        else:
            # If fastjet is needed but not found, make sure to stop
            if (not PYJET_AVAILABLE) and n_jets_allowed_to_be_clustered>0:
                raise MadEvent7Error("Fast-jet python bindings are necessary for integrating"+
                             " real-emission type of contributions. Please install pyjet.")

        # The PS point in input is sometimes provided as a dictionary or a flat list, but
        # we need it as a flat list here, so we force the conversion
        if isinstance(PS_point, LorentzVectorDict):
            PS_point = PS_point.to_list()
        elif not isinstance(PS_point, LorentzVectorList):
            PS_point = LorentzVectorList(LorentzVector(v) for v in PS_point)    

        if PYJET_AVAILABLE and drjj_cut > 0.:

            # First identify all partonic jets
            jets_list = []
            for i, p in enumerate(PS_point[self.n_initial:]):
                if is_a_jet(process_pdgs[1][i]):
                    jets_list.append(tuple(list(p)+[i+self.n_initial+1,]))
            # Count partonic jets
            starting_n_jets = len(jets_list)

            # Cluster them with fastjet
            this_event = np.array(jets_list,dtype=np.dtype(
                    [('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8'), ('id', 'i8')]) )
            sequence = pyjet.cluster(this_event, R=drjj_cut, p=-1, ep=True)
            jets = sequence.inclusive_jets()

#            misc.sprint('Process flavors: %s'%str(process_pdgs))
#            misc.sprint('Starting from configuration:\n%s'%str(PS_point))
#            misc.sprint('The following %d jets were found:'%len(jets))
#            misc.sprint("{0: <5} {1: >10} {2: >10} {3: >10} {4: >10} {5: >10} {6: >10} {7: >10}".format(
#                                      "jet#", "E", "Px", "Py", "P.z", "mass", "pT", "#constit."))
#            for i, jet in enumerate(jets):
#                misc.sprint("{0: <5} {1: >10} {2: >10} {3: >10} {4: >10} {5: >10} {6: >10} {7: >10}".format(
#                    i + 1, jet.e, jet.px, jet.py, jet.pz, jet.mass, jet.pt, 
#                    str(tuple(constituent.id for constituent in jet))
#                    ))
#                misc.sprint('\n'+'\n'.join('  | %s'%constituent for
#                                           constituent in jet.constituents_array(ep=True)))

            # Remove jets whose pT is below the user defined cut:
            if ptj_cut > 0.:
                jets = [jet for jet in jets if jet.pt >= ptj_cut]
            
#            # A useless (in production) consistency check for deltaR:
#            for i, j1 in enumerate(jets):
#                for j, j2 in enumerate(jets):
#                    if i==j: continue
#                    pj1 = LorentzVector([j1.e, j1.px, j1.py, j1.pz])
#                    pj2 = LorentzVector([j2.e, j2.px, j2.py, j2.pz])
#                    if pj1.deltaR(pj2) < drjj_cut:
#                        raise MadGraph5Error("Inconsistency with fastjet in pass_cuts : "+
#                                                  " %.5f < %.5f"%(pj1.deltaR(pj2), drjj_cut))
#                    misc.sprint(pj1.deltaR(pj2), drjj_cut)                
            
            # Make sure that the number of clustered jets is at least larger or equal to the
            # starting list of jets minus the number of particles that are allowed to go
            # unresolved in this contribution.
            if debug_cuts: logger.debug("Number of identified jets: %d (min %d)"%
                           ( len(jets), (starting_n_jets-n_jets_allowed_to_be_clustered) ))
            if len(jets) < (starting_n_jets-n_jets_allowed_to_be_clustered):

                return False
            
            all_jets = LorentzVectorList([LorentzVector(
                               [a_jet.e, a_jet.px, a_jet.py, a_jet.pz]) for a_jet in jets])
            
        else:
            all_jets = LorentzVectorList([p for i, p in enumerate(PS_point[self.n_initial:])
                                                          if is_a_jet(process_pdgs[1][i])])
            if ptj_cut > 0.:
                # Apply the Ptj cut first
                for i, p in enumerate(all_jets):
                    if debug_cuts: logger.debug('pj_%i.pt()=%.5e'%((i+1),p.pt()))
                    if p.pt() < ptj_cut:
                        return False
    
            if drjj_cut > 0.:
                # And then the drjj cut
                for i, p1 in enumerate(all_jets):
                    for j, p2 in enumerate(all_jets):
                        if j <= i:
                            continue
                        if debug_cuts: logger.debug('deltaR(pj_%i,pj_%i)=%.5e'%(
                             i+1, j+1, p1.deltaR(p2)))
                        if p1.deltaR(p2) < drjj_cut:
                            return False
    
        # Now handle all other cuts
        if etaj_cut > 0.:
            for i, p_jet in enumerate(all_jets):
                if debug_cuts: logger.debug('eta(pj_%i)=%.5e'%(i+1,p_jet.pseudoRap()))
                if abs(p_jet.pseudoRap()) > etaj_cut:
                    return False

        ###################################################################################
        # LEPTON AND PHOTON CUTS
        ###################################################################################
        
        for i, p in enumerate(PS_point[self.n_initial:]):
            # photons
            if is_a_photon(process_pdgs[1][i]):
                if debug_cuts: logger.debug('pta_%i.pt()=%.5e'%((i+1),p.pt()))
                if self.run_card['pta'] > 0.0 and p.pt() < self.run_card['pta']:
                    return False
                if debug_cuts: logger.debug('eta(pa_%i)=%.5e'%(i+1,p.pseudoRap()))
                if self.run_card['etaa'] > 0.0 and abs(p.pseudoRap()) > self.run_card['etaa']:
                    return False
                for j, p2 in enumerate(PS_point[self.n_initial:]):
                    if j > i and is_a_photon(process_pdgs[1][j]):
                        if debug_cuts: logger.debug('deltaR(pa_%i,pa_%i)=%.5e'%(i+1, j+1, p.deltaR(p2)))
                        if self.run_card['draa'] > 0.0 and p.deltaR(p2) < self.run_card['draa']:
                            return False
                    if is_a_lepton(process_pdgs[1][j]):
                        if debug_cuts: logger.debug('deltaR(pa_%i,pl_%i)=%.5e'%(i+1, j+1, p.deltaR(p2)))
                        if self.run_card['dral'] > 0.0 and p.deltaR(p2) < self.run_card['dral']:
                            return False
                for j, p_jet in enumerate(all_jets):
                    if debug_cuts: logger.debug('deltaR(pa_%i,pj_%i)=%.5e'%(i+1, j+1, p.deltaR(p_jet)))
                    if self.run_card['draj'] > 0.0 and p.deltaR(p_jet) < self.run_card['draj']:
                        return False
            
            # leptons
            if is_a_lepton(process_pdgs[1][i]):
                if debug_cuts: logger.debug('ptl_%i.pt()=%.5e'%((i+1),p.pt()))
                if self.run_card['ptl'] > 0.0 and p.pt() < self.run_card['ptl']:
                    return False
                if debug_cuts: logger.debug('eta(pl_%i)=%.5e'%(i+1,p.pseudoRap()))
                if self.run_card['etal'] > 0.0 and abs(p.pseudoRap()) > self.run_card['etal']:
                    return False
                for j, p2 in enumerate(PS_point[self.n_initial:]):
                    if j <= i:
                        continue
                    if is_a_lepton(process_pdgs[1][j]):
                        if debug_cuts: logger.debug('deltaR(pl_%i,pl_%i)=%.5e'%(i+1, j+1, p.deltaR(p2)))
                        if self.run_card['drll'] > 0.0 and p.deltaR(p2) < self.run_card['drll']:
                            return False
                for j, p_jet in enumerate(all_jets):
                    if debug_cuts: logger.debug('deltaR(pl_%i,pj_%i)=%.5e'%(i+1, j+1, p.deltaR(p_jet)))
                    if self.run_card['drjl'] > 0.0 and p.deltaR(p_jet) < self.run_card['drjl']:
                        return False  
            

        # All cuts pass, therefore return True
        return True

    def pass_flavor_sensitive_cuts(self, PS_point, flavors, xb_1 = None, xb_2 = None):
        """ Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and 
        will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
        We consider here a two-level cuts system, this second one of which is flavour sensitive.
        The xb_1 and xb_2 Bjorken x's can be specified, which can be useful in cases where
        the cut depends on quantities not invariant under boosts along the beam axis. """

        # None implemented by default
        return True

        # If the cuts depend on the boost to the lab frame in case of hadronic collision
        # then the quantity below can be used:
#        boost_vector_to_lab_frame = None
#        if (xb_1 is not None) and (xb_2 is not None) and (xb_1, xb_2)!=(1.,1.):
#            boost_vector_to_lab_frame = PS_point.get_boost_vector_to_lab_frame(xb_1, xb_2)

        debug_cuts = False
        from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList, LorentzVector

        # The PS point in input is sometimes provided as a dictionary or a flat list, but
        # we need it as a flat list here, so we force the conversion
        if isinstance(PS_point, LorentzVectorDict):
            PS_point = PS_point.to_list()
        elif not isinstance(PS_point, LorentzVectorList):
            PS_point = LorentzVectorList(LorentzVector(v) for v in PS_point)   

        if debug_cuts: logger.debug( "Processing flavor-sensitive cuts for flavors %s and PS point:\n%s"%(
            str(flavors), LorentzVectorList(PS_point).__str__(n_initial=self.phase_space_generator.n_initial) ))
        
        # All cuts pass, therefore return True
        return True

    @staticmethod
    def Lambda(s,sqrMA,sqrMB):
        """ Kahlen function."""

        return s**2 + sqrMA**2 + sqrMB**2 - 2.*s*sqrMA - 2.*sqrMB*sqrMA - 2.*s*sqrMB

    def get_scales(self, PS_point):
        """ Returns mu_r, mu_f1, mu_f2 for that PS point."""
        
        if not self.run_card['fixed_ren_scale'] or not self.run_card['fixed_fac_scale']:
            raise InvalidCmd("MadEvent7 only supports fixed mu_f and mu_r scales for now.")
        
        return self.run_card['scale'], self.run_card['dsqrt_q2fact1'], self.run_card['dsqrt_q2fact2'] 

    def get_pdfQ2(self, pdf, pdg, x, scale2):
        """ Call the PDF and return the corresponding density."""

        if pdf is None:
            return 1.
       
        if pdg not in [21,22] and abs(pdg) not in range(1,7):
            return 1.
                
        if (pdf, pdg, x, scale2) in self.PDF_cache:
            return self.PDF_cache[(pdf, pdg, x, scale2)]
        
        # Call to lhapdf API
        f = pdf.xfxQ2(pdg, x, scale2)/x
        
        # Update the PDF cache
        self.PDF_cache[(pdf, pdg,x,scale2)] = f
        self.PDF_cache_entries.append((pdf, pdg,x,scale2)) 
        if len(self.PDF_cache_entries) > self.PDF_cache_max_size:
            del self.PDF_cache[self.PDF_cache_entries.pop(0)]

        return f

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Main function of the integrand, returning the weight to be passed to the integrator."""

        # A unique float must be returned
        wgt = 1.0
        # And the conversion from GeV^-2 to picobarns
        wgt *= 0.389379304e9
        
        if __debug__: logger.debug("="*80)       
        if __debug__: logger.debug('Starting a new evaluation of the integrand from contribution:\n%s',
                                                    self.contribution_definition.nice_string())
        
        # Random variables sent
        random_variables    = list(continuous_inputs)
        if __debug__: logger.debug('Random variables received: %s',str(random_variables))        
    
        # Now assign the variables pertaining to PS generations
        PS_random_variables = [random_variables[self.dim_name_to_position[name]] for name in 
                                                self.phase_space_generator.dim_ordered_names]
        
        PS_point, PS_weight, x1s, x2s = self.phase_space_generator.get_PS_point(PS_random_variables)
        # Unpack the initial momenta rescalings (if present) so as to access both Bjorken
        # rescalings xb_<i> and the ISR factorization convolution rescalings chsi<i>.
        xb_1, chsi1 = x1s
        xb_2, chsi2 = x2s

        if PS_point is None:
            if __debug__:
                if xb_1 > 1. or xb_2 > 1.:
                    logger.debug('Unphysical configuration: x1, x2 = %.5e, %.5e'%(xb_1, xb_2))
                    logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
                elif (chsi1 is not None) and xb_1 > chsi1:
                    logger.debug('Configuration above chsi1 rescaling: x1 > chsi1 : %.5e > %.5e'%(xb_1, chsi1))
                elif (chsi2 is not None) and xb_2 > chsi2:
                    logger.debug('Configuration above chsi1 rescaling: x2 > chsi2 : %.5e > %.5e'%(xb_2, chsi2))
                else:
                    logger.debug('Phase-space generation failed.')
                logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
            return 0.0
                
        if __debug__: logger.debug("Considering the following PS point:\n%s"%(PS_point.__str__(
                                            n_initial=self.phase_space_generator.n_initial) ))

        # Account for PS weight
        wgt *= PS_weight
        if __debug__: logger.debug("PS_weight: %.5e"%PS_weight)

        # The E_cm entering the flux factor is computed *without* including the chsi<i> rescalings
        # This is necessary and must be kept so.
        E_cm = math.sqrt(xb_1*xb_2)*self.collider_energy
        
        # Include the flux factor
        flux = 1.
        if self.n_initial == 2:
            flux = 1. / (2.*math.sqrt(self.Lambda(E_cm**2, self.masses[0][0]**2, self.masses[0][1]**2)))
        elif self.n_initial == 1:
            flux = 1. / (2.*E_cm)
        flux /= math.pow(2.*math.pi, 3*self.n_final - 4)
        wgt *= flux
        if __debug__: logger.debug("Flux factor: %.5e"%flux)

        # Recover scales to be used
        mu_r, mu_f1, mu_f2 = self.get_scales(PS_point)

        # Apply the recomputed alpha_s from the PDFs or a dedicated model_reader.Alphas_Runner
        alpha_s = self.alpha_s_runner(mu_r)
        
        # Notice here that we do *not* update all the couplings/parameters dependent on this new value of mu_r / alpha_s
        # We reset only mu_r and alpha_s since this is the only thing our integrands directly depend on so far (along with
        # the masses and widths which are of course alpha_s independent.)
        # Indeed the matrix element take their couplings and others directly from the fortran exported version (i.e via f2py).
        # Also it can be quite time consuming to update the whole set of dependent parameters that are python expression, so
        # it is best to avoid it if possible.
        model_param_dict = self.model.get('parameter_dict')
        model_param_dict['aS'] = alpha_s
        if 'MU_R' in model_param_dict:
            model_param_dict['MU_R'] = mu_r

        # Now loop over processes
        total_wgt = 0.
        for process_key, (process, mapped_processes) in self.processes_map.items():
            if __debug__: logger.debug('Now considering the process group from %s.'%process.nice_string())
            
            process_pdgs = process.get_cached_initial_final_pdgs()
            
            all_processes = [process,]+mapped_processes
            all_flavor_configurations = []
            
            # The process mirroring is accounted for at the very end only            
            for proc in all_processes:
                initial_final_pdgs = proc.get_cached_initial_final_pdgs()
                all_flavor_configurations.append(initial_final_pdgs)

            # Compute the short distance cross-section. The 'events' returned is an instance
            # of EventList, specifying all the contributing kinematic configurations, 
            # and for each all the weights of the potentially contributing flavors.
            events = self.sigma(
                PS_point.to_dict(), process_key, process, all_flavor_configurations, 
                                        wgt, mu_r, mu_f1, mu_f2, chsi1, chsi2, xb_1, xb_2)
            # misc.sprint(events)
            # Apply flavor blind cuts
            if not events.filter_with_flavor_blind_cuts(self.pass_flavor_blind_cuts, process_pdgs,
                n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles,
                xb_1 = xb_1, xb_2 = xb_2):
                if __debug__: logger.debug('All events failed the flavour_blind generation-level cuts.')
                if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)            
                return 0.0
            
            # Select particular terms of the EpsilonExpansion terms stored as weigts
            events.select_epsilon_expansion_term('finite')
            
            # Apply PDF convolution
            if self.run_card['lpp1']==self.run_card['lpp2']==1:
                events.apply_PDF_convolution( self.get_pdfQ2, 
                                 (self.pdf, self.pdf), (xb_1, xb_2), (mu_f1**2, mu_f2**2) )
            # Make sure Bjorken-x rescalings don't matter anymore
            for event in events:
                event.set_Bjorken_rescalings(None, None)
                
            # Now that the PDF convolution has been performed and that the Bjorken rescaling
            # attributed of the events no longer matter, it may be possible here to combine
            # several events having the same kinematic support. This is a further
            # optimization that we can easily implement when proven needed.
            #   all_events_generated.combine_events_with_identical_kinematics()
            # (PS: the function above is a place-holder and is not implemented yet.
            
            # Apply flavor sensitive cuts
            if not events.filter_flavor_configurations(self.pass_flavor_sensitive_cuts,
                                                                 xb_1 = xb_1, xb_2 = xb_2):
                if __debug__: logger.debug('Events failed the flavour-sensitive generation-level cuts.')
                if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)            
                return 0.0

            # Aggregate the total weight to be returned to the integrator
            total_wgt += events.get_total_weight()
            
            # Select a unique flavor configuration for each event
            events.select_a_flavor_configuration()
            #misc.sprint(events)
            if __debug__: logger.debug('Short-distance events for this subprocess:\n%s\n'%str(events))
            
            # Finally apply observables
            if self.apply_observables:
                events.apply_observables(self.observable_list, xb_1, xb_2)

        # Now finally return the total weight for this contribution
        if __debug__: logger.debug(misc.bcolors.GREEN + "Final weight returned: %.5e"%total_wgt + misc.bcolors.ENDC)
        if __debug__: logger.debug("="*80)
        
        return total_wgt
    
    def sigma(self, PS_point, process_key, process, all_flavor_configurations, 
                                            base_weight, mu_r, mu_f1, mu_f2, *args, **opts):
        """ 
        This is the core function of the integrand where the short-distance objects like the matrix elements,
        the counterterms, the mappings, etc.. will be evaluated.
        The process has a lot of extra meta information that sigma may consider, but flavors provides what flavour 
        assignment should be used for some of the flavour dependent steps (like calling the cuts for instance).
        Notice that this function may be called several times with the exact same arguments but different 
        processes / flavour, so some caching may be in order. 
        The output of sigma, sigma_wgt, should be the weight summed over the various contributions building sigma
        (ME's and counterterms).
        Finally, this function is also responsible for calling the observable, using process_wgt*sigma_wgt.
        It can potentially call it several times, for the various counterterms.
        """

        sigma_wgt = base_weight

        alpha_s = self.model.get('parameter_dict')['aS']

        # Access to the matrix element. ME_result is an instance of a subclassed dictionary which includes
        # all independent results available as of now for the particular arguments specified.
        # We specify pdgs to None her to avoid the need of any permutation since we follow the order of
        # the defining process here which is the one that was exported.
        # For the reduced matrix elements however, this cannot be done.
        ME_evaluation, all_results = self.all_MEAccessors(
                        process, PS_point, alpha_s, mu_r, pdgs=all_flavor_configurations[0])
        
        ## Some code to test the color correlated MEs
        ##color_correlation_to_consider = ( ((3, -2, 3), (-2, -2, -1)), ((3, -1, 3), (-1, -1, -2)) )
        ##color_correlation_to_consider = ( -1, ((3, -1, 3), (-1, -1, -2)) )
        ##color_correlation_to_consider = ( -1, ((3, -1, 3), ) )        
        ##ME_evaluation, all_results = self.all_MEAccessors(
        ##    process, PS_point, alpha_s, mu_r, pdgs=flavors, 
        ##    color_correlation=( color_correlation_to_consider, ),
        ##    return_all_res = True)
        ##misc.sprint(str(ME_evaluation))
        ##misc.sprint(str(all_results))
        ##ME_evaluation, all_results = self.all_MEAccessors(
        ##    process, PS_point, alpha_s, mu_r, pdgs=flavors, 
        ##    color_correlation=( ( ((4, -1, 4), ), ((3, -1, 3), ) ), ),
        ##    return_all_res = True)
        ##misc.sprint(str(ME_evaluation))
        ##misc.sprint(str(all_results))
        
        sigma_wgt *= ME_evaluation['finite']

        # Return the lone LO event
        return ME7EventList([
            ME7Event( PS_point, {fc : sigma_wgt for fc in all_flavor_configurations},
                is_mirrored = process.get('has_mirror_process'),
                host_contribution_definition = self.contribution_definition,
                counterterm_structure  = None
            )
        ])

    def convolve_event_with_beam_factorization_currents(self, input_event_to_convolve, 
            beam_factorization_currents, PS_point, process, mu_r, mu_f1, mu_f2, chsi1, chsi2):
        """ Calls the currents specified in the 'beam_factorization_currents' argument, of the form:
            [{ 'beam_one': <BeamCurrent>,
               'beam_two': <BeamCurrent> }]
        and convolve their evaluations with the ME7Event specified in input so to construct
        a new convolved ME7Event which is then returned. This function is placed in the mother
        class because it is useful both in ME7_integrand_R and ME7_integrand_V."""

        # Adjust the format of the input:
        beam_factorization_currents = [
            (  bc['beam_one'],  bc['beam_two'] ) for bc in beam_factorization_currents ]

        convolved_event = None
        for bc1, bc2 in beam_factorization_currents:
            event_to_convolve = input_event_to_convolve.get_copy()
            if bc1 is not None:
                assert(isinstance(bc1, subtraction.BeamCurrent) and bc1['distribution_type']=='bulk')
                current_evaluation, all_current_results = self.all_MEAccessors(
                    bc1, lower_PS_point=PS_point, reduced_process=process, chsi=chsi1, mu_r=mu_r, mu_f=mu_f1)
                assert(current_evaluation['spin_correlations']==[None,])
                assert(current_evaluation['color_correlations']==[None,])
                event_to_convolve.convolve_flavors( 
                    base_objects.EpsilonExpansion(current_evaluation['values'][(0,0)]), leg_index=0 )
            if bc2 is not None:
                assert(isinstance(bc2, subtraction.BeamCurrent) and bc2['distribution_type']=='bulk')
                current_evaluation, all_current_results = self.all_MEAccessors(
                    bc2, lower_PS_point=PS_point, reduced_process=process, chsi=chsi2, mu_r=mu_r, mu_f=mu_f2)
                assert(current_evaluation['spin_correlations']==[None,])
                assert(current_evaluation['color_correlations']==[None,])
                event_to_convolve.convolve_flavors( 
                       base_objects.EpsilonExpansion(current_evaluation['values'][(0,0)]), leg_index=1 )
            if convolved_event is None:
                convolved_event = event_to_convolve
            else:
                # Combine the weights and flavor of the event obtained from this convolution
                # with the ones obtained from the previous convolutions in this loop.
                convolved_event += event_to_convolve
        
        # Determine Bjorken scalings.
        assert( all(bc[0] is None for bc in beam_factorization_currents) or
                all(bc[0] is not None for bc in beam_factorization_currents) )
        Bjorken_scaling_beam_one = 1. if beam_factorization_currents[0][0] is None else 1./chsi1
        assert( all(bc[1] is None for bc in beam_factorization_currents) or
                all(bc[1] is not None for bc in beam_factorization_currents) )
        Bjorken_scaling_beam_two = 1. if beam_factorization_currents[0][1] is None else 1./chsi2
        # Assign the Bjorken scaling found
        convolved_event.set_Bjorken_rescalings(Bjorken_scaling_beam_one, Bjorken_scaling_beam_two)
        return convolved_event

##########################################################################################
# Daughter classes for the various type of ME7Integrands which require a redefinition 
# of sigma() and possibly other functions.
##########################################################################################
class ME7Integrand_B(ME7Integrand):
    """ME7Integrand for the computation of a Born type of contribution."""
    
    def sigma(self, PS_point, process_key, process, all_flavor_configurations, 
                                            base_weight, mu_r, mu_f1, mu_f2, *args, **opts):

        return super(ME7Integrand_B, self).sigma(
            PS_point, process_key, process, all_flavor_configurations, base_weight,
            mu_r, mu_f1, mu_f2, *args, **opts
        )

class ME7Integrand_LIB(ME7Integrand):
    """ ME7Integrand for the computation of a Loop-Induced Born type of contribution."""
    def sigma(self, PS_point, process_key, process, all_flavor_configurations, 
                                            base_weight, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_LIB, self).sigma(
                  PS_point, process_key, process, all_flavor_configurations, 
                                            base_weight, mu_r, mu_f1, mu_f2, *args, **opts)
        
class ME7Integrand_V(ME7Integrand):
    """ ME7Integrand for the computation of a one-loop virtual type of contribution."""
    
    def __init__(self, *args, **opts):
        """Initialize a virtual type of integrand, adding additional relevant attributes."""
        
        try:  
            self.integrated_counterterms = opts.pop('integrated_counterterms')
        except KeyError:
            raise MadEvent7Error(
                "Constructor of class ME7Integrand_V requires the option "
                "'integrated_counterterms' to be specified."
            )

        super(ME7Integrand_V, self).__init__(*args, **opts)
        self.initialization_inputs['options']['integrated_counterterms'] = \
                                                              self.integrated_counterterms

    def has_integrated_counterterms(self):
        """ A virtual type of integrand can contain integrated counterterms."""
        return True

    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this contribution."""
        res = []
        if self.integrated_counterterms:
            res.append('%-30s:   %d'%('Nb. of integrated counterterms', 
                                      len(sum(self.integrated_counterterms.values(),[]))))
        return res
    
    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """Return a nicely formatted process line for the function nice_string of this
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
                    long_res.append( '   | %s'%CT.__str__(
                                        print_n=True, print_pdg=False, print_state=False )  )
                elif format==3:
                    long_res.append( '   | %s'%CT.__str__(
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

    def evaluate_integrated_counterterm(self, integrated_CT_characteristics, PS_point, 
        base_weight, mu_f1, mu_f2, chsi1, chsi2, is_process_mirrored, input_mapping, 
        all_virtual_ME_flavor_configurations, hel_config=None, compute_poles=True):
        """ Evaluates the specified integrated counterterm, provided along with its other
        characteristics, like for example the list of flavors assignments that the resolved
        process it corresponds to can take. This function returns an ME7Event specifying the
        counterevent for the specified input_mapping considered (i.e. an integrated CT like
        g > d d~ must be "attached" to all final-state gluons of the virtual process definition).
        """

        # Access the various characteristics of the integrated counterterm passed to this
        # function.
        counterterm = integrated_CT_characteristics['integrated_counterterm']
        # The resolved flavor configurations dictionary (mapping resolved to reduced) is
        # typically no longer useful here.
        resolved_flavors = integrated_CT_characteristics['resolved_flavors_combinations']
        reduced_flavors = integrated_CT_characteristics['reduced_flavors_combinations']
        symmetry_factor = integrated_CT_characteristics['symmetry_factor']


        # And the multiplicity prefactor coming from the several *resolved* flavor assignment
        # that this counterterm can lead to. Typically an integrated counterterm for g > q qbar
        # splitting will have the same reduced flavors, but with the gluon coming from
        # n_f different massless flavors. So that its multiplcitiy factor is n_f.        
        # Also, if you think of the resolved flavors e+ e- > c c~ u u~, the two counterterms
        # C(3,4) and C(5,6) will lead to the reduced flavors g u u~ and c c~ g respectively.
        # These two are however mapped in the same virtual group. In this case, the flavor 
        # configuration 'g u u~' in the ME7Event would therefore not receive any contribution
        # from the integrated CT C(5,6), but the configuration 'g c c~' will.
        n_initial = len(all_virtual_ME_flavor_configurations[0][0])
        n_final   = len(all_virtual_ME_flavor_configurations[0][1])
        all_mapped_flavors = {}
        for flavors in all_virtual_ME_flavor_configurations:
            mapped_flavors = ( 
                tuple( flavors[0][input_mapping[i]] for i in range(n_initial) ),
                tuple( flavors[1][input_mapping[i]-n_initial] for i in 
                                                      range(n_initial, n_initial+n_final) )
            )
            if mapped_flavors in reduced_flavors:
                assert ((mapped_flavors not in all_mapped_flavors))
                # Also account for the multiplicity from overall symmetry factors S_t
                all_mapped_flavors[mapped_flavors] = base_objects.EpsilonExpansion({0:
                                   float(symmetry_factor*reduced_flavors[mapped_flavors])})

        # Now map the momenta
        if isinstance(PS_point,dict):
            # Dictionary format LorentzVectorDict starts at 1
            mapped_PS_point = phase_space_generators.LorentzVectorDict(   
                (i+1, PS_point[input_mapping[i]+1]) for i in range(n_initial+n_final) )
        else:
            # List formatLorentzVectorList starts at 0
            mapped_PS_point = phase_space_generators.LorentzVectorDict(   
                (i+1, PS_point[input_mapping[i]]) for i in range(n_initial+n_final) )
        
        # Retrieve some possibly relevant model parameters
        alpha_s = self.model.get('parameter_dict')['aS']
        mu_r = self.model.get('parameter_dict')['MU_R']   
        
        # Now compute the reduced quantities which will be necessary for evaluating the
        # integrated current
        reduced_PS = counterterm.get_reduced_kinematics(mapped_PS_point)
        
        # Make sure no helicity configuration is specified since this is not supported yet.
        assert ((hel_config is None))
        # all_necessary_ME_calls is a list inputs to call the Matrix Element and the weights
        # that we must multiply/convolve them with. We start with empty entries.
        all_necessary_ME_calls = [ 
            {   'spin_correlations'             : [ ],
                'color_correlations'            : [ ],
                'main_weights'                  : [ ],
                'flavor_matrices_beam_one'      : [ ],
                'flavor_matrices_beam_two'      : [ ],
                'Bjorken_rescalings_beam_one'   : [ ],
                'Bjorken_rescalings_beam_two'   : [ ],
            },
        ]
        total_jacobian = 1.
        disconnected_currents_weight = base_objects.EpsilonExpansion({'finite': 1.0})
        
        # First call the non-beam factorization currents
        for integrated_current in counterterm.get_all_currents():
            if isinstance(integrated_current, (subtraction.BeamCurrent,subtraction.IntegratedBeamCurrent)):
                continue
            # /!\ Warnings the flavors of the reduced process as well as the ones of the current
            # are tokens that will apply to all possible flavor configuration in this contribution
            # This should however be irrelevant for the evaluation of the counterterm.
            current_evaluation, all_results = self.all_MEAccessors(
                integrated_current, lower_PS_point=reduced_PS,
                higher_PS_point=None,
                reduced_process = counterterm.process,
                leg_numbers_map = counterterm.momenta_dict,
                hel_config      = hel_config,
                compute_poles   = compute_poles )

            # Now loop over all spin- and color- correlators required for this current
            # and update the necessary calls to the ME
            if not integrated_current['resolve_mother_spin_and_color']:
                # Make sure no spin- or color-correlations were produced by the current
                assert(current_evaluation['spin_correlations']==[None,])
                assert(current_evaluation['color_correlations']==[None,])
                assert(current_evaluation['values'].keys()==[(0,0),])
                disconnected_currents_weight *= \
                         base_objects.EpsilonExpansion(current_evaluation['values'][(0,0)])
            else:
                all_necessary_ME_calls = ME7Integrand_R.update_all_necessary_ME_calls(
                     all_necessary_ME_calls, current_evaluation, weight_type='main_weight')
                
        # Then evaluate the beam factorization currents
        all_necessary_ME_calls = ME7Integrand_R.process_beam_factorization_currents(
            all_necessary_ME_calls, counterterm.get_beam_currents(), self.all_MEAccessors,
            reduced_PS, counterterm.process, chsi1, chsi2, mu_r, mu_f1, mu_f2)

        # Now perform the combination of the list of spin- and color- correlators to be merged
        # for each necessary ME call identified
        all_necessary_ME_calls = ME7Integrand_R.merge_correlators_in_necessary_ME_calls(
                                                                    all_necessary_ME_calls)

        # Finally treat the call to the reduced connected matrix elements
        alpha_s = self.model.get('parameter_dict')['aS']
        mu_r = self.model.get('parameter_dict')['MU_R']
        return ME7Integrand_R.generate_event_for_counterterm(
            ME7Event( mapped_PS_point, 
                {fc : base_weight for fc in all_mapped_flavors},
                is_mirrored                     = is_process_mirrored,
                host_contribution_definition    = self.contribution_definition,
                counterterm_structure           = counterterm),
            disconnected_currents_weight,
            ( counterterm.prefactor / total_jacobian ),
            all_necessary_ME_calls,
            counterterm.process,
            reduced_PS,
            alpha_s, mu_r,
            self.all_MEAccessors
        )
        # We can now create the ME7 Event that will store this integrated CT
        
    def test_IR_poles(self, test_options):
        """ Compare the IR poles residues in dimensional regularization from the virtual
        contribution and from the integrated counterterm. """
        
        if test_options['seed']:
            random.seed(test_options['seed']) 

        # First generate a kinematic point
        # Specifying None forces to use uniformly random generating variables.
        # Make sure to generate a point within the cuts if necessary:
        max_attempts = 10000
        n_attempts   = 0
        while True:
            n_attempts += 1
            if n_attempts > max_attempts:
                break
            virtual_PS_point = None
            # Phase-space generation can fail when beam convolution is active, because of the condition
            # Bjorken_x_i < chsi_i
            while virtual_PS_point is None:
                virtual_PS_point, jac, x1s, x2s = self.phase_space_generator.get_PS_point(None)
            xb_1, chsi1 = x1s
            xb_2, chsi2 = x2s
            if test_options['apply_higher_multiplicity_cuts'] or test_options['apply_lower_multiplicity_cuts']:
                if not self.pass_flavor_blind_cuts(
                    virtual_PS_point,
                    self.processes_map.values()[0][0].get_cached_initial_final_pdgs(),
                    n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles,
                    xb_1 = xb_1, xb_2 = xb_2):
                    continue
            break
        if n_attempts > max_attempts:
            raise MadEvent7Error(
                "Could not generate a random kinematic configuration that passes " +
                "the flavour blind cuts in less than %d attempts." % max_attempts )
        n_attempts = 0

        # Now keep track of the results from each process and poles checked
        all_evaluations = {}
        for process_key, (defining_process, mapped_processes) in self.processes_map.items():
            misc.sprint('\nConsidering %s'%
                              defining_process.nice_string().replace('Process','process'))
            # Make sure that the selected process satisfies the selected process
            if not self.is_part_of_process_selection(
                [defining_process,]+mapped_processes, selection = test_options['process'] ):
                continue

            all_processes = [defining_process,]+mapped_processes
            all_flavor_configurations = []
            # The process mirroring is accounted for at the very end only            
            for proc in all_processes:
                initial_final_pdgs = proc.get_cached_initial_final_pdgs()
                all_flavor_configurations.append(initial_final_pdgs)
            
            a_virtual_PS_point = virtual_PS_point.get_copy()
            a_xb_1, a_chsi1 = x1s
            a_xb_2, a_chsi2 = x2s
            while ( ( test_options['apply_higher_multiplicity_cuts'] or 
                      test_options['apply_higher_multiplicity_cuts']    ) and
                   not self.pass_flavor_sensitive_cuts(
                       a_virtual_PS_point,
                       defining_process.get_cached_initial_final_pdgs(),
                       xb_1 = a_xb_1, xb_2 = a_xb_2 ) ):
                n_attempts += 1
                if n_attempts > max_attempts:
                    break
                a_virtual_PS_point = None
                # Phase-space generation can fail when beam convolution is active, because of the condition
                # Bjorken_x_i < chsi_i
                while a_virtual_PS_point is None:
                    a_virtual_PS_point, a_jac, a_x1s, a_x2s = self.phase_space_generator.get_PS_point(None)
                a_xb_1, a_chsi1 = a_x1s
                a_xb_2, a_chsi2 = a_x2s
            if n_attempts > max_attempts:
                raise MadEvent7Error(
                    "Could not generate a random kinematic configuration that passes " +
                    "the flavour blind cuts in less than %d attempts." % max_attempts )
            n_attempts = 0

            process_evaluation = self.test_IR_poles_for_process(
                    test_options, defining_process, process_key, all_flavor_configurations, 
                    a_virtual_PS_point, a_chsi1, a_chsi2, a_xb_1, a_xb_2 )

            all_evaluations[process_key] = process_evaluation

        # Now produce a nice output of the evaluations and assess whether this test passed or not.
        return self.analyze_IR_poles_check(all_evaluations, test_options['acceptance_threshold'])    

    def test_IR_poles_for_process(self, test_options, defining_process, process_key, 
        all_flavor_configurations, a_virtual_PS_point, a_chsi1, a_chsi2, a_xb_1, a_xb_2):
        """ Given a test_options, a process and all flavor configurations mapped to that 
        particular process, a PS-point onto which to evaluate the various contributions, 
        as well as its corresponding convolution variables chsi1/2; perform the IR poles 
        test that consists in making sure all poles in the epsilon Laurent expansion have
        residues of zero, exhibiting the various pieces (virtual ME, integrated counterterms)
        part of the KLN cancellation."""
        
        # Select which integrated counterterms are selected for the test
        if self.has_integrated_counterterms():
            integrated_counterterms_to_consider = [ ct for ct in self.integrated_counterterms[process_key] 
                if ct['integrated_counterterm'].get_number_of_couplings() <= test_options['correction_order'].count('N') ]
            if test_options['counterterms']:
                counterterm_pattern = test_options['counterterms']
                logger.debug("Integrated counterterms before selection")
                for ct in integrated_counterterms_to_consider:
                    logger.debug("    "+str(ct['integrated_counterterm']))
                integrated_counterterms_to_consider = self.find_counterterms_matching_regexp(
                    integrated_counterterms_to_consider, counterterm_pattern, integrated_CT=True )
            logger.debug("Selected integrated counterterms")
            for ct in integrated_counterterms_to_consider:
                logger.debug("    "+str(ct['integrated_counterterm']))

        # We must also include local counterterms as they too can have poles
        if self.has_local_counterterms():
            local_counterterms_to_consider = [ ct for ct in self.counterterms[process_key]
                if ct.get_number_of_couplings() <= test_options['correction_order'].count('N') ]
            if test_options['counterterms']:
                counterterm_pattern = test_options['counterterms']
                logger.debug("Local counterterms before selection")
                for ct in local_counterterms_to_consider:
                    logger.debug("    "+str(ct))
                local_counterterms_to_consider = self.find_counterterms_matching_regexp(
                    local_counterterms_to_consider, counterterm_pattern )
            logger.debug("Selected local counterterms")
            for ct in local_counterterms_to_consider:
                logger.debug("    "+str(ct))

        mu_r, mu_f1, mu_f2 = self.get_scales(a_virtual_PS_point)

        # Specify the counterterms to be considered and backup the original values
        if self.has_local_counterterms():
            local_CT_backup = self.counterterms[process_key]
            self.counterterms[process_key] = local_counterterms_to_consider
        if self.has_integrated_counterterms():
            integrated_CT_backup = self.integrated_counterterms[process_key]
            self.integrated_counterterms[process_key] = integrated_counterterms_to_consider
        try:
            # Now call sigma in order to gather all events
            events = self.sigma(
                a_virtual_PS_point.to_dict(), process_key, defining_process, all_flavor_configurations, 
                1.0, mu_r, mu_f1, mu_f2, a_chsi1, a_chsi2, a_xb_1, a_xb_2, 
                compute_poles            = True,
                apply_flavour_blind_cuts = (test_options['apply_lower_multiplicity_cuts'] or
                                            test_options['apply_higher_multiplicity_cuts'])
            )
        except Exception as e:
            logger.critical("The following exception occurred when generating events for this integrand %s:\n%s"%(
                                                  self.contribution_short_name, str(e)))
            if self.has_local_counterterms():
                self.counterterms[process_key] = local_CT_backup
            if self.has_integrated_counterterms():
                self.integrated_counterterms[process_key] = integrated_CT_backup
            # Forward the exception while not resetting the stack trace.
            raise

#        misc.sprint('Events generated before post-processing.')
#        misc.sprint(events)
        
        # Now post-process the events as done in __call__ of the integrand.
        # Apply flavor blind cuts
        events.filter_with_flavor_blind_cuts(self.pass_flavor_blind_cuts, 
            defining_process.get_cached_initial_final_pdgs(),
            n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles,
            xb_1 = a_xb_1, xb_2 = a_xb_2 )
        # Apply PDF convolution. It is important to convolve with the PDFs, so as to 
        # correctly get the factors 1/chsi_i multiplying the weights in beam 
        # factorization counterterms. But either the PDF density can be set to 1.
        if self.run_card['lpp1']==self.run_card['lpp2']==1:
            pdf = None if test_options['set_PDFs_to_unity'] else self.pdf
            events.apply_PDF_convolution( self.get_pdfQ2, (pdf, pdf), (a_xb_1, a_xb_2), 
                                                                 (mu_f1**2, mu_f2**2) )
        # Make sure Bjorken-x rescalings chsi_i don't matter anymore
        for event in events:
            event.set_Bjorken_rescalings(None, None)
        # Apply flavor sensitive cuts
        events.filter_flavor_configurations(self.pass_flavor_sensitive_cuts,
                                                          xb_1 = a_xb_1, xb_2 = a_xb_2)

#        misc.sprint('Events generated after post-processing:')
#        misc.sprint(events)

        # Now collect the results necessary to test that the poles cancel in the dictionary below
        evaluation = {
             'virtual_ME'           : None,
             'integrated_CTs'       : None,
             'defining_process'     : defining_process,
             'PS_point'             : a_virtual_PS_point }
        events_sum = None
        for event in events:
            # TODO: instead of looking at the total weight, we could consider doing the
            # check for one particular flavor configuration which could be specified
            # as an extra test parameter.
            # An alternative way is to look at the overall summed event which should have
            # zero poles for all its flavour configurations
            if events_sum is None:
                events_sum = event.get_copy()
            else:
                events_sum += event
            event_wgt = event.get_total_weight()
            if event.counterterm_structure is None:
                if evaluation['virtual_ME'] is None:
                    evaluation['virtual_ME'] = event_wgt
                else:
                    evaluation['virtual_ME'] += event_wgt
            else:
                if str(event.counterterm_structure) in evaluation:
                    evaluation[str(event.counterterm_structure)] += event_wgt           
                else:
                    evaluation[str(event.counterterm_structure)] = event_wgt
                if evaluation['integrated_CTs'] is None:
                    evaluation['integrated_CTs'] = event_wgt
                else:
                    evaluation['integrated_CTs'] += event_wgt

        # Small monitoring of the various contributions:
        logger.debug('-'*50)
        for entry, value in evaluation.items():
            if entry not in ['defining_process','PS_point','integrated_CTs','virtual_ME']:
                logger.debug('%-20s : %s'%(entry, value.__str__(format='.16e')))
        logger.debug('-'*50)
        logger.debug('%-20s : %s'%('integrated_CTs', evaluation['integrated_CTs'].__str__(format='.16e')))
        logger.debug('%-20s : %s'%('virtual_ME', evaluation['virtual_ME'].__str__(format='.16e')))        
        logger.debug('-'*50)
        logger.debug('\nAll events summed:\n%s\n'%str(events_sum))
        logger.debug('-'*50)

        relative_diff = evaluation['virtual_ME'].relative_diff(evaluation['integrated_CTs']*-1.)
        # Finite parts are of course expected to differ, so let's not show them
        relative_diff.truncate(min_power = -2, max_power = -1)
        # To be commented out when we will have a full-fledged analysis coded up 
        # in analyze_IR_poles_check()
        misc.sprint('Summary for that PS point:\n%-20s : %s\n%-20s : %s\n%-20s : %s'%(
             'virtual contrib.',
             evaluation['virtual_ME'].__str__(format='.16e'),
             'integrated CTs.',
             evaluation['integrated_CTs'].__str__(format='.16e'),
             'relative diff.',
             relative_diff.__str__(format='.16e')
        ))
        
        return evaluation

    def analyze_IR_poles_check(self, all_evaluations, acceptance_threshold):
        """ Analyze the results of the check_IR_pole_residues command. """
        
        #TODO
#        misc.sprint("----- SUMMARY -----")
#        for key, evaluation in all_evaluations.items():
#            misc.sprint("Result for test: %s | %s"%(str(dict(key[0])['PDGs']),key[1]))
#            virtual_ME_expansion = evaluation['virtual_ME']
#            all_integrated_CT_summed_res = evaluation['integrated_CTs']
#            relative_diff = virtual_ME_expansion.relative_diff(all_integrated_CT_summed_res)
#            misc.sprint('\n%-30s : %s\n%-30s : %s\n%-30s : %s'%(
#                 'virtual contrib.',
#                 virtual_ME_expansion.__str__(format='.16e'),
#                 'integrated CTs.',
#                 all_integrated_CT_summed_res.__str__(format='.16e'),
#                 'relative diff.',
#                 relative_diff.__str__(format='.16e')
#            ))
#            This is to be compared with the acceptance_threshold
#            difference_norm = relative_diff.norm()
        
        return True

    def sigma(self, PS_point, process_key, process, all_flavor_configurations, 
                             base_weight, mu_r, mu_f1, mu_f2, chsi1, chsi2, *args, **opts):
        """ Overloading of the sigma function from ME7Integrand to include necessary 
        additional contributions. """
        
        all_events_generated = ME7EventList()

        alpha_s = self.model.get('parameter_dict')['aS']

        compute_poles = opts.get('compute_poles', False)

        ME_evaluation, all_results = self.all_MEAccessors(
                        process, PS_point, alpha_s, mu_r, pdgs=all_flavor_configurations[0])

        event_weight = base_objects.EpsilonExpansion(ME_evaluation) * base_weight

        # Notice that the Bjorken rescalings for this even will be set during the convolution
        # performed in the next line, and not directly here in the constructor.
        event_to_convolve = ME7Event( PS_point, 
                {fc : event_weight for fc in all_flavor_configurations},
                is_mirrored                     = process.get('has_mirror_process'),
                host_contribution_definition    = self.contribution_definition,
                counterterm_structure           = None
        )
        convolved_event = self.convolve_event_with_beam_factorization_currents(event_to_convolve, 
            process['beam_factorization'], PS_point, process, mu_r, mu_f1, mu_f2, chsi1, chsi2 )
        
        # Now register the convoluted event
        all_events_generated.append(convolved_event)
        
        # Now loop over all integrated counterterms
        for counterterm_characteristics in self.integrated_counterterms[process_key]:

            # And over all the ways in which this current PS point must be remapped to
            # account for all contributions of the integrated CT. (e.g. the integrated
            # splitting g > d d~ must be "attached" to all final state gluons appearing
            # in this virtual ME process definition.)
            for input_mapping in counterterm_characteristics['input_mappings']:
                # At NLO at least, it is OK to save a bit of time by enforcing 'compute_poles=False').
                # This will need to be re-assessed at NNLO for RV contributions.
                CT_event = self.evaluate_integrated_counterterm( 
                    counterterm_characteristics, PS_point, base_weight, mu_f1, mu_f2, chsi1, chsi2,
                    process.get('has_mirror_process'),
                    input_mapping, all_flavor_configurations,
                    hel_config      = None, 
                    compute_poles   = compute_poles)
            
            if CT_event is not None:
                all_events_generated.append( CT_event )
            else:
                logger.warning('The evaluation of integrated counterterms '+
                                                    'should never yield a None ME7 Event.')

        # Notice that it may be possible in some subtraction scheme to combine
        # several events having the same kinematics support. This is a further
        # optimization that we can easily implement when proven needed.
        # all_events_generated.combine_events_with_identical_kinematics()
        # (PS: the function above is a place-holder and is not implemented yet.
        
        return all_events_generated

class ME7Integrand_R(ME7Integrand):
    """ME7Integrand for the computation of a single real-emission type of contribution."""

    divide_by_jacobian = True
    
    def __init__(self, *args, **opts):
        """Initialize a real-emission type of integrand,
        adding additional relevant attributes.
        """
        
        # Initialize the (counter)terms that make up this integrand
        requires = "Constructor of class ME7Integrand_R requires the option '%s'."
        try:
            self.counterterms = opts.pop('counterterms')
        except KeyError:
            raise MadEvent7Error(requires % 'counterterms')
        # Initialize a mapping walker to handle the limits of this integrand
        try:
            self.subtraction_mappings_scheme = opts.pop('subtraction_mappings_scheme')
        except KeyError:
            raise MadEvent7Error(requires % 'subtraction_mappings_scheme')
        try:
            self.walker = walkers.VirtualWalker(self.subtraction_mappings_scheme)
        except KeyError:
            raise MadEvent7Error(
                "Invalid subtraction_mappings_scheme '%s'." %
                self.subtraction_mappings_scheme)
        # Initialize the ME7Integrand
        super(ME7Integrand_R, self).__init__(*args, **opts)
        # Update the initialization inputs
        self.initialization_inputs['options']['counterterms'] = self.counterterms
        self.initialization_inputs['options']['subtraction_mappings_scheme'] = \
                                                          self.subtraction_mappings_scheme

    def has_local_counterterms(self):
        """ A real-emission type of integrand can contain local counterterms."""
        return True

    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this integrand."""
        res = []
        if self.counterterms:
            res.append('%-30s:   %d'%('Number of local counterterms', 
               len([1 for CT in sum(self.counterterms.values(),[]) if CT.is_singular()]) ))
        return res
        
    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        integrand."""
        
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        process_string = defining_process.nice_string(print_weighted=False)
        res = GREEN + '  ' + process_string.replace('Process: ','') + ENDC

        if not self.counterterms:
            return res                                                               

        if format < 2:
            if process_key in self.counterterms:
                res += ' | %d local counterterms'%len([
                    1 for CT in self.counterterms[process_key] if CT.is_singular() ])
            else:
                res += ' | 0 local counterterm'
                
        else:
            long_res = [' | with the following local counterterms:']
            for CT in self.counterterms[process_key]:
                if CT.is_singular():
                    if format==2:
                        long_res.append( '   | %s'%CT.__str__(
                            print_n=True, print_pdg=False, print_state=False )  )
                    elif format==3:
                        long_res.append( '   | %s'%CT.__str__(
                            print_n=True, print_pdg=True, print_state=True )  )
                    elif format>3:
                        long_res.append( '   | %s'%str(CT))
            res += '\n'.join(long_res)

        return res

    @classmethod
    def combine_color_correlators(cls, color_correlators):
        """ This function takes a list of color_correlators in argument, each specified as:
               color_correlator = ( connection_left, connection_right )
            where connection_<x> is given as:
               connection = ( (motherA, emittedA, daughterA), (motherB, emittedB, daughterB), etc...)
            and returns two lists: 
                combined_color_correlators: the list of of new color correlator specifiers 
                                            that arose from the combination of all of them
                multiplier_factors: a list of float factors to multiply the contribution
                                    of each of the combined color correlator returned above.
        """

        # Trivial combination if there is a single one:
        if len(color_correlators)==1:
            return [color_correlators[0],], [1.0,]

        # First weed out the trivial color correlators set to None.
        non_trivial_color_correlators = [ccs for ccs in color_correlators if ccs is not None]
        
        # We don't support terms that factorize a *sum* of color correlator
        if any(len(ccs)!=1 for ccs in non_trivial_color_correlators):
            raise NotImplementedError(
""" The combination of color correlators only the supports the case where each term returned by
    the current factorizes a *single* correlator not a sum of several ones. So re-implement as follows:
    from:
        evaluation = {
              'spin_correlations'  = [ None, ]
              'color_correlations' = [ (correlatorA, correlatorB) ]
              'values'             = [ (0,0) : my_weight ]
        }
    into:
        evaluation = {
              'spin_correlations'  = [ None, ]
              'color_correlations' = [ (correlatorA,), (correlatorB,) ]
              'values'             = [ (0,0) : my_weight,
                                       (0,1) : my_weight ]
        }
""")
        non_trivial_color_correlators = [ccs[0] for ccs in non_trivial_color_correlators]

        # Convert the NLO short-hand convention if present to the general one
        normalized_non_trivial_color_correlators = []
        for cc in non_trivial_color_correlators:
            normalized_non_trivial_color_correlators.append(
                ( cc[0] if not isinstance(cc[0],int) else ((cc[0],-1,cc[0]),),
                  cc[1] if not isinstance(cc[1],int) else ((cc[1],-1,cc[1]),)  )
            )
        non_trivial_color_correlators = normalized_non_trivial_color_correlators
        
        # Place some limitation of the combinations we understand and can support now
        # Only combinations of two color correlators for now
        if len(non_trivial_color_correlators)>2:
            raise NotImplementedError("Only combinations of at most TWO color correlators are supported for now.")
        # Only NLO-type of correlators for now
        if any( len(cc[0])>1 for cc in non_trivial_color_correlators):
            raise NotImplementedError("Only combinations of at most NLO-type of color correlators are supported for now.")
        
        # Treat the easy case first where there is only one or zero color correlator
        if len(non_trivial_color_correlators)==0:
            return [ None,], [1.0,]
        elif len(non_trivial_color_correlators)==1:
            return [ (non_trivial_color_correlators[0],) ], [1.0,]
        
        # Now treat the non-trivial case of the combination of two NLO-type of color correlators
        # For details of this implementation, I refer the reader to the corresponding documentation
        # produced on this topic
        assert(len(non_trivial_color_correlators)==2)
        ccA = non_trivial_color_correlators[0]
        ccB = non_trivial_color_correlators[1]
        # Identify the terms of the compound soft operator S^{i1,j1} \otimes S^{i2,j2} 
        # Each correlator ccA/B will be given in the form of:
        #    (((i, -1, i),), ((j, -1, j),))
        i1_, j1_ = ccA[0][0][0], ccA[1][0][0]
        i2_, j2_ = ccB[0][0][0], ccB[1][0][0]
        
        # Below is my original way of combining the color structure of the two single soft
        def VH_two_iterated_soft_combination((i1,j1),(i2,j2)):
            # now construct all four ways in which the two gluons could have been connected
            # between these four lines.
            if i1==i2:
                connections_A = [
                    ( (i1,-1,i1),(i2,-2,i2) ),
                    ( (i2,-2,i2),(i1,-1,i1) )
                ]
            else:
                connections_A = [
                    ( (i1,-1,i1),(i2,-2,i2) ) if i1>i2 else
                    ( (i2,-2,i2),(i1,-1,i1) )
                ]
            if j1==j2:
                connections_B = [
                    ( (j1,-1,j1),(j2,-2,j2) ),
                    ( (j2,-2,j2),(j1,-1,j1) )
                ]
            else:
                connections_B = [
                    ( (j1,-1,j1),(j2,-2,j2) ) if j1>j2 else
                    ( (j2,-2,j2),(j1,-1,j1) )
                ]
            combined_color_correlators = []
            for connection_A in connections_A:
                for connection_B in connections_B:
                    combined_color_correlators.append(
                        ( connection_A, connection_B )
                    )
            
            multiplier = 1.0/float(len(combined_color_correlators))
            
             # Implementing them as separate calls to the ME, as done in the commented line below,
            # return [ tuple(combined_color_correlators), ], [ multiplier, ]
            # would only yield the same result once the accessor correctly sums over the
            # different color-correlated MEs, so avoid for now.
            
            # Now return the combined color correlators identified
            return [ (ccc,) for ccc in combined_color_correlators ],\
                   [ multiplier, ]*len(combined_color_correlators)
       
        def double_correlator((i,j),(k,l)):
            """ Returns the double correlator of Catani-Grazzini (Eq.113 of hep-ph/9908523v1)
                <M| ( T^-1_i \dot T^-1_j ) * ( T^-1_k \dot T^-1_l ) | M > 
                
            converted into MadGraph's conventions:
            
              ( (a,-1,a),(b,-2,b) ) , ( (c,-1,c),(d,-2,d) ) --> T^-1_a T^-2_b T^-1_c T^-2_d
            """

            # It is important to never commute two color operators acting on the same index, so we must chose carefully which
            # index to pick to carry the gluon index '-2' of the first connection. This can be either 'k' or 'l'.
            if j!=k and j!=l:
                # If all indices are different, we can pick either k or l, it is irrelevant
                index1, index2, index3, index4 = i, k, j, l
            elif j==k and j!=l:
                # If j is equal to k, we must pick l
                index1, index2, index3, index4 = i, l, j, k
            elif j==l and j!=k:
                # If j is equal to l, we must pick k
                index1, index2, index3, index4 = i, k, j, l
            elif j==l and j==k:
                # If j is equal to both l and k, then agin it doesn't matter and we can pick k
                index1, index2, index3, index4 = i, k, j, l
    
            # The sorting according to the first index of each tuple of each of the two convention is to match
            # Madgraph's convention for sorting color connection in the color correlators definition
            return (
                tuple(sorted([ (index1,-1,index1), (index2,-2,index2) ], key = lambda el: el[0], reverse=True)),
                tuple(sorted([ (index3,-1,index3), (index4,-2,index4) ], key = lambda el: el[0], reverse=True))
            )
       
        # Below is a purely abelian combination of the two single softs:
        #  -> second line of Eq. (6.13) of Gabor's hep-ph/0502226v2
        
        def abelian_combination((i1,j1),(i2,j2)):
            # The two terms of the symmetrised sum
            correlator_A   = double_correlator((i1,j1),(i2,j2))
            correlator_B   = double_correlator((i2,j2),(i1,j1))
            # There is an extra factor 2 because in the IR subtraction module, only one combination
            # S(r) S(s) is considered and not the symmetric version S(s) S(r)
            overall_factor = (1.0/4.0)*2.0

            # Group them if equal:
            if correlator_A==correlator_B:
                return [ (correlator_A,), ], [ 2.0*overall_factor ]
            else:
                # Implementing them as separate calls to the ME, as done in the commented line below,
                # return [ (correlator_A,),  (correlator_B,) ], [ overall_factor, overall_factor ]
                # would only yield the same result once the accessor correctly sums over the
                # different color-correlated MEs, so avoid for now.
                return [ (correlator_A,), (correlator_B,) ], [ overall_factor, overall_factor ]

        #return VH_two_iterated_soft_combination((i1_,j1_),(i2_,j2_))
        return abelian_combination((i1_,j1_),(i2_,j2_))

    @classmethod
    def combine_spin_correlators(cls, spin_correlators):
        """ This function takes several spin-correlators specified int the form
              spin_correlator = ( (legIndex, ( vector_A, vector_B, ...) ),
                                 (anotherLegIndex, ( vector_C, vector_D, ...) ),
                                 etc... ) 
            and returns the list of new correlator specifiers that arises from the
            combination of all of them.
        """
        
        # Trivial combination if there is a single one:
        if len(spin_correlators):
            return spin_correlators[0]

        # This combination is done with a simple concatenation of the lists. Example:
        # Let's say correlator A only correlates to leg #1 with two four-vectors v_A and v_B
        # (these can be seen as a replacement of polarization vectors to consider and summed over)
        #     correlators_A[0]  = ( (1, (v_A, v_B)), )
        # Now correlator B correlates with the two legs #4 and #7 (they must be different 
        #  by construction!), each defined with a single vector v_C and v_B
        #     correlators_B[0]  = ( (4, (v_C,)), (7, (v_D,)) )
        # Then it is clear tha the combined spin correlation should be:
        #     combined_spin_correlator = ( (1, (v_A, v_B)), (4, (v_C,)), (7, (v_D,)) )
        # Notice that this implies that both combinations :
        #    pol_vec_1 = v_A, pol_vec_4 = v_C,  pol_vec_7 = v_D
        # as well as:
        #    pol_vec_1 = v_B, pol_vec_4 = v_C,  pol_vec_7 = v_D
        # will be computed and summed in the resulting spin-correlated matrix element call.
        
        # Make sure the spin correlators don't share common legs
        for i, sc_a in enumerate(spin_correlators):
            for sc_b in spin_correlators[i+1:]:
                assert (len( set(c[0] for c in sc_a)&set(c[0] for c in sc_b))==0 )
        
        return tuple(sum([ (list(sc) if not sc is None else []) for sc in spin_correlators],[]))

    @classmethod
    def merge_correlators_in_necessary_ME_calls(cls, all_necessary_ME_calls):
        """ Merges the correlators defined the list of inputs to provide to the various calls
        to the reduced matrix elements specified by the argument all_necessary_ME_calls. """
        
        new_all_necessary_ME_calls = []
        for ME_call in all_necessary_ME_calls:
            combined_spin_correlator  = cls.combine_spin_correlators(ME_call['spin_correlations'])
            # The combination of color correlators can give rise to several ones, each with its multiplier
            combined_color_correlators, multipliers = cls.combine_color_correlators(ME_call['color_correlations'])
            # For now only multiply the finite part of currents
            combined_weight = reduce(lambda x,y: x*y, ME_call['main_weights'])
            flavor_matrices_beam_one = [fv for fv in ME_call['flavor_matrices_beam_one'] if fv is not None]
            flavor_matrices_beam_two = [fv for fv in ME_call['flavor_matrices_beam_two'] if fv is not None]
            if len(flavor_matrices_beam_one)>1 or len(flavor_matrices_beam_two)>1:
                raise MadGraph5Error("MadNkLO currently does not support the combination of multiple beam convolution flavor matrices.")
            flavor_matrix_beam_one = flavor_matrices_beam_one[0] if len(flavor_matrices_beam_one)==1 else None
            flavor_matrix_beam_two = flavor_matrices_beam_two[0] if len(flavor_matrices_beam_two)==1 else None
            Bjorken_rescalings_beam_one = [br for br in ME_call['Bjorken_rescalings_beam_one'] if br is not None]
            Bjorken_rescalings_beam_two = [br for br in ME_call['Bjorken_rescalings_beam_two'] if br is not None]
            if len(Bjorken_rescalings_beam_one)>1 or len(Bjorken_rescalings_beam_two)>1:
                raise MadGraph5Error("MadNkLO currently does not support the combination of multiple beam Bjorken x's rescalings.")
            Bjorken_rescaling_beam_one = Bjorken_rescalings_beam_one[0] if len(Bjorken_rescalings_beam_one)==1 else 1.0
            Bjorken_rescaling_beam_two = Bjorken_rescalings_beam_two[0] if len(Bjorken_rescalings_beam_two)==1 else 1.0
            for combined_color_correlator, multiplier in zip(combined_color_correlators,multipliers):
                # Finally add the processed combination as a new ME call
                new_all_necessary_ME_calls.append(
                    {   'spin_correlation'            : combined_spin_correlator,
                        'color_correlation'           : combined_color_correlator,
                        'main_weight'                 : combined_weight*multiplier,
                        'flavor_matrix_beam_one'      : flavor_matrix_beam_one,
                        'flavor_matrix_beam_two'      : flavor_matrix_beam_two,
                        'Bjorken_rescaling_beam_one'  : Bjorken_rescaling_beam_one,
                        'Bjorken_rescaling_beam_two'  : Bjorken_rescaling_beam_two,
                    },
                )
        return new_all_necessary_ME_calls

    @classmethod
    def update_all_necessary_ME_calls(cls, all_necessary_ME_calls, new_evaluation, 
            weight_type='main_weight', Bjorken_rescaling_beam_one=None, Bjorken_rescaling_beam_two=None):
        """ Combined previous current evaluation with the new one in argument so as to setup
        the input of the matrix elements to be computed. The current type can be either:
            'main_weight', 'flavor_matrices_beam_one', 'flavor_matrices_beam_two'
        which indicates if this is a weight that necessitate flavor convolutions (i.e. from
        beam factorization).
        """ 
        new_all_necessary_ME_calls = []
        for ((spin_index, color_index), current_wgt) in new_evaluation['values'].items():
            # Now combine the correlators necessary for this current
            # with those already specified in 'all_necessary_ME_calls'
            for ME_call in all_necessary_ME_calls:
                new_all_necessary_ME_calls.append({
                    # Append this spin correlation to those already present for that call
                    'spin_correlations' : ME_call['spin_correlations'] + [new_evaluation['spin_correlations'][spin_index], ],
                    # Append this color correlation to those already present for that call
                    'color_correlations' : ME_call['color_correlations'] + [new_evaluation['color_correlations'][color_index], ],
                    # Append this weight to those already present for that call
                    'main_weights' : ME_call['main_weights'] + [ base_objects.EpsilonExpansion(current_wgt) if 
                        weight_type=='main_weight' else base_objects.EpsilonExpansion({'finite':1.0}) ],
                    'flavor_matrices_beam_one' : ME_call['flavor_matrices_beam_one'] + [ current_wgt if weight_type=='flavor_matrix_beam_one' else None ],
                    'flavor_matrices_beam_two' : ME_call['flavor_matrices_beam_two'] + [ current_wgt if weight_type=='flavor_matrix_beam_two' else None ],
                    'Bjorken_rescalings_beam_one' : ME_call['Bjorken_rescalings_beam_one'] + [ Bjorken_rescaling_beam_one,],
                    'Bjorken_rescalings_beam_two' : ME_call['Bjorken_rescalings_beam_two'] + [ Bjorken_rescaling_beam_two,],
                })
        # Return the new list of necessary ME calls
        return new_all_necessary_ME_calls

    @classmethod
    def process_beam_factorization_currents(cls, all_necessary_ME_calls, all_beam_currents, 
                     all_MEAccessors, PS_point, process, chsi1, chsi2, mu_r, mu_f1, mu_f2):
        """ Calls the beam currents specified in the argument all_beam_currents with format:
              [{ 'beam_one': <BeamCurrent>,
               'beam_two': <BeamCurrent> }]
        and combine their evaluation with the list of already existing inputs provided in
        the list of all_necessary_ME_calls."""

        new_all_necessary_ME_calls = []
        for beam_currents in all_beam_currents:
            new_necessary_ME_calls = list(all_necessary_ME_calls)
            if beam_currents['beam_one'] is not None:
                if isinstance(beam_currents['beam_one'], subtraction.BeamCurrent) and beam_currents['beam_one']['distribution_type']=='bulk':
                    rescaling = 1./chsi1
                else:
                    rescaling = 1.
                current_evaluation, all_current_results = all_MEAccessors(beam_currents['beam_one'], 
                    lower_PS_point=PS_point, reduced_process=process, chsi=chsi1, mu_r=mu_r, mu_f=mu_f1)
                new_necessary_ME_calls = ME7Integrand_R.update_all_necessary_ME_calls(
                    new_necessary_ME_calls, current_evaluation,
                    weight_type='flavor_matrix_beam_one', Bjorken_rescaling_beam_one = rescaling)
            if beam_currents['beam_two'] is not None:
                if isinstance(beam_currents['beam_two'], subtraction.BeamCurrent) and beam_currents['beam_two']['distribution_type']=='bulk':
                    rescaling = 1./chsi2
                else:
                    rescaling = 1.
                current_evaluation, all_current_results = all_MEAccessors(beam_currents['beam_two'],
                    lower_PS_point=PS_point, reduced_process=process, chsi=chsi2, mu_r=mu_r, mu_f=mu_f2) 
                new_necessary_ME_calls = ME7Integrand_R.update_all_necessary_ME_calls(
                    new_necessary_ME_calls, current_evaluation,  
                    weight_type='flavor_matrix_beam_two', Bjorken_rescaling_beam_two = rescaling)
            new_all_necessary_ME_calls.extend(new_necessary_ME_calls)
        return new_all_necessary_ME_calls

    def evaluate_counterterm( self, counterterm, PS_point, base_weight, mu_r, mu_f1, mu_f2, 
        xb_1, xb_2, chsi1, chsi2, is_process_mirrored, all_resolved_flavors, hel_config=None, 
        apply_flavour_blind_cuts=True):
        """Evaluate a counterterm for a given PS point and flavors."""

        # Now call the walker which hikes through the counterterm structure
        # to return the list of currents and PS points for their evaluation
        # hike_output = {
        #         'currents' : [stroll_output1, stroll_output2, ...],
        #         'matrix_element': (ME_process, ME_PS),
        #         'kinematic_variables' : kinematic_variables (a dictionary) }

        hike_output = self.walker.walk_to_lower_multiplicity(
            PS_point, counterterm, compute_jacobian=self.divide_by_jacobian )

        # Access the matrix element characteristics
        ME_process, ME_PS = hike_output['matrix_element']
        
        # Generate what is the kinematics (reduced_PS) returned as a list
        # and the reduced_flavors for this counterterm by using the default reduced flavors
        # originating from the defining process and the real-emission kinematics dictionary
        reduced_PS, reduced_flavors = counterterm.get_reduced_quantities(
                                                              ME_PS, defining_flavors=None)
        
        n_unresolved_left = self.contribution_definition.n_unresolved_particles
        n_unresolved_left -= counterterm.count_unresolved()
        # Apply cuts if requested and return immediately if they do not pass
        if apply_flavour_blind_cuts and not self.pass_flavor_blind_cuts(reduced_PS,
                reduced_flavors, xb_1 = xb_1, xb_2 = xb_2, 
                n_jets_allowed_to_be_clustered=n_unresolved_left):
            # Return None to indicate that no counter-event was generated.
            return None

        # Compute all the reduced flavor configurations for this counterterm
        all_reduced_flavors = [ counterterm.get_reduced_flavors(resolved_flavors) 
                                             for resolved_flavors in all_resolved_flavors ]

        # The above "hike" can be used to evaluate the currents first and the ME last.
        # Note that the code below can become more complicated when tracking helicities,
        # but let's forget this for now.
        assert ((hel_config is None))
        # all_necessary_ME_calls is a list inputs to call the Matrix Element and the weights
        # that we must multiply/convolve them with. We start with empty entries.
        all_necessary_ME_calls = [ 
            {   'spin_correlations'             : [ ],
                'color_correlations'            : [ ],
                'main_weights'                  : [ ],
                'flavor_matrices_beam_one'      : [ ],
                'flavor_matrices_beam_two'      : [ ],
                'Bjorken_rescalings_beam_one'   : [ ],
                'Bjorken_rescalings_beam_two'   : [ ],
            },
        ]
        total_jacobian = 1.
        disconnected_currents_weight = base_objects.EpsilonExpansion({'finite': 1.0})

        # First evaluate the currents building the counterterm
        for stroll_output in hike_output['currents']:
            stroll_currents = stroll_output['stroll_currents']
            higher_PS_point = stroll_output['higher_PS_point']
            lower_PS_point  = stroll_output['lower_PS_point']
            stroll_vars     = stroll_output['stroll_vars']
            if stroll_vars.has_key('jacobian'):
                total_jacobian *= stroll_vars['jacobian']
            for current in stroll_currents:
                # WARNING The use of reduced_process here is fishy (for all but the last)

                current_evaluation, all_current_results = self.all_MEAccessors(
                    current,
                    higher_PS_point=higher_PS_point, lower_PS_point=lower_PS_point,
                    leg_numbers_map=counterterm.momenta_dict,
                    reduced_process=ME_process, hel_config=None,
                    **stroll_vars )

                # Now loop over all spin- and color- correlators required for this current
                # and update the necessary calls to the ME
                if not current['resolve_mother_spin_and_color']:
                    # Make sure no spin- or color-correlations were produced by the current
                    assert(current_evaluation['spin_correlations']==[None,])
                    assert(current_evaluation['color_correlations']==[None,])
                    assert(current_evaluation['values'].keys()==[(0,0),])
                    # WARNING:: this can only work for local 4D subtraction counterterms!
                    # For the integrated ones it is very likely that we cannot use a nested structure,
                    # and there will be only one level anyway,
                    # so there is not need of fancy combination of Laurent series.
                    disconnected_currents_weight *= \
                         base_objects.EpsilonExpansion(current_evaluation['values'][(0,0)])
                else:
                    all_necessary_ME_calls = ME7Integrand_R.update_all_necessary_ME_calls(
                        all_necessary_ME_calls, current_evaluation, weight_type='main_weight')

        # Then evaluate the beam factorization currents
        all_necessary_ME_calls = ME7Integrand_R.process_beam_factorization_currents(
            all_necessary_ME_calls, counterterm.get_beam_currents(), self.all_MEAccessors, 
            ME_PS, ME_process, chsi1, chsi2, mu_r, mu_f1, mu_f2)

        # Now perform the combination of the list of spin- and color- correlators to be merged
        # for each necessary ME call identified
        all_necessary_ME_calls = ME7Integrand_R.merge_correlators_in_necessary_ME_calls(
                                                                    all_necessary_ME_calls)

        # Finally treat the call to the reduced connected matrix elements
#        misc.sprint('I got for %s:'%str(counterterm.nice_string()))
        alpha_s = self.model.get('parameter_dict')['aS']
        mu_r = self.model.get('parameter_dict')['MU_R']

        return ME7Integrand_R.generate_event_for_counterterm(
            ME7Event( ME_PS, 
                {fc : base_weight for fc in all_reduced_flavors},
                is_mirrored                     = is_process_mirrored,
                host_contribution_definition    = self.contribution_definition,
                counterterm_structure           = counterterm),
            disconnected_currents_weight,
            ( counterterm.prefactor / total_jacobian ), 
            all_necessary_ME_calls,
            ME_process,
            ME_PS, 
            alpha_s, mu_r,
            self.all_MEAccessors
        )

    @classmethod
    def generate_event_for_counterterm(cls, template_MEEvent, disconnected_currents_weight,
            overall_prefactor, all_necessary_ME_calls, ME_process, ME_PS, alpha_s, mu_r, all_MEAccessors):

        # Initialize the placeholder which stores the event that we will progressively build
        # below.
        ME7_event_to_return = None
        for ME_call in all_necessary_ME_calls:
            color_correlators = tuple(ME_call['color_correlation']) if ME_call['color_correlation'] else None
            spin_correlators = tuple(ME_call['spin_correlation']) if ME_call['spin_correlation'] else None
            connected_currents_weight = ME_call['main_weight']
#            misc.sprint(ME_PS,ME_process.nice_string(),ME_process.get_cached_initial_final_numbers())
#            misc.sprint(spin_correlators)
#            misc.sprint(color_correlators, connected_currents_weight)
            try:
                ME_evaluation, all_ME_results = all_MEAccessors(
                   ME_process, ME_PS, alpha_s, mu_r,
                   # Let's worry about the squared orders later, we will probably directly fish
                   # them out from the ME_process, since they should be set to a unique combination
                   # at this stage.
                   squared_orders    = None,
                   color_correlation = color_correlators,
                   spin_correlation  = spin_correlators, 
                   hel_config        = None 
                )
#                misc.sprint(ME_evaluation)

            except MadGraph5Error as e:
                logger.critical("""
A reduced matrix element is missing in the library of automatically generated matrix elements.
This is typically what can happen when your process definition is not inclusive over all IR sensitive particles.
Make sure that your process definition is specified using the relevant multiparticle labels (typically 'p' and 'j').
Also make sure that there is no coupling order specification which receives corrections.
The missing process is: %s"""%ME_process.nice_string())
                raise e
            # for i in ME_PS.keys():
            #     if i in (1, 2): continue
            #     for j in ME_PS.keys():
            #         if j in (1, 2) or j <= i: continue
            #         misc.sprint(i, j, (ME_PS[i]+ME_PS[j]).square())
            # misc.sprint(current_weight.get_term('finite'), ME_evaluation['finite'])
            # misc.sprint('reduced process = %s' % (
            #     ' '.join('%d(%d)' % (l.get('number'), l.get('id')) for l in
            #              counterterm.process.get_initial_legs()) + ' > ' +
            #     ' '.join('%d(%d)' % (l.get('number'), l.get('id')) for l in
            #              counterterm.process.get_final_legs())
            # ))
            # misc.sprint(counterterm.prefactor)
            # misc.sprint(
            #     'color corr. = %-20s | current = %-20.16f | ME = %-20.16f | Prefactor = %-3f  |  Final = %-20.16f ' % (
            #         str(color_correlators),
            #         current_weight.get_term('finite'),
            #         ME_evaluation['finite'],
            #         counterterm.prefactor,
            #         current_weight.get_term('finite') * ME_evaluation['finite'] * counterterm.prefactor
            #     ))
            
            # Multiply the various pieces building the event weight (most being Epsilon expansions):
            # The ME evaluation, disconnected current weights and connected current weight
            event_weight = base_objects.EpsilonExpansion(ME_evaluation) * disconnected_currents_weight * connected_currents_weight                
            # And the prefactor
            event_weight *= overall_prefactor

            # Skip an event with no contribution (some dipoles of the eikonal for example)
            if event_weight.norm() == 0.:
                continue

            # Now build the event from the template provided
            event_to_convolve  = template_MEEvent.get_copy()
            event_to_convolve *= event_weight
            
            # The PDF Bjorken x's arguments will need a 1/z rescaling due to the change
            # of variable making the + distribution act on the PDF only.
            event_to_convolve.set_Bjorken_rescalings(ME_call['Bjorken_rescaling_beam_one'],
                                                     ME_call['Bjorken_rescaling_beam_two'] )

            # And convolve it if necessary
            if ME_call['flavor_matrix_beam_one'] is not None:
                event_to_convolve.convolve_flavors( ME_call['flavor_matrix_beam_one'], leg_index=0 )
            if ME_call['flavor_matrix_beam_two'] is not None:    
                event_to_convolve.convolve_flavors( ME_call['flavor_matrix_beam_two'], leg_index=1 )
            
            # Aggregate this event with the previous one, they should always be compatible
            if ME7_event_to_return is None:
                ME7_event_to_return = event_to_convolve
            else:
                ME7_event_to_return += event_to_convolve

        # Finally return the ME7 event constructed for this counterterm
        return ME7_event_to_return
        
    def sigma(self, PS_point, process_key, process, all_flavor_configurations, 
                  base_weight, mu_r, mu_f1, mu_f2, chsi1, chsi2, xb_1, xb_2, *args, **opts):
        """ Implementation of the short-distance cross-section for the real-emission integrand.
        Counterterms will be computed on top of the actual real-emission integrand.
        Note that the PS_point specified corresponds to a point already multiplied by both
        the bjorken x's and boosted back to the c.o.m frame."""

        all_events_generated = ME7EventList()

        alpha_s = self.model.get('parameter_dict')['aS']

        apply_flavour_blind_cuts = opts.get('apply_flavour_blind_cuts', True)

        ME_evaluation, all_results = self.all_MEAccessors(
                        process, PS_point, alpha_s, mu_r, pdgs=all_flavor_configurations[0])

        event_weight = base_objects.EpsilonExpansion(ME_evaluation) * base_weight
        
        # Notice that the Bjorken rescalings for this even will be set during the convolution
        # performed in the next line, and not directly here in the constructor.
        event_to_convolve = ME7Event( PS_point, 
                {fc : event_weight for fc in all_flavor_configurations},
                is_mirrored                     = process.get('has_mirror_process'),
                host_contribution_definition    = self.contribution_definition,
                counterterm_structure           = None
        )
        convolved_event = self.convolve_event_with_beam_factorization_currents(event_to_convolve, 
            process['beam_factorization'], PS_point, process, mu_r, mu_f1, mu_f2, chsi1, chsi2 )
        
        # Now register the convoluted event
        all_events_generated.append(convolved_event)
        
        for counterterm in self.counterterms[process_key]:
            if not counterterm.is_singular():
                continue
            CT_event = self.evaluate_counterterm(
                  counterterm, PS_point, base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, chsi1, chsi2, 
                  process.get('has_mirror_process'), all_flavor_configurations, 
                  hel_config = None, apply_flavour_blind_cuts = apply_flavour_blind_cuts)
            # The function above returns None if the counterterms is removed because of
            # flavour_blind_cuts.
            if CT_event is not None:
                all_events_generated.append( CT_event )

        # Notice that it may be possible in some subtraction scheme to combine
        # several events having the same kinematics support. This is a further
        # optimization that we can easily implement when proven needed.
        # all_events_generated.combine_events_with_identical_kinematics()
        # (PS: the function above is a place-holder and is not implemented yet.
        
        return all_events_generated

    def test_IR_limits(self, test_options):
        """Test how well local counterterms approximate a real-emission matrix element."""

        # Apply the passed options
        seed = test_options.get('seed', None)
        if seed: random.seed(seed)
        walker_name = test_options.get('walker', None)
        if walker_name is None:
            walker = self.walker
        else:
            walker = walkers.VirtualWalker(walker_name)

        # First generate an underlying Born
        # Specifying None forces to use uniformly random generating variables.
        # Make sure to generate a point within the cuts if necessary:
        max_attempts = 10000
        n_attempts   = 0
        while True:
            n_attempts += 1
            if n_attempts > max_attempts:
                break
            real_emission_PS_point = None
            # Phase-space generation can fail when beam convolution is active, because of the condition
            # Bjorken_x_i < chsi_i
            while real_emission_PS_point is None:
                real_emission_PS_point, jac, x1s, x2s = self.phase_space_generator.get_PS_point(None)
            xb_1, chsi1 = x1s
            xb_2, chsi2 = x2s
            if test_options['apply_higher_multiplicity_cuts']:
                if not self.pass_flavor_blind_cuts(
                    real_emission_PS_point,
                    self.processes_map.values()[0][0].get_cached_initial_final_pdgs(),
                    n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles,
                    xb_1 = xb_1, xb_2 = xb_2):
                    continue
            break
        if n_attempts > max_attempts:
            raise MadEvent7Error(
                "Could not generate a random kinematic configuration that passes " +
                "the flavour blind cuts in less than %d attempts." % max_attempts )
        n_attempts = 0

        # Loop over processes
        all_evaluations = {}
        for process_key, (defining_process, mapped_processes) in self.processes_map.items():
            logger.debug("Considering %s"%defining_process.nice_string())
            # Make sure that the selected process satisfies the selection requirements
            if not self.is_part_of_process_selection(
                [defining_process, ]+mapped_processes,
                selection=test_options['process'] ):
                continue
            
            all_processes = [defining_process,]+mapped_processes
            all_flavor_configurations = []
            # The process mirroring is accounted for at the very end only            
            for proc in all_processes:
                initial_final_pdgs = proc.get_cached_initial_final_pdgs()
                all_flavor_configurations.append(initial_final_pdgs)
            
            a_real_emission_PS_point = real_emission_PS_point.get_copy()
            a_xb_1, a_chsi1 = x1s
            a_xb_2, a_chsi2 = x2s
            while (test_options['apply_higher_multiplicity_cuts'] and
                   not self.pass_flavor_sensitive_cuts(
                       a_real_emission_PS_point,
                       defining_process.get_cached_initial_final_pdgs(),
                       xb_1 = xb_1, xb_2 = xb_2 ) ):
                n_attempts += 1
                if n_attempts > max_attempts:
                    break
                a_real_emission_PS_point = None
                # Phase-space generation can fail when beam convolution is active, because of the condition
                # Bjorken_x_i < chsi_i
                while a_real_emission_PS_point is None:
                    a_real_emission_PS_point, a_jac, a_x1s, a_x2s = self.phase_space_generator.get_PS_point(None)
                a_xb_1, a_chsi1 = a_x1s
                a_xb_2, a_chsi2 = a_x2s
            if n_attempts > max_attempts:
                raise MadEvent7Error(
                    "Could not generate a random kinematic configuration that passes " +
                    "the flavour blind cuts in less than %d attempts." % max_attempts )
            n_attempts = 0

            # Make sure to have the PS point provided as LorentzVectorDict
            a_real_emission_PS_point = phase_space_generators.LorentzVectorDict(
                (i+1, mom) for i, mom in enumerate(a_real_emission_PS_point) )

            # Use correction_order to select CT subset
            counterterms_to_consider = [
                ct for ct in self.counterterms[process_key]
                if ct.get_number_of_couplings() <= test_options['correction_order'].count('N') ]

            selected_singular_structures = []

            for limit_specifier in test_options['limits']:
                # Select the limits to be probed interpreting limits as a regex pattern.
                # If no match is found, then reconstruct the singular structure from the limits
                # provided
                selected_counterterms = self.find_counterterms_matching_regexp(
                    counterterms_to_consider, limit_specifier )
                if selected_counterterms:
                    selected_singular_structures.extend([
                        ct.reconstruct_complete_singular_structure()
                        for ct in selected_counterterms])
                else:
                    ss = subtraction.SingularStructure.from_string(
                        limit_specifier, defining_process)
                    if ss is None:
                        logger.info("No limit matching %s for process %s." % 
                            (limit_specifier, defining_process.nice_string()) )
                        continue
                    selected_singular_structures.append(ss)

            # Filter the singular structures to consider so as to remove those involving
            # a beam factorisation not present in this integrand
            selected_singular_structures = [
                ss for ss in selected_singular_structures if not (
                    (1 in ss.get_beam_factorization_legs() and a_chsi1 is None) or 
                    (2 in ss.get_beam_factorization_legs() and a_chsi2 is None))
            ]
            
            if len(selected_singular_structures)==0:
                logger.warning('Empty selection of limits when investigating process %s'%defining_process.nice_string())
                continue

            logger.debug('Reconstructed complete singular structure: \n'+'\n'.join(
                str(ss) for ss in selected_singular_structures ))

            # Loop over approached limits
            process_evaluations = {}
            for limit in selected_singular_structures:
                limit_evaluations = self.test_IR_limits_for_limit_and_process(
                    test_options, walker, limit, defining_process, process_key, all_flavor_configurations, 
                    a_real_emission_PS_point, a_chsi1, a_chsi2, a_xb_1, a_xb_2 )
                process_evaluations[str(limit)] = limit_evaluations

            process_string = defining_process.base_string()
            if defining_process.has_key('n_loops'):
                process_string += " @ " + str(defining_process['n_loops']) + " loops"
            all_evaluations[process_string] = process_evaluations

        # Now produce a nice matplotlib of the evaluations
        # and assess whether this test passed or not
        return self.analyze_IR_limits_test(
            all_evaluations, 
            test_options['acceptance_threshold'],
            seed               = seed,
            show               = test_options['show_plots'],
            save_plots         = test_options['save_plots'],
            save_results_path  = test_options['save_results_to_path'],
            plots_suffix       = test_options['plots_suffix'],
        )

    def test_IR_limits_for_limit_and_process(self, test_options, walker, limit, defining_process, process_key, 
        all_flavor_configurations, a_real_emission_PS_point, a_chsi1, a_chsi2, a_xb_1, a_xb_2):
        """ Given a test_options, a process, a specific limit, a walker to approach it,
        all flavor configurations mapped to that particular process, a 'starting' 
        real-emission type of PS-point from which to scale towards the limit, as well as
        its corresponding convolution variables chsi1/2, perform the IR limit test that 
        consists in approaching the limit and computing both the real-emission matrix 
        element and all counterterms that should subtract its singularities."""
                
        logger.debug("Approaching limit %s " % str(limit) )
        
        # First select the counterterms to evaluate and temporarily assign them to this
        # integrand instance so that the sigma function will run on them.
        if self.has_local_counterterms():
            local_counterterms_to_consider = [ ct for ct in self.counterterms[process_key]
                if ct.get_number_of_couplings() <= test_options['correction_order'].count('N') ]
            if test_options['counterterms']:
                counterterm_pattern = test_options['counterterms']
                if counterterm_pattern.startswith('def'):
                    counterterm_pattern = str(limit)
                logger.debug("Local counterterms before selection")
                for ct in local_counterterms_to_consider:
                    logger.debug("    "+str(ct))
                local_counterterms_to_consider = self.find_counterterms_matching_regexp(
                    local_counterterms_to_consider, counterterm_pattern )
            logger.debug("Selected local counterterms")
            for ct in local_counterterms_to_consider:
                logger.debug("    "+str(ct))
        
        # We must also include the integrated counterterms as they may have a chsi dependence
        if self.has_integrated_counterterms():
            integrated_counterterms_to_consider = [ ct for ct in self.integrated_counterterms[process_key] 
                if ct['integrated_counterterm'].get_number_of_couplings() <= test_options['correction_order'].count('N') ]
            if test_options['counterterms']:
                counterterm_pattern = test_options['counterterms']
                if counterterm_pattern.startswith('def'):
                    counterterm_pattern = str(limit)
                logger.debug("Integrated counterterms before selection")
                for ct in integrated_counterterms_to_consider:
                    logger.debug("    "+str(ct['integrated_counterterm']))
                integrated_counterterms_to_consider = self.find_counterterms_matching_regexp(
                    integrated_counterterms_to_consider, counterterm_pattern, integrated_CT=True )
            logger.debug("Selected integrated counterterms")
            for ct in integrated_counterterms_to_consider:
                logger.debug("    "+str(ct['integrated_counterterm']))

        # Now progressively approach the limit, using a log scale
        limit_evaluations = {}
        n_steps = test_options['n_steps']
        min_value = test_options['min_scaling_variable']
        base = min_value ** (1./n_steps)
        for step in range(n_steps+1):
            # Determine the new phase-space point
            scaling_parameter = base ** step
            scaled_real_PS_point = walker.approach_limit(
                     a_real_emission_PS_point, limit, scaling_parameter, defining_process )
            # Also scale the chsi initial-state convolution parameters if the limit
            # specifies a beam factorization structure for that initial state:
            beam_factorisation_legs = limit.get_beam_factorization_legs()
            if 1 in beam_factorisation_legs:
                scaled_chsi1 = 1.-(1.-a_chsi1)*scaling_parameter
            else:
                scaled_chsi1 = a_chsi1
            if 2 in beam_factorisation_legs:
                scaled_chsi2 = 1.-(1.-a_chsi2)*scaling_parameter
            else:
                scaled_chsi2 = a_chsi2
            if test_options['apply_higher_multiplicity_cuts']:
                if not self.pass_all_cuts( scaled_real_PS_point,
                        self.processes_map.values()[0][0].get_cached_initial_final_pdgs(),
                        n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles,
                        xb_1 = a_xb_1, xb_2 = a_xb_2 ):
                    logger.warning('Aborting prematurely since the following scaled real-emission point'+
                                  ' does not pass higher multiplicity cuts.')
                    break

#            misc.sprint('Scaled PS point: %s'%str(a_real_emission_PS_point))
            mu_r, mu_f1, mu_f2 = self.get_scales(scaled_real_PS_point)

            # Specify the counterterms to be considered and backup the original values
            if self.has_local_counterterms():
                local_CT_backup = self.counterterms[process_key]
                self.counterterms[process_key] = local_counterterms_to_consider
            if self.has_integrated_counterterms():
                integrated_CT_backup = self.integrated_counterterms[process_key]
                self.integrated_counterterms[process_key] = integrated_counterterms_to_consider
            try:
                # Now call sigma in order to gather all events
                events = self.sigma(
                    scaled_real_PS_point.to_dict(), process_key, defining_process, all_flavor_configurations, 
                    1.0, mu_r, mu_f1, mu_f2, scaled_chsi1, scaled_chsi2, a_xb_1, a_xb_2, 
                    compute_poles            = False,
                    apply_flavour_blind_cuts = test_options['apply_lower_multiplicity_cuts']
                )
            except Exception as e:
                logger.critical("The following exception occurred when generating events for this integrand %s:\n%s"%(
                                                      self.contribution_short_name, str(e)))
                if self.has_local_counterterms():
                    self.counterterms[process_key] = local_CT_backup
                if self.has_integrated_counterterms():
                    self.integrated_counterterms[process_key] = integrated_CT_backup
                # Forward the exception while not resetting the stack trace.
                raise

            # Restore the original list of counterterms
            if self.has_local_counterterms():
                self.counterterms[process_key] = local_CT_backup
            if self.has_integrated_counterterms():
                self.integrated_counterterms[process_key] = integrated_CT_backup

#            misc.sprint('Events generated before post-processing.')
#            misc.sprint(events)
            
            # Now post-process the events as done in __call__ of the integrand.
            # Apply flavor blind cuts
            events.filter_with_flavor_blind_cuts(self.pass_flavor_blind_cuts, 
                defining_process.get_cached_initial_final_pdgs(),
                n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles,
                xb_1 = a_xb_1, xb_2 = a_xb_2 )
            # Select particular terms of the EpsilonExpansion terms stored as weights.
            # For the test_IR_limits, we are only interested in the finite part.
            events.select_epsilon_expansion_term(test_options['epsilon_expansion_term'])
            # Apply PDF convolution. It is important to convolve with the PDFs, so as to 
            # correctly get the factors 1/chsi_i multiplying the weights in beam 
            # factorization counterterms. But either the PDF density can be set to 1.
            if self.run_card['lpp1']==self.run_card['lpp2']==1:
                pdf = None if test_options['set_PDFs_to_unity'] else self.pdf
                events.apply_PDF_convolution( self.get_pdfQ2, (pdf, pdf), (a_xb_1, a_xb_2), 
                                                                     (mu_f1**2, mu_f2**2) )
            # Make sure Bjorken-x rescalings chsi_i don't matter anymore

            for event in events:
                event.set_Bjorken_rescalings(None, None)
            # Apply flavor sensitive cuts
            if test_options['apply_lower_multiplicity_cuts']:
                events.filter_flavor_configurations(self.pass_flavor_sensitive_cuts,
                                                              xb_1 = a_xb_1, xb_2 = a_xb_2)
    
#            misc.sprint('Events generated after post-processing:')
#            misc.sprint(events)
            
            # Now store the results to be returned which will eventually be passed to 
            # the IR test analyzer.

            # Initialize result
            this_eval = {}
            # Loop over all events produced and distribute them in several entries of what
            # will be returned to the IR_analyzer
            for event in events:
                # TODO: instead of looking at the total weight, we could consider doing the
                # check for one particular flavor configuration which could be specified
                # as an extra test parameter.
                event_wgt = event.get_total_weight()
                if event.counterterm_structure is None:
                    if 'ME' in this_eval:
                        this_eval['ME'] += event_wgt
                    else:
                        this_eval['ME'] = event_wgt
                else:
                    if str(event.counterterm_structure) in this_eval:
                        this_eval[str(event.counterterm_structure)] += event_wgt           
                    else:
                        this_eval[str(event.counterterm_structure)] = event_wgt                      
                        
            logger.debug('For scaling variable %.3e, weight from ME = %.16f' %(
                                      scaling_parameter, this_eval['ME'] ))
            total_CTs_wgt = 0.0
            for CT_str, CT_weight in this_eval.items():
                if CT_str=='ME': continue
                total_CTs_wgt += CT_weight
                logger.debug('Weight from CT %s = %.16f' % (CT_str, CT_weight) )
                logger.debug('Ratio: %.16f'%( CT_weight/float(this_eval['ME']) ))
            if this_eval['ME'] > 0.:
                logger.debug('Ratio sum(CTs)/ME: %.16f'%(total_CTs_wgt/float(this_eval['ME'])))
            limit_evaluations[scaling_parameter] = this_eval
            
        # Now return all evaluations performed for each value of the scale
        return limit_evaluations

    @staticmethod
    def analyze_IR_limit(
        evaluations, acceptance_threshold,
        title=None, def_ct=None, plot_all=True, show=True,
        filename=None, plots_suffix=None, number_of_FS_legs=None ):

        import matplotlib.pyplot as plt

        plot_title = True
        plot_size = (6,6)
        plot_extension = ".pdf"
        if plots_suffix:
            plot_extension = '_' + plots_suffix + plot_extension
        test_failed = False
        test_ratio  = -1.0

        # Produce a plot of all counterterms
        x_values = sorted(evaluations.keys())
        lines = evaluations[x_values[0]].keys()
        lines.sort(key=len)
        # Skip ME-def line if there is no defining ct
        plot_def = def_ct and def_ct in lines
        plot_total = len(lines) > 2

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', size=16)
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        TOTAL_color = colors.pop(3)
        MEdef_color = colors.pop(2)

        plt.figure(1, figsize=plot_size)
        plt.gca().set_prop_cycle(color=colors)
        if plot_title and title: plt.title(title)
        GeV_pow = -2*(number_of_FS_legs-2)
        units = '[GeV${}^{' + str(GeV_pow) + '}$]'
        plt.xlabel('$\lambda$')
        plt.ylabel('Integrands ' + units)
        plt.subplots_adjust(left=0.15)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        total           = [0., ] * len(x_values)
        ME_minus_def_ct = [0., ] * len(x_values)
        line_labels = {}
        for line in lines:
            y_values = [abs(evaluations[x][line]) for x in x_values]
            for i in range(len(x_values)):
                total[i] += evaluations[x_values[i]][line]
            if plot_def and (line == "ME" or line == def_ct):
                def_ct_sign = 1
                if line == def_ct:
                    def_ct_sign = (-1) ** def_ct.count("(")
                for i in range(len(x_values)):
                    ME_minus_def_ct[i] += def_ct_sign * evaluations[x_values[i]][line]
            line_label = copy.copy(line)
            if line_label.startswith("("):
                line_label = line_label[1:]
            if line_label.endswith(",)"):
                line_label = line_label[:-2]
            line_labels[line] = line_label
            if plot_all:
                if '(' in line:
                    style = '--'
                else:
                    style = '-'
                plt.plot(x_values, y_values, style, label=line_label)
        if plot_def:
            abs_ME_minus_def_ct = [abs(y) for y in ME_minus_def_ct]
            plt.plot(x_values, abs_ME_minus_def_ct, color=MEdef_color, label='ME-def')
        if plot_total:
            abs_total = [abs(y) for y in total]
            plt.plot(x_values, abs_total, color=TOTAL_color, label='TOTAL')
        plt.legend()
        if filename:
            plt.savefig(filename + '_integrands' + plot_extension)

        plt.figure(2, figsize=plot_size)
        plt.gca().set_prop_cycle(color=colors)
        if plot_title and title: plt.title(title)
        plt.xlabel('$\lambda$')
        plt.ylabel('Ratio to ME')
        plt.subplots_adjust(left=0.15)
        plt.xscale('log')
        plt.grid(True)
        for line in lines:
            y_values = [abs(evaluations[x][line]/evaluations[x]["ME"]) for x in x_values]
            if '(' in line:
                style = '--'
            else:
                style = '-'
            plt.plot(x_values, y_values, style, label=line_labels[line])
        plt.legend()
        if filename:
            plt.savefig(filename + '_ratios' + plot_extension)

        # Check that the ratio of def_ct to the ME is close to -1
        if plot_def and not test_failed:
            def_ct_2_ME_ratio = evaluations[x_values[0]][def_ct]
            def_ct_2_ME_ratio /= evaluations[x_values[0]]["ME"]
            foo_str = "The ratio of the defining CT to the ME at lambda = %s is: %s."
            logger.info(foo_str % (x_values[0], def_ct_2_ME_ratio))
            test_ratio = abs(def_ct_2_ME_ratio)-1
            test_failed = test_ratio > acceptance_threshold
        # Check that the ratio between total and ME is close to 0
        if plot_total and not test_failed:
            total_2_ME_ratio = total[0]
            total_2_ME_ratio /= evaluations[x_values[0]]["ME"]
            foo_str = "The ratio of the total to the ME at lambda = %s is: %s."
            logger.info(foo_str % (x_values[0], total_2_ME_ratio))
            test_ratio  = abs(total_2_ME_ratio)
            test_failed = test_ratio > acceptance_threshold

        plt.figure(3, figsize=plot_size)
        plt.gca().set_prop_cycle(color=colors)
        if plot_title and title: plt.title(title)
        plt.xlabel('$\lambda$')
        plt.ylabel('Weighted integrands ' + units)
        plt.subplots_adjust(left=0.15)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        if plot_def:
            wgt_ME_minus_def_ct = [abs(x_values[i] * ME_minus_def_ct[i])
                                   for i in range(len(x_values))]
            plt.plot(x_values, wgt_ME_minus_def_ct, color=MEdef_color, label='ME-def')
        if plot_total:
            wgt_total = [abs(x_values[i] * total[i]) for i in range(len(x_values))]
            plt.plot(x_values, wgt_total, color=TOTAL_color, label='TOTAL')
        plt.legend()
        if filename:
            plt.savefig(filename + '_weighted' + plot_extension)

        if show:
            plt.show()
        else:
            plt.close('all')

        return (not test_failed, test_ratio)

    def analyze_IR_limits_test(
        self, all_evaluations, acceptance_threshold,
        seed=None, show=True,
        save_plots=False, save_results_path=None, plots_suffix=None ):
        """Analyze the results of the test_IR_limits command."""

        test_failed = False
        results = dict()
        for (process, process_evaluations) in all_evaluations.items():
            results[process] = dict()
            for (limit, limit_evaluations) in process_evaluations.items():
                proc, loops = process.split("@")
                title = "$" + proc + "$"
                initial_state, final_state = proc.split('>')
                number_of_FS_legs = final_state.count(' ') - 1
                title = title.replace('~','x').replace('>','\\to').replace(' ','\\;')
                title = title.replace('+','^+').replace('-','^-')
                title += "@" + loops + " approaching " + limit
                if seed: title += " (seed %d)" % seed
                filename = None
                if save_plots:
                    filename = copy.copy(process)
                    filename = filename.replace(">", "_")
                    filename = filename.replace(" ", "").replace("~", "x")
                    filename = filename.replace("@", "_") + "_" + limit
                    filename = filename.replace(",", "").replace("(","").replace(")","")
                    if seed: filename += "_"+str(seed)
                results[process][limit] = self.analyze_IR_limit(
                    limit_evaluations, acceptance_threshold=acceptance_threshold,
                    title=title, def_ct=limit, show=show,
                    filename=filename, plots_suffix=plots_suffix,
                    number_of_FS_legs=number_of_FS_legs)
                if not results[process][limit][0]:
                    test_failed = True
        if save_results_path:
            output_lines = []
            for process in results:
                for limit, (test_outcome, limiting_ratio) in results[process].items():
                    output_lines.append('%-30s | %-20s | %-10s | %-30s'%(
                        process,
                        str(limit),
                        'PASSED' if test_outcome else 'FAILED',
                        '%.16e'%limiting_ratio
                    ))
            with open(save_results_path,'w') as results_output:     
                results_output.write('\n'.join(output_lines))

        if test_failed:
            logger.info("analyse_IR_limits_test results:")
            for process in results:
                logger.info("    " + str(process))
                for limit, (test_outcome, limiting_ratio) in results[process].items():
                    if test_outcome:
                        logger.info("    " * 2 + str(limit) + ": PASSED")
                    else:
                        logger.info("    " * 2 + str(limit) + ": FAILED")
            return False
        return True
    

class ME7Integrand_RV(ME7Integrand_R, ME7Integrand_V):
    """ ME7Integrand for the computation of integrands featuring both local and integrated
    counterterms."""
    
    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this contribution."""
        res = ME7Integrand_R.get_additional_nice_string_printout_lines(self, format=format)
        res.extend(ME7Integrand_V.get_additional_nice_string_printout_lines(self, format=format))
        return res
    
    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        contribution."""
        GREEN = '\033[92m'
        ENDC = '\033[0m'        
        res = GREEN+'  %s'%(defining_process.nice_string(print_weighted=False)
                                                             .replace('Process: ',''))+ENDC

        if not self.integrated_counterterms or not self.counterterms:
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
    
    def sigma(self, PS_point, process_key, process, all_flavor_configurations, base_weight, mu_r, 
                                                mu_f1, mu_f2, chsi1, chsi2, xb_1, xb_2, *args, **opts):
        # Start with the events from the process itself and its local counterterms.
        all_events_generated = ME7Integrand_R.sigma(self, PS_point, process_key, process, 
            all_flavor_configurations, base_weight, mu_r, mu_f1, mu_f2, chsi1, chsi2, xb_1, xb_2, *args, **opts)

        compute_poles = opts.get('compute_poles', False)

        # Add the integrated counterterms
        # Now loop over all integrated counterterms
        for counterterm_characteristics in self.integrated_counterterms[process_key]:
            # And over all the ways in which this current PS point must be remapped to
            # account for all contributions of the integrated CT. (e.g. the integrated
            # splitting g > d d~ must be "attached" to all final state gluons appearing
            # in this virtual ME process definition.)
            for input_mapping in counterterm_characteristics['input_mappings']:
                # At NLO at least, it is OK to save a bit of time by enforcing 'compute_poles=False').
                # This will need to be re-assessed at NNLO for RV contributions.
                CT_event = self.evaluate_integrated_counterterm( 
                    counterterm_characteristics, PS_point, base_weight, mu_f1, mu_f2, chsi1, chsi2,
                    process.get('has_mirror_process'),
                    input_mapping, all_flavor_configurations,
                    hel_config    = None, 
                    compute_poles = compute_poles)
            if CT_event is not None:
                all_events_generated.append( CT_event )
            else:
                logger.warning('The evaluation of integrated counterterms '+
                                                    'should never yield a None ME7 Event.')

        # Notice that it may be possible in some subtraction scheme to combine
        # several events having the same kinematics support. This is a further
        # optimization that we can easily implement when proven needed.
        # all_events_generated.combine_events_with_identical_kinematics()
        # (PS: the function above is a place-holder and is not implemented yet.

        return all_events_generated

        
class ME7Integrand_RR(ME7Integrand_R):
    """ ME7Integrand for the computation of a double real-emission type of contribution."""
    def sigma(self, PS_point, process_key, process, all_flavor_configurations, base_weight, mu_r, 
                                                mu_f1, mu_f2, chsi1, chsi2, xb_1, xb_2, *args, **opts):
        return super(ME7Integrand_RR, self).sigma(PS_point, process_key, process, 
            all_flavor_configurations, base_weight, mu_r, mu_f1, mu_f2, chsi1, chsi2, xb_1, xb_2, *args, **opts)

class ME7Integrand_RRR(ME7Integrand_R):
    """ ME7Integrand for the computation of a double real-emission type of contribution."""
    def sigma(self, PS_point, process_key, process, all_flavor_configurations, base_weight, 
                                          mu_r, mu_f1, mu_f2, chsi1, chsi2, xb_1, xb_2, *args, **opts):
        return super(ME7Integrand_RRR, self).sigma(PS_point, process_key, process,
            all_flavor_configurations, base_weight, mu_r, mu_f1, mu_f2, chsi1, chsi2, xb_1, xb_2, *args, **opts)

#===============================================================================
# ME7IntegrandList
#===============================================================================
class ME7IntegrandList(base_objects.PhysicsObjectList):
    """ Container for ME7Integrnds."""
    
        
    integrands_natural_order = [
        ('LO',    (ME7Integrand_B, ME7Integrand_LIB) ),
        ('NLO',   (ME7Integrand_R, ME7Integrand_V) ),
        ('NNLO',  (ME7Integrand_RR, ) ),
        ('NNNLO', (ME7Integrand_RRR, ) )
    ]
    
    def is_valid_element(self, obj):
        """Test if object obj is a valid instance of ME7Integrand."""
        return isinstance(obj, ME7Integrand)
    
    def get_integrands_of_order(self, correction_order):
        """ Returns a list of all contributions of a certain correction_order in argument."""
        return ME7IntegrandList([integrand for integrand in self if
                integrand.contribution_definition.correction_order==correction_order])

    def get_integrands_of_type(self, correction_classes):
        """ Returns a list of all contributions that are direct instances of certain classes."""
        if not isinstance(correction_classes, tuple):
            if isinstance(correction_classes, list):
                correction_classes = tuple(correction_classes)
            else:
                correction_classes = (correction_classes,)                
        return ME7IntegrandList([integrand for integrand in self if 
                                                isinstance(integrand, correction_classes)])

    def nice_string(self, format=0):
        """ A nice representation of a list of contributions. 
        We can reuse the function from ContributionDefinitions."""
        return base_objects.ContributionDefinitionList.contrib_list_string(self, 
                                                                            format=format)

    def sort_integrands(self):
        """ Sort integrands according to the order dictated by the class attribute
        'integrands_natural_order'"""
        
        new_order = []
        for correction_order, integrand_types in self.integrands_natural_order:
            for integrand_type in integrand_types:
                selected_integrands = self.get_integrands_of_order(correction_order).\
                                        get_integrands_of_type(integrand_type)
                new_order.extend(selected_integrands)
                for integrand in selected_integrands:
                    self.pop(self.index(integrand))
    
        # Finally all remaining contributions of unknown types
        new_order.extend(self)
        self[:] = new_order


# Integrand classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Integrand.
# Notice that this must be placed after all the Integrand daughter classes in this module have been declared.
ME7Integrand_classes_map = {'Born': ME7Integrand_B,
                            'LoopInduced_Born': ME7Integrand_LIB,
                            'Virtual': ME7Integrand_V,
                            'SingleReals': ME7Integrand_R,
                            'RealVirtual': ME7Integrand_RV,
                            'DoubleReals': ME7Integrand_RR,
                            'TripleReals': ME7Integrand_RRR,
                            'Unknown': None}
