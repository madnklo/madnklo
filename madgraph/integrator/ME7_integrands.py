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
from pprint import pformat

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
logger1 = logging.getLogger('madgraph')
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
import madgraph.integrator.vectors as vectors
#    import madgraph.various.histograms as histograms  # imported later to not slow down the loading of the code
import models.check_param_card as check_param_card
import models.model_reader as model_reader
import models.import_ufo as import_ufo

from madgraph.various.math_tools.su_n_group_constants import SU3

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
                 requires_mirroring           = False,
                 host_contribution_definition = 'N/A',
                 # This entry is mostly for debugging purposes, but it is useful to be able
                 # to know what is the list of defining PDGs of the resolved process from 
                 # which this counterterm originates. For this reason the counterterm_structure
                 # is (possible a list of) tuple(s) of the form:
                 #   (
                 #      integrated_CT_instance,
                 #      list_of_all combination_of_PDGs_of_the_resolved_process_this_integrated_CT_can_comes_from,
                 #      identifier_of_the_reduced_kinematics (None if unique for that counterterm or a string otherwise)
                 #   )
                 counterterm_structure        = None,
                 Bjorken_xs                   = (1.0,1.0),
                 Bjorken_x_rescalings         = (1.0,1.0),
                 is_a_mirrored_event          = False,
                 # The user is free to add below any additional information that 
                 # plays no structural role in in the MadNkLO construction. Example:
                 # test_IR_poles will store here information such as the input_mapping used
                 # and the convolutional mask employed for nicer printouts.
                 additional_information       = {}
                ):
        """Initialize this particular event with a unique PS point
        and a dictionary of the form:
             ( (initial_state_flavors,), (final_state_flavors,) ) : weight
        to specify all possible flavor configurations applying to this event.
        For PDF counterterms and integrated collinear CT counterterms, we must
        provide the possibility of applying a rescaling to the Bjorken x that will be
        used to compute each of the two beams. By default this rescaling is of course one."""
        
        self.PS_point = PS_point.to_list()
        self.weights_per_flavor_configurations = weights_per_flavor_configurations
        self.requires_mirroring = requires_mirroring
        self.is_a_mirrored_event = is_a_mirrored_event
        self.host_contribution_definition = host_contribution_definition
        self.counterterm_structure = counterterm_structure
        # This can either be None or a dictionary of the format:
        # {  'beam_one' : 'ALL' or (PDG_allowed_1, PDG_allowed_2, PDG_allowed_3, ...)
        #    'beam_two' : 'ALL' or (PDG_allowed_1, PDG_allowed_2, PDG_allowed_3, ...)
        # }
        self.Bjorken_x_rescalings = Bjorken_x_rescalings
        self.Bjorken_xs = Bjorken_xs
        self.additional_information = additional_information

        # This entry specifies which flavor has been selected for this event.
        # None implies that it has not been selected yet.
        self.selected_flavors = None
        
    def get_copy(self):
        """ Returns a shallow copy of self, except for the
        'weights_per_flavor_configurations' attribute. """
        
        return ME7Event(
            self.PS_point, dict(self.weights_per_flavor_configurations),
            requires_mirroring           = self.requires_mirroring,
            host_contribution_definition = self.host_contribution_definition,
            counterterm_structure        = self.counterterm_structure,
            Bjorken_x_rescalings         = self.Bjorken_x_rescalings,
            Bjorken_xs                   = self.Bjorken_xs,
            is_a_mirrored_event          = self.is_a_mirrored_event,
            additional_information       = copy.deepcopy(self.additional_information)
        )
        
    def set_Bjorken_rescalings(self, *args):
        """ Assigns the specified Bjorken x rescaling to this event."""
        self.Bjorken_x_rescalings = tuple(args)

    def set_Bjorken_xs(self, *args):
        """ Assigns the specified Bjorken x's."""
        self.Bjorken_xs = tuple(args)

    def apply_PDF_convolution(self, pdf_accessor, all_pdf, all_mu_f_squared):
        """ Convolute the weights_per_flavor_configurations."""

        # In order for the + distributions of the PDF counterterms and integrated
        # collinear ISR counterterms to hit the PDF only (and not the matrix elements or
        # observables functions), a change of variable is necessary: xb_1' = xb_1 * xi1
        # In this case, the xi1 to factor to apply is given by the Bjorken_x_rescalings
        # factor, which must be used for rescaling the argument of the PDF.
        # The 1/xi from the jacobian of this change of variable will typically be applied right after
        # this call to apply_PDF_convolution, because this must be done independently of
        # the user selection for the beam type (i.e. lpp in the run card).
# #gl
#         for flavors in self.weights_per_flavor_configurations:
#             PDFs = 1.
#             for i, flavor in enumerate(flavors[0]):
#                 if self.Bjorken_x_rescalings[0] != 1.0 and self.Bjorken_x_rescalings[1] != 1.0:
#                     Bjorken_x_rescalings_temp = (self.Bjorken_x_rescalings[0], 1.0)
#                     PDFs *= pdf_accessor(all_pdf[i], flavor,
#                          self.Bjorken_xs[i]*Bjorken_x_rescalings_temp[i], all_mu_f_squared[i])
#                 else:
#                     PDFs *= pdf_accessor(all_pdf[i], flavor,
#                          self.Bjorken_xs[i]*self.Bjorken_x_rescalings[i], all_mu_f_squared[i])
#             self.weights_per_flavor_configurations[flavors] *= PDFs

            #     if self.Bjorken_x_rescalings[0] != 1.0 and self.Bjorken_x_rescalings[1] != 1.0:
            #         Bjorken_x_rescalings_temp = (1.0, self.Bjorken_x_rescalings[1])
            #         PDFs *= pdf_accessor(all_pdf[i], flavor,
            #              self.Bjorken_xs[i]*Bjorken_x_rescalings_temp[i], all_mu_f_squared[i])
            #     else:
            #         PDFs *= pdf_accessor(all_pdf[i], flavor,
            #              self.Bjorken_xs[i]*self.Bjorken_x_rescalings[i], all_mu_f_squared[i])
            # self.weights_per_flavor_configurations[flavors] *= PDFs

        for flavors in self.weights_per_flavor_configurations:
            PDFs = 1.
            for i, flavor in enumerate(flavors[0]):
                PDFs *= pdf_accessor(all_pdf[i], flavor,
                     self.Bjorken_xs[i]*self.Bjorken_x_rescalings[i], all_mu_f_squared[i])
            self.weights_per_flavor_configurations[flavors] *= PDFs

    def store_information(self, key, value):
        """ This function registers a new additional information that will be attached to this
        event."""
        self.additional_information[key] = value

    def get_information(self, key):
        """ Tries to fetch a particular additional information registered in this event.
        If not present, it returns None."""
        return self.additional_information.get(key, None)

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
        
        return summed_weight

    def generate_mirrored_event(self):
        """ Creates a mirror event of self if necessary."""
        if not self.requires_mirroring or self.is_empty():
            return None
        

        # Now that the mirroring was applied to this event, it no longer needs any
        self.requires_mirroring = False
        # Neither does the mirrored event we instantiate below
        mirrored_event = self.get_copy()
        # Set the mirrored event flag "is_a_mirrored_event" appropriately
        mirrored_event.is_a_mirrored_event = True

        # First swap initial state flavors
        mirrored_event.weights_per_flavor_configurations = {}
        for flavor_config, wgt in self.weights_per_flavor_configurations.items():
            swapped_flavor_config = ( (flavor_config[0][1],flavor_config[0][0]),
                                           flavor_config[1] )
            mirrored_event.weights_per_flavor_configurations[swapped_flavor_config] = wgt
        
        # Then swap the kinematic configurations:
        # a) Initial-state momenta must be swapped
        # b) Final-state momenta must have their z-component swapped (for this to be correct
        # however, first make sure that initial momenta are indeed on the z-axis in the c.o.m frame)
        # Normally all events generated should be in the c.o.m frame already, but it may be
        # that in test_IR_limits, for debugging purposes only, the option 'boost_back_to_com' was
        # set to False, in which case we must here first *temporarily boost* it back to the c.o.m.
        if __debug__:
            # Two initial states
            assert(len(self.weights_per_flavor_configurations.keys()[0][0])==2)
            sqrts = math.sqrt((self.PS_point[0]+self.PS_point[1]).square())
            # Assert initial states along the z axis
            assert(abs(self.PS_point[0][1]/sqrts)<1.0e-9)
            assert(abs(self.PS_point[1][1]/sqrts)<1.0e-9)
            assert(abs(self.PS_point[0][2]/sqrts)<1.0e-9)
            assert(abs(self.PS_point[1][2]/sqrts)<1.0e-9)
        # If initial states are back to back we can directly proceed with a simple swap of the
        # z-axis, otherwise we must first boost to the c.o.m
        PS_point_to_swap = self.PS_point.get_copy()
        boost_vector = None
        if abs((self.PS_point[0]+self.PS_point[1])[3]/sqrts)>1.0e-9:
            boost_vector = (PS_point_to_swap[0]+PS_point_to_swap[1]).boostVector()
            for vector in PS_point_to_swap:
                vector.boost(-boost_vector)
        # debug: now make sure the event is back to back
        if __debug__:
            sqrts = math.sqrt((self.PS_point[0]+self.PS_point[1]).square())
            assert(abs((PS_point_to_swap[0]+PS_point_to_swap[1])[3]/sqrts)<1.0e-9)
        mirrored_event.PS_point = vectors.LorentzVectorList([ PS_point_to_swap[1],PS_point_to_swap[0] ])
        for vector in PS_point_to_swap[2:]:
            mirrored_event.PS_point.append(vectors.LorentzVector(
                                               [vector[0],vector[1],vector[2],-vector[3]]))
        # And now if we had boosted the event we must now boost it back
        if boost_vector is not None:
            for vector in PS_point_to_swap:
                vector.boost(boost_vector)
        
        # Then swap Bjorken x's and rescaling.
        mirrored_event.Bjorken_x_rescalings = (self.Bjorken_x_rescalings[1], self.Bjorken_x_rescalings[0])
        mirrored_event.Bjorken_xs = (self.Bjorken_xs[1], self.Bjorken_xs[0])
        
        return mirrored_event

    def filter_flavor_configurations(self, flavor_cut, **opts):
        """ Apply the flavor cut on this event for each flavor configuration, removing
        all those failings. If there is None left, this returns False."""
        
        new_weights_per_flavor_configurations = {}
        for flavors, weight in self.weights_per_flavor_configurations.items():
            if flavor_cut(self.PS_point, flavors, xb_1=self.Bjorken_xs[0],
                                                        xb_2=self.Bjorken_xs[1], **opts):
                new_weights_per_flavor_configurations[flavors] = weight
        
        self.weights_per_flavor_configurations = new_weights_per_flavor_configurations
        
        return (not self.is_empty())

    def apply_flavor_blind_cuts(self, flavor_blind_cut, *args, **opts):
        """ Apply the flavor-blind cut to this event, returning False if it failed."""
        
        return flavor_blind_cut(self.PS_point, *args, xb_1=self.Bjorken_xs[0],
                                                        xb_2=self.Bjorken_xs[1], **opts)

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

        if not isinstance(term_name, int):
            is_finite_specified = (term_name.lower()=='finite')
        for flavors, weight in self.weights_per_flavor_configurations.items():
            if isinstance(weight, float):
                if not is_finite_specified:
                    del self.weights_per_flavor_configurations[flavors]
            else:
                self.weights_per_flavor_configurations[flavors] = weight.get_term(term_name)

    def convolve_initial_states(self, flavor_matrix):
        """ Convolves the flavor matrix in argument (of the form:
            {   start_flavor_a :  {  end_flavors_1 : wgt_x,
                                     end_flavors_2 : wgt_y,
                                     [...]
                                 },
                start_flavor_b :  {  end_flavors_1 : wgt_z,
                                     end_flavors_2 : wgt_w,
                                     [...]
                                 },
                [...]
            }
            where 'wgt_x' can be EpsilonExpansion instances.
            and modifies self.weights_per_flavor_configurations with the convolved result.

            Note that end_flavors is a tuple of flavor with identical matrix elements (such as in q -> q g, which is flavor independent)

            Finally each element of these flavor configurations is a 2-tuple indicating the flavor of each of the two beams,
            however the flavor can also remain unspecified and set to None.

            Example:
            {   (21, None) :     {  ( (2, None),(-2, None),(1, None),(-1, None) ): wgt_x,
                                    ( (21, None),) : wgt_y,
                                    [...]
                                 },
                (21,2)     :     {  ( (2, 21),(-2, 21),(1, 21),(-1, 21) ) : wgt_z,
                                     [...]
                                 },
                [...]
            }
        """

        def substitute_initial_state_flavors(flavor_config_to_substitute, new_flavors):
            """ Substitute in the flavor config the convoluted flavor with the one in
            argument and returns the corresponding new flavor configuration as a two tuples,
            for the initial and final states respectively."""
            return (
                tuple([
                    flavor_config_to_substitute[0][0] if new_flavors[0] is None else new_flavors[0],
                    flavor_config_to_substitute[0][1] if new_flavors[1] is None else new_flavors[1],
                ]),
                flavor_config_to_substitute[1]
            )

        if len(flavor_matrix) == 0:
            # shortcut in case an empty convolutional flavor matrix is specified which can happen
            # if the beam currents are zero but don't have the attribute is_zero = True
            self.weights_per_flavor_configurations = {}
            return

        all_new_flavor_configurations = {}
        for flavor_config, wgt in self.weights_per_flavor_configurations.items():
            start_flavors = flavor_config[0]
            for matrix_start_flavors, matrix_all_end_flavors in flavor_matrix.items():
                if (matrix_start_flavors[0] in (None, start_flavors[0])) and (matrix_start_flavors[1] in (None, start_flavors[1])):
                    new_flavor_configs = []
                    for end_flavors_groups, matrix_wgt in matrix_all_end_flavors.items():
                        new_flavor_configs.extend([
                            ( substitute_initial_state_flavors( flavor_config, (end_flavor_beam_one, end_flavor_beam_two) ),
                              base_objects.EpsilonExpansion(matrix_wgt) * wgt )
                            for (end_flavor_beam_one, end_flavor_beam_two) in end_flavors_groups ])

                    for new_flavor_config, new_wgt in new_flavor_configs:
                        try:
                            all_new_flavor_configurations[new_flavor_config] += new_wgt
                        except KeyError:
                            all_new_flavor_configurations[new_flavor_config] = new_wgt

        # Now assign the newly create flavor configurations
        self.weights_per_flavor_configurations = all_new_flavor_configurations

    def convolve_flavors(self, flavor_matrix, leg_index, initial_state=True):
        """ Convolves the flavor matrix in argument (of the form:
            {   start_flavor_a :  {  end_flavors_1 : wgt_x,
                                     end_flavors_2 : wgt_y,
                                     [...]
                                 },
                start_flavor_b :  {  end_flavors_1 : wgt_z,
                                     end_flavors_2 : wgt_w,
                                     [...]
                                 },
                [...]
            }
            where 'wgt_x' can be EpsilonExpansion instances.
            and modifies self.weights_per_flavor_configurations with the convolved result.
            leg index specifies which leg must be convolved.
               (e.g. beam #1 correpsonds to leg_index=0 and initial_state = True.
            Note that end_flavors is a tuple of flavor with identical matrix elements (such as in q -> q g, which is flavor independent)

        """
            
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
            
        if len(flavor_matrix)==0:
            # shortcut in case an empty convolutional flavor matrix is specified which can happen
            # if the beam currents are zero but don't have the attribute is_zero = True
            self.weights_per_flavor_configurations = {}
            return

        new_flavor_configurations = {}
        for flavor_config, wgt in self.weights_per_flavor_configurations.items():
            start_flavor = flavor_config[0][leg_index] if initial_state else flavor_config[1][leg_index]
            if start_flavor in flavor_matrix:
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

    def is_empty(self):
        """
        Check if the weight vector (weights_per_flavor_configurations) is an empty dictionary.
        For now does not check if all weights are zero
        :return: bool
        """
        return len(self.weights_per_flavor_configurations) == 0

    def __add__(self, other):
        """ overload the '+' operator."""
        new_event = self.get_copy()
        return new_event.__iadd__(other)
            
    def __iadd__(self, other):
        """ overload the '+=' operator."""
        if __debug__: assert(self.requires_mirroring == other.requires_mirroring)
        if __debug__: assert(self.PS_point == other.PS_point)
        if __debug__: assert(self.Bjorken_x_rescalings  == other.Bjorken_x_rescalings)
        if __debug__: assert(self.Bjorken_xs  == other.Bjorken_xs)
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
            res.append("%s%svent from %s contribution '%s'%s%s:"%(
                GREEN,
                'E' if not self.is_a_mirrored_event else 'Mirrored e',
                self.host_contribution_definition.correction_order,
                self.host_contribution_definition.short_name(),
                ' (requires mirroring)' if self.requires_mirroring else '',
                ENDC))
        else:
            if isinstance(self.counterterm_structure, list):
                all_ct_strucs = self.counterterm_structure
            else:
                all_ct_strucs = [self.counterterm_structure,]
            
            counterterm_structure_str_elems = []
            for ct_struct, all_resolved_PDGs, kinematics_identifier in all_ct_strucs:
                # We show only one representative resolved_PDGs structure: the first one.
                if isinstance(all_resolved_PDGs, dict):
                    representative_resolved_PDGs = all_resolved_PDGs.keys()[0]
                else:
                    representative_resolved_PDGs = all_resolved_PDGs[0]
                counterterm_structure_str_elems.append('%s%s@%s'%(str(ct_struct),
                    '' if kinematics_identifier is None else '|%s'%kinematics_identifier,
                    str(representative_resolved_PDGs).replace(' ','')))
            
            if len(all_ct_strucs)>1:
                counterterm_structure_str = '( %s )'%(' + '.join(counterterm_structure_str_elems))
            else:
                counterterm_structure_str = counterterm_structure_str_elems[0]

            res.append('%s%sounterterm event from limit%s %s of %s contribution %s%s'%(
                GREEN,
                'C' if not self.is_a_mirrored_event else 'Mirrored c',
                's' if isinstance(self.counterterm_structure, list) and len(self.counterterm_structure)>1 else '',
                counterterm_structure_str,
                self.host_contribution_definition.correction_order,
                self.host_contribution_definition.short_name(),
                ENDC))
            
        res.append('%s  Kinematic configuration:%s'%(BLUE, ENDC))
        res.extend('     %s'%line for line in str(self.PS_point).split('\n'))
        if not any(bs is None for bs in self.Bjorken_x_rescalings):
            res.append("     Bjorken x's scalings: %s"%(', '.join('%.16e'%bs for bs in self.Bjorken_x_rescalings)))
        if not any(bs is None for bs in self.Bjorken_xs):
            res.append("     Bjorken x's         : %s"%(', '.join('%.16e'%bs for bs in self.Bjorken_xs)))
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
        
        if len(self.additional_information)>0:
            res.append('%s  Additional information:%s'%(BLUE,  ENDC))
            for key, value in self.additional_information.items():
                res.append('    %s : %s'%(str(key), str(value)))

        return '\n'.join(res)

    def counterterm_structure_short_string(self):
        """ Return a short string specifying the counterterm structure(s) of this event."""

        if self.counterterm_structure is None:
            if not self.is_a_mirrored_event:
                return 'Event'
            else:
                return '<->Event'

        event_str = str(self.counterterm_structure)
        if isinstance(self.counterterm_structure, list):
            event_ct_structs = self.counterterm_structure
        else:
            event_ct_structs = [self.counterterm_structure,]
        
        event_str_elems = []
        for ct_struct, all_resolved_PDGs, kinematics_identifier in event_ct_structs:
            # We show only one representative resolved_PDGs structure: the first one.
            if isinstance(all_resolved_PDGs, dict):
                representative_resolved_PDGs = all_resolved_PDGs.keys()[0]
            else:
                representative_resolved_PDGs = all_resolved_PDGs[0]
            str_ct_struct = str(ct_struct)
            # Clean-up of the embedding overall parenthesis for the short string
            while str_ct_struct.startswith("(") and str_ct_struct.endswith(",)"):
                str_ct_struct = str_ct_struct[1:]
                str_ct_struct = str_ct_struct[:-2]
            event_str_elems.append('%s%s@%s'%(str_ct_struct,
                '' if kinematics_identifier is None else '|%s' % kinematics_identifier,
                str(representative_resolved_PDGs).replace(' ','')))
        
        if len(event_ct_structs)>1:
            event_str = '( %s )'%(' + '.join(event_str_elems))
        else:
            event_str = event_str_elems[0]
        input_mapping = self.get_information('input_mapping')
        if input_mapping is not None:
             event_str += 'x{%s}'%(','.join('%d:%d'%(k,input_mapping[k]) for k in
                                                        sorted(input_mapping.keys()) ))
        if self.is_a_mirrored_event:
            event_str = '<->%s'%event_str
        beam_convolution_masks = self.get_information('beam_convolution_masks')
        if beam_convolution_masks is not None:
            for (beam_name, beam_short_name) in [('beam_one','F1'),('beam_two','F2')]:
                if beam_convolution_masks[beam_name] != 'ALL':
                    event_str += '@%s<-(%s)'%(beam_short_name,'|'.join('%d'%pdg for pdg in
                                                    beam_convolution_masks[beam_name]))
        return event_str

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
    
    def apply_observables(self, observable_list, integrator_weight=1.):
        """ Applies the observable list to the entire set of events of this list."""

        for event in self:
            observable_list.apply_observables(event, integrator_weight)

    def generate_mirrored_events(self):
        """ Creates the mirror configurations if necessary."""
        for event in list(self):
            mirrored_event = event.generate_mirrored_event()
            if mirrored_event:
                self.append(mirrored_event)

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
    
    FROZEN_DIMENSIONS = {}#{'x_1': 0.13, 'x_2': 0.26, 'x_3': 0.35}

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
                if contribution_definition.correlated_beam_convolution:
                    # print('ME7Integrands.py - correlated_beam_convolution : ' + str(contribution_definition.correlated_beam_convolution))
                    # print('ME7Integrands.py - torino_sub_BS : ' + str(contribution_definition.torino_sub_BS))
                    target_type = 'BeamSoft'
                else:
                    target_type = 'RealVirtual'
            elif n_loops == 0 and n_unresolved_particles == 0:
                target_type = 'Born'
            elif n_loops == 1 and n_unresolved_particles == 0:
                target_type = 'Virtual'
            elif n_loops == 2 and n_unresolved_particles == 0:
                target_type = 'DoubleVirtual'
            elif n_loops == 0 and n_unresolved_particles == 1:
                target_type = 'SingleReals'
            elif n_loops == 0 and n_unresolved_particles == 2:
                target_type = 'DoubleReals'
            elif n_loops == 0 and n_unresolved_particles == 3:
                target_type = 'TripleReals'
            elif n_loops > 0 and n_unresolved_particles > 0:
                target_type = 'RealVirtual'
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
                       ME7_configuration,
                       sectors = None,
        ):
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
                                       'options'                    : {
                                           'sectors': sectors
                                        }
                                    }
        
        # The original ContributionDefinition instance at the origin this integrand
        self.contribution_definition    = contribution_definition

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
        # When accessing currents we must know whether leg numbers must be tracked.
        # By default, this is not the case and in higher order ME7Integrands, this will depend on the
        # particular subtraction scheme loaded.
        self.accessor_tracks_leg_numbers = False

        # Save identified sectors if any
        self.sectors = sectors
        # For now specify a specific sector by hardcoding a selector function
        # None is equivalent to `lambda defining_process, sector: True`
        self.is_sector_selected = None

        # Update and define many properties of self based on the provided run-card and model.
        self.synchronize(model, run_card, ME7_configuration)

        # Initialize the call counter
        # This counter is incremented for each time self.__call__ is called and reinitialized in self.synchronize
        self.n_observable_calls = 0

        # Store a flavor cut definition, which the ME7 interface can overwrite. This will be used
        # by the function pass_flavor_sensitive cut and it is always None by default.
        self.flavor_cut_function = None

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

    def get_nice_string_sector_lines(self, process_key, format=0):
        """ Returns a nicely formatted summary of the sector information for this integrand."""

        if self.sectors is None or self.sectors[process_key] is None:
            return []

        all_sectors_info = self.sectors[process_key]

        res = []
        if format <2:
            res.append('   | %d sectors: %s'%(len(all_sectors_info),
                                            ' | '.join('%s'%str(s['sector']) for s in all_sectors_info)))
        else:
            res.append('   | with %d sectors:'%len(all_sectors_info))
            for sector_info in all_sectors_info:
                line_elems = [ '   | > sector %s'%str(sector_info['sector']) ]
                if sector_info['counterterms'] is not None:
                    line_elems.append('local counterterms [%s]'%(','.join('%d'%i_ct for i_ct in sector_info['counterterms'])))
                else:
                    if self.has_local_counterterms():
                        line_elems.append('all local counterterms')
                if sector_info['integrated_counterterms'] is not None:
                    i_cts= sector_info['integrated_counterterms']
                    line_elems.append('integrated counterterms [%s].'%(','.join(
                        '%d@%s'%(i_ct, 'ALL' if i_cts[i_ct] is None else '{%s}'%(','.join(
                            '%d'%i_mapping for i_mapping in i_cts[i_ct])) ) for i_ct in sorted(i_cts.keys()))))
                else:
                    if self.has_integrated_counterterms():
                        line_elems.append('all integrated counterterms')

                if len(line_elems)==1:
                    res.append('%s'%line_elems[0])
                elif len(line_elems)==2:
                    res.append('%s with %s'%tuple(line_elems))
                else:
                    res.append('%s with %s and %s'%tuple(line_elems))

        return res

    def get_short_name(self):
        """ Returns the short-name for this integrand, typically extracted from the one
        of the contribution definition."""
        return self.contribution_definition.short_name()

    def nice_string(self, format=0):
        """ Nice string representation of self."""
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        res = ['< %s%s%s >'%(BLUE,self.get_short_name(),ENDC)]
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
                res.extend(self.get_nice_string_sector_lines(process_key, format=format))
                for mapped_process in mapped_processes:
                    res.append(BLUE+u'   \u21b3  '+mapped_process.nice_string(print_weighted=False)\
                                                                        .replace('Process: ','')+ENDC)

        return '\n'.join(res).encode('utf-8')

    def __str__(self):
        return self.nice_string()

    @staticmethod
    def build_flavor_cut_function(flavor_cut_string, contrib_name):
        """ Given a string defining a flavor cut function specification (taken from the run
        card), build here the flavor cut function to apply at run-time. Example of a complicated
        cut function:
          
          "{
            'R':'* 21 > (21,2)x3 * | * 21 > (21,3)x3 *',
            'V':'21 (-2,2) > 21x2 (range(1,7)+range(-1,-7,-1))x5
           }"

          Means that the real must have:
            a) anything as first incoming particle and a gluon the second incoming particle
            b) Three 'gluons or up quarks' in the final state  + anything
               *or* Three 'gluons or strange quarks' in the final state  + anything
          
          Means that the virtual must have:
            a) a gluon and a 'up quark or anti-up quark' as first and second particle
            b) 2 gluons and 5 quarks or antiquarks *exactly* (all particles must fit in one category at least)
            
          Notice that one can also speficy just a string, in which case the cut will apply
          to all contributions.
    
        """

        # First check if the format is a dictionary that specifies a different function
        # for different integrand or if it applies to all integrands.
        try:
            flavor_func_dic = eval(flavor_cut_string)
            is_a_dict = isinstance(flavor_func_dic, dict)
        except:
            is_a_dict = False
        if is_a_dict:
            if contrib_name in flavor_func_dic:
                flavor_cut_string = flavor_func_dic[contrib_name]
            else:
                return None
        
        if flavor_cut_string.lower() == 'none':
            return None
        if flavor_cut_string.lower() == 'hardcoded':
            return 'hardcoded'
        
        patterns = [p.strip() for p in flavor_cut_string.split('|')]
        processed_patterns = []
        for pattern in patterns:
            processed_pattern = {}
            try:
                initial_state_pattern, final_state_pattern = pattern.split('>')
                initial_state_pattern = initial_state_pattern.strip()
                final_state_pattern = final_state_pattern.strip()
            except:
                raise MadEvent7Error('Malformed flavor-cut pattern: %s.'%pattern)
            initial_state_patterns = [isp.strip() for isp in re.split(r'\s+',initial_state_pattern)]
            processed_pattern['initial_states'] = []
            for isp in initial_state_patterns:
                if isp=='*':
                    processed_pattern['initial_states'].append(None)
                else:
                    try:
                        specifier = eval(isp)
                    except:
                        raise MadEvent7Error('Malformed flavor-cut pattern: %s.'%isp)
                    if isinstance(specifier,int):
                        processed_pattern['initial_states'].append((specifier,))
                    elif isinstance(specifier, (tuple, list)):
                        processed_pattern['initial_states'].append(tuple(specifier))
                    else:
                        raise MadEvent7Error('Malformed flavor-cut pattern: %s.'%isp)
            processed_pattern['final_states'] = []
            processed_pattern['needs_exact_final_states'] = True
            final_state_patterns = [fsp.strip() for fsp in re.split(r'\s+',final_state_pattern)]
            for fsp in final_state_patterns:
                try:
                    matcher, multiplier = fsp.split('x')
                except:
                    multiplier = '1'
                    matcher = fsp
                try:
                    multiplier = int(multiplier.strip())
                except:
                    raise MadEvent7Error('Malformed flavor-cut pattern: %s.'%multiplier)
                matcher = matcher.strip()
                if matcher=='*':
                    processed_pattern['needs_exact_final_states'] = False
                else:
                    try:
                        specifier = eval(matcher)
                    except:
                        raise MadEvent7Error('Malformed flavor-cut pattern: %s.'%matcher)
                    if isinstance(specifier,int):
                        processed_pattern['final_states'].append(((specifier,),multiplier))
                    elif isinstance(specifier, (tuple, list)):
                        processed_pattern['initial_states'].append((tuple(specifier),multiplier))
                    else:
                        raise MadEvent7Error('Malformed flavor-cut pattern: %s.'%matcher)
            processed_patterns.append(processed_pattern)
            
        # Ok, now that we've processed the patterns to consider, we can create the function
        # that checks it:
        def pass_flavor_cut_func(flavor_config):
            initial_pdgs, final_pdgs = flavor_config[0], flavor_config[1]
            # Loop over all capturing patterns
            for pattern in processed_patterns:
                # First check the initial states
                if len(initial_pdgs)!=len(pattern['initial_states']):
                    continue
                must_continue = False
                for i, pdg in enumerate(initial_pdgs):
                    if (pattern['initial_states'][i] is not None) and \
                                    (pdg not in pattern['initial_states'][i]):
                        must_continue = True
                        break
                if must_continue:
                    continue
                
                # Then check the final states.
                # First make sure all constraints specified by the user apply
                must_continue = False
                for (matcher, multiplier) in pattern['final_states']:
                    if len([1 for pdg in final_pdgs if pdg in matcher])!=multiplier:
                        must_continue = True
                        break
                if must_continue:
                    continue
                # And if 'needs_exact_final_states' then also make sure that all final
                # states belong in at least one category
                if pattern['needs_exact_final_states']:
                    must_continue = False
                    for pdg in final_pdgs:
                        if not any(pdg in matcher for (matcher, multiplier) in pattern['final_states']):
                            must_continue = True
                            break
                    if must_continue:
                        continue
                
                # If we haven't "continued" by now, then it means that the pattern was
                # capturing and we can return True
                return True
                
            # No pattern matched upt to this point so we must return False
            return False
            
        # We can now return the flavor cut function created
        return pass_flavor_cut_func
                
    def synchronize(self, model, run_card, ME7_configuration):
        """ Synchronize this integrand with the most recent run_card and model."""

        # The option dictionary of ME7
        self.ME7_configuration          = ME7_configuration

        # Load or access the already loaded subtraction scheme module
        if self.contribution_definition.overall_correction_order.count('N')>0:
            self.subtraction_scheme = subtraction.SubtractionCurrentExporter.get_subtraction_scheme_module(
                            self.ME7_configuration['subtraction_scheme'], root_path = self.ME7_configuration['me_dir'])
            self.accessor_tracks_leg_numbers = self.subtraction_scheme.are_current_instances_for_specific_leg_numbers
        else:
            self.subtraction_scheme = None
            self.accessor_tracks_leg_numbers = False

        
        # A ModelReader instance, initialized with the values of the param_card.dat of this run
        self.model                      = model
        if not isinstance(self.model, model_reader.ModelReader):
            raise MadGraph5Error("The ME7Integrand must be initialized with a ModelReader instance.")

        # A RunCardME7 instance, properly initialized with the values of the run_card.dat of this run
        self.run_card                   = run_card
        
        self.flavor_cut_function = ME7Integrand.build_flavor_cut_function(
                                                self.run_card['flavor_cuts'], self.get_short_name())
        
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
                          self.contribution_definition.is_beam_active('beam_two') ),
            correlated_beam_convolution = self.contribution_definition.correlated_beam_convolution,
            torino_sub_BS = self.contribution_definition.torino_sub_BS  #gl
        )

        # Add a copy of the PS generator dimensions here.
        # Notice however that we could add more dimensions pertaining to this integrand only, and PS generation.
        # This is in particular true for discrete integration dimension like sectors, helicities, etc...
        integrand_dimensions = integrands.DimensionList(self.phase_space_generator.dimensions)

        # Now remove the frozen dimensions from the list of continuous ones.
        integrand_dimensions = integrands.DimensionList(
                            d for d in integrand_dimensions if d.name not in self.FROZEN_DIMENSIONS)
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

        #Import the observables from the FO_analysis folder
        if run_card['FO_analysis'].upper() == "NONE":
            self.apply_observables = False
        else:
            self.apply_observable = True

            #The FO_analysis parameter points to a python module located in
            #me_dir/FO_analysis. It must contain observable_list, defined as a
            #madgraph.integrators.observables.HwUObservableList

            sys.path.append(self.ME7_configuration['me_dir'])   # Load the process dir in the path
            try:
                analysis_module = __import__('FO_analysis.'+run_card['FO_analysis'], fromlist=[''])
                reload(analysis_module)    # Make sure to reinitialize the plots (empty the HwUs)
            except ImportError:
                sys.path.pop()  # Clean the path after importing
                raise ImportError("The FO_analysis specified in the run_card could not be imported")
            try:
                self.observable_list = analysis_module.observable_list
            except NameError:
                sys.path.pop()  # Clean the path after importing
                raise NameError("Failed to access specified FO_analysis.observable_list")
            sys.path.pop()  # Clean the path after importing
        self.n_observable_calls = 0    # Re-initialize the counting tools



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
        
        def pdg_list_match(target_list, selection_list, exact=False):
            if len(target_list) != len(selection_list):
                return False
            if exact:
                for target_pdg, selected_pdgs in zip(target_list, selection_list):
                    if target_pdg not in selected_pdgs:
                        return False
                return True
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

            # Ask for exact match for the initial states
            if (not selection['in_pdgs'] is None) and \
               (not pdg_list_match(process.get_initial_ids(), selection['in_pdgs'], exact=True)):
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
        it can be useful for both the ME7Integrnd_V and ME7_integrand_R.
        Note here that 'counterterms' is not a simple list of counterterms but rather a list
        of two-tuples (ID, ct) where ID is the "identifier" of the CT (i.e. its position in
        the list of counterterms for that particular process.
        """

        # First select only the counterterms which are not pure matrix elements
        # (i.e. they have singular structures) and also exclude here soft-collinear
        # counterterms since they do not have an approach limit function.
        if not integrated_CT:
            selected_counterterms = [ ct for ct in counterterms if ct[1].is_singular() ]
        else:
            selected_counterterms = [ ct for ct in counterterms if ct[1]['integrated_counterterm'].is_singular() ]

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
                        limit_pattern = '((%s,),)' % limit_pattern
                    elif integrated_CT and not limit_pattern.startswith('['):
                        # If not specified as a raw string, we take the liberty of adding
                        # the enclosing parenthesis.
                        limit_pattern = '[(%s,),]' % limit_pattern
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
                ct = counterterm[1] if not integrated_CT else counterterm[1]['integrated_counterterm']
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

        # The PS point in input is sometimes provided as a dictionary or a flat list, but
        # we need it as a flat list here, so we force the conversion
        if isinstance(PS_point, LorentzVectorDict):
            PS_point = PS_point.to_list()
        elif not isinstance(PS_point, LorentzVectorList):
            PS_point = LorentzVectorList(LorentzVector(v) for v in PS_point)

        # If the cuts depend on the boost to the lab frame in case of hadronic collision
        # then the boost below can be used. Notice that the jet cuts is typically formulated in terms of
        # *pseudo* rapidities and not rapidities, so that when jets are clustered into a massive combined
        # momentum, this makes it important to boost to the lab frame. We therefore turn on this boost here
        # by default.
        if (xb_1 is not None) and (xb_2 is not None) and (xb_1, xb_2) != (1., 1.):
            # Prevent border effects
            PS_point = PS_point.get_copy()
            PS_point.boost_from_com_to_lab_frame(xb_1, xb_2, self.run_card['ebeam1'], self.run_card['ebeam2'])

        for i, p in enumerate(PS_point[self.n_initial:]):
            if is_a_photon(process_pdgs[1][i]):
                if p.pt() < 100.0:
                    return False

        for i, p in enumerate(PS_point[self.n_initial:]):
            if is_a_photon(process_pdgs[1][i]):
                for j, p2 in enumerate(PS_point[self.n_initial:]):
                    if (j != i) and process_pdgs[1][j]!=21:
                        if p.deltaR(p2) < 0.4:
                            return False

        ###################################################################################
        # JET CLUSTERING AND CUTS
        ###################################################################################
        
        ptj_cut = self.run_card['ptj']
        drjj_cut = self.run_card['drjj']
        etaj_cut = self.run_card['etaj']

        # #gl
        # p = PS_point[2]
        # #print('ME7 - p : ' + str(p))
        # if p[1]**2 + p[2]**2 < 10**2:
        #     return False
        
        if ptj_cut <= 0. and    \
           drjj_cut <= 0. and   \
           etaj_cut <= 0. :
            return True
        else:
            # If fastjet is needed but not found, make sure to stop
            if (not PYJET_AVAILABLE) and n_jets_allowed_to_be_clustered>0:
                raise MadEvent7Error("Fast-jet python bindings are necessary for integrating"+
                             " real-emission type of contributions. Please install pyjet.")

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
            sequence = pyjet.cluster(this_event, R=drjj_cut, p=-1, ep=True, algo='genkt')
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
        if self.flavor_cut_function is None:
            return True
        elif callable(self.flavor_cut_function):
            return self.flavor_cut_function(flavors)
        elif self.flavor_cut_function.lower() != 'hardcoded':
            raise MadEvent7Error("Flavor cut function '%s' not reckognized."%self.flavor_cut_function)

        debug_cuts = False
        from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList, LorentzVector

        # The PS point in input is sometimes provided as a dictionary or a flat list, but
        # we need it as a flat list here, so we force the conversion
        if isinstance(PS_point, LorentzVectorDict):
            PS_point = PS_point.to_list()
        elif not isinstance(PS_point, LorentzVectorList):
            PS_point = LorentzVectorList(LorentzVector(v) for v in PS_point)

        # If the cuts depend on the boost to the lab frame in case of hadronic collision
        # then the boost below can be used. Notice that the jet cuts is typically formulated in terms of
        # *pseudo* rapidities and not rapidities, so that when jets are clustered into a massive combined
        # momentum, this makes it important to boost to the lab frame. We therefore turn on this boost here
        # by default.
        if (xb_1 is not None) and (xb_2 is not None) and (xb_1, xb_2) != (1., 1.):
            # Prevent border effects
            PS_point = PS_point.get_copy()
            PS_point.boost_from_com_to_lab_frame(xb_1, xb_2, self.run_card['ebeam1'], self.run_card['ebeam2'])

        if debug_cuts: logger.debug( "Processing flavor-sensitive cuts for flavors %s and PS point:\n%s"%(
            str(flavors), LorentzVectorList(PS_point).__str__(n_initial=self.phase_space_generator.n_initial) ))
        
        ###########################################
        # User can define his own flavor cut below
        ###########################################
        
        #--> USER-DEFINED CODE

        ###########################################
        # User can define his own flavor cut above
        ###########################################
        
        # Apply the flavor_cut_definition
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

        # For testing, one can easily substitute PDFs with a 1./x distribution by uncommenting the
        # line below
        #return 1./x

        if pdf is None:
            return 1.

        if pdg not in [21,22] and abs(pdg) not in range(1,7):
            return 1.
                
        if (pdf, pdg, x, scale2) in self.PDF_cache:
            return self.PDF_cache[(pdf, pdg, x, scale2)]

        # Call to lhapdf API
        f = pdf.xfxQ2(pdg, x, scale2)/x
#gl
        # if pdg == 21:
        #     f *=10000.

        # Update the PDF cache
        self.PDF_cache[(pdf, pdg,x,scale2)] = f
        self.PDF_cache_entries.append((pdf, pdg,x,scale2))
        if len(self.PDF_cache_entries) > self.PDF_cache_max_size:
            del self.PDF_cache[self.PDF_cache_entries.pop(0)]

        return f

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Main function of the integrand, returning the weight to be passed to the integrator."""
        
        # Check if an integrator_jacobian is specified in the options to be applied to
        # the weight when registering observables (i.e. filling observables)
        if 'integrator_jacobian' in opts:
            integrator_jacobian = opts['integrator_jacobian']
        else:
            integrator_jacobian = 1.0
        
        # The integrator may have passed an observable thread lock to make sure to avoid
        # the concurrency issues during parallel runs
        if 'observables_lock' in opts:
            observables_lock = opts['observables_lock']
        else:
            observables_lock = misc.dummy_lock()
        
        # Increment the number of calls in a thread-safe manner
        with observables_lock:
            if self.apply_observables:
                self.n_observable_calls += 1

        if __debug__: logger.debug("="*80)
        if __debug__: logger.debug('Starting a new evaluation of the integrand from contribution:\n%s',
                                                    self.contribution_definition.nice_string())
        
        # Random variables sent
        random_variables    = list(continuous_inputs)
        if __debug__: logger.debug('Random variables received: %s',str(random_variables))
    
        # Now assign the variables pertaining to PS generations
        PS_random_variables = [
                ( self.FROZEN_DIMENSIONS[name] if name in self.FROZEN_DIMENSIONS else
                  random_variables[self.dim_name_to_position[name]] )
                                          for name in self.phase_space_generator.dim_ordered_names ]

        if self.sectors is None:
            weight = self.evaluate(PS_random_variables, integrator_jacobian, observables_lock,
                                                                            selected_process_key=None, sector_info=None)
            return weight
        else:
            total_weight = 0.
            for process_key, (process, mapped_processes) in self.processes_map.items():
                if self.sectors is None or self.sectors[process_key] is None:
                    all_sectors = [None, ]
                else:
                    all_sectors = self.sectors[process_key]
                for sector_info in all_sectors:
                    if (sector_info is not None) and (self.is_sector_selected is not None) and \
                                                    not self.is_sector_selected(process, sector_info['sector']):
                        continue
                    total_weight += self.evaluate(PS_random_variables, integrator_jacobian, observables_lock,
                                                  selected_process_key=process_key, sector_info=sector_info)
            return total_weight


    def evaluate(self,PS_random_variables, integrator_jacobian, observables_lock, selected_process_key=None, sector_info=None):
        """ Evaluate this integrand given the PS generating random variables and
        possibily for a particular defining process and sector. If no process_key is specified, this funcion will
        loop over and aggregate all of them"""

        # A unique float must be returned
        wgt = 1.0
        # And the conversion from GeV^-2 to picobarns
        wgt *= 0.389379304e9

        # Note: if the process_key is defined as well as the sector info, one may consider using a dedicated
        # PS generator for those (i.e. g g > X is mapped with d d~ > X but the bjorken x's pre-sampling could be
        # improved and for the sector, the parametrisation of the unresolved degrees of freedom could be chosen
        # accordingly to the singularities present in this sector.)
        PS_point, PS_weight, x1s, x2s = self.phase_space_generator.get_PS_point(PS_random_variables)

        # Unpack the initial momenta rescalings (if present) so as to access both Bjorken
        # rescalings xb_<i> and the ISR factorization convolution rescalings xi<i>.
        xb_1, xi1 = x1s
        xb_2, xi2 = x2s

        if xb_1 > 1. or xb_2 > 1. or math.isnan(xb_1) or math.isnan(xb_2):
            raise MadEvent7Error('Unphysical configuration: x1, x2 = %.5e, %.5e'%(xb_1, xb_2))
            #logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
            #return 0.0
        
        if PS_point is None:
            if __debug__:
                if xb_1 > 1. or xb_2 > 1.:
                    logger.debug('Unphysical configuration: x1, x2 = %.5e, %.5e'%(xb_1, xb_2))
                else:
                    logger.debug('Phase-space generation failed.')
                #raise MadEvent7Error('Unphysical configuration encountered.')
                logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
            return 0.0
                
        if __debug__: logger.debug("Considering the following PS point:\n%s"%(PS_point.__str__(
                                            n_initial=self.phase_space_generator.n_initial) ))

        # Account for PS weight
        wgt *= PS_weight
        if __debug__: logger.debug("PS_weight: %.5e"%PS_weight)

        # The E_cm entering the flux factor is computed *without* including the xi<i> rescalings
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

            # Make sure to skip processes that should not be considered if the process_key is specified
            if selected_process_key is not None and selected_process_key != process_key:
                continue

            # If one wishes to integrate only one particular subprocess, it can be done by uncommenting
            # and modifying the lines below.
#            if process.get_cached_initial_final_pdgs() in [((2,-2),(23,1,-1)), ((2,-2),(23,1,-1))] :
#               continue
#            else:
#                misc.sprint("Now doing: %s"%str(process.get_cached_initial_final_pdgs()))
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
                wgt, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, sector_info=sector_info)
            
            if events is None or len(events)==0:
                continue
            # Now handle the process mirroring, by adding, for each ME7Event defined with
            # requires_mirroring = True, a new mirrored event with the initial states swapped,
            # as well as the Bjorken x's and rescalings and with p_z -> -p_z on all final
            # state momenta
            events.generate_mirrored_events()

            # misc.sprint(events)
            # Apply flavor blind cuts
            if not events.filter_with_flavor_blind_cuts(self.pass_flavor_blind_cuts, process_pdgs,
                n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles):
                if __debug__: logger.debug('All events failed the flavour_blind generation-level cuts.')
                if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
                continue

            # Select particular terms of the EpsilonExpansion terms stored as weigts
            events.select_epsilon_expansion_term('finite')

            # Apply PDF convolution
            if self.run_card['lpp1']==self.run_card['lpp2']==1:
                events.apply_PDF_convolution( self.get_pdfQ2,
                                               (self.pdf, self.pdf), (mu_f1**2, mu_f2**2) )
            # Apply the 1/xi**2 factor from the change of variable xi -> xi' / xi
            # Remember that the Bjorken rescalings are defined as 1/xi at this stage
# #gl
#             for event in events:
#                 if self.contribution_definition.torino_sub_BS:
#                     event *= event.Bjorken_x_rescalings[0]*1.0
#                 else:
#                     event *= event.Bjorken_x_rescalings[0]*event.Bjorken_x_rescalings[1]

            for event in events:
                event *= event.Bjorken_x_rescalings[0]*event.Bjorken_x_rescalings[1]

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
            if not events.filter_flavor_configurations(self.pass_flavor_sensitive_cuts):
                if __debug__: logger.debug('Events failed the flavour-sensitive generation-level cuts.')
                if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
                continue

            # Aggregate the total weight to be returned to the integrator
            total_wgt += events.get_total_weight()
            
            # Select a unique flavor configuration for each event
            events.select_a_flavor_configuration()

            #misc.sprint(events)
            if __debug__: logger.debug('Short-distance events for this subprocess:\n%s\n'%str(events))
            
            # Finally apply observables in a thread-safe manner
            with observables_lock:
                if self.apply_observables:
                    events.apply_observables(self.observable_list, integrator_jacobian)


        # Now finally return the total weight for this contribution
        if __debug__: logger.debug(misc.bcolors.GREEN + "Final weight returned: %.5e"%total_wgt + misc.bcolors.ENDC)
        if __debug__: logger.debug("="*80)

        return total_wgt

    def compute_matrix_element_event_weight(self,PS_point, process_key, process, all_flavor_configurations,
                  base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts):
        """ Returns the weight for the matrix element event, which may need further processing in daughter classes, for
          example in the virtual where one may use a dummy dressed with the I operator."""

        alpha_s = self.model.get('parameter_dict')['aS']

        ME_evaluation, all_results = self.all_MEAccessors(
                        process, PS_point, alpha_s, mu_r, pdgs=all_flavor_configurations[0])

        return base_objects.EpsilonExpansion(ME_evaluation) * base_weight


    def generate_matrix_element_event(self, PS_point, process_key, process, all_flavor_configurations,
                  base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts):
        """ Generates one event which corresponds to the physical contribution (i.e. not the
        counterterms) hosted by this contribution."""

        event_weight = self.compute_matrix_element_event_weight(
            PS_point, process_key, process, all_flavor_configurations,
                  base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts)

        sector_info = opts.get('sector_info', None)
        if sector_info is not None and sector_info['sector'] is not None:
            event_weight *= sector_info['sector'](PS_point,all_flavor_configurations[0],
                                                  counterterm_index=-1, input_mapping_index=-1)

        # Notice that the Bjorken rescalings for this even will be set during the convolution
        # performed in the next line, and not directly here in the constructor.
        event_to_convolve = ME7Event( PS_point,
                {fc : event_weight for fc in all_flavor_configurations},
                requires_mirroring              = process.get('has_mirror_process'),
                host_contribution_definition    = self.contribution_definition,
                counterterm_structure           = None,
                Bjorken_xs = (xb_1, xb_2)
        )
        convolved_event = self.convolve_event_with_beam_factorization_currents(event_to_convolve,
            process['beam_factorization'], PS_point, process, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2 )

        return convolved_event
    
    def sigma(self, PS_point, process_key, process, all_flavor_configurations,
                                base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, *args, **opts):
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

        # Apply flavor blind cuts
        if not self.pass_flavor_blind_cuts(PS_point, all_flavor_configurations[0]):
            if __debug__: logger.debug('Event failed the flavour_blind generation-level cuts.')
            if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
            return None

        sigma_wgt = self.compute_matrix_element_event_weight(PS_point, process_key, process, all_flavor_configurations,
                                                             base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, *args, **opts)

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

        # Return the lone LO event
        return ME7EventList([
            ME7Event( PS_point, {fc : sigma_wgt for fc in all_flavor_configurations},
                requires_mirroring = process.get('has_mirror_process'),
                host_contribution_definition = self.contribution_definition,
                counterterm_structure = None,
                Bjorken_xs = (xb_1, xb_2)
            )
        ])

    def convolve_event_with_beam_factorization_currents(self, input_event_to_convolve,
            beam_factorization_currents, PS_point, process, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2):
        """ Calls the currents specified in the 'beam_factorization_currents' argument, of the form:
            [{ 'beam_one': <BeamCurrent>,
               'beam_two': <BeamCurrent> }]
        and convolve their evaluations with the ME7Event specified in input so to construct
        a new convolved ME7Event which is then returned. This function is placed in the mother
        class because it is useful both in ME7_integrand_R and ME7_integrand_V."""

        # Adjust the format of the input:
        assert('correlated_convolution' not in beam_factorization_currents)
        assert('non_factorisable_convolution' not in beam_factorization_currents)
        beam_factorization_currents = [
            (  bc['beam_one'],  bc['beam_two'] ) for bc in beam_factorization_currents ]

        # Compute the 4-vector Q characterizing this PS point, defined as the sum of all
        # initial_state momenta. This is not expected to be used or relevant here, but we
        # specify it nonetheless so as to conform with the prototype of the call to
        # BeamFactorization currents.
        Q = sum(p for i, p in enumerate(PS_point.to_list()) if i<self.n_initial)

        xis = (xi1, xi2)
        mu_fs = (mu_f1, mu_f2)
        convolved_event = None
        Bjorken_scalings = [None, None]

        for bcs in beam_factorization_currents:
            event_to_convolve = input_event_to_convolve.get_copy()
            for bc in bcs:
                if bc is None:
                    continue
                # This should be a physical contribution with therefore no distribution
                assert(isinstance(bc, subtraction.BeamCurrent) and bc['distribution_type']=='bulk')
                current_evaluation, all_current_results = self.all_MEAccessors(
                    bc, track_leg_numbers=self.accessor_tracks_leg_numbers,
                    higher_PS_point=PS_point, reduced_process=process,
                    xis=xis, mu_r=mu_r, mu_fs=mu_fs, Q=Q, allowed_backward_evolved_flavors=('ALL','ALL') )
                assert(current_evaluation['spin_correlations']==[None,])
                assert(current_evaluation['color_correlations']==[None,])
                assert(len(current_evaluation['Bjorken_rescalings'])==1)
                if current_evaluation['values'].keys()!=[(0,0,0,0),]:
                    raise MadGraph5Error("The beam factorisation currents attached to the physical 'matrix element' event"+
                        " of a contribution must only involve a single value under key (0,0,0,0) in the evaluation they return.")
                event_to_convolve.convolve_initial_states( current_evaluation['values'][(0,0,0,0)] )
                for beam_id in [0,1]:
                    if current_evaluation['Bjorken_rescalings'][0][beam_id] is not None:
                        if Bjorken_scalings[beam_id] is not None:
                            if abs(Bjorken_scalings[beam_id]-current_evaluation['Bjorken_rescalings'][0][beam_id])\
                                                                                /abs(Bjorken_scalings[beam_id]) > 1.0e-10:
                                raise MadGraph5Error("MadNkLO does not support the combination of multiple different "+
                                                                                                "beam Bjorken x's rescalings.")
                        Bjorken_scalings[beam_id] = current_evaluation['Bjorken_rescalings'][0][beam_id]

                # Check if the flavor matrix that multiplied this event was zero
                if event_to_convolve.is_empty():
                    break
            if event_to_convolve.is_empty():
                continue

            if not convolved_event:
                convolved_event = event_to_convolve
            else:
                # Combine the weights and flavor of the event obtained from this convolution
                # with the ones obtained from the previous convolutions in this loop.
                convolved_event += event_to_convolve

        if convolved_event is None or convolved_event.is_empty():
            return None

        Bjorken_scalings = [ (1. if bs is None else 1./bs) for bs in Bjorken_scalings]

        # Now apply the bounds on the xb_i convolution due to the change of variable
        if ((xb_1 is not None) and (xi1 is not None) and xb_1 > xi1) or \
            ((xb_2 is not None) and (xi2 is not None) and  xb_2 > xi2):
            return None

        # Assign the Bjorken scaling found
        convolved_event.set_Bjorken_rescalings(Bjorken_scalings[0], Bjorken_scalings[1])

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

class ME7Integrand_V(ME7Integrand):
    """ ME7Integrand for the computation of a virtual type of contribution."""
    
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
            if len(self.integrated_counterterms[process_key]) > 0:
                long_res = [' | with the following integrated counterterms:']
            else:
                long_res = [' | with no integrated counterterm']
            for i_CT, CT_properties in enumerate(self.integrated_counterterms[process_key]):
                CT = CT_properties['integrated_counterterm']
                if format==2:
                    long_res.append( '   | I%d : %s' % (i_CT, str(CT)))
                elif format==3:
                    long_res.append( '   | I%d : %s' % (i_CT, CT.__str__(
                        print_n=True, print_pdg=True, print_state=True )))
                elif format==4:
                    long_res.append(CT.nice_string("   | I%-3d : "%i_CT))
                elif format>4:
                    long_res.append(CT.nice_string("   | I%-3d : "%i_CT))
                    for key, value in CT_properties.items():
                        if not key in ['integrated_counterterm', 'matching_process_key']:
                            long_res.append( '     + %s : %s'%(key, str(value)))
            res += '\n'.join(long_res)

        return res

    def compute_matrix_element_event_weight_with_I_operator(self,
            PS_point, process_key, process, all_flavor_configurations, base_weight,
                                        mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts):
        """ Computes the virtual matrix element weight as a DUMMY with correct poles (but not finite part of course)
        using Catani's I operator."""

        # Now this function below is for now a copy of the `generate_matrix_element_event` function of the
        # base class but it should be modified in order to include the I1^(0) operator

        is_loop_induced = self.contribution_definition.process_definition.get('NLO_mode').startswith('sqrvirt')

        #TODO DEV NICO

        model = process.get('model')

        # Prepare the overall result
        evaluation = base_objects.EpsilonExpansion({0:0})

        # Prepare the overall prefactor
        ## Retrieve alpha_s and mu_r
        model_param_dict = model.get('parameter_dict')
        alpha_s = model_param_dict['aS']

        mu_r = model_param_dict['MU_R']
        # ## Build the total momentum
        # Q_square = (
        #     sum([PS_point[l['number']] for l in process.get_initial_legs()])
        # ).square()

        ## The prefactor is alpha_s/(2pi) * (mu_r^2/Q^2)^epsilon* (4pi)^epsilon/Gamma(1-epsilon) )

        prefactor = base_objects.EpsilonExpansion({0:- alpha_s / (2 * math.pi)})

        # ### We build the expansion of (mu_r^2/Q^2)^epsilon
        # logMuQ = math.log(mu_r ** 2 / Q_square)
        # prefactor *= base_objects.EpsilonExpansion({0: 1., 1: logMuQ, 2: 0.5 * logMuQ ** 2})

        ## we divide by an overall Sepsilon = (4pi)^epsilon exp(-epsilon EulerGamma) which cancels the
        ## factor (4pi)^epsilon/Gamma(1-epsilon) up to order epsilon^2. Here's the ratio
        prefactor *= base_objects.EpsilonExpansion({0:1. , 1:0. , 2:math.pi**2/12.})


        # Get the squared matrix element without color correlations:
        ME_evaluation, all_results = self.all_MEAccessors(
            process, PS_point, alpha_s, mu_r, pdgs=all_flavor_configurations[0])

        ME_diagonal = base_objects.EpsilonExpansion(ME_evaluation)
        diagonal_terms = 0

        # Obtain the number of light flavors
        Nf = model.get_nflav()

        anomalous_dimension = {
            3 : 1.5*SU3.CF,
           -3 : 1.5*SU3.CF,
            8 : 11./6.*SU3.CA - 2./3.*SU3.TR*Nf
        }
        # Get the list of colored particles in the process
        colored_legs = process.get_legs_with_color()
        # Build color correlators and call each color-correlated matrix element
        for i, a in enumerate(colored_legs):
            # Collect info about parton a
            part_a = model.get_particle(a.get('id'))
            pa = PS_point[a.get('number')]

            # Build the diagonal part of the Catani operator
            # 0908.4272 eq B.2 line 1
            if part_a.get('mass') == 'ZERO':
                term = base_objects.EpsilonExpansion({
                    -2: SU3.casimir(part_a.get('color')),
                    -1: anomalous_dimension[part_a.get('color')]

                })
            else:
                term = base_objects.EpsilonExpansion({
                    -1: SU3.casimir(part_a.get('color'))
                })
            diagonal_terms += term

            # Build the dipole part, using the symmetry (a,b) <-> (b,a)
            for b in colored_legs[i+1:]:
                # collect info about b
                part_b = model.get_particle(a.get('id'))
                pb = PS_point[b.get('number')]


                factor_ab = base_objects.EpsilonExpansion(0)
                # 0908.4272 eq B.2 line 2: massless a, any b
                if part_a.get('mass') == 'ZERO':
                    factor_ab += base_objects.EpsilonExpansion({
                        -1: math.log(2.*pa.dot(pb)/mu_r**2)
                    })
                # 0908.4272 eq B.2 line 3: massive a, massive b
                elif part_b.get('mass') != 'ZERO':
                    vkl = math.sqrt(
                        1 - pa.square()*pb.square()/ (pa.dot(pb))**2
                    )
                    factor_ab += base_objects.EpsilonExpansion({
                        -1 : 0.5/vkl * math.log( (1+vkl)/(1-vkl) )
                    })
                # 0908.4272 eq B.2 line 3: massive a, massless b
                else:
                    m_a = model_param_dict[part_a.get('mass')]
                    factor_ab += base_objects.EpsilonExpansion({
                        -1: -0.5 * math.log(m_a)
                    })


                correlator = [tuple([a.get('number'),b.get('number')])]
                ME_evaluation, all_results = self.all_MEAccessors(
                    process, PS_point, alpha_s, mu_r, pdgs=all_flavor_configurations[0],
                    color_correlation=correlator)
                ME_ab = base_objects.EpsilonExpansion(ME_evaluation)
                evaluation+= 2.*ME_ab*factor_ab


        evaluation += diagonal_terms * ME_diagonal

        evaluation *= prefactor
        evaluation.truncate(-2,0)

        event_weight = evaluation * base_weight

        return event_weight

    def compute_matrix_element_event_weight(self,PS_point, process_key, process, all_flavor_configurations,
                  base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts):
        """ Returns the weight for the matrix element event, which may need further processing in daughter classes, for
          example in the virtual where one may use a dummy dressed with the I operator."""

        # NLO_mode for the virtuals can be of the form 'tree_DUMMY%dLOOP' with %d being the loop count.
        # It can also be sqrtvirt_DUMMY%dLOOP for multi-loop-induced virtual.
        NLO_modes = self.contribution_definition.process_definition.get('NLO_mode').split('_')

        if len(NLO_modes)==2 and NLO_modes[1].startswith('DUMMY'):
            return self.compute_matrix_element_event_weight_with_I_operator(
                            PS_point, process_key, process, all_flavor_configurations,
                                      base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts)
        else:
            return super(ME7Integrand_V, self).compute_matrix_element_event_weight(
                            PS_point, process_key, process, all_flavor_configurations,
                                      base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts)

    def evaluate_integrated_counterterm(self, integrated_CT_characteristics, PS_point,
                                            base_weight, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, input_mapping,
                                            all_virtual_ME_flavor_configurations, hel_config=None, compute_poles=True,
                                            sector=(None, -1, -1), **opts):
        """ Evaluates the specified integrated counterterm, provided along with its other
        characteristics, like for example the list of flavors assignments that the resolved
        process it corresponds to can take. This function returns an ME7Event specifying the
        counterevent for the specified input_mapping considered (i.e. an integrated CT like
        g > d d~ must be "attached" to all final-state gluons of the virtual process definition).
        """
        #logger.info('Entering evaluate_integrated_counterterm') 

        # Make sure no helicity configuration is specified since this is not supported yet.
        assert ((hel_config is None))

        # Access the various characteristics of the integrated counterterm passed to this
        # function.
        counterterm = integrated_CT_characteristics['integrated_counterterm']

        # The resolved flavor configurations dictionary (mapping resolved to reduced) is
        # typically no longer useful here.
        resolved_flavors = integrated_CT_characteristics['resolved_flavors_combinations']
        reduced_flavors = integrated_CT_characteristics['reduced_flavors_combinations']
        reduced_flavors_with_resolved_initial_states = \
            integrated_CT_characteristics['reduced_flavors_with_resolved_initial_states_combinations']
        symmetry_factor = integrated_CT_characteristics['symmetry_factor']
        beam_convolution_masks = integrated_CT_characteristics['allowed_backward_evolved_flavors']
        allowed_backward_evolved_flavors1 = beam_convolution_masks['beam_one']
        allowed_backward_evolved_flavors2 = beam_convolution_masks['beam_two']
        is_reduced_process_mirrored = counterterm.process.get('has_mirror_process')

        # Typical debug lines below
        # misc.sprint(counterterm)
        # misc.sprint(beam_convolution_masks)
        # misc.sprint(resolved_flavors)
        # misc.sprint(reduced_flavors)
        # misc.sprint(reduced_flavors_with_resolved_initial_states)
        # import pdb
        # pdb.set_trace()

        # And the multiplicity prefactor coming from the several *resolved* flavor assignment
        # that this counterterm can lead to. Typically an integrated counterterm for g > q qbar
        # splitting will have the same reduced flavors, but with the gluon coming from
        # n_f different massless flavors. So that its multiplicitiy factor is n_f.
        # --
        # Also, if you think of the resolved flavors e+ e- > c c~ u u~, the two counterterms
        # C(3,4) and C(5,6) will lead to the reduced flavors g u u~ and c c~ g respectively.
        # These two are however mapped in the same virtual group. In this case, the flavor
        # configuration 'g u u~' in the ME7Event would therefore not receive any contribution
        # from the integrated CT C(5,6), but the configuration 'g c c~' will.
        # --
        # One other subtlety occurs because the backward flavor evolution in the initial
        # states will be registered in the event since it must be convoluted with the correct PDFs.
        # Building upon the previous example, we can consider the following process:
        #             g(1) s(2) > s(3) c(4) c~(5) u(6) u~(7)
        #   + mapped  g(1) q(2) > q(3) qprime(4) qprime~(5) u(6) u~(7)
        # And the counterterm C(2,3)C(4,5) which yields the following reduced mapped flavor:
        #             g g > g u u~
        # With a multiplicity of naively n_f**2. However, the actual events that will be generated
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
        # which will allow to apply the corresponding multiplicity factor 'n_f' at the end.
        # --
        # Finally, one last issue when considering the grouping of flavors is affecting cases like:
        #
        #    C(2,5)  of              u(1) u~(2) > z(3) u(4) u~(5)
        #                 + mapped   c(1) c~(2) > z(3) c(4) c~(5)
        #
        # whose reduced flavors are:
        #                             u(1) g > z(3) u(4)
        #                 + mapped    c(1) g > z(3) c(4)
        # unfortunately, the backward evolution matrix returned by the current C(2,5) with
        # beam_convolution_masks of (-2,-4) will yield contributions like:
        #
        #                             u(1) c~(2) > z(3) u(4)
        #
        # which do not belong here (and would be double-counted if kept here). For this reason,
        # the dictionary 'reduced_flavors_with_resolved_initial_states' discussed above is used
        # not only so as to include the correct multiplicity factor to each flavor, but also
        # as a post-flavor-backward-evolution mask which will remove all spurious contributions like
        # the ones above.
        n_initial = len(all_virtual_ME_flavor_configurations[0][0])
        n_final = len(all_virtual_ME_flavor_configurations[0][1])
        all_mapped_flavors = []
        for flavors in all_virtual_ME_flavor_configurations:
            mapped_flavors = (
                tuple(flavors[0][input_mapping[i]] for i in range(n_initial)),
                tuple(flavors[1][input_mapping[i] - n_initial] for i in
                      range(n_initial, n_initial + n_final))
            )
            if mapped_flavors in reduced_flavors:
                assert ((mapped_flavors not in all_mapped_flavors))
                all_mapped_flavors.append(mapped_flavors)

        # Now map the momenta
        if isinstance(PS_point, dict):
            # Dictionary format LorentzVectorDict starts at 1
            mapped_PS_point = phase_space_generators.LorentzVectorDict(
                (i + 1, PS_point[input_mapping[i] + 1]) for i in range(n_initial + n_final))
        else:
            # List formatLorentzVectorList starts at 0
            mapped_PS_point = phase_space_generators.LorentzVectorDict(
                (i + 1, PS_point[input_mapping[i]]) for i in range(n_initial + n_final))

        # We must also map the Bjorken x's and the xi rescalings
        xi1, xi2 = [xi1, xi2][input_mapping[0]], [xi1, xi2][input_mapping[1]]
        xb_1, xb_2 = [xb_1, xb_2][input_mapping[0]], [xb_1, xb_2][input_mapping[1]]
        # print('xi1 : ' + str(xi1) + ', ' + 'xi2 : ' + str(xi2) + ', ' + 'xb_1 : ' + str(xb_1) + ', ' + 'xb_2 : ' + str(xb_2))

        # Now compute the reduced quantities which will be necessary for evaluating the
        # integrated current
        reduced_PS = counterterm.get_reduced_kinematics(mapped_PS_point)
        # print('mapped_PS_point : ' + str(mapped_PS_point))
        # print('reduced_PS : ' + str(reduced_PS))


        # Retrieve some possibly relevant model parameters
        alpha_s = self.model.get('parameter_dict')['aS']
        mu_r = self.model.get('parameter_dict')['MU_R']

        # Compute the 4-vector Q characterizing this PS point, defined as the sum of all
        # initial_state momenta, before any mapping is applied.
        # Notice that the initial state momenta in mapped_PS_point include *both* the
        # Bjorken x's rescalings *and* the xi rescalings. However, what is factorised in the
        # initial-state PS factorization is the sum of initial-state momenta *without* the
        # xi rescalings (see Eq.5.20 of https://arxiv.org/pdf/0903.1218.pdf),
        # So we must divide here each momenta by the corresponding xi rescalings.
        rescalings = (
            1. if xi1 is None else 1. / xi1,
            1. if xi2 is None else 1. / xi2,
        )
# #gl
#         total_incoming_momentum = vectors.LorentzVector()
#         for i, p in enumerate(PS_point.to_list()[:self.n_initial]):
#             if xi1 != None and xi2 != None:
#                 rescalings_temp = (
#                     1. / xi1,
#                     1. ,
#                 )
#                 total_incoming_momentum += p * rescalings_temp[i]
#             else:
#                 total_incoming_momentum += p * rescalings[i]

        # print('mapped_PS_point : ' + str(mapped_PS_point))
        # print('ME7 - rescalings : ' + str(rescalings))
        total_incoming_momentum = vectors.LorentzVector()
        for i, p in enumerate(mapped_PS_point.to_list()[:self.n_initial]):
            total_incoming_momentum += p * rescalings[i]
            # print('p * rescalings[i] : ' + str(p * rescalings[i]))

        # With the current design we only consider and support the case where there is only
        # *one* regular (i.e. non-beam) "mapping currents" in the counterterm.
        # Notice that exactly *one* of such currents must return a specific reduced kinematics
        # as it does not make sense to be combine several together
        non_beam_factorization_currents = [
            block for block in counterterm.get_currents_blocks()
            if not any(type(c) in (subtraction.BeamCurrent, subtraction.IntegratedBeamCurrent) for c in block) ]

        all_necessary_ME_calls, disconnected_currents_weight = ME7Integrand_R.generate_all_necessary_ME_calls(
            non_beam_factorization_currents, counterterm.process, reduced_PS,
            self.all_MEAccessors, self.accessor_tracks_leg_numbers,
            compute_poles=compute_poles, momenta_dict=counterterm.momenta_dict,
            xis=(xi1, xi2), mu_fs=(mu_f1, mu_f2), mu_r=mu_r, Q=total_incoming_momentum, sector_info=sector,
            n_initial_legs = self.n_initial
        )

        # For now integrated counterterms are supposed to return a single reduced kinematics with None as identifier
        assert(len(all_necessary_ME_calls)==1)
        assert(all_necessary_ME_calls.keys()[0]==None)
        assert(all_necessary_ME_calls.values()[0] is not None)

        # We can now loop over the reduced kinematics produced by the currents:
        all_events = ME7EventList()

        for beam_currents in counterterm.get_beam_currents():

            # Then evaluate the beam factorization currents this must be done for each reduced kinematics configuration
            # so that it is important that whatever that can be cached in these currents is cached.
            all_necessary_ME_calls_for_these_beam_currents = ME7Integrand_R.process_beam_factorization_currents(
                all_necessary_ME_calls, beam_currents, self.all_MEAccessors,
                self.accessor_tracks_leg_numbers, reduced_PS, counterterm.process, xb_1, xb_2,
                xi1, xi2, mu_r, mu_f1, mu_f2, total_incoming_momentum,
                allowed_backward_evolved_flavors1=allowed_backward_evolved_flavors1,
                allowed_backward_evolved_flavors2=allowed_backward_evolved_flavors2,
                momenta_dict=counterterm.momenta_dict, sector_info = sector,
                n_initial_legs = self.n_initial
            )

            for reduced_kinematics_identifier, (reduced_kinematics, necessary_ME_calls) in \
                                                            all_necessary_ME_calls_for_these_beam_currents.items():

                this_base_weight = base_weight
                cut_weight = 1.0

                # Now perform the combination of the list of spin- and color- correlators to be merged
                # for each necessary ME call identified
                necessary_ME_calls = ME7Integrand_R.merge_correlators_in_necessary_ME_calls(necessary_ME_calls)

                # If there is no necessary ME call left, it is likely because the xi upper bound of the
                # Bjorken x's convolution were not respected. We must now abort the event.
                if len(necessary_ME_calls) == 0:
                    continue

                if sector[0] is not None:
                    this_base_weight *= sector[0](reduced_PS, all_mapped_flavors[0],
                                             counterterm_index=sector[1], input_mapping_index=sector[2])

                # Finally treat the call to the reduced connected matrix elements
                alpha_s = self.model.get('parameter_dict')['aS']
                mu_r = self.model.get('parameter_dict')['MU_R']

                #event_PS = reduced_PS.to_list(ordered_keys=[l.get('number') for l in counterterm.process.get('legs')])

                template_event = ME7Event(
                    mapped_PS_point,
                    {fc: base_objects.EpsilonExpansion({0: cut_weight * this_base_weight}) for fc in all_mapped_flavors},
                    requires_mirroring=is_reduced_process_mirrored,
                    host_contribution_definition=self.contribution_definition,
                    counterterm_structure=(counterterm, resolved_flavors, reduced_kinematics_identifier),
                    Bjorken_xs=(xb_1, xb_2)
                )

                integrated_CT_event = ME7Integrand_R.generate_event_for_counterterm(
                    template_event,
                    disconnected_currents_weight,
                    counterterm.prefactor,
                    necessary_ME_calls,
                    counterterm.process,
                    reduced_PS,
                    alpha_s, mu_r,
                    self.all_MEAccessors
                )

                # Immediately return the integrated CTevent if it is None (which can happen during the flavor
                # convolution for instance:
                if integrated_CT_event is None:
                    continue

                # Finally account for the integrated counterterm multiplicity
                new_weights_per_flavor_configurations = {}
                for fc, wgt in integrated_CT_event.weights_per_flavor_configurations.items():
                    # Account for both the overall symmetry factors S_t and the flavor symmetry factor
                    # identified when the contribution class built this counterterm
                    if fc in reduced_flavors_with_resolved_initial_states:
                        # See comments at the beginning of this function as to how it may be that this
                        # flavor configuration is not part of the 'reduced_flavors_with_resolved_initial_states'
                        # and why it should then not be included.
                        new_weights_per_flavor_configurations[fc] = \
                            wgt * float(symmetry_factor) * reduced_flavors_with_resolved_initial_states[fc]

                integrated_CT_event.weights_per_flavor_configurations = new_weights_per_flavor_configurations
                # Make sure to crash if the Event is now empty since the lines above are not
                # supposed to kill all flavor contributions of the event.
                if integrated_CT_event.is_empty():
                    raise MadEvent7Error(
                        "The post-flavor-convolution masking step disabled all flavor configuration for the " +
                        " following integrated counterterm event. This should never happen.\n%s" % str(integrated_CT_event))

                all_events.append(integrated_CT_event)
        #print('all_events : ' + str(all_events))
        return all_events
        
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
            while virtual_PS_point is None:
                virtual_PS_point, jac, x1s, x2s = self.phase_space_generator.get_PS_point(None)
            xb_1, xi1 = x1s
            xb_2, xi2 = x2s
            # Make sure that xb_k < xik for k=1,2 if xik!=None because we want to make sure that
            # we capture all the contributions
            if not xi1 is None and xi1 < xb_1:
                continue
            if not xi2 is None and xi2 < xb_1:
                continue
            # Then make sure it passes the generation cuts if present
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

            if self.sectors is None or self.sectors[process_key] is None:
                all_sectors = [None, ]
            else:
                all_sectors = self.sectors[process_key]

            for i_sector, sector_info in enumerate(all_sectors):
                if (sector_info is not None) and (self.is_sector_selected is not None) and \
                                    not self.is_sector_selected(defining_process, sector_info['sector']):
                    continue

                if sector_info is not None:
                    misc.sprint('\nConsidering sector %s' % sector_info['sector'])

                all_processes = [defining_process,]+mapped_processes
                all_flavor_configurations = []
                # The process mirroring is accounted for at the very end only
                for proc in all_processes:
                    initial_final_pdgs = proc.get_cached_initial_final_pdgs()
                    all_flavor_configurations.append(initial_final_pdgs)

                a_virtual_PS_point = virtual_PS_point.get_copy()
                a_xb_1, a_xi1 = x1s
                a_xb_2, a_xi2 = x2s
                # We should not test here the flavour sensitive cuts because counterterms
                # act on the flavor space and make an event pass the cut even if the resolved one
                # does not.
                while ( ( test_options['apply_higher_multiplicity_cuts'] or
                          test_options['apply_higher_multiplicity_cuts']    ) and
                       not self.pass_flavor_blind_cuts(
                           a_virtual_PS_point,
                           defining_process.get_cached_initial_final_pdgs(),
                           n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles,
                           xb_1 = a_xb_1, xb_2 = a_xb_2 ) ):
                    n_attempts += 1
                    if n_attempts > max_attempts:
                        break
                    a_virtual_PS_point = None
                    # Phase-space generation can fail when beam convolution is active, because of the condition
                    # Bjorken_x_i < xi_i
                    while a_virtual_PS_point is None:
                        a_virtual_PS_point, a_jac, a_x1s, a_x2s = self.phase_space_generator.get_PS_point(None)
                    a_xb_1, a_xi1 = a_x1s
                    a_xb_2, a_xi2 = a_x2s
                if n_attempts > max_attempts:
                    raise MadEvent7Error(
                        "Could not generate a random kinematic configuration that passes " +
                        "the flavour blind cuts in less than %d attempts." % max_attempts )
                n_attempts = 0

                process_evaluation = self.test_IR_poles_for_process(
                        test_options, defining_process, process_key, all_flavor_configurations,
                        a_virtual_PS_point, a_xi1, a_xi2, a_xb_1, a_xb_2, sector_info=sector_info )

                all_evaluations[(process_key, -1 if sector_info is None else i_sector)] = process_evaluation

        # Now produce a nice output of the evaluations and assess whether this test passed or not.
        return self.analyze_IR_poles_check(all_evaluations, test_options['acceptance_threshold'])

    def test_IR_poles_for_process(self, test_options, defining_process, process_key,
        all_flavor_configurations, a_virtual_PS_point, a_xi1, a_xi2, a_xb_1, a_xb_2, sector_info=None):
        """ Given a test_options, a process and all flavor configurations mapped to that
        particular process, a PS-point onto which to evaluate the various contributions,
        as well as its corresponding convolution variables xi1/2; perform the IR poles
        test that consists in making sure all poles in the epsilon Laurent expansion have
        residues of zero, exhibiting the various pieces (virtual ME, integrated counterterms)
        part of the KLN cancellation."""
        
        # Select which integrated counterterms are selected for the test
        if self.has_integrated_counterterms():
            integrated_counterterms_to_consider = [ (i_ct, ct) for i_ct, ct in enumerate(self.integrated_counterterms[process_key])
                if ct['integrated_counterterm'].get_number_of_couplings() <= test_options['correction_order'].count('N') ]
            if test_options['counterterms']:
                counterterm_pattern = test_options['counterterms']
                logger.debug("Integrated counterterms before selection")
                for i_ct, ct in integrated_counterterms_to_consider:
                    logger.debug("    "+str(ct['integrated_counterterm']))
                integrated_counterterms_to_consider = self.find_counterterms_matching_regexp(
                    integrated_counterterms_to_consider, counterterm_pattern, integrated_CT=True )
            logger.debug("Selected integrated counterterms")
            for i_ct, ct in integrated_counterterms_to_consider:
                logger.debug("    "+str(ct['integrated_counterterm']))
            integrated_counterterms_to_consider = [i_ct for i_ct, ct in integrated_counterterms_to_consider]

        # We must also include local counterterms as they too can have poles
        if self.has_local_counterterms():
            local_counterterms_to_consider = [ (i_ct, ct) for i_ct, ct in enumerate(self.counterterms[process_key])
                if ct.get_number_of_couplings() <= test_options['correction_order'].count('N') ]
            if test_options['counterterms']:
                counterterm_pattern = test_options['counterterms']
                logger.debug("Local counterterms before selection")
                for i_ct, ct in local_counterterms_to_consider:
                    logger.debug("    "+str(ct))
                local_counterterms_to_consider = self.find_counterterms_matching_regexp(
                    local_counterterms_to_consider, counterterm_pattern )
            logger.debug("Selected local counterterms")
            for i_ct, ct in local_counterterms_to_consider:
                logger.debug("    "+str(ct))
            local_counterterms_to_consider = [i_ct for i_ct, ct in local_counterterms_to_consider]

        mu_r, mu_f1, mu_f2 = self.get_scales(a_virtual_PS_point)

        # Specify the counterterms to be considered and backup the original values
        if sector_info is None:
            input_sector_info = {'sector': None,
                                 'counterterms': None,
                                 'integrated_counterterms': None}
        else:
            input_sector_info = copy.deepcopy(sector_info)
        if self.has_local_counterterms():
            if input_sector_info['counterterms'] is None:
                input_sector_info['counterterms'] = local_counterterms_to_consider
            else:
                input_sector_info['counterterms'] = filter(lambda i_ct: i_ct in local_counterterms_to_consider,
                                                           input_sector_info['counterterms'])
        if self.has_integrated_counterterms():
            if input_sector_info['integrated_counterterms'] is None:
                input_sector_info['integrated_counterterms'] = {i_ct: None for i_ct in
                                                                integrated_counterterms_to_consider}
            else:
                for i_ct in list(input_sector_info['integrated_counterterms'].keys()):
                    if i_ct not in integrated_counterterms_to_consider:
                        del input_sector_info['integrated_counterterms'][i_ct]

        try:
            # Now call sigma in order to gather all events
            events = self.sigma(
                a_virtual_PS_point.to_dict(), process_key, defining_process, all_flavor_configurations,
                1.0, mu_r, mu_f1, mu_f2, a_xb_1, a_xb_2, a_xi1, a_xi2,
                compute_poles            = True,
                apply_flavour_blind_cuts = (test_options['apply_lower_multiplicity_cuts'] or
                                            test_options['apply_higher_multiplicity_cuts']),
                sector_info = input_sector_info
            )
        except Exception as e:
            logger.critical("The following exception occurred when generating events for this integrand %s:\n%s"%(
                                                  self.get_short_name(), str(e)))
            # Forward the exception while not resetting the stack trace.
            raise

#        misc.sprint('Events generated before post-processing.')
#        misc.sprint(events)
#        misc.sprint('Before post-processing:',len(events))

        # Now post-process the events as done in __call__ of the integrand.
        
        # Now handle the process mirroring, by adding, for each ME7Event defined with
        # requires_mirroring = True, a new mirrored event with the initial states swapped,
        # as well as the Bjorken x's and rescalings and with p_z -> -p_z on all final
        # state momenta
        events.generate_mirrored_events()

        # Apply flavor blind cuts
        if test_options['apply_higher_multiplicity_cuts'] or test_options['apply_lower_multiplicity_cuts']:
            events.filter_with_flavor_blind_cuts(self.pass_flavor_blind_cuts,
                defining_process.get_cached_initial_final_pdgs(),
                n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles)
        # Apply PDF convolution.
        if self.run_card['lpp1']==self.run_card['lpp2']==1:
            pdf = None if test_options['set_PDFs_to_unity'] else self.pdf
            events.apply_PDF_convolution( self.get_pdfQ2, (pdf, pdf), (mu_f1**2, mu_f2**2) )
        # Apply the 1/xi**2 factor from the change of variable xi -> xi' / xi
        # Remember that the Bjorken rescalings are defined as 1/xi at this stage
# #gl
#         for event in events:
#             if self.contribution_definition.torino_sub_BS:
#                 event *= event.Bjorken_x_rescalings[0]*1.0
#             else:
#                 event *= event.Bjorken_x_rescalings[0]*event.Bjorken_x_rescalings[1]
        for event in events:
            event *= event.Bjorken_x_rescalings[0] * event.Bjorken_x_rescalings[1]

        # Make sure Bjorken-x rescalings xi_i don't matter anymore
        for event in events:
            event.set_Bjorken_rescalings(None, None)

        # Apply flavor sensitive cuts
        if test_options['apply_higher_multiplicity_cuts'] or test_options['apply_lower_multiplicity_cuts']:
            events.filter_flavor_configurations(self.pass_flavor_sensitive_cuts)

#        misc.sprint('Events generated after post-processing:')
#        misc.sprint(events)

        # Now collect the results necessary to test that the poles cancel in the dictionary below
        # Some contributions like the beam-soft ones may not have any physical contributions
        # but just are a receptacle for integrated softs. We therefore initialize it at zero.
        evaluation = {
             'virtual_ME'           : base_objects.EpsilonExpansion({0: 0.0}),
             'integrated_CTs'       : base_objects.EpsilonExpansion({0: 0.0}),
             'defining_process'     : defining_process,
             'PS_point'             : a_virtual_PS_point }
        event_weight_sum = base_objects.EpsilonExpansion({0: 0.0})

        for event in events:
            # TODO: instead of looking at the total weight, we could consider doing the
            # check for one particular flavor configuration which could be specified
            # as an extra test parameter.
            # An alternative way is to look at the overall summed event which should have
            # zero poles for all its flavor configurations
            event_wgt = event.get_total_weight()
            event_weight_sum += event_wgt
            event_str = event.counterterm_structure_short_string()
            if event_str in evaluation:
                evaluation[event_str] += event_wgt
            else:
                evaluation[event_str] = event_wgt
            # Also add the weight from this event in the two aggregation keys
            if event.counterterm_structure is None:
                evaluation['virtual_ME'] += event_wgt
            else:
                evaluation['integrated_CTs'] += event_wgt

        # Small monitoring of the various contributions:
        max_length = max(len(k) for k in evaluation)+1
        separator_length = max_length+2
        logger.debug('-'*separator_length)
        # Fancy ordering key to get a nice printout order of the event weights.
        for entry in sorted(evaluation.keys(), key=lambda e:
            (e[3:].replace('Event','0_')+'1' if e.startswith('<->') else e.replace('Event','0_')) ):
            value = evaluation[entry]
            if entry not in ['defining_process','PS_point','integrated_CTs','virtual_ME']:
                logger.debug(('%%-%ds : %%s'%max_length)%(entry, value.__str__(format='.16e')))
        logger.debug('-'*separator_length)
        logger.debug(('%%-%ds : %%s'%max_length)%('integrated_CTs', evaluation['integrated_CTs'].__str__(format='.16e')))
        logger.debug(('%%-%ds : %%s'%max_length)%('virtual_ME', evaluation['virtual_ME'].__str__(format='.16e')))
        logger.debug('-'*separator_length)
        logger.debug('\nAll events summed:\n%s'%str(event_weight_sum))
        logger.debug('-'*separator_length)

        diff = evaluation['virtual_ME']+evaluation['integrated_CTs']
        relative_diff = evaluation['virtual_ME'].relative_diff(evaluation['integrated_CTs']*-1.)
        # Finite parts are of course expected to differ, so let's not show them
        diff.truncate(max_power = -1)
        relative_diff.truncate(max_power = -1)
        # To be commented out when we will have a full-fledged analysis coded up
        # in analyze_IR_poles_check()
        logger.info('Summary for that PS point:\n%-20s : %s\n%-20s : %s\n%-20s : %s'%(
             'virtual contrib.',
             evaluation['virtual_ME'].__str__(format='.16e'),
             'integrated CTs.',
             evaluation['integrated_CTs'].__str__(format='.16e'),
             'relative diff.',
             relative_diff.__str__(format='.16e')
        ))
        normalization = max(
            max(abs(v) for v in evaluation['virtual_ME'].values()) if len(evaluation['virtual_ME'])>0 else 0.,
            max(abs(v) for v in evaluation['integrated_CTs'].values()) if len(evaluation['integrated_CTs'])>0 else 0.)
        if normalization == 0.:
            logger.info('%sResults identically zero.%s'%(misc.bcolors.GREEN, misc.bcolors.ENDC))
        else:
            max_diff = max( abs(v/normalization) for v in diff.values()) if len(diff)>0 else 0.
            if max_diff > test_options['acceptance_threshold']:
                logger.info('%sFAILED%s (%.2e > %.2e)'%(misc.bcolors.RED, misc.bcolors.ENDC,
                                                max_diff, test_options['acceptance_threshold']))
            else:
                logger.info('%sPASSED%s (%.2e < %.2e)'%(misc.bcolors.GREEN, misc.bcolors.ENDC,
                                                max_diff, test_options['acceptance_threshold']))
            
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
                             base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts):
        """ Overloading of the sigma function from ME7Integrand to include necessary
        additional contributions. """

        all_events_generated = ME7EventList()

        compute_poles = opts.get('compute_poles', False)
        sector_info = opts.get('sector_info', None)

        matrix_element_event = self.generate_matrix_element_event(
                PS_point, process_key, process, all_flavor_configurations,
                  base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts)

        # Some contributions might have not physical contributions and overloaded the above
        # so as to return None
        if matrix_element_event is not None:
            all_events_generated.append(matrix_element_event)

        # Now loop over all integrated counterterms
        for i_ct, counterterm_characteristics in enumerate(self.integrated_counterterms[process_key]):

            # Skip integrated counterterms that do not belong to this sector
            selected_input_mappings = None
            if (sector_info is not None) and (sector_info['integrated_counterterms'] is not None):
                if i_ct not in sector_info['integrated_counterterms']:
                    continue
                selected_input_mappings = sector_info['integrated_counterterms'][i_ct]

            # Example of a hack below to include only soft integrated CT. Uncomment to enable.
            #ss = counterterm_characteristics['integrated_counterterm'].reconstruct_complete_singular_structure()
            #ss = ss.substructures[0].substructures[0].substructures[0]
            #integrated_current = counterterm_characteristics['integrated_counterterm'].get_integrated_current()
            #if not (ss.name()=='C' and len(ss.substructures)==0 and isinstance(integrated_current, subtraction.BeamCurrent) and \
            #                                                    integrated_current['distribution_type']=='counterterm'):
            #    continue

            # And over all the ways in which this current PS point must be remapped to
            # account for all contributions of the integrated CT. (e.g. the integrated
            # splitting g > d d~ must be "attached" to all final state gluons appearing
            # in this virtual ME process definition.)
            for i_mapping, input_mapping in enumerate(counterterm_characteristics['input_mappings']):

                # Skip integrated counterterms that do not belong to this sector
                if (selected_input_mappings is not None) and i_mapping not in selected_input_mappings:
                    continue

                # At NLO at least, it is OK to save a bit of time by enforcing 'compute_poles=False').
                # This will need to be re-assessed at NNLO for RV contributions.
                CT_events = self.evaluate_integrated_counterterm(
                    counterterm_characteristics, PS_point, base_weight, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2,
                    input_mapping, all_flavor_configurations,
                    hel_config      = None,
                    compute_poles   = compute_poles,
                    sector          = (sector_info['sector'] if sector_info else None, i_ct, i_mapping)
                )
            
                for CT_event in CT_events:
                    # Attach additional information to this CT_event which plays no role in the
                    # MadNkLO construction but which may be used, in test_IR_poles for instance,
                    # for improving the printout of the event record.
                    if ('allowed_backward_evolved_flavors' in counterterm_characteristics) and \
                      not all(aef=='ALL' for aef in counterterm_characteristics['allowed_backward_evolved_flavors'].values()):
                        CT_event.store_information('beam_convolution_masks',
                                counterterm_characteristics['allowed_backward_evolved_flavors'])
                    CT_event.store_information('input_mapping',input_mapping)
                    all_events_generated.append( CT_event )

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
        # Initialize the ME7Integrand
        super(ME7Integrand_R, self).__init__(*args, **opts)
        # Update the initialization inputs
        self.initialization_inputs['options']['counterterms'] = self.counterterms

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
            if len(self.counterterms[process_key])>0:
                long_res = [' | with the following local counterterms:']
            else:
                long_res = [' | with no local counterterm']
            for i_CT, CT in enumerate(self.counterterms[process_key]):
                if CT.is_singular():
                    if format==2:
                        long_res.append( '   | L%d : %s' % (i_CT, str(CT)))
                    elif format==3:
                        long_res.append( '   | L%d : %s' % (i_CT, CT.__str__(
                            print_n=True, print_pdg=True, print_state=True ) ))
                    elif format>3:
                        long_res.append(CT.nice_string("   | L%-3d : "%i_CT))
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
        """This function takes several spin correlators specified in the form
        spin_correlator = ( (legIndex_1, ( vector_A, vector_B, ...) ),
                            (legIndex_2, ( vector_C, vector_D, ...) ),
                            etc... )
        and returns the new correlator specifier
        that arises from the combination of all of them.
        """
        
        # This combination is done with a simple concatenation of the lists.
        # Example:
        # Say correlator A only correlates to leg #1 with two four-vectors v_A and v_B
        # (these are replacement of polarization vectors to consider and summed over)
        #     correlators_A[0]  = ( (1, (v_A, v_B)), )
        # Now correlator B correlates with the two legs #4 and #7
        # (they must be different by construction!)
        # each defined with a single vector v_C and v_B
        #     correlators_B[0]  = ( (4, (v_C,)), (7, (v_D,)) )
        # Then it is clear tha the combined spin correlation should be:
        #     combined_spin_correlator = ( (1, (v_A, v_B)), (4, (v_C,)), (7, (v_D,)) )
        # Notice that this implies that both combinations
        #    pol_vec_1 = v_A,  pol_vec_4 = v_C,  pol_vec_7 = v_D
        # and
        #    pol_vec_1 = v_B,  pol_vec_4 = v_C,  pol_vec_7 = v_D
        # will be computed and summed in the resulting
        # spin-correlated matrix element call.
        
        # Trivial combination if there is a single one or none (all are None)
        if len(spin_correlators) == 1:
            return spin_correlators[0]
        if all(sc is None for sc in spin_correlators):
            return None

        # Make sure the spin correlators don't share common legs
        for i, sc_a in enumerate(spin_correlators[:-1]):
            if sc_a is None:
                continue
            set_a = set(c[0] for c in sc_a)
            for sc_b in spin_correlators[i+1:]:
                if sc_b is None:
                    continue
                set_b = set(c[0] for c in sc_b)
                assert set_a.isdisjoint(set_b)

        # Combine spin correlators
        combined = []
        for sc in spin_correlators:
            if not sc is None:
                combined += sc

        return tuple(combined)

    @classmethod
    def flavor_matrix_sanity_check(cls, flavor_matrix):
        """ Sanity check of a single flavor matrix to be applied. """

        beam_one_active = False
        beam_two_active = False
        msg = None
        for start_flavors, all_end_flavors_groups in flavor_matrix.items():
            beam_one_active = beam_one_active or (start_flavors[0] is not None)
            beam_two_active = beam_two_active or (start_flavors[1] is not None)
            for beam_id in [0,1]:
                if ((start_flavors[beam_id] is None) and (not all(
                        all(end_flavors[beam_id] is None for end_flavors in end_flavors_group)
                        for end_flavors_group in all_end_flavors_groups.keys() ))) or \
                    ((start_flavors[beam_id] is not None) and (not all(
                        all(end_flavors[beam_id] is not None for end_flavors in end_flavors_group)
                        for end_flavors_group in all_end_flavors_groups.keys() ))):
                    return beam_one_active, beam_two_active, "Inconsistent flavor matrix: %s"%pformat(flavor_matrix)

        return beam_one_active, beam_two_active, msg

    @classmethod
    def flavor_matrices_sanity_check(cls, flavor_matrices):
        """ Sanity check of multiple flavor matrices to be combined. """

        # Combining more than 2 flavor matrices should never make sense
        if len(flavor_matrices)>2:
            return 'MadNkLO does not support the combination of more than two beam convolution flavor matrices.'

        # First assess the sanity of each matrix
        found_beam_one_active = False
        found_beam_two_active = False
        for flavor_matrix in flavor_matrices:
            beam_one_active, beam_two_active, msg = cls.flavor_matrix_sanity_check(flavor_matrix)
            if msg is not None:
                return msg
            if (found_beam_one_active and beam_one_active) or (found_beam_two_active and beam_two_active):
                return 'MadNkLO does not support the iterative application of flavor matrices affecting the same beam.'
            found_beam_one_active = found_beam_one_active or beam_one_active
            found_beam_two_active = found_beam_one_active or beam_two_active

        return None

    @classmethod
    def merge_correlators_in_necessary_ME_calls(cls, all_necessary_ME_calls):
        """ Merges the correlators defined the list of inputs to provide to the various calls
        to the reduced matrix elements specified by the argument all_necessary_ME_calls. """
        
        new_all_necessary_ME_calls = []
        for ME_call in all_necessary_ME_calls:
            combined_spin_correlator  = cls.combine_spin_correlators(ME_call['spin_correlations'])
            # The combination of color correlators can give rise to several ones, each with its multiplier
            combined_color_correlators, multipliers = cls.combine_color_correlators(ME_call['color_correlations'])

            combined_weight = reduce(lambda x,y: x*y, ME_call['main_weights'])
            flavor_matrices = [fv for fv in ME_call['flavor_matrix'] if fv is not None]
            # Now perform a sanity check of the flavor matrices
            msg = cls.flavor_matrices_sanity_check(flavor_matrices)
            if msg is not None:
                raise MadGraph5Error(msg)
            Bjorken_rescalings_beam_one = [br for br in ME_call['Bjorken_rescalings_beam_one'] if br is not None]
            Bjorken_rescalings_beam_two = [br for br in ME_call['Bjorken_rescalings_beam_two'] if br is not None]
            # Make sure beam rescalings are either unique or equal.
            for beam_rescalings in [Bjorken_rescalings_beam_one, Bjorken_rescalings_beam_two]:
                if len(beam_rescalings)>1:
                    if not all( (abs(rescaling-beam_rescalings[0])/abs(beam_rescalings[0])) < 1.0e-10 for rescaling in beam_rescalings ):
                        raise MadGraph5Error("MadNkLO does not support the combination of multiple different "+
                                                                                            "beam Bjorken x's rescalings.")
            Bjorken_rescaling_beam_one = Bjorken_rescalings_beam_one[0] if len(Bjorken_rescalings_beam_one)>=1 else 1.0
            Bjorken_rescaling_beam_two = Bjorken_rescalings_beam_two[0] if len(Bjorken_rescalings_beam_two)>=1 else 1.0

            for combined_color_correlator, multiplier in zip(combined_color_correlators,multipliers):
                # Finally add the processed combination as a new ME call
                new_all_necessary_ME_calls.append(
                    {   'spin_correlation'            : combined_spin_correlator,
                        'color_correlation'           : combined_color_correlator,
                        'main_weight'                 : combined_weight*multiplier,
                        'flavor_matrices'             : flavor_matrices,
                        'Bjorken_rescaling_beam_one'  : Bjorken_rescaling_beam_one,
                        'Bjorken_rescaling_beam_two'  : Bjorken_rescaling_beam_two,
                    },
                )
            # print('LINEA 3392 - combined_weight : ' + str(combined_weight))
        return new_all_necessary_ME_calls


    @classmethod
    def update_all_necessary_ME_calls_for_specific_reduced_kinematics(cls,
            necessary_ME_calls_fixed_kinematics, new_evaluation_fixed_kinematics, weight_type='main_weight'):
        """ Combined previous current evaluation with the new one in argument so as to setup
        the input of the matrix elements to be computed. The current type can be either:
            'main_weight', 'flavor_matrices_beam_one', 'flavor_matrices_beam_two'
        which indicates if this is a weight that necessitate flavor convolutions (i.e. from
        beam factorization).
        This combination is done for a single entry specifying a particular reduced kinematics.
        """

        new_necessary_ME_calls = []
        for ((spin_index, color_index, bjorken_rescalings_index), current_wgt) in new_evaluation_fixed_kinematics['values'].items():
            # Now combine the correlators necessary for this current
            # with those already specified in 'all_necessary_ME_calls'
            for ME_call in necessary_ME_calls_fixed_kinematics:
                rescaling_one = new_evaluation_fixed_kinematics['Bjorken_rescalings'][bjorken_rescalings_index][0]
                rescaling_two = new_evaluation_fixed_kinematics['Bjorken_rescalings'][bjorken_rescalings_index][1]
                new_necessary_ME_calls.append({
                    # Append this spin correlation to those already present for that call
                    'spin_correlations': ME_call['spin_correlations'] + [new_evaluation_fixed_kinematics['spin_correlations'][spin_index], ],
                    # Append this color correlation to those already present for that call
                    'color_correlations': ME_call['color_correlations'] + [
                        new_evaluation_fixed_kinematics['color_correlations'][color_index], ],
                    # Append this weight to those already present for that call
                    'main_weights': ME_call['main_weights'] + [ base_objects.EpsilonExpansion(current_wgt) if
                            weight_type == 'main_weight' else base_objects.EpsilonExpansion({'finite': 1.0}) ],
                    'flavor_matrix': ME_call['flavor_matrix'] + [ current_wgt if weight_type == 'flavor_matrix' else None ],
                    'Bjorken_rescalings_beam_one' : ME_call['Bjorken_rescalings_beam_one'] + [
                                None if rescaling_one is None else 1./rescaling_one,],
                    'Bjorken_rescalings_beam_two' : ME_call['Bjorken_rescalings_beam_two'] + [
                                None if rescaling_two is None else 1. / rescaling_two,],
                })
        return new_necessary_ME_calls

    @classmethod
    def update_all_necessary_ME_calls(cls, all_necessary_ME_calls, new_evaluation, **opts):
        """ Combined previous current evaluation with the new one in argument so as to setup
        the input of the matrix elements to be computed. The current type can be either:
            'main_weight', 'flavor_matrices_beam_one', 'flavor_matrices_beam_two'
        which indicates if this is a weight that necessitate flavor convolutions (i.e. from
        beam factorization).
        This combination is done across all specified reduced kinematics.
        """

        # First breakdown the current evaluation into chunks for each different reduced kinematics
        current_evaluation_per_reduced_kinematics = {}
        for ((spin_index, color_index, bjorken_rescalings_index, reduced_kinematics_index), current_wgt) in new_evaluation['values'].items():
            if reduced_kinematics_index not in current_evaluation_per_reduced_kinematics:
                current_evaluation_per_reduced_kinematics[reduced_kinematics_index] = dict(new_evaluation)
                # We want to overwrite the values as the reduced_kinematics index is treated externally
                current_evaluation_per_reduced_kinematics[reduced_kinematics_index]['values'] = {}
            current_evaluation_per_reduced_kinematics[reduced_kinematics_index]['values']\
                                                [(spin_index, color_index, bjorken_rescalings_index)] = current_wgt

        # Now make sure that we will not be combining reduced kinematics from different currents
        # This check is redundant with he one done below
        #if new_all_necessary_ME_calls.keys() != [None,] and \
        #                any(new_evaluation['reduced_kinematics'][ri]!=None for ri in reduced_kinematics_chunks):
        #    raise MadEvent7Error("MadNkLO cannot combine several subtraction currents that each return mapped kinematics.")

        new_all_necessary_ME_calls = {}
        for reduced_kinematics_identifier, (reduced_kinematics, necessary_ME_calls) in all_necessary_ME_calls.items():

            for reduced_kinematics_index, current_evaluation in current_evaluation_per_reduced_kinematics.items():

                # "combine" reduced kinematics
                if new_evaluation['reduced_kinematics'][reduced_kinematics_index] is None:
                    new_reduced_kinematics_identifier = reduced_kinematics_identifier
                    new_reduced_kinematics = reduced_kinematics
                elif reduced_kinematics_identifier is None:
                    new_reduced_kinematics_identifier = new_evaluation['reduced_kinematics'][reduced_kinematics_index][0]
                    new_reduced_kinematics = new_evaluation['reduced_kinematics'][reduced_kinematics_index][1]
                else:
                    # Note: in subtraction scheme where one may want to combine the rescalings
                    raise MadEvent7Error(
                        "MadNkLO cannot combine several subtraction currents evaluation that each different return mapped kinematics.")

                new_necessary_ME_calls = cls.update_all_necessary_ME_calls_for_specific_reduced_kinematics(
                                                                    necessary_ME_calls , current_evaluation, **opts)

                new_all_necessary_ME_calls[new_reduced_kinematics_identifier] = (new_reduced_kinematics, new_necessary_ME_calls)

        # Return the new list of necessary ME calls
        return new_all_necessary_ME_calls

    @classmethod
    def process_beam_factorization_currents(cls, all_necessary_ME_calls, beam_currents,
                     all_MEAccessors, track_leg_numbers, PS_point, process, xb_1, xb_2, xi1, xi2, mu_r, mu_f1, mu_f2, Q,
                     allowed_backward_evolved_flavors1='ALL',
                     allowed_backward_evolved_flavors2='ALL', **opts):
        """ Calls the beam currents specified in the argument all_beam_currents with format:
              [{       'beam_one': <BeamCurrent>,
                       'beam_two': <BeamCurrent>,
           <optional>  'correlated_convolution'       : <beam_current>,
           <optional>  'non_factorisable_convolution' : <beam_current>
               }]
        and combine their evaluation with the list of already existing inputs provided in
        the list of all_necessary_ME_calls. The two options allowed_backward_evolved_flavors<i>
        allow to specify to the BeamCurrents a mask to apply, selecting only particular
        flavors that the current can backward evolve onto."""

        xis = (xi1, xi2)
        mu_fs = (mu_f1, mu_f2)
        allowed_backward_evolved_flavors = (allowed_backward_evolved_flavors1, allowed_backward_evolved_flavors2)

        new_all_necessary_ME_calls = dict(all_necessary_ME_calls)

        beam_currents_to_call = []

        if beam_currents['beam_one'] is not None:
            beam_currents_to_call.append(beam_currents['beam_one'])

        if beam_currents['beam_two'] is not None:
            beam_currents_to_call.append(beam_currents['beam_two'])

        if 'correlated_convolution' in beam_currents and beam_currents['correlated_convolution'] is not None:

            if beam_currents['beam_one'] is not None or beam_currents['beam_two'] is not None:
                raise MadGraph5Error('MadNkLO does not support the specification of both a current'+
                    ' requiring correlated beam convolution with one-sided convolution currents.')

            if (not isinstance(beam_currents['correlated_convolution'],
                        subtraction.IntegratedBeamCurrent)) and (xi1 != xi2 or xi1 is None):
                if 'torino_sub_BS' in beam_currents:    #gl
                    # print('PASSING')
                    pass
                else:
                    raise MadGraph5Error('Currents requiring correlated beam convolutions must be'
                        ' evaluated within a contribution in which xi1 == xi2 and different than None.')

            # if (not isinstance(beam_currents['correlated_convolution'],
            #         subtraction.IntegratedBeamCurrent)) and (xi1 != xi2 or xi1 is None):
                
            #     raise MadGraph5Error('Currents requiring correlated beam convolutions must be'
            #         ' evaluated within a contribution in which xi1 == xi2 and different than None.')

            beam_currents_to_call.append(beam_currents['correlated_convolution'])

        if 'non_factorisable_convolution' in beam_currents and beam_currents['non_factorisable_convolution'] is not None:

            if beam_currents['beam_one'] is not None or beam_currents['beam_two'] is not None:
                raise MadGraph5Error('MadNkLO does not support the specification of both a current'+
                    ' requiring a non-factorisable beam convolution with one-sided convolution currents.')

            beam_currents_to_call.append(beam_currents['non_factorisable_convolution'])

        for beam_current in beam_currents_to_call:

            current_evaluation, all_current_results = all_MEAccessors(beam_current,
                track_leg_numbers=track_leg_numbers, higher_PS_point=PS_point, reduced_process=process,
                mu_r=mu_r, mu_fs=mu_fs, xis=xis, Q=Q,
                allowed_backward_evolved_flavors = allowed_backward_evolved_flavors, **opts)

            # weight_type is a class attribute of the current_evaluation which is an instance of either the class
            # utils.SubtractionCurrentEvaluation or utils.BeamFactorizationCurrentEvaluation
            # and will be 'main_weight' in the former case and 'flavor_matrix' in the latter
            new_all_necessary_ME_calls = ME7Integrand_R.update_all_necessary_ME_calls(
                        new_all_necessary_ME_calls, current_evaluation, weight_type=current_evaluation.weight_type )

        return new_all_necessary_ME_calls

    def evaluate_local_counterterm(
        self, counterterm, PS_point, base_weight,
        mu_r, mu_f1, mu_f2,
        xb_1, xb_2, xi1, xi2,
        is_process_mirrored, all_resolved_flavors,
        hel_config=None, apply_flavour_blind_cuts=True,
        boost_back_to_com=True, always_generate_event=False, sector=(None,-1), **opts ):
        """Evaluate a counterterm for a given PS point and flavors.
        The option 'boost_back_to_com' allows to prevent this function to boost back the
        PS point of the Event generated to the c.o.m frame. This would yields *incorrect*
        distributions but it is useful for diagnostic purposes when called from within
        test_IR_limits. It however *must* be kept to True by default."""

        # Note that the code below can become more complicated when tracking helicities,
        # but let's forget this for now.
        assert ((hel_config is None))

        # Ez
        # logger.debug("Beginning of evaluate_local_counterterm")
        # logger.debug("counterterm: " + str(counterterm))
        # logger.debug("boost_back_to_com: " + str(boost_back_to_com))
        # logger.debug(str(PS_point))
        # logger.debug("+" * 100)

        # Compute the 4-vector Q characterizing this PS point, defined as the sum of all
        # initial_state momenta, before any mapping is applied.
        # Notice that the initial state momenta in mapped_PS_point include *both* the
        # Bjorken x's rescalings *and* the xi rescalings. However, what is factorised in the
        # initial-state PS factorization is the sum of initial-state moment *without* the
        # xi rescalings (see Eq.5.20 of https://arxiv.org/pdf/0903.1218.pdf),
        # So we must divide here each momentum by the xi rescalings.
        rescalings = (
            1. if xi1 is None else 1. / xi1,
            1. if xi2 is None else 1. / xi2,
        )
        total_incoming_momentum = vectors.LorentzVector()
        for i, p in enumerate(PS_point.to_list()[:self.n_initial]):
            total_incoming_momentum += p * rescalings[i]


        #logger1.info('options : ' + str(opts))
# Gl
        #In the following we'll work with a copy of PS_point in order to preserve the original one from any border effects

#        copied_PS_point = PS_point.get_copy()

        #logger1.info('NEW EVALUATION: PS_point evaluate= ' + str(PS_point))



        # With the current design we only consider and support the case where there is only
        # *one* regular (i.e. non-beam) "mapping currents" in the counterterm.
        # Notice that exactly *one* of such currents must return a specific reduced kinematics
        # as it does not make sense to be combine several together
        non_beam_factorization_currents = [
            block for block in counterterm.get_currents_blocks()
            if not any(type(c) in (subtraction.BeamCurrent, subtraction.IntegratedBeamCurrent) for c in block) ]

        # Access the matrix element characteristics
        ME_process = counterterm.current

        n_unresolved_left = self.contribution_definition.n_unresolved_particles
        n_unresolved_left -= counterterm.count_unresolved()

        all_necessary_ME_calls, disconnected_currents_weight = ME7Integrand_R.generate_all_necessary_ME_calls(
            non_beam_factorization_currents, ME_process, PS_point,
            self.all_MEAccessors, self.accessor_tracks_leg_numbers,
            momenta_dict = counterterm.momenta_dict,
            xis=(xi1, xi2), mu_fs=(mu_f1, mu_f2), mu_r=mu_r, Q=total_incoming_momentum,
            sector_info = sector,
            n_initial_legs = self.n_initial
        )

        # We can now loop over the reduced kinematics produced by the currents:
        all_events = ME7EventList()

        # Local counterterms can return several reduced kinematics (because of the soft treatment in Catani's dipoles
        # for instance). We therefore now iterate over all of them and apply the beam currents to each separately
        for reduced_kinematics_identifier, (reduced_kinematics, necessary_ME_calls) in all_necessary_ME_calls.items():

            all_necessary_ME_calls_for_this_reduced_kinematics = {
                                            reduced_kinematics_identifier: (reduced_kinematics, necessary_ME_calls) }

            # And now we apply the beam currents
            for beam_currents in counterterm.get_beam_currents():

                # Evaluating the beam factorization currents this must be done for each reduced kinematics configuration
                # so that it is important that whatever that can be cached in these currents is cached.
                # Notice that we pass here the PS_point to the currents and not any reduced kinematics
                all_necessary_ME_calls_for_these_beam_currents = ME7Integrand_R.process_beam_factorization_currents(
                    all_necessary_ME_calls_for_this_reduced_kinematics, beam_currents, self.all_MEAccessors,
                    self.accessor_tracks_leg_numbers, reduced_kinematics, ME_process, xb_1, xb_2,
                    xi1, xi2, mu_r, mu_f1, mu_f2, total_incoming_momentum, momenta_dict = counterterm.momenta_dict,
                    sector_info = sector, n_initial_legs = self.n_initial)
                #logger1.info('all_necessary_ME_calls_for_these_beam_currents : ' + str(all_necessary_ME_calls_for_these_beam_currents))    

                for reduced_kinematics_identifier, (reduced_kinematics, necessary_ME_calls) in \
                                                                all_necessary_ME_calls_for_these_beam_currents.items():

                    cut_weight = 1.

                    # Now perform the combination of the list of spin- and color- correlators to be merged
                    # for each necessary ME call identified
                    necessary_ME_calls = ME7Integrand_R.merge_correlators_in_necessary_ME_calls(necessary_ME_calls)

                    if len(necessary_ME_calls) == 0:
                        reduced_kinematics_identifier = 'IS_CUT'

                    # Make sure to skip this configuration if zero because counterterm is cut, unless one must always generate
                    # an event
                    if reduced_kinematics_identifier=='IS_CUT':
                        if not always_generate_event:
                            continue
                        else:
                            # Force the even to have zero weight then but still generate it
                            cut_weight = 0.
                            reduced_kinematics_identifier = None

                    this_base_weight = base_weight

# Gl
                    # logger1.info('reduced_kinematics_1= ' + str(reduced_kinematics))


                    # Now that all currents have been evaluated using the PS points with initial-state
                    # momenta that can be boosted because of initial-collinear mappings, we must boost
                    # back in the c.o.m frame the PS point that will be used for generating the
                    # ME7Event as well as for calling the reduced ME.
                    # First avoid possible border effects by making a copy (should be removed for performance
                    # gain, after it is checked to be safe).
                    # And now boost it back in the c.o.m frame.
                    if boost_back_to_com:
                        reduced_kinematics.boost_to_com(tuple([l.get('number') for l in counterterm.process.get_initial_legs()]))

# Gl
                    # logger1.info('reduced_kinematics_2= ' + str(reduced_kinematics))
                    #logger1.info('AFTER REDUCING_2: PS_point evaluate= ' + str(PS_point))


                    # Generate what is the kinematics (reduced_PS) returned as a list
                    # and the reduced_flavors for this counterterm by using the default reduced flavors
                    # originating from the defining process and the real-emission kinematics dictionary

                    reduced_kinematics_as_list, reduced_flavors = counterterm.get_reduced_quantities(reduced_kinematics, defining_flavors=None)
                    #logger1.info('reduced_flavors= ' + str(reduced_flavors))
                    #logger1.info('post boost sbc = ' + str((reduced_kinematics_as_list[1]+reduced_kinematics_as_list[2])))

 
                    if apply_flavour_blind_cuts and not self.pass_flavor_blind_cuts(
                            reduced_kinematics_as_list, reduced_flavors, xb_1=xb_1, xb_2=xb_2,
                            n_jets_allowed_to_be_clustered=n_unresolved_left):
                        # this configuration must be skipped
                        if not always_generate_event:
                            continue
                        else:
                            cut_weight = 0.

                    # Now apply the sectoring function if specified
                    if sector[0] is not None:
                        #print('ME7Int. - sector : ' + str(sector[0](PS_point, all_resolved_flavors[0],counterterm_index=sector[1], input_mapping_index=-1)))
                        this_base_weight *= sector[0](PS_point, all_resolved_flavors[0],
                                                      counterterm_index=sector[1], input_mapping_index=-1)

                    # Finally treat the call to the reduced connected matrix elements
                    #misc.sprint('I got for %s:'%str(counterterm.nice_string()))
                    alpha_s = self.model.get('parameter_dict')['aS']
                    mu_r = self.model.get('parameter_dict')['MU_R']

                    # Compute all the reduced flavor configurations for this counterterm
                    all_reduced_flavors = [counterterm.get_reduced_flavors(resolved_flavors)
                                           for resolved_flavors in all_resolved_flavors]

                    #logger1.info('1) all_reduced_flavors= ' + str(all_reduced_flavors))
                    # VERY IMPORTANT: We must convolve the counter-event with the initial state PDFs
                    # corresponding to the *RESOLVED* flavors, and the flavour_sensitive_cuts also
                    # applied on the resolved flavors (since for the initial-state the cuts are basically
                    # use to select particular PDF components). The correct thing to do is to overwrite
                    # the initial states flavors to *always* match their resolve counterpart.
                    # The counterpart of this is that the integrated collinear initial-state counterterms
                    # will get their initial state flavor backward evolved to also match the initial
                    # states of this local counterterm, so that the cancellation between local and integrated
                    # collinear counterterm is maintained. Of coure this plays no role for purely final-state
                    # counterterms.
                    all_reduced_flavored_with_initial_states_substituted = []
                    for i_config, reduced_flavors in enumerate(all_reduced_flavors):
                        all_reduced_flavored_with_initial_states_substituted.append(
                            (all_resolved_flavors[i_config][0], reduced_flavors[1]))


                    # Certain identical reduced flavor combinations can arise multiple times from different resolved ones,
                    # for instance the process e+ e- > u d u~ d~ has four mapped configurations:
                    #    -11 11 -> 2 1 -2 -1
                    #    -11 11 -> 2 3 -2 -3
                    #    -11 11 -> 4 1 -4 -1
                    #    -11 11 -> 4 3 -4 -3
                    # Which however yields only two mapped flavour configurations, each appearing twice:
                    #    -11 11 -> 1 -1 21
                    #    -11 11 -> 3 -3 21
                    # It is important then to keep this multiplicity factor
                    all_unique_reduced_flavored_with_initial_states_substituted = {}
                    for fc in all_reduced_flavored_with_initial_states_substituted:
                        if fc in all_unique_reduced_flavored_with_initial_states_substituted:
                            all_unique_reduced_flavored_with_initial_states_substituted[fc] += 1.0
                        else:
                            all_unique_reduced_flavored_with_initial_states_substituted[fc] = 1.0

                    # Now the phase-space point stored in the event generated is not a dictionary but
                    # a LorentzVectorList which must be ordered exactly like the flavor configurations
                    # in it are.
                    event_PS = reduced_kinematics.to_list(ordered_keys=[l.get('number') for l in counterterm.process.get('legs')])
                    #logger1.info('event_PS ordered according flavour configuration : ' + str(event_PS))

                    template_event = ME7Event(
                        event_PS,
                        {fc: cut_weight * this_base_weight * multiplicity
                         for fc, multiplicity in all_unique_reduced_flavored_with_initial_states_substituted.items()},
                        requires_mirroring=is_process_mirrored,
                        host_contribution_definition=self.contribution_definition,
                        counterterm_structure=(counterterm, all_resolved_flavors, reduced_kinematics_identifier),
                        Bjorken_xs=(xb_1, xb_2)
                    )

                    CT_event = ME7Integrand_R.generate_event_for_counterterm(
                        template_event,
                        disconnected_currents_weight,
                        counterterm.prefactor,
                        necessary_ME_calls,
                        ME_process,
                        reduced_kinematics_as_list,
                        alpha_s, mu_r,
                        self.all_MEAccessors,
                        always_generate_event = always_generate_event
                    )
                    if CT_event is not None:
                        all_events.append(CT_event)
        # print('ME7 - all events: ' + str(all_events))
        return all_events

    @classmethod
    def generate_all_necessary_ME_calls(cls, currents, ME_process, PS_point, all_MEAccessors, track_leg_numbers, **opts):
        """ Generates a list of counterterm events using the following specified inputs:
                > List of non-beam-factorisation currents
                > Reduced process
                > Input kinematics for the current
        """

        # all_necessary_ME_calls is a list inputs to call the Matrix Element and the weights
        # that we must multiply/convolve them with. We start with empty entries.
        all_necessary_ME_calls = {
            None:  # This key corresponds to the reduced_kinematics identifier (only esthetics for tagging the event)
                (None,  # This first tuple entry will be the specific reduced kinematics to consider
                 [
                     {'spin_correlations': [],
                      'color_correlations': [],
                      'main_weights': [],
                      'flavor_matrix': [],
                      'Bjorken_rescalings_beam_one': [],
                      'Bjorken_rescalings_beam_two': [],
                      },
                 ]
                 )
            }
        disconnected_currents_weight = base_objects.EpsilonExpansion({'finite': 1.0})

        # Now evaluate the mapping currents identified
        for current in currents:
            current_evaluation, all_current_results = all_MEAccessors(
                current,
                track_leg_numbers=track_leg_numbers,
                higher_PS_point=PS_point,
                reduced_process=ME_process, hel_config=None,
                **opts)

            # Now loop over all spin- and color- correlators required for this current
            # and update the necessary calls to the ME
            if not current.does_resolve_mother_spin_and_color():
                # Make sure no spin- or color-correlations were produced by the current
                assert (current_evaluation['spin_correlations'] == [None, ])
                assert (current_evaluation['color_correlations'] == [None, ])
                assert (current_evaluation['reduced_kinematics'] == [None, ])
                assert (current_evaluation['Bjorken_rescalings'] == [(None, None),])
                assert (current_evaluation['values'].keys() == [(0, 0, 0, 0), ])
                # Note: this can only work for local 4D subtraction counterterms!
                # For the integrated ones it is very likely that we cannot use a nested structure,
                # and there will be only one counterterm node level in this case anyway,
                disconnected_currents_weight *= \
                    base_objects.EpsilonExpansion(current_evaluation['values'][(0, 0, 0, 0)])
            else:
                all_necessary_ME_calls = ME7Integrand_R.update_all_necessary_ME_calls(
                    all_necessary_ME_calls, current_evaluation, weight_type='main_weight')

        # Check if any current produced a particular reduced kinematics.
        # Whenever a specific reduced kinematic identifier is present, a corresponding kinematic configuration
        # should be present, so this needs not be checked.
        if None in all_necessary_ME_calls and all_necessary_ME_calls[None][0] is None:
            all_necessary_ME_calls[None] = (PS_point , all_necessary_ME_calls[None][1])


        return all_necessary_ME_calls, disconnected_currents_weight

    @classmethod
    def generate_event_for_counterterm(cls, template_MEEvent, disconnected_currents_weight,
            overall_prefactor, all_necessary_ME_calls, ME_process, ME_PS, alpha_s, mu_r, all_MEAccessors,
            always_generate_event=False):


        # Initialize the placeholder which stores the event that we will progressively build
        # below.
        ME7_event_to_return = None
        
        for ME_call in all_necessary_ME_calls:

            # print('ME_call : ' + str(ME_call))
            # Immediately skip this contribution if it does not pass the boundary check
            # This comes from the change of variable xb_i' = xb_i * xi_i
            # Remember that Bjorken_rescaling_beam_X stored in ME_call refers to 1/xi_i
            xb_1 = template_MEEvent.Bjorken_xs[0]
            xb_2 = template_MEEvent.Bjorken_xs[1]
            if ((xb_1 is not None) and xb_1 > 1./ME_call['Bjorken_rescaling_beam_one']) or \
                ((xb_2 is not None) and xb_2 > 1./ME_call['Bjorken_rescaling_beam_two']):
                if not always_generate_event:
                    continue
                else:
                    # Set disconnected currents weight to zero in order to disable the actual
                    # contribution from this even while still producing one
                    disconnected_currents_weight = 0.

            color_correlators = tuple(ME_call['color_correlation']) if ME_call['color_correlation'] else None
            spin_correlators = tuple(ME_call['spin_correlation']) if ME_call['spin_correlation'] else None
            connected_currents_weight = ME_call['main_weight']
            #misc.sprint(ME_PS, ME_process.nice_string())
            #misc.sprint(ME_process.get_cached_initial_final_numbers())
            try:
                ME_evaluation, all_ME_results = all_MEAccessors(
                    ME_process, ME_PS, alpha_s, mu_r,
                    # Let's worry about the squared orders later,
                    # we will probably directly fish them out from the ME_process,
                    # since they should be set to a unique combination at this stage.
                    squared_orders    = None,
                    color_correlation = color_correlators,
                    spin_correlation  = spin_correlators,
                    hel_config        = None
                )

            except MadGraph5Error as e:
                logger.critical("""
A reduced matrix element is missing in the library of automatically generated matrix elements.
This is typically what can happen when your process definition is not inclusive over all IR sensitive particles.
Make sure that your process definition is specified using the relevant multiparticle labels (typically 'p' and 'j').
Also make sure that there is no coupling order specification which receives corrections.
The missing process is: %s"""%ME_process.nice_string())
                raise e

            # Multiply the various pieces building the event weight
            # (most are EpsilonExpansion's)
            event_weight = base_objects.EpsilonExpansion(ME_evaluation)
            event_weight *= disconnected_currents_weight
            event_weight *= connected_currents_weight
            event_weight *= overall_prefactor
            # misc.sprint("")
            # misc.sprint(spin_correlators, color_correlators)
            # misc.sprint(base_objects.EpsilonExpansion(ME_evaluation))
            # misc.sprint(disconnected_currents_weight)
            # misc.sprint(connected_currents_weight)
            # misc.sprint(event_weight)

            # Skip an event with no contribution (some dipoles of the eikonal for example)
            if not always_generate_event and event_weight.norm() == 0.:
                continue

            # Now build the event from the template provided
            event_to_convolve  = template_MEEvent.get_copy()

            event_to_convolve *= event_weight

            # The PDF Bjorken x's arguments will need a 1/z rescaling due to the change
            # of variable making the + distribution act on the PDF only and on the boundary of
            # the convolution over the Bjorken x's.
            event_to_convolve.set_Bjorken_rescalings(ME_call['Bjorken_rescaling_beam_one'],
                                                     ME_call['Bjorken_rescaling_beam_two'] )

            # And convolve if necessary
            for flavor_matrix_convolution in ME_call['flavor_matrices']:
                event_to_convolve.convolve_initial_states( flavor_matrix_convolution )

            # If the convolution is with a zero flavor matrix, this event does not contribute
            if event_to_convolve.is_empty():
                continue

            # Aggregate this event with the previous one, they should always be compatible
            if ME7_event_to_return is None:
                ME7_event_to_return = event_to_convolve
            else:
                ME7_event_to_return += event_to_convolve

        # Finally return the ME7 event constructed for this counterterm
        # Notice that this can be None if all terms were zero (for example during the flavor
        # convolution).
        # logger1.info('ME7_event_to_return : ' + str(ME7_event_to_return))
        return ME7_event_to_return
    
    def sigma(self, PS_point, process_key, process, all_flavor_configurations,
                  base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts):
        """ Implementation of the short-distance cross-section for the real-emission integrand.
        Counterterms will be computed on top of the actual real-emission integrand.
        Note that the PS_point specified corresponds to a point already multiplied by both
        the bjorken x's and boosted back to the c.o.m frame."""

        all_events_generated = ME7EventList()
        sector_info = opts.get('sector_info', None)

        matrix_element_event = self.generate_matrix_element_event(
            PS_point, process_key, process, all_flavor_configurations,
            base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts)
        #print('ME7 - SIGMA - matrix_element_event : ' + str(matrix_element_event))

        # Some contributions might have not physical contributions and overloaded the above
        # so as to return None
        if matrix_element_event is not None:
            all_events_generated.append(matrix_element_event)

        for i_ct, counterterm in enumerate(self.counterterms[process_key]):

            # Skip counterterms that do not belong to this sector
            if (sector_info is not None) and (sector_info['counterterms'] is not None) and i_ct not in sector_info['counterterms']:
                continue

            if not counterterm.is_singular():
                continue

            # Select a particular counterterm as follows for debugging:
            #singular_structure = counterterm.reconstruct_complete_singular_structure().substructures[0].substructures[0]
            #if not (singular_structure.name() == 'C' and len(singular_structure.substructures) == 0):
            #    continue

            CT_events = self.evaluate_local_counterterm(
                counterterm, PS_point, base_weight, mu_r, mu_f1, mu_f2,
                xb_1, xb_2, xi1, xi2,
                process.get('has_mirror_process'), all_flavor_configurations,
                sector = (sector_info['sector'] if sector_info else None, i_ct),
                **opts)

            # The function above returns None if the counterterms is removed because of
            # flavour_blind_cuts.
            if CT_events is not None:
                all_events_generated.extend(CT_events)

        # Notice that it may be possible in some subtraction scheme to combine
        # several events having the same kinematics support. This is a further
        # optimization that we can easily implement when proven needed.
        # all_events_generated.combine_events_with_identical_kinematics()
        # (PS: the function above is a place-holder and is not implemented yet.
        return all_events_generated

    def counterterms_to_consider(self, process_key, test_options):
        """Function called by the function test_IR_limits so that it can be overloaded in
        daughter classes like the beam-soft (BS) one."""

        # Use correction_order to select CT subset
        return [ ct for ct in self.counterterms[process_key]
                 if ct.get_number_of_couplings() <= test_options['correction_order'].count('N') ]

    def test_IR_limits(self, test_options):
        """Test how well local counterterms approximate a real-emission matrix element."""

        # Apply the passed options
        seed = test_options['seed']
        if seed: random.seed(seed)

        walker_name = test_options['walker']
        walker = walkers.VirtualWalker(walker_name)

#Gl
        #logger1.info("walker name = " + str(walker))

        # First generate an underlying Born
        # Specifying None forces to use uniformly random generating variables.
        # Make sure to generate a point within the cuts if necessary:
        max_attempts = 10000
        n_attempts   = 0
        while n_attempts < max_attempts:
            n_attempts += 1
            real_emission_PS_point = None
            # Phase-space generation can fail when beam convolution is active,
            # because of the condition Bjorken_x_i < xi_i
            while real_emission_PS_point is None:
                real_emission_PS_point, jac, x1s, x2s = self.phase_space_generator.get_PS_point(None)
            xb_1, xi1 = x1s
            xb_2, xi2 = x2s
            if test_options['apply_higher_multiplicity_cuts']:
                if not self.pass_flavor_blind_cuts(
                    real_emission_PS_point,
                    self.processes_map.values()[0][0].get_cached_initial_final_pdgs(),
                    n_jets_allowed_to_be_clustered=self.contribution_definition.n_unresolved_particles,
                    xb_1=xb_1, xb_2=xb_2):
                    continue
            break
        if not n_attempts < max_attempts:
            raise MadEvent7Error(
                "Could not generate a random kinematic configuration that passes " +
                "the flavour blind cuts in less than %d attempts." % max_attempts )
        n_attempts = 0

        # Loop over processes
        all_evaluations = {}
        for process_key, (defining_process, mapped_processes) in self.processes_map.items():

            # Make sure that the selected process satisfies the selection requirements
            if not self.is_part_of_process_selection(
                    [defining_process, ] + mapped_processes,
                    selection=test_options['process']):
                continue
            logger.info("Considering %s" % defining_process.nice_string())

            if self.sectors is None or self.sectors[process_key] is None:
                all_sectors = [None, ]
            else:
                all_sectors = self.sectors[process_key]
# Ez
#            logger.info("test_IR_limits: Sector list: ")
#            for sector_info in all_sectors:
#                logger.info(" %s" %str(sector_info['sector']))

            for sector_info in all_sectors:
                if (sector_info is not None) and (self.is_sector_selected is not None) and \
                                    not self.is_sector_selected(defining_process, sector_info['sector']):
                    continue
                if sector_info is not None:
                    logger.info("Considering sector: %s" %str(sector_info['sector']))

                all_processes = [defining_process,]+mapped_processes
                all_flavor_configurations = []
                # The process mirroring is accounted for at the very end only
                for proc in all_processes:
                    initial_final_pdgs = proc.get_cached_initial_final_pdgs()
                    all_flavor_configurations.append(initial_final_pdgs)

                # Reset the seed before each process so that results from a global scan of test_IR_limits
                # can easily be reproduced for a specific limit.
                seed = test_options['seed']
                if seed: random.seed(seed)

                a_real_emission_PS_point = real_emission_PS_point.get_copy()
                a_xb_1, a_xi1 = x1s
                a_xb_2, a_xi2 = x2s
                # We should not test here the flavour sensitive cuts because counterterms
                # act on the flavor space and make an event pass the cut even if the resolved one
                # does not.
                while (test_options['apply_higher_multiplicity_cuts'] and
                       not self.pass_flavor_blind_cuts(
                           a_real_emission_PS_point,
                           defining_process.get_cached_initial_final_pdgs(),
                           n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles,
                           xb_1 = xb_1, xb_2 = xb_2 ) ):
                    n_attempts += 1
                    if n_attempts > max_attempts:
                        break
                    a_real_emission_PS_point = None
                    # Phase-space generation can fail when beam convolution is active, because of the condition
                    # Bjorken_x_i < xi_i
                    while a_real_emission_PS_point is None:
                        a_real_emission_PS_point, a_jac, a_x1s, a_x2s = self.phase_space_generator.get_PS_point(None)
                    a_xb_1, a_xi1 = a_x1s
                    a_xb_2, a_xi2 = a_x2s
                if n_attempts > max_attempts:
                    raise MadEvent7Error(
                        "Could not generate a random kinematic configuration that passes " +
                        "the flavour blind cuts in less than %d attempts." % max_attempts )
                n_attempts = 0

                # Make sure to have the PS point provided as LorentzVectorDict
                a_real_emission_PS_point = phase_space_generators.LorentzVectorDict(
                                (i+1, mom) for i, mom in enumerate(a_real_emission_PS_point) )

                counterterms_to_consider = self.counterterms_to_consider(process_key, test_options)

                if len(counterterms_to_consider)==0:
                    logger.info('No counterterms to investigate found.')
                    return True

                selected_singular_structures = []

                for limit_specifier in test_options['limits']:
                    # Select the limits to be probed interpreting limits as a regex pattern.
                    # If no match is found, then reconstruct the singular structure from the limits
                    # provided
                    # The CT identifier is irrelevant here, but we must supply the correct format for the function below
                    selected_counterterms = self.find_counterterms_matching_regexp(
                                                    [(None, ct) for ct in counterterms_to_consider], limit_specifier )
                    selected_counterterms = [ct for i_ct, ct in selected_counterterms]
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
# Ez
                logger.debug('test_IR_limits: complete singular structure before filtering: \n'+'\n'.join(
                    str(ss) for ss in selected_singular_structures ))

                # Filter the singular structures to consider so as to remove those involving
                # a beam factorization not present in this integrand
                selected_singular_structures = [
                    ss for ss in selected_singular_structures if not (
                        (1 in ss.get_beam_factorization_legs() and a_xi1 is None) or
                        (2 in ss.get_beam_factorization_legs() and a_xi2 is None))
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
                        a_real_emission_PS_point, a_xi1, a_xi2, a_xb_1, a_xb_2, sector_info=sector_info )
                    process_evaluations[str(limit)] = limit_evaluations

                process_string = defining_process.base_string()
                if defining_process.has_key('n_loops'):
                    process_string += " @ " + str(defining_process['n_loops']) + " loops"

                if sector_info is not None:
                    process_string += "|sector:%s"%str(sector_info['sector'])

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

    def scale_configuration(self, xi1, xi2, real_emission_PS_point, limit,
                          scaling_parameter, defining_process, walker, boost_back_to_com):
        """ Function used by test_IR_limits_for_limit_and_process to scale a given
        configuration and return the scaled xi's and PS_point. The definition of this
        function is useful so that it can be overloaded in the beam-soft integrands
        (BS, VS, etc...)"""

        scaled_real_PS_point = walker.approach_limit(
                       real_emission_PS_point, limit, scaling_parameter, defining_process )

        if boost_back_to_com:
            scaled_real_PS_point.boost_to_com((1,2))

        # Also scale the xi initial-state convolution parameters if the limit
        # specifies a beam factorization structure for that initial state:
        beam_factorization_legs = limit.get_beam_factorization_legs()

        # Notice that if we wanted a particular overall scaling of the counterterms,
        # we would need to correctly adjust the power of the multiplying scaling parameter.
        if 1 in beam_factorization_legs:
            scaled_xi1 = 1.-(1.-xi1)*scaling_parameter
        else:
            scaled_xi1 = xi1
        if 2 in beam_factorization_legs:
            scaled_xi2 = 1.-(1.-xi2)*scaling_parameter
        else:
            scaled_xi2 = xi2

        return scaled_real_PS_point, scaled_xi1, scaled_xi2

    def test_IR_limits_for_limit_and_process(self, test_options, walker, limit, defining_process, process_key,
        all_flavor_configurations, a_real_emission_PS_point, a_xi1, a_xi2, a_xb_1, a_xb_2, sector_info=None):
        """ Given a test_options, a process, a specific limit, a walker to approach it,
        all flavor configurations mapped to that particular process, a 'starting'
        real-emission type of PS-point from which to scale towards the limit, as well as
        its corresponding convolution variables xi1/2, perform the IR limit test that
        consists in approaching the limit and computing both the real-emission matrix
        element and all counterterms that should subtract its singularities."""
                
        logger.info("Approaching limit %s " % str(limit) )
        
        # First select the counterterms to evaluate and temporarily assign them to this
        # integrand instance so that the sigma function will run on them.
        if self.has_local_counterterms():
            local_counterterms_to_consider = [ (i_ct, ct) for i_ct, ct in enumerate(self.counterterms[process_key])
                if ct.get_number_of_couplings() <= test_options['correction_order'].count('N') ]
            if test_options['counterterms']:
                counterterm_pattern = test_options['counterterms']
                if counterterm_pattern.startswith('def'):
                    counterterm_pattern = str(limit)
                logger.debug("Local counterterms before selection")
                for i_ct, ct in local_counterterms_to_consider:
                    logger.debug("    "+str(ct))
                local_counterterms_to_consider = self.find_counterterms_matching_regexp(
                    local_counterterms_to_consider, counterterm_pattern )
            logger.debug("Selected local counterterms")
            for i_ct, ct in local_counterterms_to_consider:
                logger.debug("    "+str(ct))
            local_counterterms_to_consider = [i_ct for i_ct, ct in local_counterterms_to_consider]
        
        # We must also include the integrated counterterms as they may have a xi dependence
        if self.has_integrated_counterterms():
            integrated_counterterms_to_consider = [ (i_ct, ct) for i_ct, ct in enumerate(self.integrated_counterterms[process_key])
                if ct['integrated_counterterm'].get_number_of_couplings() <= test_options['correction_order'].count('N') ]
            if test_options['counterterms']:
                counterterm_pattern = test_options['counterterms']
                if counterterm_pattern.startswith('def'):
                    counterterm_pattern = str(limit)
                logger.debug("Integrated counterterms before selection")
                for i_ct, ct in integrated_counterterms_to_consider:
                    logger.debug("    "+str(ct['integrated_counterterm']))
                integrated_counterterms_to_consider = self.find_counterterms_matching_regexp(
                    integrated_counterterms_to_consider, counterterm_pattern, integrated_CT=True )
            logger.debug("Selected integrated counterterms")
            for i_ct, ct in integrated_counterterms_to_consider:
                logger.debug("    "+str(ct['integrated_counterterm']))
            integrated_counterterms_to_consider = [i_ct for i_ct, ct in integrated_counterterms_to_consider]
        
        # Now progressively approach the limit, using a log scale
        limit_evaluations = {}
        #n_steps = test_options['n_steps']
        #min_value = test_options['min_scaling_variable']
        #max_value = test_options['max_scaling_variable']
#gl
        n_steps = 10
        min_value = 1.e-10
        max_value = 1
        base = (min_value/max_value) ** (1./n_steps)

        ME_step_array = []  
        CTs_step_array = []
        ratio_CTs_array = []
        ratio_ME_CTs_array = []

        for step in range(n_steps+1):

            # Determine the new configuration
            #scaling_parameter = max_value * (base ** step)
#gl            
            scaling_parameter = min_value * 10**(n_steps - step)
            scaled_real_PS_point, scaled_xi1, scaled_xi2 = self.scale_configuration(
                a_xi1, a_xi2, a_real_emission_PS_point, limit, scaling_parameter,
                defining_process, walker, test_options['boost_back_to_com'])
            if test_options['apply_higher_multiplicity_cuts']:
                if not self.pass_flavor_blind_cuts( scaled_real_PS_point,
                        self.processes_map.values()[0][0].get_cached_initial_final_pdgs(),
                        n_jets_allowed_to_be_clustered = self.contribution_definition.n_unresolved_particles,
                        xb_1 = a_xb_1, xb_2 = a_xb_2 ):
                    logger.warning('Aborting prematurely since the following scaled real-emission point'+
                                  ' does not pass higher multiplicity cuts.')
                    logger.warning(str(scaled_real_PS_point))
# Ez
                    logger.warning(str(self.processes_map.values()[0][0].get_cached_initial_final_pdgs()))
                    break

            # misc.sprint('Scaled PS point: %s'%str(scaled_real_PS_point))
            # misc.sprint('Scaled Bjorken rescalings: %s %s'%(scaled_xi1, scaled_xi2))
            mu_r, mu_f1, mu_f2 = self.get_scales(scaled_real_PS_point)

            # Specify the counterterms to be considered and backup the original values
            if sector_info is None:
                input_sector_info = {'sector': None,
                                     'counterterms':None,
                                     'integrated_counterterms':None}
            else:
                input_sector_info = copy.deepcopy(sector_info)
            if self.has_local_counterterms():
                if input_sector_info['counterterms'] is None:
                    input_sector_info['counterterms'] = local_counterterms_to_consider
                else:
                    input_sector_info['counterterms'] = filter(lambda i_ct: i_ct in local_counterterms_to_consider,
                                                                                input_sector_info['counterterms'])
            if self.has_integrated_counterterms():
                if input_sector_info['integrated_counterterms'] is None:
                    input_sector_info['integrated_counterterms'] = {i_ct: None for i_ct in integrated_counterterms_to_consider}
                else:
                    for i_ct in list(input_sector_info['integrated_counterterms'].keys()):
                        if i_ct not in integrated_counterterms_to_consider:
                            del input_sector_info['integrated_counterterms'][i_ct]

            try:
                # Now call sigma in order to gather all events
                events = self.sigma(
                    scaled_real_PS_point.to_dict(), process_key, defining_process, all_flavor_configurations,
                    1.0, mu_r, mu_f1, mu_f2, a_xb_1, a_xb_2, scaled_xi1, scaled_xi2,
                    compute_poles            = False,
                    apply_flavour_blind_cuts = test_options['apply_lower_multiplicity_cuts'],
                    boost_back_to_com        = test_options['boost_back_to_com'],
                    always_generate_event    = True,
                    sector_info = input_sector_info,
                )
            except Exception as e:
                logger.critical("The following exception occurred when generating events for this integrand %s:\n%s"%(
                                                      self.get_short_name(), str(e)))
                # Forward the exception while not resetting the stack trace.
                raise

            # misc.sprint('Events generated before post-processing.')
            # misc.sprint(events)
            
            # Now post-process the events as done in __call__ of the integrand.
            
            # Now handle the process mirroring, by adding, for each ME7Event defined with
            # requires_mirroring = True, a new mirrored event with the initial states swapped,
            # as well as the Bjorken x's and rescalings and with p_z -> -p_z on all final
            # state momenta
            events.generate_mirrored_events()

            # Apply flavor blind cuts
            if test_options['apply_higher_multiplicity_cuts']:
                events.filter_with_flavor_blind_cuts(self.pass_flavor_blind_cuts,
                    defining_process.get_cached_initial_final_pdgs(),
                    n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles)
            # Select particular terms of the EpsilonExpansion terms stored as weights.
            # For the test_IR_limits, we are only interested in the finite part.
            events.select_epsilon_expansion_term(test_options['epsilon_expansion_term'])

            if self.run_card['lpp1']==self.run_card['lpp2']==1:
                pdf = None if test_options['set_PDFs_to_unity'] else self.pdf
                events.apply_PDF_convolution( self.get_pdfQ2, (pdf, pdf), (mu_f1**2, mu_f2**2) )
            # Apply the 1/xi**2 factor from the change of variable xi -> xi' / xi
            # Remember that the Bjorken rescalings are defined as 1/xi at this stage
# #gl
#             for event in events:
#                 if self.contribution_definition.torino_sub_BS:
#                     event *= event.Bjorken_x_rescalings[0]*1.0
#                 else:
#                     event *= event.Bjorken_x_rescalings[0]*event.Bjorken_x_rescalings[1]

            for event in events:
                event *= event.Bjorken_x_rescalings[0] * event.Bjorken_x_rescalings[1]

            # Make sure Bjorken-x rescalings xi_i don't matter anymore
            for event in events:
                event.set_Bjorken_rescalings(None, None)
            # Apply flavor sensitive cuts
            if test_options['apply_lower_multiplicity_cuts']:
                events.filter_flavor_configurations(self.pass_flavor_sensitive_cuts)
    
            # misc.sprint('Events generated after post-processing:')
            # misc.sprint(events)
            
            # Now store the results to be returned which will eventually be passed to
            # the IR test analyzer.

            # Initialize result
            # Contributions such as the beam soft ones (BS, VS, etc...) do not have
            # a physical ME contributions, so we initialize it here to zero.
            this_eval = { 'ME': 0. }

            # Loop over all events produced and distribute them in several entries of what
            # will be returned to the IR_analyzer
            for event in events:
                # TODO: instead of looking at the total weight, we could consider doing the
                # check for one particular flavor configuration which could be specified
                # as an extra test parameter.
                event_wgt = event.get_total_weight()
                if event.counterterm_structure is None:
                    this_eval['ME'] += event_wgt
                else:
                    event_str = event.counterterm_structure_short_string()
                    if test_options['minimal_label']:
                        event_str = event_str.split('@')[0]
                    if event_str in this_eval:
                        this_eval[event_str] += event_wgt
                    else:
                        this_eval[event_str] = event_wgt

            logger.debug('For scaling variable %.3e, weight from ME = %.16e' %( scaling_parameter, this_eval['ME'] ))
            total_CTs_wgt = 0.0
            total_absCTs_wgt = 0.0
            for CT_str, CT_weight in this_eval.items():
                if CT_str=='ME': continue
                total_CTs_wgt += CT_weight
                total_absCTs_wgt += abs(CT_weight)
                logger.debug('Weight from CT %s = %.16e' % (CT_str, CT_weight) )
                if this_eval['ME'] != 0.:
                    logger.debug('Ratio: %.16f'%( CT_weight/float(this_eval['ME']) ))

#gl
            CTs_step_array.append(total_CTs_wgt)
            ME_step_array.append(this_eval['ME'])
#gl
            if step > 0:
                if CTs_step_array[step -1] != 0.0:
                    ratio_CTs_array.append(CTs_step_array[step] / CTs_step_array[step-1])
                    logger.debug('CTs_n/CTs_n-1: %.16f'%( CTs_step_array[step] / CTs_step_array[step-1] ))
                if ME_step_array[step -1] != 0.0:
                    logger.debug('ME_n/ME_n-1: %.16f'%( ME_step_array[step] / ME_step_array[step-1] ))
                if CTs_step_array[step] != 0.0 and ME_step_array[step] != 0.0:
                    ratio_ME_CTs_array.append(CTs_step_array[step] / ME_step_array[step] )


            if step==n_steps:
                printout_func = logger.info
            else:
                printout_func = logger.debug
            if this_eval['ME'] != 0.:
                test_result = total_CTs_wgt/float(this_eval['ME'])
                logger.debug('total_CTs_wgt: %.16f'%( total_CTs_wgt )) #gl
                printout_func('%sRatio sum(CTs)/ME: %.16e%s'%(
                    misc.bcolors.RED if abs(test_result+1.0) > test_options['acceptance_threshold'] else misc.bcolors.GREEN
                    , test_result,misc.bcolors.ENDC) )
            else:
                if total_absCTs_wgt != 0.:
                    test_result = total_CTs_wgt/total_absCTs_wgt
                    logger.debug('total_CTs_wgt: %.16f'%( total_CTs_wgt )) #gl
                    printout_func('Ratio sum(CTs)/sum(absCTs): %s%.16e%s'%(
                        misc.bcolors.RED if abs(test_result) > test_options['acceptance_threshold'] else misc.bcolors.GREEN
                        , test_result, misc.bcolors.ENDC))
                else:
                    test_result = total_CTs_wgt
                    logger.debug('total_CTs_wgt: %.16f'%( total_CTs_wgt )) #gl
                    printout_func('sum(CTs): %s%.16e%s'%(
                        misc.bcolors.RED if abs(test_result) > test_options['acceptance_threshold'] else misc.bcolors.GREEN
                        , test_result, misc.bcolors.ENDC))
            limit_evaluations[scaling_parameter] = this_eval

#gl
        #logger.info('ME/CTs :' + str(ratio_ME_CTs_array))
        
        # Pad missing evaluations (typically counterterms that were evaluated outside of their active range)
        # by zeros so that it can be uniformly treated by analyze_IR_limits
        all_keys = set([])
        for entry, evaluations in limit_evaluations.items():
            all_keys |= set(evaluations.keys())
        for entry, evaluations in limit_evaluations.items():
            for key in all_keys:
                if key not in evaluations:
                    evaluations[key] = 0.0
        
        # Now return all evaluations performed for each value of the scale
        return limit_evaluations

    @staticmethod
    def analyze_IR_limit(
        evaluations, acceptance_threshold,
        title=None, def_ct=None, plot_all=True, show=True,
        filename=None, plots_suffix=None, number_of_FS_legs=None,
        string_idenfier=None):

        # chose whether to use a grid display or a figure display
        display_mode = 'grid' # any value in ['figure','grid'] is legal.
        if display_mode not in ['figure','grid']:
            raise MadEvent7Error('Display mode %s not recognized in analyze_IR_limit.'%display_mode)
        
        import matplotlib
        # If possible use Agg as it allows display-less systems to use pyplot (i.e. work over ssh)
        #try:
        #    matplotlib.use('Agg')
        #except ValueError:
        #    pass
        import matplotlib.pyplot as plt

        plot_title = True
        plot_extension = ".pdf"
        if plots_suffix:
            plot_extension = '_' + plots_suffix + plot_extension
        test_failed = False
        test_ratio  = -1.0

        # Produce a plot of all counterterms
        x_values = sorted(evaluations.keys())
        lines = evaluations[x_values[0]].keys()
        lines.sort(key=len)
##        from pprint import pprint
##        pprint(evaluations)
##        misc.sprint(lines)
        # Skip ME-def line if there is no defining ct
        plot_def = def_ct and def_ct in lines
        plot_def = False
        plot_total = len(lines) > 2 or (not plot_def)

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', size=16)
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        TOTAL_color = colors.pop(3)
        MEdef_color = colors.pop(2)

        if display_mode == 'figure':
            plot_size = (6,6)
            plt.figure(1, figsize=plot_size)
            plt.subplots_adjust(left=0.15)
        elif display_mode == 'grid':
            # We will use a 4x4 grid
            plot_size = (18,10)
            plt.figure(figsize=plot_size)
            plt.subplots_adjust(top=0.95)
            plt.subplots_adjust(bottom=0.075)
            plt.subplot(2,2,1)

        plt.gca().set_prop_cycle(color=colors)
        if plot_title and title: plt.title(title)
        GeV_pow = -2*(number_of_FS_legs-2)
        units = '[GeV${}^{' + str(GeV_pow) + '}$]'
        plt.xlabel('$\lambda$')
        plt.ylabel('Integrands ' + units)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        total               = [0., ] * len(x_values)
        total_of_abs_values = [0., ] * len(x_values)
        ME_minus_def_ct = [0., ] * len(x_values)
        line_labels = {}
        for line in lines:
            y_values = [abs(evaluations[x][line]) for x in x_values]
            for i in range(len(x_values)):
                total[i] += evaluations[x_values[i]][line]
                total_of_abs_values[i] += abs(evaluations[x_values[i]][line])
            if plot_def and (line == "ME" or line == def_ct):
                def_ct_sign = 1
                if line == def_ct:
                    def_ct_sign = (-1) ** def_ct.count("(")
                for i in range(len(x_values)):
                    ME_minus_def_ct[i] += def_ct_sign * evaluations[x_values[i]][line]
            line_label = copy.copy(line)
            # The matplotlib display of <-> is corrupted on some system
            line_label = line_label.replace('<->','[M] ')
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
        if display_mode == 'figure' and filename:
            #pass
            plt.savefig(filename + '_integrands' + plot_extension)

        if display_mode == 'figure':
            plt.figure(2, figsize=plot_size)
        elif display_mode == 'grid':
            #pass
            plt.subplot(2,2,2)

        plt.gca().set_prop_cycle(color=colors)
        if plot_title and title: plt.title(title)
        plt.xlabel('$\lambda$')
        plt.xscale('log')
        plt.grid(True)
        found_zero_ME = False
        for line in lines:
            if any(evaluations[x]["ME"]==0. for x in x_values):
                if any(total_of_abs_values[i]==0. for i in range(len(x_values))):
                    y_values = [0.0 for _ in range(len(x_values))]
                else:
                    y_values = [abs(evaluations[x][line]/total_of_abs_values[i]) for i, x in enumerate(x_values)]
                found_zero_ME = True
            else:
                y_values = [abs(evaluations[x][line]/evaluations[x]["ME"]) for x in x_values]
            if '(' in line:
                style = '--'
            else:
                style = '-'
            plt.plot(x_values, y_values, style, label=line_labels[line])
        if found_zero_ME:
            plt.ylabel('Ratio to sum of abs CT')
        else:
            plt.ylabel('Ratio to ME')
        plt.legend()
        if display_mode == 'figure' and filename:
            plt.savefig(filename + '_ratios' + plot_extension)

        # Check that the ratio of def_ct to the ME is close to -1
        string_idenfier_with_limit = string_idenfier+'<>'+def_ct
        if plot_def and not test_failed:
            def_ct_2_ME_ratio = evaluations[x_values[0]][def_ct]
            if evaluations[x_values[0]]["ME"] != 0.:
                def_ct_2_ME_ratio /= evaluations[x_values[0]]["ME"]
                foo_str = "{}: One minus the ratio of the defining CT to the ME at lambda = %s is: %s.".format(string_idenfier_with_limit)
            else:
                def_ct_2_ME_ratio /= total_of_abs_values[0]
                foo_str = "{}: One minus the ratio of the defining CT to the total abs sum at lambda = %s is: %s.".format(string_idenfier_with_limit)
            # If there is no counterterms at all then the test must pass trivially
            if def_ct_2_ME_ratio == 0.:
                def_ct_2_ME_ratio = 1.
            test_ratio = abs(def_ct_2_ME_ratio)-1.
            logger.debug(foo_str % (x_values[0], test_ratio))
            test_failed = test_ratio > acceptance_threshold
        # Check that the ratio between total and ME is close to 0
        if plot_total and not test_failed:
            total_2_ME_ratio = total[0]
            if evaluations[x_values[0]]["ME"] == 0.:
                if total_of_abs_values[0] != 0.:
                    total_2_ME_ratio /= total_of_abs_values[0]
                else:
                    total_2_ME_ratio = 0.
                foo_str = "{}: Ratio of the total to the sum of abs CT at lambda = %s is: %s.".format(string_idenfier_with_limit)
            else:
                total_2_ME_ratio /= evaluations[x_values[0]]["ME"]
                foo_str = "{}: Ratio of the total to the ME at lambda = %s is: %s.".format(string_idenfier_with_limit)
            # If there is no counterterms at all then the test must pass trivially
            if total_2_ME_ratio==1.0:
                total_2_ME_ratio = 0.
            logger.debug(foo_str % (x_values[0], total_2_ME_ratio))
            test_ratio  = abs(total_2_ME_ratio)
            test_failed = test_ratio > acceptance_threshold

        if display_mode == 'figure':
            plt.figure(3, figsize=plot_size)
        elif display_mode == 'grid':
            plt.subplot(2,2,3)
            #plt.subplot(2,2,2)

        plt.gca().set_prop_cycle(color=colors)
        if plot_title and title: plt.title(title)
        plt.xlabel('$\lambda$')
        plt.ylabel('Weighted integrands ' + units)
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
        if display_mode == 'figure' and filename:
            plt.savefig(filename + '_weighted' + plot_extension)

        if display_mode == 'grid' and filename:
            plt.savefig(filename + plot_extension)

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
                # Clean-up of the embedding overall parenthesis for the title label
                limit_str = limit
                while limit_str.startswith("(") and limit_str.endswith(",)"):
                    limit_str = limit_str[1:]
                    limit_str = limit_str[:-2]
                proc, loops = process.split("@")
                title = "$" + proc + "$"
                initial_state, final_state = proc.split('>')
                number_of_FS_legs = final_state.count(' ') - 1
                title = title.replace('~','x').replace('>','\\to').replace(' ','\\;')
                title = title.replace('+','^+').replace('-','^-')
                title += "@" + loops + " approaching " + limit_str
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
                    number_of_FS_legs=number_of_FS_legs, string_idenfier=process)
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
                        logger.info('%s%s%s'%(misc.bcolors.GREEN,"    " * 2 + str(limit) + ": PASSED",misc.bcolors.ENDC))
                    else:
                        logger.info('%s%s%s'%(misc.bcolors.RED,"    " * 2 + str(limit) + ": FAILED",misc.bcolors.ENDC))
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
            if len(self.counterterms[process_key])+len(self.integrated_counterterms[process_key])>0:
                long_res = [' | with the following local and integrated counterterms:']
            else:
                long_res = [' | with no local or integrated counterterm']
            for i_CT, CT in enumerate(self.counterterms[process_key]):
                if CT.is_singular():
                    if format==2:
                        long_res.append( '   | L%d : %s' % (i_CT, str(CT)))
                    elif format==3:
                        long_res.append( '   | L%d : %s' % (i_CT, CT.__str__(
                            print_n=True, print_pdg=True, print_state=True ) ))
                    elif format>3:
                        long_res.append(CT.nice_string("   | L%-3d : "%i_CT))
            for i_CT, CT_properties in enumerate(self.integrated_counterterms[process_key]):
                CT = CT_properties['integrated_counterterm']
                if format==2:
                    long_res.append( '   | I%d : %s' % (i_CT, str(CT)))
                elif format==3:
                    long_res.append( '   | I%d : %s' % (i_CT, CT.__str__(
                        print_n=True, print_pdg=True, print_state=True )))
                elif format==4:
                    long_res.append(CT.nice_string("   | I%-3d : "%i_CT))
                elif format>4:
                    long_res.append(CT.nice_string("   | I%-3d : "%i_CT))
                    for key, value in CT_properties.items():
                        if not key in ['integrated_counterterm', 'matching_process_key']:
                            long_res.append( '     + %s : %s'%(key, str(value)))
            res += '\n'.join(long_res)

        return res
    
    def sigma(self, PS_point, process_key, process, all_flavor_configurations, base_weight, mu_r,
                                                mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts):

        # Start with the events from the process itself and its local counterterms.
        all_events_generated = ME7Integrand_R.sigma(self, PS_point, process_key, process,
            all_flavor_configurations, base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts)

        compute_poles = opts.get('compute_poles', False)
        sector_info = opts.get('sector_info', None)

        # Add the integrated counterterms
        # Now loop over all integrated counterterms
        for i_ct, counterterm_characteristics in enumerate(self.integrated_counterterms[process_key]):

            # Skip integrated counterterms that do not belong to this sector
            selected_input_mappings = None
            if (sector_info is not None) and (sector_info['integrated_counterterms'] is not None):
                if i_ct not in sector_info['integrated_counterterms']:
                    continue
                selected_input_mappings = sector_info['integrated_counterterms'][i_ct]

            # And over all the ways in which this current PS point must be remapped to
            # account for all contributions of the integrated CT. (e.g. the integrated
            # splitting g > d d~ must be "attached" to all final state gluons appearing
            # in this virtual ME process definition.)
            for i_mapping, input_mapping in enumerate(counterterm_characteristics['input_mappings']):

                # Skip integrated counterterms that do not belong to this sector
                if (selected_input_mappings is not None) and i_mapping not in selected_input_mappings:
                    continue

                # Example of a hack below to include only soft integrated CT. Uncomment to enable.
                #ss = counterterm_characteristics['integrated_counterterm'].reconstruct_complete_singular_structure()
                #ss = ss.substructures[0].substructures[0].substructures[0]
                #integrated_current = counterterm_characteristics['integrated_counterterm'].get_integrated_current()
                #if not (ss.name()=='C' and len(ss.substructures)==0 and integrated_current['distribution_type']=='bulk'):
                #    continue

                # At NLO at least, it is OK to save a bit of time by enforcing 'compute_poles=False').
                # This will need to be re-assessed at NNLO for RV contributions.
                CT_events = self.evaluate_integrated_counterterm(
                    counterterm_characteristics, PS_point, base_weight, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2,
                    input_mapping, all_flavor_configurations,
                    hel_config    = None,
                    compute_poles = compute_poles,
                    sector        = (sector_info['sector'] if sector_info else None, i_ct, i_mapping)
                )
                for CT_event in CT_events:
                    # Attach additional information to this CT_event which plays no role in the
                    # MadNkLO construction but which may be used, in test_IR_poles for instance,
                    # for improving the printout of the event record.
                    if ('allowed_backward_evolved_flavors' in counterterm_characteristics) and \
                          not all(aef=='ALL' for aef in counterterm_characteristics['allowed_backward_evolved_flavors'].values()):
                        CT_event.store_information('beam_convolution_masks',
                                counterterm_characteristics['allowed_backward_evolved_flavors'])
                    CT_event.store_information('input_mapping',input_mapping)
                    all_events_generated.append( CT_event )

        # Notice that it may be possible in some subtraction scheme to combine
        # several events having the same kinematics support. This is a further
        # optimization that we can easily implement when proven needed.
        # all_events_generated.combine_events_with_identical_kinematics()
        # (PS: the function above is a place-holder and is not implemented yet.
        return all_events_generated


class ME7Integrand_BS(ME7Integrand_RV):
    """ Specialisation of the RV integrand for handling contribution featuring
        *correlated( convolution of the beams (i.e. 'BS', 'VS', etc... contributions).
        This is the case only for integrated counterterms with possibly disjoint structures
        having each at least one soft *BeamCurrents* that originated from a colorful
        subtraction currents scheme with mappings such as ppToOneWalker that recoils the
        soft momentum equally against both initial-state beam"""
        
    def generate_matrix_element_event(self, *args, **opts):
        """ This particular type of integrand should have no physical contribution. It is only
        a receptacle for the integrated soft counterterms that recoil equally against the
        initial-state beams."""
        return None

    def counterterms_to_consider(self, process_key, test_options):
        """Overload this function so as to return integrated counterterms since it is the only
        thing present for now."""

        return [
            ct['integrated_counterterm'] for ct in self.integrated_counterterms[process_key]
            if ct['integrated_counterterm'].get_number_of_couplings() <= test_options['correction_order'].count('N') ]
        

    def scale_configuration(self, xi1, xi2, real_emission_PS_point, limit,
                           scaling_parameter, defining_process, walker, boost_back_to_com):
        """ Function used by test_IR_limits_for_limit_and_process to scale a given
        configuration and return the scaled xi's and PS_point. We must specialize it here
        so as to leave the real_emission_PS untouched since the limit anyway correspond to
        rescaling of xi only."""

        # Leave the PS point untouched.
        scaled_real_PS_point = real_emission_PS_point.get_copy()

        # For BS convolutions the limit always involve rescaling both initial state symetrically.
        scaled_xi1 = 1.-(1.-xi1)*scaling_parameter
        scaled_xi2 = 1.-(1.-xi2)*scaling_parameter

        return scaled_real_PS_point, scaled_xi1, scaled_xi2


class ME7Integrand_VV(ME7Integrand_V):
    """ ME7Integrand for the computation of a double virtual type of contribution."""

    def compute_matrix_element_event_weight_with_I_operator(self,
            PS_point, process_key, process, all_flavor_configurations, base_weight,
                                        mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts):
        """ Computes the virtual matrix element weight as a DUMMY with correct poles (but not finite part of course)
        using Catani's I operator."""

        # We do not yet support the implementation of the I2 operator for loop-induced as this would
        # necessitate two-loop matrix elements.
        if self.contribution_definition.process_definition.get('NLO_mode').startswith('sqrvirt'):
            return ME7Integrand.compute_matrix_element_event_weight(self,
                            PS_point, process_key, process, all_flavor_configurations,
                                      base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, *args, **opts)

        # Now this function below is for now a copy of the `generate_matrix_element_event` function of the
        # base class but it should be modified in order to include the I2^(0), I2^(1) and I1**2 operators

        #TODO

        alpha_s = self.model.get('parameter_dict')['aS']

        ME_evaluation, all_results = self.all_MEAccessors(
            process, PS_point, alpha_s, mu_r, pdgs=all_flavor_configurations[0])

        event_weight = base_objects.EpsilonExpansion(ME_evaluation) * base_weight

        # As a test for now, set the poles to 2 to 5
        event_weight[-1] = 2.0
        event_weight[-2] = 3.0
        event_weight[-3] = 4.0
        event_weight[-4] = 5.0

        return event_weight

class ME7Integrand_RR(ME7Integrand_R):
    """ ME7Integrand for the computation of a double real-emission type of contribution."""
    def sigma(self, PS_point, process_key, process, all_flavor_configurations,
                                base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, *args, **opts):
        return super(ME7Integrand_RR, self).sigma(PS_point, process_key, process, all_flavor_configurations,
                                base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, *args, **opts)

class ME7Integrand_RRR(ME7Integrand_R):
    """ ME7Integrand for the computation of a double real-emission type of contribution."""
    def sigma(self, PS_point, process_key, process, all_flavor_configurations,
                                base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, *args, **opts):
        return super(ME7Integrand_RRR, self).sigma(PS_point, process_key, process, all_flavor_configurations,
                                base_weight, mu_r, mu_f1, mu_f2, xb_1, xb_2, *args, **opts)

#===============================================================================
# ME7IntegrandList
#===============================================================================
class ME7IntegrandList(base_objects.PhysicsObjectList):
    """ Container for ME7Integrnds."""
    
        
    integrands_natural_order = [
        ('LO',    (ME7Integrand_B,) ),
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
        return base_objects.ContributionDefinitionList.contrib_list_string(self, format=format)

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
                            'Virtual': ME7Integrand_V,
                            'DoubleVirtual': ME7Integrand_VV,
                            'SingleReals': ME7Integrand_R,
                            'RealVirtual': ME7Integrand_RV,
                            'BeamSoft': ME7Integrand_BS,
                            'DoubleReals': ME7Integrand_RR,
                            'TripleReals': ME7Integrand_RRR,
                            'Unknown': None}

