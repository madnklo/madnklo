################################################################################
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
################################################################################
from __future__ import division

import madgraph.iolibs.template_files.subtraction.\
                      subtraction_current_implementations_utils as current_utils

import os
import unittest
import copy
import sys
import tests.IOTests as IOTests
import madgraph.core.subtraction as subtraction
import madgraph.core.base_objects as base_objects
import madgraph.various.misc as misc
import models.import_ufo as import_ufo
import madgraph.core.contributions as contributions
import madgraph.core.accessors as accessors
import madgraph.integrator.walkers as walkers
import madgraph.integrator.phase_space_generators as phase_space_generators
import madgraph.iolibs.template_files.subtraction.QCD_local_currents as QCD_local_currents

from madgraph import MG5DIR
from bidict import bidict

FINAL   = subtraction.SubtractionLeg.FINAL
INITIAL = subtraction.SubtractionLeg.INITIAL

pjoin = os.path.join

class SubtractionCurrentTest(unittest.TestCase):
    """Test and run currents, comparing against target results."""
   
    def setUp(self):
        """ Setup the model and the overhead common to all tests  """
        
        model_with_params_set = import_ufo.import_model(
                pjoin(MG5DIR,'models','loop_sm'), prefix=True,
                complex_mass_scheme = False )
        model_with_params_set.pass_particles_name_in_mg_default()
        model_with_params_set.set_parameters_and_couplings(
                param_card = pjoin(MG5DIR,'models','loop_sm','restrict_default.dat'),
                complex_mass_scheme=False )

        self.model = model_with_params_set
        self.current_exporter = subtraction.SubtractionCurrentExporter(
            self.model, export_dir=None, current_set='colorful')
        
        self.walker = walkers.FinalRescalingNLOWalker
        
        legs = base_objects.LegList([
            base_objects.Leg(
                    {'id': 1, 'state': base_objects.Leg.INITIAL, 'number': 1}),
            base_objects.Leg(
                    {'id': -1, 'state': base_objects.Leg.INITIAL, 'number': 2}),
            base_objects.Leg(
                    {'id': 22, 'state': base_objects.Leg.FINAL, 'number': 3}),
            base_objects.Leg(
                    {'id': 1,  'state': base_objects.Leg.FINAL, 'number': 4}),
            base_objects.Leg(
                    {'id': -1, 'state': base_objects.Leg.FINAL, 'number': 5}),
            base_objects.Leg(
                    {'id': 21, 'state': base_objects.Leg.FINAL, 'number': 6}),
            base_objects.Leg(
                    {'id': 21, 'state': base_objects.Leg.FINAL, 'number': 7}),
        ])
        
        self.reduced_process = base_objects.Process({
            'legs': legs,
            'model': self.model
        })
        
    def add_currents_to_accessor(self, currents, ME_accessor_dict):
        """ Add the currents in argument to the specified ME_accessor_dict."""

        mapped_currents = self.current_exporter.export(currents)

        for (module_path, class_name, _), current_properties in mapped_currents.items():
            ME_accessor_dict.add_MEAccessor(accessors.VirtualMEAccessor(
                current_properties['defining_current'],
                module_path,
                class_name,
                'subtraction.subtraction_current_implementations_utils', 
                current_properties['instantiation_options'], 
                mapped_process_keys=current_properties['mapped_process_keys'], 
                root_path=pjoin(MG5DIR,'madgraph','iolibs','template_files'),
                model=self.model
            ))
    
    def generate_PS_point(self,n_initial, n_final, initial_masses=None, final_masses=None):
        """ Generate a PS point with specified number of initial and final legs, with
        possibly a masses specified."""
        
        if initial_masses is None:
            initial_masses = tuple([0.0]*n_initial)
        if final_masses is None:
            final_masses = tuple([0.0]*n_final)
        
        PS_generator = phase_space_generators.FlatInvertiblePhasespace(
            initial_masses,
            final_masses,
            (500.0, 500.0),
            (0, 0) )

        PS_point, wgt, xb_1, xb_2 = PS_generator.get_PS_point(None)
        
        misc.sprint('PS point:\n\n%s\n\n'%PS_point.__str__(n_initial=len(initial_masses)))

        return dict( (i, momentum) for i, momentum in enumerate(PS_point) )

    def test_NLO_FF_currents(self):
        """Test various NLO FF currents."""
        
        base_PS = self.generate_PS_point(2,6)
        Q = sum(pi for (i, pi) in base_PS.items() if i > 1)

        accessors_dict = accessors.MEAccessorDict()
 
        currents = [
            subtraction.Current({
                'n_loops'                       : 0,
                'squared_orders'                : {'QCD':2, 'QED':0},
                'resolve_mother_spin_and_color' : True,
                'singular_structure'            : subtraction.CollStructure(
                    subtraction.SubtractionLeg(4, 1,FINAL),
                    subtraction.SubtractionLeg(5,-1,FINAL)
                )
            }),
            subtraction.Current({
                'n_loops'                       : 0,
                'squared_orders'                : {'QCD':2, 'QED':0},
                'resolve_mother_spin_and_color' : True,
                'singular_structure'            : subtraction.CollStructure(
                    subtraction.SubtractionLeg(5, -1,FINAL),
                    subtraction.SoftStructure(
                        subtraction.SubtractionLeg(7, 21,FINAL)
                    )
                )
            }),
            subtraction.Current({
                'n_loops'                       : 0,
                'squared_orders'                : {'QCD':2, 'QED':0},
                'resolve_mother_spin_and_color' : True,
                'singular_structure'            : subtraction.SoftStructure(
                    subtraction.SubtractionLeg(7, 21,FINAL)
                )
            }),
        ]
        self.add_currents_to_accessor(currents, accessors_dict)

#       ------------------------------------------
#       ---- Testing g > q q~ hard collinear
#       ------------------------------------------
        misc.sprint('Testing hard collinear g > q q~:')
        n_parent = len(base_PS)
        momenta_map = bidict({i: frozenset((i, )) for i in base_PS.keys()})
        momenta_map[n_parent] = frozenset((4, 5))
        reduced_PS = copy.copy(base_PS)
        pC = reduced_PS.pop(4) + reduced_PS.pop(5)
        reduced_PS[n_parent] = pC
        Q = pC + base_PS[6]
        # Put the mapped momentum on-shell (this is not a well-defined mapping,
        # but it is sufficient for now to test this current)
        reduced_PS[n_parent].set_square(0)

#       This would be a more generic way of doing this, but it would involve
#       instantiating a counterterm, which I would like to avoid for now.
#        self.walker.walk_to_lower_multiplicity(
#            a_PS, counterterm, kinematic_variables = False)
        current_evaluation, all_current_results = accessors_dict(
            currents[0], higher_PS_point=base_PS, lower_PS_point=reduced_PS,
            leg_numbers_map=momenta_map, reduced_process=None,
            hel_config=None, Q=Q)

        misc.sprint(current_evaluation)
        misc.sprint(all_current_results)

#       ------------------------------------------
#       ---- Testing q~ > q~ g  soft collinear
#       ------------------------------------------
        misc.sprint('Testing soft collinear q~ > q~ g:')
        n_parent = len(base_PS)
        momenta_map = bidict({i: frozenset((i, )) for i in base_PS.keys()})
        momenta_map[n_parent] = frozenset((5, 7))
        reduced_PS = copy.copy(base_PS)
        pC = reduced_PS.pop(5) + reduced_PS.pop(7)
        reduced_PS[n_parent] = pC
        Q = pC + base_PS[4] + base_PS[6]
        # Put the mapped momentum on-shell (this is not a well-defined mapping,
        # but it is sufficient for now to test this current)
        reduced_PS[n_parent].set_square(0)

        current_evaluation, all_current_results = accessors_dict(
            currents[1], higher_PS_point=base_PS, lower_PS_point=reduced_PS,
            leg_numbers_map=momenta_map, reduced_process=None,
            hel_config=None, Q=Q)

        misc.sprint(current_evaluation)
        misc.sprint(all_current_results)
        
#       ------------------------------------------
#       ---- Testing soft gluon current
#       ------------------------------------------
        misc.sprint('Testing soft gluon current:')
        momenta_map = bidict({i: frozenset((i, )) for i in base_PS.keys()})
        momenta_map[7] = frozenset((7, ))
        reduced_process = self.reduced_process.get_copy()
        # Remove the last leg that is going soft
        reduced_process.get('legs').pop(-1)
        reduced_PS = copy.copy(base_PS)
        pS = reduced_PS.pop(7)
        Q = pS + base_PS[5] + base_PS[6]

        current_evaluation, all_current_results = accessors_dict(
            currents[2], higher_PS_point=base_PS, lower_PS_point=reduced_PS,
            leg_numbers_map=momenta_map, reduced_process=reduced_process,
            hel_config=None, Q=Q)

        misc.sprint(current_evaluation)
        misc.sprint(all_current_results)
