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
import madgraph.various.misc as misc
import models.import_ufo as import_ufo
import madgraph.core.contributions as contributions
import madgraph.integrator.phase_space_generators as phase_space_generators

from madgraph import MG5DIR
from bidict import bidict

FINAL   = subtraction.SubtractionLeg.FINAL
INITIAL = subtraction.SubtractionLeg.INITIAL

pjoin = os.path.join
nice_momenta_string = phase_space_generators.VirtualPhaseSpaceGenerator.nice_momenta_string

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
                complex_mass_scheme=False)

        self.model = model_with_params_set
        self.current_exporter = subtraction.SubtractionCurrentExporter(
                                                    self.model, export_dir=None)
        
        self.mapper = phase_space_generators.VirtualWalker(
                                            map_type='FlatCollinear', model = self.model)
        
    def add_currents_to_accessor(self, currents, ME_accessor_dict):
        """ Add the currents in argument to the specified ME_accessor_dict."""

        mapped_currents = self.current_exporter.export(currents)

        for (module_path, class_name, _), current_properties in mapped_currents.items():
            ME_accessor_dict.add_MEAccessor(contributions.VirtualMEAccessor(
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
            (0,0))

        PS_point, wgt, xb_1, xb_2 = PS_generator.get_PS_point(None)
        
        misc.sprint('PS point:\n\n%s\n\n'%nice_momenta_string(PS_point))

        return dict( (i, momentum) for i, momentum in enumerate(PS_point) )

    def test_NLO_FF_currents(self):
        """Test various NLO FF currents."""
       
        accessors_dict = contributions.MEAccessorDict()
 
        currents = [
            subtraction.Current({
                'n_loops'                       : 0,
                'squared_orders'                : {'QCD':2, 'QED':0},
                'resolve_mother_spin_and_color' : True,
                'singular_structure'            : subtraction.CollStructure(
                    subtraction.SubtractionLeg(3, 1,FINAL),
                    subtraction.SubtractionLeg(4,-1,FINAL)
                )
            }),
        ]
        self.add_currents_to_accessor(currents, accessors_dict)

        

        a_PS = self.generate_PS_point(2,3)
        a_PS[7] = a_PS[3] + a_PS[4]
        
        # Put the mapped momentum onshell, this is not a well-defined
        # mapping, but it is sufficient for now to test this current.
        a_PS[7].rescaleEnergy()
        momenta_map = bidict( { 7 : frozenset((3,4)) } )        

#       This would be a more generic way of doing this, but it would involve
#       instantiating a counterterm, which I would like to avoid for now.
#        self.mapper.walk_to_lower_multiplicity(
#            a_PS, counterterm, kinematic_variables = False)

        current_evaluation, all_current_results = accessors_dict(
                currents[0], a_PS, hel_config=None, 
                reduced_process = None,
                leg_numbers_map = momenta_map,
                mapping_variables={})

        misc.sprint(current_evaluation)
        misc.sprint(all_current_results)
