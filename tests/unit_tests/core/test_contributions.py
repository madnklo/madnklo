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

"""Unit test library for the routines of the core library related to writing
color information for diagrams.
"""

import copy
import fractions
import os

import madgraph.core.base_objects as base_objects
import madgraph.loop.loop_base_objects as loop_base_objects
import madgraph.core.subtraction as sub
import madgraph.core.color_algebra as color
import madgraph.various.misc as misc
import madgraph.core.contributions as contributions
import madgraph.interface.master_interface as cmd
import madgraph.various.misc as misc
import madgraph.iolibs.export_ME7 as export_ME7
import madgraph.iolibs.save_load_object as save_load_object
import madgraph.integrator.phase_space_generators as phase_space_generators
import models.import_ufo as import_ufo
from madgraph import MG4DIR, MG5DIR

import tests.unit_tests as unittest

pjoin = os.path.join

root_path = os.path.dirname(os.path.realpath( __file__ ))

class ME7ContributionTest(unittest.TestCase):
    """Test class for functionalities related to contributions."""
    
    mymodel = loop_base_objects.LoopModel()
    mylegs = base_objects.LegList()
    myprocess = base_objects.Process()
    # To be set during setUp
    LO_contributions  = None
    NLO_contributions = None
 
    @misc.mute_logger()
    def setUp(self):

        self.mymodel = import_ufo.import_model(
            pjoin(MG5DIR,'tests','input_files','LoopSMTest'), prefix=True,
                                            complex_mass_scheme = False )
        self.mymodel.pass_particles_name_in_mg_default()
        self.mymodel.set_parameters_and_couplings( param_card = 
                pjoin(MG5DIR,'tests','input_files','LoopSMTest','restrict_default.dat'),
                scale=91.188,
                complex_mass_scheme=False)

        # Setting up the process p p > h j j and its subtraction
        self.mylegs = base_objects.MultiLegList([
            base_objects.MultiLeg(
                    {'ids': [1,2,21], 'state': base_objects.Leg.INITIAL}),
            base_objects.MultiLeg(
                    {'ids': [1,2,21], 'state': base_objects.Leg.INITIAL}),
            base_objects.MultiLeg(
                    {'ids': [22], 'state': base_objects.Leg.FINAL}),
            base_objects.MultiLeg(
                    {'ids': [21,1,-1], 'state': base_objects.Leg.FINAL})
        ])

        self.myprocdef = base_objects.ProcessDefinition({
            'legs': self.mylegs,
            'model': self.mymodel
        })

        # Generate only Born LO contributions
        with misc.TMP_directory(debug=False) as tmp_path:
            
            # Generate the output for this.
            self.madgraph_cmd = cmd.MasterCmd(main='MadGraph')
            self.madgraph_cmd._curr_model = self.mymodel
            self.madgraph_cmd._export_dir = pjoin(tmp_path,'ME7ContributionTest_LO')

            # Generate contributions
            generation_options = {'ME7_definition': True, 
                                  'diagram_filter': False, 
                                  'LO': True, 
                                  'NNLO': [], 
                                  'NNNLO': [],
                                  'optimize': False, 
                                  'NLO': [], 
                                  'loop_induced': [],
                                  'ignore_contributions' : []}

            self.madgraph_cmd.add_contributions(self.myprocdef, generation_options)
            LO_contributions = self.madgraph_cmd._curr_contribs
            LO_contributions.apply_method_to_all_contribs(
                    'generate_amplitudes', log='Generate diagrams for')

            self.exporter = export_ME7.ME7Exporter(self.madgraph_cmd, False, group_subprocesses=True)
            self.exporter.pass_information_from_cmd(self.madgraph_cmd)
            self.exporter.copy_template(self.madgraph_cmd._curr_model)
            self.exporter.export(True, args=[])
            # We don't want to finalize and output the model which is not complete.
            # The export above was enough to have fully functional contributions to test
#            self.exporter.finalize(['nojpeg'], self.madgraph_cmd.history)
            self.LO_contributions = self.madgraph_cmd._curr_contribs         


        # Generate all NLO contributions 
        with misc.TMP_directory(debug=False) as tmp_path:

            # Generate the output for this.
            self.madgraph_cmd = cmd.MasterCmd(main='MadGraph')
            self.madgraph_cmd._curr_model = self.mymodel
            self.madgraph_cmd._export_dir = pjoin(tmp_path,'ME7ContributionTest_LO')

            # Generate contributions
            generation_options = {'ME7_definition': True, 
                                  'diagram_filter': False, 
                                  'LO': True, 
                                  'NNLO': [], 
                                  'NNNLO': [],
                                  'optimize': False, 
                                  'NLO': ['QCD'], 
                                  'loop_induced': [],
                                  'ignore_contributions' : []}

            self.madgraph_cmd.add_contributions(self.myprocdef, generation_options)  

            self.madgraph_cmd._curr_contribs.apply_method_to_all_contribs(
                    'generate_amplitudes', log='Generate diagrams for')

            self.exporter = export_ME7.ME7Exporter(self.madgraph_cmd, False, group_subprocesses=True)
            self.exporter.pass_information_from_cmd(self.madgraph_cmd)
            self.exporter.copy_template(self.madgraph_cmd._curr_model)
            self.exporter.export(True, args=[])
            # The export above was enough to have fully functional contributions to test            
#            self.exporter.finalize(['nojpeg'], self.madgraph_cmd.history)
            self.NLO_contributions = self.madgraph_cmd._curr_contribs 

    def test_current_generation_and_access(self):
        """Test the generation of counterterms in single real-emission and virtual 
        type of contributions. Also make sure that they can be exported."""
        
        real_emission_contrib = self.NLO_contributions.get_contributions_of_type(
                                                          contributions.Contribution_R)[0]
        accessors_dict = contributions.MEAccessorDict()
        all_subtraction_currents = \
              real_emission_contrib.get_all_necessary_subtraction_currents(accessors_dict)
        # Print all subtraction currents
        misc.sprint('Local subtraction currents')
        misc.sprint('\n'.join([str(current) for current in all_subtraction_currents]))
        with misc.TMP_directory(debug=False) as tmp_path:
            real_emission_contrib.add_current_accessors(
                          self.mymodel, accessors_dict, tmp_path, all_subtraction_currents)
            # Print all accessor keys
            #misc.sprint('\n'.join([str(key) for key in accessors_dict]))
            misc.sprint("A total of %d accessor keys have been generated."%len(accessors_dict))
        
        virtual_contrib = self.NLO_contributions.get_contributions_of_type(
                                                          contributions.Contribution_V)[0]
        accessors_dict = contributions.MEAccessorDict()
        all_integrated_currents = virtual_contrib.get_all_necessary_integrated_currents(
                                                                            accessors_dict)
        # Print all currents
        misc.sprint('Integrated currents')
        misc.sprint('\n'.join([str(current) for current in all_integrated_currents]))
        with misc.TMP_directory(debug=False) as tmp_path:
            virtual_contrib.add_current_accessors(
                           self.mymodel, accessors_dict, tmp_path, all_integrated_currents)
            # Print all accessor keys
            #misc.sprint('\n'.join([str(key) for key in accessors_dict]))
            misc.sprint("A total of %d accessor keys have been generated."%len(accessors_dict))

     
        
