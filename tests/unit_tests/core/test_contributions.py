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

        # Setting up the process p p > h j j and its subtraction
        self.mylegs = base_objects.MultiLegList([
            base_objects.MultiLeg(
                    {'ids': [1,2,-1,-2,21], 'state': base_objects.Leg.INITIAL}),
            base_objects.MultiLeg(
                    {'ids': [1,2,-1,-2,21], 'state': base_objects.Leg.INITIAL}),
            base_objects.MultiLeg(
                    {'ids': [22], 'state': base_objects.Leg.FINAL}),
            base_objects.MultiLeg(
                    {'ids': [21,1,-1,2,-2], 'state': base_objects.Leg.FINAL})
        ])

        self.myprocdef = base_objects.ProcessDefinition({
            'legs': self.mylegs,
            'model': self.mymodel,
            'split_orders': ['QCD','QED']
        })

        # The general accessor with the Born ME registered
        self.all_born_MEs_accessor = contributions.MEAccessorDict()
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
            # We want to finalize and output the model for the Born, because need to
            # register its MEs in the accessor.
            self.exporter.finalize(['nojpeg'], self.madgraph_cmd.history)
            self.LO_contributions = self.madgraph_cmd._curr_contribs         
            # Add the Born ME accessors to the dictionary
            self.LO_contributions[0].add_ME_accessors(self.all_born_MEs_accessor, 
                                                  pjoin(tmp_path,'ME7ContributionTest_LO'))

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
        
        verbose = False
        
        # Now fetch the real virtual contribution
        real_emission_contrib = self.NLO_contributions.get_contributions_of_type(
                                                          contributions.Contribution_R)[0]
        # Use the ME accessor dictionary with all Born MEs to filter what are the 
        # unphysical counterterms, for example those with the reduced process g g > a g
        n_CTs_before_filtering = len(sum(real_emission_contrib.counterterms.values(),[]))
        for CT_list in real_emission_contrib.counterterms.values():
            contributions.Contribution_R.remove_counterterms_with_no_reduced_process(
                                                       self.all_born_MEs_accessor, CT_list)
        n_CTs_after_filtering = len(sum(real_emission_contrib.counterterms.values(),[]))
        if verbose: misc.sprint('A total of %d local subtraction counterterms filtered according to'%\
                     (n_CTs_before_filtering-n_CTs_after_filtering)+' reduced process inexistence.')
        self.assertEqual(n_CTs_before_filtering-n_CTs_after_filtering,6)
        target=sorted("""
S((4, 21, f))
C((4, 21, f),(5, -1, f))
C((1, -1, i),(4, 21, f))
C((2, 21, i),(4, 21, f))
C((2, 21, i),(5, -1, f))
C(S((4, 21, f)),(5, -1, f))
C(S((4, 21, f)),(1, -1, i))
C(S((4, 21, f)),(2, 21, i))

C((1, 2, i),(4, 2, f))
C((2, 2, i),(4, 2, f))
C((1, 2, i),(5, 2, f))
C((2, 2, i),(5, 2, f))

C((4, 2, f),(5, -2, f))

C((1, -1, i),(4, -1, f))
C((2, -1, i),(4, -1, f))
C((1, -1, i),(5, -1, f))
C((2, -1, i),(5, -1, f))

S((4, 21, f))
S((5, 21, f))
C((4, 21, f),(5, 21, f))
C((1, 2, i),(4, 21, f))
C((2, -2, i),(4, 21, f))
C((1, 2, i),(5, 21, f))
C((2, -2, i),(5, 21, f))
C(S((4, 21, f)),(5, 21, f))
C(S((5, 21, f)),(4, 21, f))
C(S((4, 21, f)),(1, 2, i))
C(S((4, 21, f)),(2, -2, i))
C(S((5, 21, f)),(1, 2, i))
C(S((5, 21, f)),(2, -2, i))

C((1, 21, i),(4, 1, f))
C((2, 21, i),(4, 1, f))
C((1, 21, i),(5, -1, f))
C((2, 21, i),(5, -1, f))

C((1, 1, i),(4, 1, f))
C((2, 1, i),(4, 1, f))
C((1, 1, i),(5, 1, f))
C((2, 1, i),(5, 1, f))

C((1, 1, i),(4, 1, f))
C((2, 2, i),(5, 2, f))

C((1, -1, i),(4, -1, f))
C((2, -2, i),(5, -2, f))

S((4, 21, f))
C((4, 21, f),(5, 1, f))
C((1, 1, i),(4, 21, f))
C((2, 21, i),(4, 21, f))
C((2, 21, i),(5, 1, f))
C(S((4, 21, f)),(5, 1, f))
C(S((4, 21, f)),(1, 1, i))
C(S((4, 21, f)),(2, 21, i))

C((4, 2, f),(5, -2, f))
C((1, 2, i),(4, 2, f))
C((2, -2, i),(5, -2, f))

C((4, 1, f),(5, -1, f))
C((1, 1, i),(4, 1, f))
C((2, -1, i),(5, -1, f))

S((4, 21, f))
C((4, 21, f),(5, -2, f))
C((1, -2, i),(4, 21, f))
C((2, 21, i),(4, 21, f))
C((2, 21, i),(5, -2, f))
C(S((4, 21, f)),(5, -2, f))
C(S((4, 21, f)),(1, -2, i))
C(S((4, 21, f)),(2, 21, i))

C((4, 1, f),(5, -1, f))

C((1, 21, i),(4, 2, f))
C((2, 21, i),(4, 2, f))
C((1, 21, i),(5, -2, f))
C((2, 21, i),(5, -2, f))

C((2, -1, i),(4, -1, f))
C((1, 2, i),(5, 2, f))

S((4, 21, f))
C((4, 21, f),(5, 2, f))
C((1, 2, i),(4, 21, f))
C((2, 21, i),(4, 21, f))
C((2, 21, i),(5, 2, f))
C(S((4, 21, f)),(5, 2, f))
C(S((4, 21, f)),(1, 2, i))
C(S((4, 21, f)),(2, 21, i))

S((4, 21, f))
S((5, 21, f))
C((4, 21, f),(5, 21, f))
C((1, 1, i),(4, 21, f))
C((2, -1, i),(4, 21, f))
C((1, 1, i),(5, 21, f))
C((2, -1, i),(5, 21, f))
C(S((4, 21, f)),(5, 21, f))
C(S((5, 21, f)),(4, 21, f))
C(S((4, 21, f)),(1, 1, i))
C(S((4, 21, f)),(2, -1, i))
C(S((5, 21, f)),(1, 1, i))
C(S((5, 21, f)),(2, -1, i))

C((1, 1, i),(4, 1, f))
C((2, -2, i),(5, -2, f))

C((1, -2, i),(4, -2, f))
C((2, -2, i),(4, -2, f))
C((1, -2, i),(5, -2, f))
C((2, -2, i),(5, -2, f))""".split('\n'))
        self.assertListEqual( sorted([CT.get_singular_structure_string(print_n=True, print_pdg=True, 
                print_state=True) for CT in sum(real_emission_contrib.counterterms.values(),[])]),
            target )
#        if verbose: misc.sprint('\n'+'\n'.join(CT.get_singular_structure_string(print_n=True, print_pdg=True, 
#                   print_state=True) for CT in sum(real_emission_contrib.counterterms.values(),[])))

        # Now fetch the virtual contribution
        virtual_contrib = self.NLO_contributions.get_contributions_of_type(
                                                          contributions.Contribution_V)[0]
        # Also apply the filter to integrated counterterms
        n_CTs_before_filtering = len(sum(virtual_contrib.integrated_counterterms.values(),[]))
        for CT_list in virtual_contrib.integrated_counterterms.values():
            contributions.Contribution_V.remove_counterterms_with_no_reduced_process(
                                                       self.all_born_MEs_accessor, CT_list)
        n_CTs_after_filtering = len(sum(virtual_contrib.integrated_counterterms.values(),[]))
        if verbose: misc.sprint('A total of %d integrated subtraction counterterms filtered according to'%\
              (n_CTs_before_filtering-n_CTs_after_filtering)+' reduced process inexistence.')        
        self.assertEqual(n_CTs_before_filtering-n_CTs_after_filtering,0)
        target=sorted("""(C((2, 21, i),(5, 2, f)),)
(S((4, 21, f)),)
(S((5, 21, f)),)
(C((4, 21, f),(5, 21, f)),)
(C((1, 2, i),(4, 21, f)),)
(C((2, -2, i),(4, 21, f)),)
(C((1, 2, i),(5, 21, f)),)
(C((2, -2, i),(5, 21, f)),)
(C(S((4, 21, f)),(5, 21, f)),)
(C(S((5, 21, f)),(4, 21, f)),)
(C(S((4, 21, f)),(1, 2, i)),)
(C(S((4, 21, f)),(2, -2, i)),)
(C(S((5, 21, f)),(1, 2, i)),)
(C(S((5, 21, f)),(2, -2, i)),)
(C((4, 2, f),(5, -2, f)),)
(C((4, 1, f),(5, -1, f)),)
(C((2, 21, i),(5, -2, f)),)
(C((1, 2, i),(4, 2, f)),)
(C((1, 1, i),(4, 1, f)),)
(S((4, 21, f)),)
(C((4, 21, f),(5, -2, f)),)
(C((1, -2, i),(4, 21, f)),)
(C((2, 21, i),(4, 21, f)),)
(C(S((4, 21, f)),(5, -2, f)),)
(C(S((4, 21, f)),(1, -2, i)),)
(C(S((4, 21, f)),(2, 21, i)),)
(C((1, 21, i),(4, 2, f)),)
(C((2, 21, i),(4, 2, f)),)
(C((1, -2, i),(4, -2, f)),)
(C((2, -2, i),(4, -2, f)),)
(C((1, -2, i),(5, -2, f)),)
(C((2, -2, i),(5, -2, f)),)
(C((1, -1, i),(4, -1, f)),)
(S((4, 21, f)),)
(C((4, 21, f),(5, -1, f)),)
(C((1, -1, i),(4, 21, f)),)
(C((2, 21, i),(4, 21, f)),)
(C(S((4, 21, f)),(5, -1, f)),)
(C(S((4, 21, f)),(1, -1, i)),)
(C(S((4, 21, f)),(2, 21, i)),)
(C((1, -1, i),(4, -1, f)),)
(C((2, -1, i),(4, -1, f)),)
(C((1, -1, i),(5, -1, f)),)
(C((2, -1, i),(5, -1, f)),)
(C((1, 21, i),(4, 1, f)),)
(C((2, 21, i),(4, 1, f)),)
(C((1, 2, i),(5, 2, f)),)
(C((1, 1, i),(4, 1, f)),)
(C((2, -2, i),(5, -2, f)),)
(S((4, 21, f)),)
(C((4, 21, f),(5, 2, f)),)
(C((1, 2, i),(4, 21, f)),)
(C((2, 21, i),(4, 21, f)),)
(C(S((4, 21, f)),(5, 2, f)),)
(C(S((4, 21, f)),(1, 2, i)),)
(C(S((4, 21, f)),(2, 21, i)),)
(C((1, 1, i),(4, 1, f)),)
(C((2, -2, i),(5, -2, f)),)
(C((1, 2, i),(4, 2, f)),)
(C((2, 2, i),(4, 2, f)),)
(C((1, 2, i),(5, 2, f)),)
(C((2, 2, i),(5, 2, f)),)
(C((1, 21, i),(5, -2, f)),)
(C((2, 21, i),(5, -2, f)),)
(C((2, -1, i),(4, -1, f)),)
(C((1, 21, i),(5, -1, f)),)
(C((2, 21, i),(5, -1, f)),)
(C((1, 1, i),(4, 1, f)),)
(C((2, 1, i),(4, 1, f)),)
(C((1, 1, i),(5, 1, f)),)
(C((2, 1, i),(5, 1, f)),)
(C((2, 2, i),(5, 2, f)),)
(S((4, 21, f)),)
(C((4, 21, f),(5, 1, f)),)
(C((1, 1, i),(4, 21, f)),)
(C((2, 21, i),(4, 21, f)),)
(C(S((4, 21, f)),(5, 1, f)),)
(C(S((4, 21, f)),(1, 1, i)),)
(C(S((4, 21, f)),(2, 21, i)),)
(C((2, -2, i),(5, -2, f)),)
(C((2, -1, i),(5, -1, f)),)
(C((2, 21, i),(5, -1, f)),)
(C((4, 2, f),(5, -2, f)),)
(C((2, 21, i),(5, 1, f)),)
(S((4, 21, f)),)
(S((5, 21, f)),)
(C((4, 21, f),(5, 21, f)),)
(C((1, 1, i),(4, 21, f)),)
(C((2, -1, i),(4, 21, f)),)
(C((1, 1, i),(5, 21, f)),)
(C((2, -1, i),(5, 21, f)),)
(C(S((4, 21, f)),(5, 21, f)),)
(C(S((5, 21, f)),(4, 21, f)),)
(C(S((4, 21, f)),(1, 1, i)),)
(C(S((4, 21, f)),(2, -1, i)),)
(C(S((5, 21, f)),(1, 1, i)),)
(C(S((5, 21, f)),(2, -1, i)),)
(C((4, 1, f),(5, -1, f)),)""".split('\n'))
        self.assertListEqual( sorted([CT['integrated_counterterm'].get_singular_structure_string(print_n=True, print_pdg=True, 
                print_state=True) for CT in sum(virtual_contrib.integrated_counterterms.values(),[])]),
            target )
#        if verbose: misc.sprint('\n'+'\n'.join(CT['integrated_counterterm'].get_singular_structure_string(print_n=True, print_pdg=True, 
#              print_state=True) for CT in sum(virtual_contrib.integrated_counterterms.values(),[])))
        if verbose: misc.sprint('...but a total of %d integrated subtraction counterterms did not find a host contribution.'%
                      len(self.exporter.integrated_counterterms_refused_from_all_contribs))
        self.assertEqual(len(self.exporter.integrated_counterterms_refused_from_all_contribs),6)

        # Initialize an empty accessor dictionary for the currents. No currents
        # is ignored because of potential pre-existing ones.
        accessors_dict = contributions.MEAccessorDict()
        all_subtraction_currents = real_emission_contrib.\
                                    get_all_necessary_subtraction_currents(accessors_dict)
        # Print all subtraction currents
        if verbose: misc.sprint('Local subtraction currents')
        if verbose: misc.sprint('\n'+'\n'.join([str(current) for current in all_subtraction_currents]))
        target = sorted("""S((0, 21, f)) @ 0 loops
C((0, -1, f),(0, 21, f)) @ 0 loops
C((0, -1, i),(0, 21, f)) @ 0 loops
C((0, 21, i),(0, 21, f)) @ 0 loops
C((0, -1, f),(0, 21, i)) @ 0 loops
C(S((0, 21, f)),(0, -1, f)) @ 0 loops
C(S((0, 21, f)),(0, -1, i)) @ 0 loops
C(S((0, 21, f)),(0, 21, i)) @ 0 loops
C((0, 2, i),(0, 2, f)) @ 0 loops
C((0, -2, f),(0, 2, f)) @ 0 loops
C((0, -1, i),(0, -1, f)) @ 0 loops
C((0, 21, f),(0, 21, f)) @ 0 loops
C((0, 2, i),(0, 21, f)) @ 0 loops
C((0, -2, i),(0, 21, f)) @ 0 loops
C(S((0, 21, f)),(0, 21, f)) @ 0 loops
C(S((0, 21, f)),(0, 2, i)) @ 0 loops
C(S((0, 21, f)),(0, -2, i)) @ 0 loops
C((0, 1, f),(0, 21, i)) @ 0 loops
C((0, 1, i),(0, 1, f)) @ 0 loops
C((0, -2, i),(0, -2, f)) @ 0 loops
C((0, 1, f),(0, 21, f)) @ 0 loops
C((0, 1, i),(0, 21, f)) @ 0 loops
C(S((0, 21, f)),(0, 1, f)) @ 0 loops
C(S((0, 21, f)),(0, 1, i)) @ 0 loops
C((0, -1, f),(0, 1, f)) @ 0 loops
C((0, -2, f),(0, 21, f)) @ 0 loops
C((0, -2, f),(0, 21, i)) @ 0 loops
C(S((0, 21, f)),(0, -2, f)) @ 0 loops
C((0, 2, f),(0, 21, i)) @ 0 loops
C((0, 2, f),(0, 21, f)) @ 0 loops
C(S((0, 21, f)),(0, 2, f)) @ 0 loops""".split('\n'))
        self.assertListEqual(sorted([str(current) for current in all_subtraction_currents]),target)
        # Use the ME accessor dictionary with all Born MEs to filter what are the 
        # unphysical counterterms, for example those with the reduced process g g > a g
        all_integrated_currents = virtual_contrib.get_all_necessary_integrated_currents(
                                                                self.all_born_MEs_accessor)
        # Print all currents
        if verbose: misc.sprint('Integrated currents')
        if verbose: misc.sprint('\n'+'\n'.join([str(current) for current in all_integrated_currents]))
        target=sorted("""[integrated] (C((0, 2, f),(0, 21, i)),) @ 0 loops
[integrated] (S((0, 21, f)),) @ 0 loops
[integrated] (C((0, 21, f),(0, 21, f)),) @ 0 loops
[integrated] (C((0, 2, i),(0, 21, f)),) @ 0 loops
[integrated] (C((0, -2, i),(0, 21, f)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, 21, f)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, 2, i)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, -2, i)),) @ 0 loops
[integrated] (C((0, -2, f),(0, 2, f)),) @ 0 loops
[integrated] (C((0, -1, f),(0, 1, f)),) @ 0 loops
[integrated] (C((0, -2, f),(0, 21, i)),) @ 0 loops
[integrated] (C((0, 2, i),(0, 2, f)),) @ 0 loops
[integrated] (C((0, 1, i),(0, 1, f)),) @ 0 loops
[integrated] (C((0, -2, f),(0, 21, f)),) @ 0 loops
[integrated] (C((0, 21, i),(0, 21, f)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, -2, f)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, 21, i)),) @ 0 loops
[integrated] (C((0, -2, i),(0, -2, f)),) @ 0 loops
[integrated] (C((0, -1, i),(0, -1, f)),) @ 0 loops
[integrated] (C((0, -1, f),(0, 21, f)),) @ 0 loops
[integrated] (C((0, -1, i),(0, 21, f)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, -1, f)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, -1, i)),) @ 0 loops
[integrated] (C((0, 1, f),(0, 21, i)),) @ 0 loops
[integrated] (C((0, 2, f),(0, 21, f)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, 2, f)),) @ 0 loops
[integrated] (C((0, -1, f),(0, 21, i)),) @ 0 loops
[integrated] (C((0, 1, f),(0, 21, f)),) @ 0 loops
[integrated] (C((0, 1, i),(0, 21, f)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, 1, f)),) @ 0 loops
[integrated] (C(S((0, 21, f)),(0, 1, i)),) @ 0 loops""".split('\n'))
        self.assertListEqual(sorted([str(current) for current in all_integrated_currents]), target)
        with misc.TMP_directory(debug=False) as tmp_path:
            # Reset the accessor dictionary so as to monitor only the newly added keys
            accessors_dict = contributions.MEAccessorDict()
            real_emission_contrib.add_current_accessors(
                          self.mymodel, accessors_dict, tmp_path, all_subtraction_currents)
            # Print all accessor keys
            #misc.sprint('\n'.join([str(key) for key in accessors_dict]))
            if verbose: misc.sprint("A total of %d accessor keys have been generated for local subtraction currents."%len(accessors_dict))
            self.assertEqual(len(accessors_dict),31)
            # Reset the accessor dictionary so as to monitor only the newly added keys
            accessors_dict = contributions.MEAccessorDict()
            virtual_contrib.add_current_accessors(
                           self.mymodel, accessors_dict, tmp_path, all_integrated_currents)
            # Print all accessor keys
            #misc.sprint('\n'.join([str(key) for key in accessors_dict]))
            if verbose: misc.sprint("A total of %d accessor keys have been generated for integrated subtraction currents."%len(accessors_dict))
            self.assertEqual(len(accessors_dict),31)
     
        
