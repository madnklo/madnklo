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

"""Unit test library for the routines of the core library related to writing
color information for diagrams.
"""

import os

import madgraph.core.base_objects as base_objects
import madgraph.loop.loop_base_objects as loop_base_objects
import madgraph.core.contributions as contributions
import madgraph.interface.master_interface as cmd
import madgraph.various.misc as misc
import madgraph.iolibs.export_ME7 as export_ME7
import models.import_ufo as import_ufo
from madgraph import MG4DIR, MG5DIR

import tests.IOTests as IOTests

pjoin = os.path.join

root_path = os.path.dirname(os.path.realpath( __file__ ))

class ME7ContributionStaticTest(IOTests.IOTestManager):
    """Test class for functionalities related to contributions which do not
    need any instantiated context (which is slow to generate at the setup stage."""

    def test_basic_permutation_functions(self):
        """ Test the functions in the virtual contributions implementation related
        to figuring out the permutation to apply to the kinematics in order to match
        the flavor orderings of the real-emission contribution reduced process."""
        
        res = contributions.Contribution_V.get_basic_permutation(
            ((11,-11),(21,2,-2)),
            ((11,-11),(21,-2,2))
        )
        self.assertDictEqual(res, {0: 0, 1: 1, 2: 2, 3: 4, 4: 3})

        res = contributions.Contribution_V.get_basic_permutation(
            ((11,-11),(21,3,-3,11,-11,2,-2)),
            ((-11,11),(21,3,-3,11,-11,2,-2))
        )
        self.assertDictEqual(res, {0:1,1:0,2:2,3:3,4:4,5:5,6:6,7:7,8:8})

        res = contributions.Contribution_V.get_basic_permutation(
            ((11,-11),(3,21,-3,-2,11,2,-11)),
            ((-11,11),(21,3,-3,11,-11,2,-2))
        )
        self.assertDictEqual(res, {0: 1, 1: 0, 2: 3, 3: 2, 4: 4, 5: 6, 6: 8, 7: 7, 8: 5})

    def test_flavor_distribution_permutation_functions(self):
        """ Test the function that determines the list of all possible permutations, i.e.
        assignations of the flavors of parent "singular" flavors to a list of flavors
        of the reduced process."""
        
        # Single flavor tests
        res = contributions.Contribution_V.distribute_parent_flavors(
            ({},{21:[2,]}),
            ((-11,11),(21,21,21))
        )
        self.assertListEqual(res,
            [{2: 2, 3: 3, 4: 4}, {2: 3, 3: 2, 4: 4}, {2: 4, 3: 2, 4: 3}])

        res = contributions.Contribution_V.distribute_parent_flavors(
            ({},{2:[2,]}),
            ((-11,11),(2,-2,2,-2))
        )
        misc.sprint(res)
        self.assertListEqual(res,
            [{2: 2, 3: 3, 4: 4}, {2: 3, 3: 2, 4: 4}, {2: 4, 3: 2, 4: 3}])

        res = contributions.Contribution_V.distribute_parent_flavors(
            ({},{21:[3,]}),
            ((-11,11),(21,21,21))
        )
        self.assertListEqual(res,
            [{2: 3, 3: 2, 4: 4}, {2: 2, 3: 3, 4: 4}, {2: 2, 3: 4, 4: 3}])        

        res = contributions.Contribution_V.distribute_parent_flavors(
            ({},{21:[2,3]}),
            ((-11,11),(21,21,21))
        )
        self.assertListEqual(res,
            [ {2: 2, 3: 3, 4: 4}, 
              {2: 2, 3: 4, 4: 3}, 
              {2: 3, 3: 2, 4: 4},
              {2: 3, 3: 4, 4: 2},
              {2: 4, 3: 2, 4: 3},
              {2: 4, 3: 3, 4: 2}  ]
        )
        
        res = contributions.Contribution_V.distribute_parent_flavors(
            ({},{21:[2,3,4]}),
            ((-11,11),(21,21,21))
        )
        self.assertListEqual(res,
            [ {2: 2, 3: 3, 4: 4}, 
              {2: 2, 3: 4, 4: 3}, 
              {2: 3, 3: 2, 4: 4},
              {2: 3, 3: 4, 4: 2},
              {2: 4, 3: 2, 4: 3},
              {2: 4, 3: 3, 4: 2}  ]
        )
    
        res = contributions.Contribution_V.distribute_parent_flavors(
            ({},{21:[2,3,4]}),
            ((-11,11),(21,21,21,21))
        )
        self.assertListEqual(res,
            [{2: 2, 3: 3, 4: 4, 5: 5}, {2: 2, 3: 3, 4: 5, 5: 4}, {2: 2, 3: 4, 4: 3, 5: 5}, 
             {2: 2, 3: 4, 4: 5, 5: 3}, {2: 2, 3: 5, 4: 3, 5: 4}, {2: 2, 3: 5, 4: 4, 5: 3}, 
             {2: 3, 3: 2, 4: 4, 5: 5}, {2: 3, 3: 2, 4: 5, 5: 4}, {2: 3, 3: 4, 4: 2, 5: 5}, 
             {2: 3, 3: 4, 4: 5, 5: 2}, {2: 3, 3: 5, 4: 2, 5: 4}, {2: 3, 3: 5, 4: 4, 5: 2}, 
             {2: 4, 3: 2, 4: 3, 5: 5}, {2: 4, 3: 2, 4: 5, 5: 3}, {2: 4, 3: 3, 4: 2, 5: 5}, 
             {2: 4, 3: 3, 4: 5, 5: 2}, {2: 4, 3: 5, 4: 2, 5: 3}, {2: 4, 3: 5, 4: 3, 5: 2}, 
             {2: 5, 3: 2, 4: 3, 5: 4}, {2: 5, 3: 2, 4: 4, 5: 3}, {2: 5, 3: 3, 4: 2, 5: 4}, 
             {2: 5, 3: 3, 4: 4, 5: 2}, {2: 5, 3: 4, 4: 2, 5: 3}, {2: 5, 3: 4, 4: 3, 5: 2}]
        )
        
        res = contributions.Contribution_V.distribute_parent_flavors(
            ({},{1:[3,]}),
            ((-11,11),(21,1,-1,1,-1))
        )
        self.assertListEqual(res,[{3: 3, 5: 5}, {3: 5, 5: 3}])

        # Multiple flavor tests
        res = contributions.Contribution_V.distribute_parent_flavors(
            ({},{1:[3,],2:[9,]}),
            ((-11,11),(21,1,-1,1,-1,2,-2,2,-2))
        )
        self.assertListEqual(res, [{9: 7, 3: 3, 5: 5, 7: 9}, {9: 9, 3: 3, 5: 5, 7: 7}, 
                          {9: 7, 3: 5, 5: 3, 7: 9}, {9: 9, 3: 5, 5: 3, 7: 7}])
        
        res = contributions.Contribution_V.distribute_parent_flavors(
            ({21:[0,]},{1:[3,],2:[9,]}),
            ((21,22),(21,1,-1,1,-1,2,-2,2,-2))
        )
        self.assertListEqual(res, 
            [ {0: 0, 9: 7, 3: 3, 5: 5, 7: 9}, {0: 0, 9: 9, 3: 3, 5: 5, 7: 7}, 
              {0: 0, 9: 7, 3: 5, 5: 3, 7: 9}, {0: 0, 9: 9, 3: 5, 5: 3, 7: 7} ]
        )
        
        res = contributions.Contribution_V.distribute_parent_flavors(
            ({21:[1,]},{1:[3,],2:[9,]}),
            ((21,21),(21,1,-1,1,-1,2,-2,2,-2))
        )
        self.assertListEqual(res,
            [{0: 1, 1: 0, 3: 3, 5: 5, 7: 9, 9: 7}, {0: 1, 1: 0, 3: 3, 5: 5, 7: 7, 9: 9}, 
             {0: 1, 1: 0, 3: 5, 5: 3, 7: 9, 9: 7}, {0: 1, 1: 0, 3: 5, 5: 3, 7: 7, 9: 9}, 
             {0: 0, 1: 1, 3: 3, 5: 5, 7: 9, 9: 7}, {0: 0, 1: 1, 3: 3, 5: 5, 7: 7, 9: 9}, 
             {0: 0, 1: 1, 3: 5, 5: 3, 7: 9, 9: 7}, {0: 0, 1: 1, 3: 5, 5: 3, 7: 7, 9: 9}]
        )
        
        res = contributions.Contribution_V.distribute_parent_flavors(
            ({},{1:[2,10],2:[5,12]}),
            ((-11,11),(1,-1,1,2,-1,2,-2,-2,1,-1,2,-2))
        )
        self.assertListEqual(res,
            [ {2: 2, 4: 10, 5: 5, 7: 12, 10: 4, 12: 7}, {2: 2, 4: 10, 5: 5, 7: 7, 10: 4, 12: 12}, 
             {2: 2, 4: 10, 5: 7, 7: 12, 10: 4, 12: 5}, {2: 2, 4: 10, 5: 7, 7: 5, 10: 4, 12: 12}, 
             {2: 2, 4: 10, 5: 12, 7: 7, 10: 4, 12: 5}, {2: 2, 4: 10, 5: 12, 7: 5, 10: 4, 12: 7}, 
             {2: 2, 4: 4, 5: 5, 7: 12, 10: 10, 12: 7}, {2: 2, 4: 4, 5: 5, 7: 7, 10: 10, 12: 12}, 
             {2: 2, 4: 4, 5: 7, 7: 12, 10: 10, 12: 5}, {2: 2, 4: 4, 5: 7, 7: 5, 10: 10, 12: 12}, 
             {2: 2, 4: 4, 5: 12, 7: 7, 10: 10, 12: 5}, {2: 2, 4: 4, 5: 12, 7: 5, 10: 10, 12: 7}, 
             {2: 4, 4: 10, 5: 5, 7: 12, 10: 2, 12: 7}, {2: 4, 4: 10, 5: 5, 7: 7, 10: 2, 12: 12}, 
             {2: 4, 4: 10, 5: 7, 7: 12, 10: 2, 12: 5}, {2: 4, 4: 10, 5: 7, 7: 5, 10: 2, 12: 12}, 
             {2: 4, 4: 10, 5: 12, 7: 7, 10: 2, 12: 5}, {2: 4, 4: 10, 5: 12, 7: 5, 10: 2, 12: 7}, 
             {2: 4, 4: 2, 5: 5, 7: 12, 10: 10, 12: 7}, {2: 4, 4: 2, 5: 5, 7: 7, 10: 10, 12: 12}, 
             {2: 4, 4: 2, 5: 7, 7: 12, 10: 10, 12: 5}, {2: 4, 4: 2, 5: 7, 7: 5, 10: 10, 12: 12}, 
             {2: 4, 4: 2, 5: 12, 7: 7, 10: 10, 12: 5}, {2: 4, 4: 2, 5: 12, 7: 5, 10: 10, 12: 7},
            {2: 10, 4: 4, 5: 5, 7: 12, 10: 2, 12: 7}, {2: 10, 4: 4, 5: 5, 7: 7, 10: 2, 12: 12}, 
            {2: 10, 4: 4, 5: 7, 7: 12, 10: 2, 12: 5}, {2: 10, 4: 4, 5: 7, 7: 5, 10: 2, 12: 12}, 
            {2: 10, 4: 4, 5: 12, 7: 7, 10: 2, 12: 5}, {2: 10, 4: 4, 5: 12, 7: 5, 10: 2, 12: 7}, 
            {2: 10, 4: 2, 5: 5, 7: 12, 10: 4, 12: 7}, {2: 10, 4: 2, 5: 5, 7: 7, 10: 4, 12: 12}, 
            {2: 10, 4: 2, 5: 7, 7: 12, 10: 4, 12: 5}, {2: 10, 4: 2, 5: 7, 7: 5, 10: 4, 12: 12}, 
            {2: 10, 4: 2, 5: 12, 7: 7, 10: 4, 12: 5}, {2: 10, 4: 2, 5: 12, 7: 5, 10: 4, 12: 7}]
        )

class ME7ContributionTest(IOTests.IOTestManager):
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
            pjoin(MG5DIR,'tests','input_files','LoopSMTest'),
            prefix=True, complex_mass_scheme = False )
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

            self.exporter = export_ME7.ME7Exporter(
                self.madgraph_cmd, False, group_subprocesses=True )
            self.exporter.pass_information_from_cmd(self.madgraph_cmd)
            self.exporter.copy_template(self.madgraph_cmd._curr_model)
            self.exporter.export(True, args=[])
            # We want to finalize and output the model for the Born, because need to
            # register its MEs in the accessor.
            self.exporter.finalize(['nojpeg'], self.madgraph_cmd.history)
            self.LO_contributions = self.madgraph_cmd._curr_contribs         
            # Add the Born ME accessors to the dictionary
            self.LO_contributions[0].add_ME_accessors(
                self.all_born_MEs_accessor, pjoin(tmp_path,'ME7ContributionTest_LO') )

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

            self.exporter = export_ME7.ME7Exporter(
                self.madgraph_cmd, False, group_subprocesses=True )
            self.exporter.pass_information_from_cmd(self.madgraph_cmd)
            self.exporter.copy_template(self.madgraph_cmd._curr_model)
            self.exporter.export(True, args=[])
            # The export above was enough to have fully functional contributions to test            
            # self.exporter.finalize(['nojpeg'], self.madgraph_cmd.history)
            self.NLO_contributions = self.madgraph_cmd._curr_contribs 

    @IOTests.createIOTest()
    def testIO_current_generation_and_access(self):
        """ target: Counterterms_R.txt
            target: Counterterms_V.txt
            target: Currents_Local.txt
            target: Currents_Integ.txt
        """

        # Test the generation of counterterms in single real-emission and virtual
        # type of contributions. Also make sure that they can be exported.

        verbose = False
        
        # Real counterterms

        # Fetch the real contribution
        real_emission_contrib = self.NLO_contributions.get_contributions_of_type(
            contributions.Contribution_R )[0]
        # Use the ME accessor dictionary with all Born MEs to filter
        # the unphysical counterterms,
        # for example those with the reduced process g g > a g
        n_CTs_before_filtering = len(
            sum(real_emission_contrib.counterterms.values(), []) )
        for CT_list in real_emission_contrib.counterterms.values():
            contributions.Contribution_R.remove_counterterms_with_no_reduced_process(
                self.all_born_MEs_accessor, CT_list )
        n_CTs_after_filtering = len(
            sum(real_emission_contrib.counterterms.values(), []) )
        n_CTs_difference = n_CTs_before_filtering - n_CTs_after_filtering
        if verbose:
            print_string = 'A total of %d local counterterms were filtered'
            print_string += 'because a reduced process did not exist.'
            misc.sprint(print_string % (n_CTs_difference))
        # Check the number of filtered local counterterms
        self.assertEqual(n_CTs_difference, 6)
        # Output all local counterterms
        counterterm_strings = [
            CT.__str__(print_n=True, print_pdg=True, print_state=True)
            for CT in sum(real_emission_contrib.counterterms.values(), []) ]
        open(pjoin(self.IOpath,'Counterterms_R.txt'),'w').write(
            "\n".join(counterterm_strings) )

        # Virtual counterterms

        # Fetch the virtual contribution
        virtual_contrib = self.NLO_contributions.get_contributions_of_type(
            contributions.Contribution_V)[0]
        # Apply the filter to integrated counterterms
        n_CTs_before_filtering = len(
            sum(virtual_contrib.integrated_counterterms.values(),[]) )
        for CT_list in virtual_contrib.integrated_counterterms.values():
            contributions.Contribution_V.remove_counterterms_with_no_reduced_process(
                self.all_born_MEs_accessor, CT_list )
        n_CTs_after_filtering = len(
            sum(virtual_contrib.integrated_counterterms.values(),[]) )
        n_CTs_difference = n_CTs_before_filtering - n_CTs_after_filtering
        if verbose:
            print_string = 'A total of %d integrated counterterms were filtered'
            print_string += 'because a reduced process did not exist.'
            misc.sprint(print_string % (n_CTs_difference))
        # Check the number of filtered integrated counterterms
        self.assertEqual(n_CTs_difference, 0)
        # Output all integrated counterterms
        counterterm_strings = [
            CT['integrated_counterterm'].__str__(
                print_n=True, print_pdg=True, print_state=True )
            for CT in sum(virtual_contrib.integrated_counterterms.values(), []) ]
        open(pjoin(self.IOpath,'Counterterms_V.txt'),'w').write(
            "\n".join(counterterm_strings) )
        # Check the number of counterterms that did not find a host contribution
        refused_cts = len(self.exporter.integrated_counterterms_refused_from_all_contribs)
        if verbose:
            print_string = 'A total of %d integrated subtraction counterterms'
            print_string += 'did not find a host contribution.'
            misc.sprint(print_string % refused_cts)
        self.assertEqual(refused_cts, 6)

        # Local currents

        # Initialize an empty accessor dictionary for the currents.
        # No currents are ignored because of potential pre-existing ones.
        accessors_dict = contributions.MEAccessorDict()
        all_local_currents = \
            real_emission_contrib.get_all_necessary_local_currents(accessors_dict)
        current_strings = [str(current) for current in all_local_currents]
        # Print all local currents
        if verbose: misc.sprint('Local currents:\n' + '\n'.join(current_strings))
        # Output all local currents
        open(pjoin(self.IOpath, 'Currents_Local.txt'), 'w').write(
            "\n".join(sorted(current_strings)) )

        # Integrated currents

        # Use the ME accessor dictionary with all Born MEs to filter what are the
        # unphysical counterterms, for example those with the reduced process g g > a g
        all_integrated_currents = virtual_contrib.get_all_necessary_integrated_currents(
            self.all_born_MEs_accessor )
        current_strings = [str(current) for current in all_integrated_currents]
        # Print all currents
        if verbose: misc.sprint('Integrated currents:\n' + '\n'.join(current_strings))
        # Output all local currents
        open(pjoin(self.IOpath, 'Currents_Integ.txt'), 'w').write(
            "\n".join(sorted(current_strings)) )

        with misc.TMP_directory(debug=False) as tmp_path:

            print_string = "A total of %d accessor keys have been generated"
            print_string += "for %s subtraction currents."

            # Reset the accessor dictionary so as to monitor only the newly added keys
            accessors_dict = contributions.MEAccessorDict()
            real_emission_contrib.add_current_accessors(
                self.mymodel, accessors_dict, tmp_path, all_local_currents )
            # Print all accessor keys
            if verbose: misc.sprint(print_string % (len(accessors_dict), "local"))
            self.assertEqual(len(accessors_dict), 31)

            # Reset the accessor dictionary so as to monitor only the newly added keys
            accessors_dict = contributions.MEAccessorDict()
            virtual_contrib.add_current_accessors(
                self.mymodel, accessors_dict, tmp_path, all_integrated_currents )
            # Print all accessor keys
            if verbose: misc.sprint(print_string % (len(accessors_dict), "integrated"))
            self.assertEqual(len(accessors_dict), 31)
