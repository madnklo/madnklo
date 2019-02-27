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
""" A test suite to compare the current ME7 versus ME6."""

import itertools
import logging
import os
import shutil
import me_comparator
import madevent_comparator
import unittest
import subprocess

from madgraph import MG5DIR

pjoin = os.path.join
_file_path = os.path.dirname(os.path.realpath(__file__))
_pickle_path = pjoin(_file_path, 'input_files')


class ME6ME7Comparator(unittest.TestCase):
    """A class to compare cross-section values from ME6 and ME7."""
    
    nb_test = 0

    def compare_processes(self, my_proc_list = [], orders = {}, model = 'sm',
                        energy = 500, filename = "", pickle_file = "",
                        tolerance = 1e-06):
        """ For now it's dummy, it's not called here. """
         
        ME6_runner = me_comparator.MG5_UFO_Runner()
        ME6_runner.setup(MG5DIR,MG5DIR)
        ME6_runner.store_proc_card = True

        ME7_runner = me_comparator.MG5_UFO_Runner()
        ME7_runner.setup(MG5DIR,MG5DIR)
        ME7_runner.store_proc_card = True


        if os.path.exists(pjoin(MG5DIR,'models','paralel_test_model_%s' % model)):
            shutil.rmtree(pjoin(MG5DIR,'models','paralel_test_model_%s' % model))
        os.system('cp -rf %s %s' % (pjoin(mg5_path,'models',model) ,
                                    pjoin(MG5DIR,'models','paralel_test_model_%s' % model)))
 
        # Create and setup a comparator
        my_comp = me_comparator.MEComparator()
        my_comp.set_me_runners(ME6_runner, ME7_runner)

        # Run the actual comparison
        my_comp.run_comparison(my_proc_list,
                               ['paralel_test_model_%s' % model,  model], orders, energy)

        # Print the output
        if filename:
            my_comp.output_result(filename=filename)
        
                # Store output to a pickle file in the input_files directory
        if pickle_file:
            me_comparator.PickleRunner.store_comparison(\
                os.path.join(_pickle_path, pickle_file),
                my_comp.get_non_zero_processes(),
                my_comp.me_runners[0].model,
                my_comp.me_runners[0].name,
                my_comp.me_runners[0].orders,
                my_comp.me_runners[0].energy)

        # Assert that all process comparisons passed the tolerance cut
        my_comp.assert_processes(self, tolerance)
            
        # Do some cleanup
        my_comp.cleanup()
       
    def compare_cross_section(self, my_proc_list = [], orders = {}, model = 'sm',
                        filename = "", print_result = False, append_output=False,
                        tolerance = 1e-01, ME6_options={}, ME7_options={}):
        """ """
     
        ME6_runner = madevent_comparator.ME6Runner(**ME6_options)
        ME6_runner.setup(MG5DIR)
        ME6_runner.store_proc_card = True

        ME7_options['PS_generator'] = 'MCPS'
        
        ME7_runner_MCPS = madevent_comparator.ME7Runner(**ME7_options)
        ME7_runner_MCPS.setup(MG5DIR)
        ME7_runner_MCPS.store_proc_card = True
        
        ME7_options['PS_generator'] = 'FLATPS'
        
        ME7_runner_FLATPS = madevent_comparator.ME7Runner(**ME7_options)
        ME7_runner_FLATPS.setup(MG5DIR)
        ME7_runner_FLATPS.store_proc_card = True
        
        self.nb_test +=1      
        if os.path.exists(pjoin(MG5DIR,'models','paralel_test_model_%s' % model)):
            shutil.rmtree(pjoin(MG5DIR,'models','paralel_test_model_%s' % model))
        os.system('cp -rf %s %s' % (pjoin(MG5DIR,'models',model) ,
                                    pjoin(MG5DIR,'models','paralel_test_model_%s' % model)))
        
        # Create and setup a comparator
        my_comp = madevent_comparator.MadEventComparator(allow_no_present=True)
        my_comp.set_me_runners(ME6_runner,ME7_runner_MCPS,ME7_runner_FLATPS)

        # Run the actual comparison
        my_comp.run_comparison(my_proc_list,
                               ['paralel_test_model_%s' % model, model, model], orders)

        # Print the output
        if filename:
            if append_output:            
                mystream = open(filename, 'a')
            else:
                mystream = open(filename, 'w')
            my_comp.output_result(filename=mystream,tolerance=tolerance)
            mystream.close()
        
        if print_result:
            print my_comp.results[0]
        # Assert that all process comparisons passed the tolerance cut
        my_comp.assert_processes(self, tolerance)
            
        # Do some cleanup
        my_comp.cleanup()
        return my_comp.results
       
    def compare_cross_section_to_values( self, values, my_proc_list = [], 
                        orders = {}, model = 'sm',
                        filename = "", print_result = False, append_output=False,
                        tolerance = 1e-02, ME7_options={}):   
                
        ME7_runner = madevent_comparator.ME7Runner(**ME7_options)
        ME7_runner.setup(MG5DIR)
        ME7_runner.store_proc_card = True
        

        # Create and setup a comparator
        my_comp = madevent_comparator.MadEventComparator(allow_no_present=True)
        my_comp.set_me_runners(ME7_runner)

        # Run the actual comparison
        my_comp.run_comparison(my_proc_list,
                               [model], orders)

        # add the default value to the comparison
        my_comp.results.insert(0,values)
        my_comp.me_runners =(madevent_comparator.FakeRunner(), my_comp.me_runners[0])
       
        # Print the output
        if filename:
            if append_output:            
                mystream = open(filename, 'a')
            else:
                mystream = open(filename, 'w')
            my_comp.output_result(filename=mystream)
            mystream.close()

        # Assert that all process comparisons passed the tolerance cut
        my_comp.assert_processes(self, tolerance)
            
        # Do some cleanup
        my_comp.cleanup()
    

    ############################################################################
    # Short test for the evaluation of the cross-section
    ############################################################################

    # An example of a comparison against values recomputed live with ME6
    def test_ME7_paralel_cross_sm(self):
        """Test a short list of sm processes"""
        
        #print(opts)
        
        #if 'PS_generator' in opts:
        #    misc.sprint('hello')
        #self.PS_gen = opts['PS_generator']
        
        # Create a list of processes to check automatically                                                                                                                             
        #proc_lists = [['p p > t t~'], ['u d~ > W+ j'], ['u d~ > W+ j j']]
        proc_lists= [(['u c > h > u c e+ e- mu+ mu- $$ c u / a s d s~ d~'],0,99)]
        #proc_lists = [['u d > u d  mu+ mu-']]
        #proc_lists = [['u d~ > W+ j j']]
        #proc_lists = [['p p > t t~']]
        #proc_lists = [['u d~ > W+ j', 'u d~ > W+ j j']]
        # Store list of non-zero processes and results in file                                                                                                     
        pickle_file = os.path.join(_pickle_path, "short_ME7_parraleltest_cross_sm.pkl")
        if os.path.isfile('short_ME7_cs_sm.log'):
            os.remove('short_ME7_cs_sm.log')
        for my_proc_list, QCDorder, QEDorder in proc_lists:
            print 'Now running process(es) %s ...'%str(my_proc_list)
            self.compare_cross_section(my_proc_list,
                             orders = {'QED':QEDorder, 'QCD':QCDorder},
                             filename = "short_ME7_cs_sm.log",
                             append_output = True,
                             ME6_options={'accuracy':0.01},
                             ME7_options={'n_points':1000, 'integrator':'VEGAS3'})

    # An example of a comparison directly against hard-coded references
    def test_ME7_short_pp_ttx(self):
        """Quick test of our favourite standard candle."""
        
        my_proc_list = ['p p > t t~']
        values = {'number_of_P0': '1',
                  'cross_P0_qq_ttx': '77.98',
                  'cross_P0_gg_ttx': '539.3',
                  'total_cross_section': '6.1728e+02'}

        # Run the comparison
        self.compare_cross_section_to_values(values, my_proc_list,
                             orders = {'QED':99, 'QCD':99},
                             filename = "short_ME7_pp_ttx.log",
                             ME7_options={'n_points':1000, 'integrator':'VEGAS3'})

    # An example of a comparison directly against many hard-coded references
    def test_ME7_short_pp_wp_jets(self):
        """Quick test of our favourite standard candle."""
       
        # Notice that here only the entry 'total_cross_section' matters.
        # The rest is just additional info from ME6. It's useful to have it in
        # just so as to have the display of what process it is when outputting the results.
        my_proc_list = [(['u d~ > W+ j'], {'QED':99, 'QCD':99},
                            {'number_of_P0': '1',
                             'cross_P0_qq_wpg': '2406',
                             'total_cross_section': '2.406e+03'}),
                        (['u d~ > W+ j j'], {'QED':99, 'QCD':99},
                            {'number_of_P0': '2',
                             'cross_P0_qq_wpgg': '240.7',
                             'cross_P0_qq_wpqq': '119.9',
                             'total_cross_section': '3.606e+02'})]

        if os.path.isfile('short_ME7_pp_wp_jets.log'):
            os.remove('short_ME7_pp_wp_jets.log')
        for proc_defs, orders, reference_result in my_proc_list:
            print 'Now running process(es) %s %s ...'%(str(proc_defs), ' '.join('%s=%d'%it for it in orders.items()))            
            self.compare_cross_section_to_values(reference_result, proc_defs,
                             orders = orders,
                             filename = "short_ME7_pp_wp_jets.log",
                             append_output = True,
                             ME7_options={'n_points':5000, 'integrator':'VEGAS3'})
