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
import subprocess
import unittest
import os
import re
import shutil
import sys
import logging
import random
import timeit

pjoin = os.path.join

logger = logging.getLogger('test_ME7')

import tests.unit_tests.iolibs.test_file_writers as test_file_writers

import madgraph.interface.master_interface as Cmd
import madgraph.interface.ME7_interface as ME7_interface
import madgraph.integrator.ME7_integrands as ME7_integrands
from madgraph.core.accessors import ProcessKey
import madgraph.various.misc as misc
_file_path = os.path.dirname(os.path.realpath(__file__))
_pickle_path =os.path.join(_file_path, 'input_files')

from madgraph import MG4DIR, MG5DIR, MadGraph5Error, InvalidCmd

debugging_warning = ' /!\ USE ONLY FOR DEBUGGING /!\ '
debugging_written = 'Output for %s written at %s'
debugging_reused  = 'Reusing output for %s written at %s'

def get_test_IR_limit_cmd(options):

    options_str = ' '.join(
        ('--%s=%s' % (key, value) if value is not None else '--%s' % key)
        for key, value in options.items()
    )
    return 'test_IR_limits ' + options_str

#===============================================================================
# TestME7 colorful output for e+ e- > j j j @NLO
#===============================================================================
class TestME7_NLO_colorful_epem_jjj(unittest.TestCase):
    """This test validates the command 'test_IR_limits' of ME7 in the colorful scheme
    as well as integrand calls for the process e+ e- > j j j --NLO=QCD"""
    
    # If the debug mode is set to True, then the process output is not refreshed
    # but reused instead
    debugging = False
    is_process_generated = False

    def setUp(self):
        """ basic building of the class to test """
        
        self.tmp_process_dir = pjoin(_file_path, 'TMP_TestME7_colorful_epem_jjj_output')
        # Generate the process output if it does not exist yet or if we
        # are not in debug mode.
        if os.path.isdir(self.tmp_process_dir):
            if not self.is_process_generated and not self.debugging:
                shutil.rmtree(self.tmp_process_dir)
            else:
                TestME7_NLO_colorful_epem_jjj.is_process_generated = True
        if not self.is_process_generated:
            self.cmd = Cmd.MasterCmd()
            if os.path.isdir(self.tmp_process_dir):
                shutil.rmtree(self.tmp_process_dir)

            # Now generate and output a process, so as to run ME7 commands on it
            self.do('import model loop_sm')
            self.do('set subtraction_currents_scheme colorful')
            self.do('set subtraction_mappings_scheme LorentzNLO')            
            self.do('generate e+ e- > j j j --NLO=QCD')
            self.do('output %s' % self.tmp_process_dir)
            TestME7_NLO_colorful_epem_jjj.is_process_generated = True
            if self.debugging:
                misc.sprint(debugging_warning)
                misc.sprint(
                    debugging_written % (self.__class__.__name__, self.tmp_process_dir))
        else:
            if self.debugging:
                misc.sprint(debugging_warning)
                misc.sprint(
                    debugging_reused % (self.__class__.__name__, self.tmp_process_dir))
            
        # Now initialize an ME7 interface on the above process output
        self.cmd = ME7_interface.MadEvent7Cmd(me_dir=self.tmp_process_dir)
        self.cmd.no_notification()
 
    def __del__(self):
        if os.path.isdir(self.tmp_process_dir) and not self.debugging:
            shutil.rmtree(self.tmp_process_dir)

    def do(self, line):
        """ exec a line in the cmd under test """
        self.cmd.exec_cmd(line)
        
    def verify_ME7_test_results(self, results_file_path):
        """ Parses and verify that all tests output in 'results_file_path' are passed."""
        
        for line in open(results_file_path,'r').read().split('\n'):
            process, limit, outcome, ratio = line.split('|')[:4]
            self.assertTrue(outcome.strip()=='PASSED', line)
        
    def test_ME7_colorful_ggqqx_collinear_limits(self):
        """Check the test of collinear limits on a particular process."""
        
        options = {'correction_order'       : 'NLO',
                   'limits'                 : 'purecollinear',
                   'counterterms'           : 'all',
                   'process'                : 'e+ e- > g g d d~',
                   'show_plots'             : False,
                   'save_plots'             : False,
                   'seed'                   : 666,
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-7,
                   'acceptance_threshold'   : 5.0e-3,
                   'save_results_to_path'   : 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir,'test_IR_limit_output_for_acceptance_test.dat'))

    def test_ME7_colorful_ggqqx_soft_limits(self):
        """Check the test of collinear limits on a particular process."""
        
        options = {'correction_order'       : 'NLO',
                   'limits'                 : 'puresoft',
                   'counterterms'           : 'all',
                   'process'                : 'e+ e- > g g d d~',
                   'show_plots'             : False,
                   'save_plots'             : False,
                   'seed'                   : 666,
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-11,
                   'acceptance_threshold'   : 1.0e-5,
                   'save_results_to_path'   : 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir,'test_IR_limit_output_for_acceptance_test.dat'))

    def test_ME7_colorful_ggqqx_softcollinear_limits(self):
        """Check the test of collinear limits on a particular process."""
        
        options = {'correction_order'       : 'NLO',
                   'limits'                 : "r'^\(C\(S.*$'",
                   'counterterms'           : 'all',
                   'process'                : 'e+ e- > g g d d~',
                   'show_plots'             : False,
                   'save_plots'             : False,
                   'seed'                   : '2',
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-7,
                   'acceptance_threshold'   : 1.0e-2,
                   'save_results_to_path'   : 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir,'test_IR_limit_output_for_acceptance_test.dat'))


    def test_ME7_colorful_qqxQQx_collinear_limits(self):
        """Check the test of collinear limits on a particular process."""
        
        main_cmd = 'test_IR_limits'

        options = {'correction_order'       : 'NLO',
                   'limits'                 : 'collinear',
                   'counterterms'           : 'all',
                   'process'                : 'e+ e- > d d~ d d~',
                   'show_plots'             : False,
                   'save_plots'             : False,
                   'seed'                   : 666,
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-7,
                   'acceptance_threshold'   : 5.0e-4,
                   'save_results_to_path'   : 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir,'test_IR_limit_output_for_acceptance_test.dat'))

    def test_ME7_born_integrand_call(self):
        """ Check the result of a single call to the born integrand_call."""
        
        born_integrand = self.cmd.all_integrands.get_integrands_of_type(
                                                          ME7_integrands.ME7Integrand_B)[0]
        
        dimensions = born_integrand.get_dimensions()
  
        def call():
            born_integrand(
                dimensions.get_continuous_dimensions().random_sample(),
                dimensions.get_discrete_dimensions().random_sample(),
                cache_active = True
            )

        n_calls = 1000
        with ME7_interface.ME7RunEnvironment( silence = True, loggers = logging.CRITICAL ):
            timing = [ 1.e3*(res/float(n_calls)) for i,res in 
                                 enumerate(timeit.repeat(call, number=n_calls, repeat=2)) ]
            res = '\n'+'\n'.join('%d : %g ms'%(i+1, res) for i,res in enumerate(timing))
            self.assertTrue((timing[0] < 50.0), 'Born integrand call too slow: %g ms'%(timing[0]))

    def test_ME7_real_integrand_call(self):
        """ Check the result of a single call to the born integrand_call."""
        
        real_integrand = self.cmd.all_integrands.get_integrands_of_type(
                                                          ME7_integrands.ME7Integrand_R)[0]
        
        dimensions = real_integrand.get_dimensions()
  
        def call():
            real_integrand(
                dimensions.get_continuous_dimensions().random_sample(),
                dimensions.get_discrete_dimensions().random_sample(),
                cache_active = False
            )

        n_calls = 100
        with ME7_interface.ME7RunEnvironment( silence = True, loggers = logging.CRITICAL ):
            timing = [ 1.e3*(res/float(n_calls)) for i,res in 
                                 enumerate(timeit.repeat(call, number=n_calls, repeat=2)) ]
            res = '\n'+'\n'.join('%d : %g ms'%(i+1, res) for i,res in enumerate(timing))
            self.assertTrue((timing[0] < 500.0), 'Real integrand call too slow: %g ms'%(timing[0]))
            
    def test_ME7_virtual_integrand_call(self):
        """ Check the result of a single call to the born integrand_call."""
        
        virtual_integrand = self.cmd.all_integrands.get_integrands_of_type(
                                                          ME7_integrands.ME7Integrand_V)[0]
        
        dimensions = virtual_integrand.get_dimensions()
  
        def call():
            virtual_integrand(
                dimensions.get_continuous_dimensions().random_sample(),
                dimensions.get_discrete_dimensions().random_sample(),
                cache_active = False
            )

        n_calls = 100
        with ME7_interface.ME7RunEnvironment( silence = True, loggers = logging.CRITICAL ):
            timing = [ 1.e3*(res/float(n_calls)) for i,res in 
                                 enumerate(timeit.repeat(call, number=n_calls, repeat=2)) ]
            res = '\n'+'\n'.join('%d : %g ms'%(i+1, res) for i,res in enumerate(timing))
            self.assertTrue((timing[0] < 500.0), 'Virtual integrand call too slow: %g ms'%(timing[0]))

#===============================================================================
# TestME7 cataniseymour output for e+ e- > j j j @NLO
#===============================================================================
class TestME7_NLO_cataniseymour_epem_jjj(unittest.TestCase):
    """This test validates the command 'test_IR_limits' of ME7 in the cataniseymour scheme
    for the process e+ e- > j j j --NLO=QCD"""
    
    # If the debug mode is set to True, then the process output is not refreshed
    # but reused instead
    debugging = False 
    is_process_generated = False

    def setUp(self):
        """ basic building of the class to test """
        
        self.tmp_process_dir = pjoin(_file_path, 'TMP_TestME7_cataniseymour_epem_jjj_output')
        # Generate the process output if it does not exist yet or if we
        # are not in debug mode.
        if os.path.isdir(self.tmp_process_dir):
            if not self.is_process_generated and not self.debugging:
                shutil.rmtree(self.tmp_process_dir)
            else:
                TestME7_NLO_cataniseymour_epem_jjj.is_process_generated = True
        if not self.is_process_generated:
            self.cmd = Cmd.MasterCmd()
            if os.path.isdir(self.tmp_process_dir):
                shutil.rmtree(self.tmp_process_dir)

            # Now generate and output a process, so as to run ME7 commands on it
            self.do('import model loop_sm')
            self.do('set subtraction_currents_scheme cataniseymour')
            self.do('set subtraction_mappings_scheme LorentzNLO')            
            self.do('generate e+ e- > j j j --NLO=QCD --ignore_contributions=V')
            self.do('output %s --ignore_integrated_counterterms=R' % self.tmp_process_dir)
            TestME7_NLO_cataniseymour_epem_jjj.is_process_generated = True
            if self.debugging:
                misc.sprint(debugging_warning)
                misc.sprint(
                    debugging_written % (self.__class__.__name__, self.tmp_process_dir))
        else:
            if self.debugging:
                misc.sprint(debugging_warning)
                misc.sprint(
                    debugging_reused % (self.__class__.__name__, self.tmp_process_dir))

        # Now initialize an ME7 interface on the above process output
        self.cmd = ME7_interface.MadEvent7Cmd(me_dir=self.tmp_process_dir)
        self.cmd.no_notification()

    def __del__(self):
        if os.path.isdir(self.tmp_process_dir) and not self.debugging:
            shutil.rmtree(self.tmp_process_dir)

    def do(self, line):
        """ exec a line in the cmd under test """
        self.cmd.exec_cmd(line)
        
    def verify_ME7_test_results(self, results_file_path):
        """ Parses and verify that all tests output in 'results_file_path' are passed."""
        
        for line in open(results_file_path,'r').read().split('\n'):
            process, limit, outcome, ratio = line.split('|')[:4]
            self.assertTrue(outcome.strip()=='PASSED', line)
        
    def test_ME7_cataniseymour_ggqqx_collinear_limits(self):
        """Check the test of collinear limits on a particular process."""
        
        options = {'correction_order'       : 'NLO',
                   'limits'                 : 'purecollinear',
                   'counterterms'           : 'all',
                   'process'                : 'e+ e- > g g d d~',
                   'show_plots'             : False,
                   'save_plots'             : False,
                   'seed'                   : 666,
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-7,
                   'acceptance_threshold'   : 5.0e-4,
                   'save_results_to_path'   : 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir,'test_IR_limit_output_for_acceptance_test.dat'))

    def test_ME7_cataniseymour_ggqqx_soft_limits(self):
        """Check the test of collinear limits on a particular process."""
        
        options = {'correction_order'       : 'NLO',
                   'limits'                 : "['S(3)','S(4)']",
                   'counterterms'           : 'all',
                   'process'                : 'e+ e- > g g d d~',
                   'show_plots'             : False,
                   'save_plots'             : False,
                   'seed'                   : 666,
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-7,
                   'acceptance_threshold'   : 5.0e-4,
                   'save_results_to_path'   : 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir,'test_IR_limit_output_for_acceptance_test.dat'))
        
    def test_ME7_cataniseymour_ggqqx_softcollinear_limits(self):
        """Check the test of collinear limits on a particular process."""
        
        options = {'correction_order'       : 'NLO',
                   'limits'                 : 
        "['C(S(3),4)','C(S(3),5)','C(S(3),6)','C(S(4),3)','C(S(4),5)','C(S(4),6)']",
                   'counterterms'           : 'all',
                   'process'                : 'e+ e- > g g d d~',
                   'show_plots'             : False,
                   'save_plots'             : False,
                   'seed'                   : 666,
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-7,
                   'acceptance_threshold'   : 5.0e-4,
                   'save_results_to_path'   : 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir,'test_IR_limit_output_for_acceptance_test.dat'))

    def test_ME7_cataniseymour_qqxQQx_collinear_limits(self):
        """Check the test of collinear limits on a particular process."""
        
        options = {'correction_order'       : 'NLO',
                   'limits'                 : 'collinear',
                   'counterterms'           : 'all',
                   'process'                : 'e+ e- > d d~ d d~',
                   'show_plots'             : False,
                   'save_plots'             : False,
                   'seed'                   : 666,
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-7,
                   'acceptance_threshold'   : 5.0e-4,
                   'save_results_to_path'   : 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir,'test_IR_limit_output_for_acceptance_test.dat'))

#===============================================================================
# TestME7 colorful output for p p > j j @NLO
#===============================================================================
class TestME7_NLO_colorful_pp_jj(unittest.TestCase):
    """This test validates the command 'test_IR_limits' of ME7 in the colorful scheme
    as well as integrand calls for the process p p > j j --NLO=QCD"""

    # If the debug mode is set to True, then the process output is not refreshed
    # but reused instead
    debugging = False
    is_process_generated = False

    def setUp(self):
        """ basic building of the class to test """

        self.tmp_process_dir = pjoin(_file_path, 'TMP_TestME7_colorful_pp_jj_output')
        # Generate the process output if it does not exist yet or if we
        # are not in debug mode.
        if os.path.isdir(self.tmp_process_dir):
            if not self.is_process_generated and not self.debugging:
                shutil.rmtree(self.tmp_process_dir)
            else:
                TestME7_NLO_colorful_pp_jj.is_process_generated = True
        if not self.is_process_generated:
            self.cmd = Cmd.MasterCmd()
            if os.path.isdir(self.tmp_process_dir):
                shutil.rmtree(self.tmp_process_dir)

            # Now generate and output a process, so as to run ME7 commands on it
            self.do('import model loop_sm')
            self.do('set subtraction_currents_scheme colorful')
            self.do('set subtraction_mappings_scheme LorentzNLO')
            self.do('generate p p > j j --NLO=QCD --ignore_contributions=V')
            self.do('output %s --ignore_integrated_counterterms=R' % self.tmp_process_dir)
            TestME7_NLO_colorful_pp_jj.is_process_generated = True
            if self.debugging:
                misc.sprint(debugging_warning)
                misc.sprint(
                    debugging_written % (self.__class__.__name__, self.tmp_process_dir))
        else:
            if self.debugging:
                misc.sprint(debugging_warning)
                misc.sprint(
                    debugging_reused % (self.__class__.__name__, self.tmp_process_dir))

        # Now initialize an ME7 interface on the above process output
        self.cmd = ME7_interface.MadEvent7Cmd(me_dir=self.tmp_process_dir)
        self.cmd.no_notification()

    def __del__(self):
        if os.path.isdir(self.tmp_process_dir) and not self.debugging:
            shutil.rmtree(self.tmp_process_dir)

    def do(self, line):
        """ exec a line in the cmd under test """
        self.cmd.exec_cmd(line)

    def verify_ME7_test_results(self, results_file_path):
        """ Parses and verify that all tests output in 'results_file_path' are passed."""

        for line in open(results_file_path, 'r').read().split('\n'):
            process, limit, outcome, ratio = line.split('|')[:4]
            self.assertTrue(outcome.strip() == 'PASSED', line)

    def test_ME7_colorful_gq_ggq_collinear_limits(self):
        """Check the test of collinear limits on a particular process."""

        options = {'correction_order': 'NLO',
                   'counterterms': 'def',
                   'process': 'g u > g g u',
                   'show_plots': False,
                   'save_plots': False,
                   'seed': 666,
                   'n_steps': 10,
                   'min_scaling_variable': 1.0e-7,
                   'acceptance_threshold': 5.0e-3,
                   'save_results_to_path': 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        options['limits'] = 'C(1,4)'
        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(
            pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))
        options['limits'] = 'C(2,4)'
        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(
            pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))
        options['limits'] = 'C(1,5)'
        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(
            pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))
        options['limits'] = 'C(2,5)'
        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(
            pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))

#===============================================================================
# TestME7 cataniseymour output for e+ e- > g u u~ @NNLO
#===============================================================================
class TestME7_NNLO_colorful_epem_guux(unittest.TestCase):
    """This test validates the command 'test_IR_limits' of ME7 in the colorful scheme
    as well as integrand calls for the process e+ e- > g u u~ --NNLO=QCD"""

    # If the debug mode is set to True, then the process output is not refreshed
    # but reused instead
    debugging = False 
    is_process_generated = False

    def setUp(self):
        """ basic building of the class to test """
        
        self.tmp_process_dir = pjoin(_file_path, 'TMP_TestME7_colorful_epem_guux_NNLO_output')
        # Generate the process output if it does not exist yet or if we
        # are not in debug mode.
        if os.path.isdir(self.tmp_process_dir):
            if not self.is_process_generated and not self.debugging:
                shutil.rmtree(self.tmp_process_dir)
            else:
                TestME7_NNLO_colorful_epem_guux.is_process_generated = True
        if not self.is_process_generated:
            self.cmd = Cmd.MasterCmd()

            # Now generate and output a process, so as to run ME7 commands on it
            self.do('import model loop_sm')
            self.do('set subtraction_currents_scheme colorful')
            self.do('set subtraction_mappings_scheme LorentzNLO')
            self.do('generate e+ e- > g u u~ --NNLO=QCD --ignore_contributions=V,VV')
            self.do('output %s --ignore_integrated_counterterms=all' % self.tmp_process_dir)
            if self.debugging:
                misc.sprint('/!\ USE ONLY FOR DEBUGGING /!\ Output for %s written at %s'
                            % (self.__class__.__name__, self.tmp_process_dir))
            TestME7_NNLO_colorful_epem_guux.is_process_generated = True
        else:
            if self.debugging:
                misc.sprint('/!\ USE ONLY FOR DEBUGGING /!\ Reusing output for %s written at %s' %
                            (self.__class__.__name__, self.tmp_process_dir))

        # Now initialize an ME7 interface on the above process output
        self.cmd = ME7_interface.MadEvent7Cmd(me_dir=self.tmp_process_dir)
        self.cmd.no_notification()

    def __del__(self):
        if os.path.isdir(self.tmp_process_dir) and not self.debugging:
            shutil.rmtree(self.tmp_process_dir)

    def do(self, line):
        """ exec a line in the cmd under test """
        self.cmd.exec_cmd(line)

    def verify_ME7_test_results(self, results_file_path):
        """ Parses and verify that all tests output in 'results_file_path' are passed."""

        for line in open(results_file_path, 'r').read().split('\n'):
            process, limit, outcome, ratio = line.split('|')[:4]
            self.assertTrue(outcome.strip() == 'PASSED', line)

    def test_ME7_g_gqqx_triple_collinear(self):
        """Check the test of collinear limits on a particular process."""

        options = {'correction_order': 'NNLO',
                   'limits': 'C(3,4,5)',
                   'counterterms': 'C(3,4,5)',
                   'process': 'e+ e- > g u u~ u~ u',
                   'show_plots': False,
                   'save_plots': False,
                   'seed': 666,
                   'n_steps': 10,
                   'min_scaling_variable': 1.0e-16,
                   'acceptance_threshold': 5.0e-4,
                   'save_results_to_path': 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))

    def test_ME7_q_qqqx_triple_collinear(self):
        """Check the test of collinear limits on a particular process."""

        options = {'correction_order': 'NNLO',
                   'limits': 'C(4,5,7)',
                   'counterterms': 'C(4,5,7)',
                   'process': 'e+ e- > g u u~ u~ u',
                   'show_plots': False,
                   'save_plots': False,
                   'seed': 666,
                   'n_steps': 10,
                   'min_scaling_variable': 1.0e-16,
                   'acceptance_threshold': 5.0e-4,
                   'save_results_to_path': 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))


    def test_ME7_q_qq2q2x_triple_collinear(self):
        """Check the test of collinear limits on a particular process."""

        options = {'correction_order': 'NNLO',
                   'limits': 'C(4,6,7)',
                   'counterterms': 'C(4,6,7)',
                   'process': 'e+ e- > g u u~ s~ s',
                   'show_plots': False,
                   'save_plots': False,
                   'seed': 666,
                   'n_steps': 10,
                   'min_scaling_variable': 1.0e-16,
                   'acceptance_threshold': 5.0e-4,
                   'save_results_to_path': 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))

    def test_ME7_q_qgg_triple_collinear(self):
        """Check the test of collinear limits on a particular process."""

        options = {'correction_order': 'NNLO',
                   'limits': 'C(3,4,6)',
                   'counterterms': 'C(3,4,6)',
                   'process': 'e+ e- > g u u~ g g',
                   'show_plots': False,
                   'save_plots': False,
                   'seed': 666,
                   'n_steps': 10,
                   'min_scaling_variable': 1.0e-16,
                   'acceptance_threshold': 5.0e-4,
                   'save_results_to_path': 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))

    def test_ME7_g_ggg_triple_collinear(self):
        """Check the test of collinear limits on a particular process."""

        options = {'correction_order': 'NNLO',
                   'limits': 'C(3,6,7)',
                   'counterterms': 'C(3,6,7)',
                   'process': 'e+ e- > g u u~ g g',
                   'show_plots': False,
                   'save_plots': False,
                   'seed': 666,
                   'n_steps': 10,
                   'min_scaling_variable': 1.0e-16,
                   'acceptance_threshold': 5.0e-4,
                   'save_results_to_path': 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))

    def test_ME7_all_triple_collinears(self):
        """Check the test of collinear limits on a particular process."""

        options = {'correction_order': 'NNLO',
                   'limits': "r'\(C\(\d,\d,\d\),\)'",
                   'counterterms': 'def',
                   'show_plots': False,
                   'save_plots': False,
                   'seed': 666,
                   'n_steps': 10,
                   'min_scaling_variable': 1.0e-16,
                   'acceptance_threshold': 8.0e-4,
                   'save_results_to_path': 'test_IR_limit_output_for_acceptance_test.dat'
                   }

        self.do(get_test_IR_limit_cmd(options))
        self.verify_ME7_test_results(pjoin(self.tmp_process_dir, 'test_IR_limit_output_for_acceptance_test.dat'))
