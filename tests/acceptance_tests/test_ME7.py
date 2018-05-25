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

#===============================================================================
# TestME7_IR_Limits
#===============================================================================
class TestME7(unittest.TestCase):
    """This test validates the command 'test_IR_limits' of ME7"""
    
    # If the debug mode is set to True, then the process output is not refreshed
    # but reused instead
    debugging = True 

    def setUp(self):
        """ basic building of the class to test """
        
        self.tmp_process_dir = pjoin(_file_path, 'TMP_TestME7_IR_Limits_output')
        # Generate the process output if it does not exist yet or if we
        # are not in debug mode.
        if not os.path.isdir(self.tmp_process_dir) or not self.debugging:
            self.cmd = Cmd.MasterCmd()
            if os.path.isdir(self.tmp_process_dir):
                shutil.rmtree(self.tmp_process_dir)

            # Now generate and output a process, so as to run ME7 commands on it
            self.do('import model loop_sm')
            self.do('generate e+ e- > j j j --NLO=QCD')
            self.do('output %s'%self.tmp_process_dir)

        # Now initialize an ME7 interface on the above process output
        self.cmd = ME7_interface.MadEvent7Cmd(me_dir=self.tmp_process_dir)
        self.cmd.no_notification()
 
    def tearDown(self):
        if os.path.isdir(self.tmp_process_dir) and not self.debugging:
            shutil.rmtree(self.tmp_process_dir)

    def do(self, line):
        """ exec a line in the cmd under test """
        self.cmd.exec_cmd(line)
        
    def test_ME7_qqxQQx_collinear_limits(self):
        """Check the test of collinear limits on a particular process."""
        
        main_cmd = 'test_IR_limits'

        options = {'correction_order'       : 'NLO',
#                   'limit_type'             : 'C(S(3),4)',
                   'limit_type'             : 'S(3)',
                   'counterterms'           : None,
                   'process'                : 'e+ e- > g g d d~',
                   'seed'                   : '666',
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-5,
                   'acceptance_threshold'   : 1.0e-4,
                   }

        self.do('%s %s'%(main_cmd, ' '.join(
            ('--%s=%s'%(key,value) if value is not None else '--%s'%key)
            for key,value in options.items()
        )))

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

        n_calls = 5000
        with ME7_interface.ME7RunEnvironment( silence = True, loggers = logging.CRITICAL ):
            misc.sprint('\n'+'\n'.join('%d : %g ms'%(i+1, 1.e3*(res/float(n_calls)))
                    for i,res in enumerate(timeit.repeat(call, number=n_calls, repeat=1))))

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

        n_calls = 5000
        with ME7_interface.ME7RunEnvironment( silence = True, loggers = logging.CRITICAL ):
            res = '\n'+'\n'.join('%d : %g ms'%(i+1, 1.e3*(res/float(n_calls)))
                    for i,res in enumerate(timeit.repeat(call, number=n_calls, repeat=1)))
        print res
            
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
            res = '\n'+'\n'.join('%d : %g ms'%(i+1, 1.e3*(res/float(n_calls)))
                    for i,res in enumerate(timeit.repeat(call, number=n_calls, repeat=1)))
        print res
            
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
            res = '\n'+'\n'.join('%d : %g ms'%(i+1, 1.e3*(res/float(n_calls)))
                    for i,res in enumerate(timeit.repeat(call, number=n_calls, repeat=1)))
        print res
