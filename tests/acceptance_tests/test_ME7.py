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

pjoin = os.path.join

logger = logging.getLogger('test_ME7')

import tests.unit_tests.iolibs.test_file_writers as test_file_writers

import madgraph.interface.master_interface as Cmd
import madgraph.interface.ME7_interface as ME7interface
import madgraph.various.misc as misc
_file_path = os.path.dirname(os.path.realpath(__file__))
_pickle_path =os.path.join(_file_path, 'input_files')

from madgraph import MG4DIR, MG5DIR, MadGraph5Error, InvalidCmd

#===============================================================================
# TestME7_IR_Limits
#===============================================================================
class TestME7_IR_Limits(unittest.TestCase):
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
            self.do('generate e+ e- > j j j --NLO=QCD --ignore_contributions=V')
            self.do('output %s'%self.tmp_process_dir)

        # Now initialize an ME7 interface on the above process output
        self.cmd = ME7interface.MadEvent7Cmd(me_dir=self.tmp_process_dir)
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
                   'limit_type'             : 'C(3,5)',
                   'process'                : 'e+ e- > u u~ s s~ ',
                   'seed'                   : '666',
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-6,
                   'acceptance_threshold'   : 1.0e-6,
                   'compute_only_limit_defining_counterterm' : True,
                   }
        
        self.do('%s %s'%(main_cmd, ' '.join( ('--%s=%s'%(key,value) if value is not None else '--%s'%key) 
                                                                   for key,value in options.items())))
        
        options = {'correction_order'       : 'NLO',
                   'limit_type'             : 'C(3,5)',
                   'process'                : 'e+ e- > d d~ g g ',
                   'seed'                   : '666',
                   'n_steps'                : 10,
                   'min_scaling_variable'   : 1.0e-6,
                   'acceptance_threshold'   : 1.0e-6,
                   'compute_only_limit_defining_counterterm' : True,
                   }

        self.do('%s %s'%(main_cmd, ' '.join( ('--%s=%s'%(key,value) if value is not None else '--%s'%key)
                                                                  for key,value in options.items())))
