"""A test suite to check that classic NLO corrections are properly calculated in MadNkLO"""


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
_input_files_path = pjoin(_file_path, 'input_files')

class MadNkLONLOCrossSectionChecker(unittest.TestCase):
    """A class to compare MadNkLO NLO cross-section in well-known processes"""

    def compare_cross_section_to_values(self, values, my_proc_list=[],
                                        orders={}, model='sm',
                                        filename="", print_result=False, append_output=False,
                                        tolerance=1e-02, ME7_options={},custom_run_card=None):

        ME7_runner = madevent_comparator.MadNkLOFinalNLORunner(**ME7_options)
        ME7_runner.setup(MG5DIR)
        ME7_runner.store_proc_card = True
        ME7_runner.custom_run_card = custom_run_card

        # Create and setup a comparator
        my_comp = madevent_comparator.MadEventComparator(allow_no_present=True)
        my_comp.set_me_runners(ME7_runner)

        # Run the actual comparison
        my_comp.run_comparison(my_proc_list,
                               [model], orders)

        # add the default value to the comparison
        my_comp.results.insert(0, values)
        my_comp.me_runners = (madevent_comparator.FakeRunner(), my_comp.me_runners[0])

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

    def test_MNK_NLO_epem_jj(self):
        """Quick test of our favourite standard candle."""

        logging.root.setLevel(logging.DEBUG)

        custom_run_card = os.path.join(_input_files_path,"test_MNK_NLO_epem_jj.Run_card.dat")

        my_proc_list = ['e+ e- > j j']
        values = {
            'LO':'0.5320',
            'NLO':'0.0198'
        }

        # Run the comparison
        self.compare_cross_section_to_values(values, my_proc_list,
                                             orders={'QED': 99, 'QCD': 99},
                                             filename="short_MNK_NLO_epem_jj.log",
                                             ME7_options={'n_points': 2,
                                                          'integrator': 'VEGAS3',
                                                          'PS_generator':'FLATPS'},
                                             custom_run_card=custom_run_card)