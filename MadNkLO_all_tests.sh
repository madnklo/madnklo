#!/bin/bash
UNIT_TEST='./tests/test_manager.py -t0 -p U'
ACCEPTANCE_TEST='./tests/test_manager.py -t0 -p A'
$UNIT_TEST $(cat MadNkLO_unit_test.txt)
$ACCEPTANCE_TEST $(cat MadNkLO_acceptance_test.txt)
$ACCEPTANCE_TEST TestME7_NLO_colorful_epem_jjj
$ACCEPTANCE_TEST TestME7_NLO_cataniseymour_epem_jjj
$ACCEPTANCE_TEST TestME7_NLO_colorful_pp_jj
$ACCEPTANCE_TEST TestME7_NNLO_colorful_epem_guux
$ACCEPTANCE_TEST SubtractionCurrentTest
