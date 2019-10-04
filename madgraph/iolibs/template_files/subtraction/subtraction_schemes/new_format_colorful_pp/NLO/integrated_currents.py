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
"""Implementation of NLO type of integrated currents for colorful_pp."""

import os
import sys
import math

import madgraph.core.subtraction as subtraction
import madgraph.integrator.phase_space_generators as PS_utils
from madgraph.core.base_objects import EpsilonExpansion
import madgraph.various.misc as misc

import commons.utils as utils
import commons.QCD_local_currents as currents
import commons.factors_and_cuts as factors_and_cuts

from integrated_current_expressions import HE

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError


