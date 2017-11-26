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
import unittest
import logging

logger = logging.getLogger('test_PS_volume')

import madgraph.integrator.functions as functions
import madgraph.integrator.integrands as integrands
import madgraph.integrator.phase_space_generators as psgen
import madgraph.integrator.vegas3_integrator as vegas
import madgraph.various.misc as misc

#=========================================================================================
# Phase space generator
#=========================================================================================

E = 1.
n_final = 2
massive = False

verbosity = 2

# This test should really use a decay,
# but FlatInvertiblePhasespace does not support it yet
IS_masses = [0.] * 2
FS_masses = [0.] * n_final
generator = psgen.FlatInvertiblePhasespace(
    IS_masses, FS_masses, (E/2., E/2.), beam_types=(0, 0))

#=========================================================================================
# Flat PS integrand
#=========================================================================================

class FlatPS(functions.VirtualFunction):

    def __call__(self, continuous_inputs, discrete_inputs, **opts):

        ps, wgt = generator.generateKinematics(E, continuous_inputs)
        return dict({'weight': wgt})

#=========================================================================================
# TestPSVolume
#=========================================================================================

class TestPSVolume(unittest.TestCase):
    """This test checks the phase-space volume of mappings
    against the one of the flat phase-space generator.
    """

    def setUp(self):

        self.function_list = functions.FunctionList()
        self.function_list.append(FlatPS())

    def test_PS_volume(self):
        """Compute the phase-space volume."""

        # Build the integrand
        my_integrand = integrands.VirtualIntegrand()
        my_integrand.set_dimensions(generator.get_dimensions())
        my_integrand.function_list = self.function_list
        misc.sprint(str(generator.get_dimensions()))
        # Setup the integrator
        my_integrator = vegas.Vegas3Integrator(
            my_integrand, verbosity=verbosity, refine_n_iterations=10)

        res = my_integrator.integrate()
        # Finally integrate
        print '\n' + '=' * 70
        print '%.4e +/- %.2e' % (res[0], res[1])
        print '=' * 70

        pass
