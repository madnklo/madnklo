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
"""Unit test library for the various phase-space generation/handling features."""

import madgraph.integrator.phase_space_generators as PS
import math
import random

import tests.unit_tests as unittest

#===============================================================================
# Test Vectors
#===============================================================================

class VectorsTest(unittest.TestCase):
    """Test class for BaseVector, Vector and LorentzVector."""

    def setUp(self):
        pass

    def test_Vector_basic(self):
        """Test the basic operations for vectors."""

        # Component-wise addition and subtraction
        v1 = PS.Vector([1., 0., 2., 3.])
        v2 = PS.Vector(4, 0.5)
        self.assertEqual(v1 + v2, PS.Vector([1.5, 0.5, 2.5, 3.5]))
        v3 = PS.Vector([1., 2., 0.5])
        self.assertEqual(v1 - v3, PS.Vector([0., -2., 1.5, 3.]))
        v3 += PS.Vector([1., 0.])
        self.assertEqual(v3, PS.Vector([2., 2., 0.5]))
        v3 -= PS.Vector([0.5, 1.5])
        self.assertEqual(v3, PS.Vector([1.5, 0.5, 0.5]))

        # Multiplication and division by scalars
        self.assertEqual(v2 * 2., PS.Vector(4, 1.))
        self.assertEqual(v1 / 4., PS.Vector([0.25, 0., 0.5, 0.75]))
        v3 *= 3
        self.assertEqual(v3, PS.Vector([4.5, 1.5, 1.5]))
        v3 /= 1.5
        self.assertEqual(v3, PS.Vector([3., 1., 1.]))
        self.assertEqual(3 * v1, v1 * 3)

        # Negative
        self.assertEqual(-v3, PS.Vector([-3., -1., -1.]))

    def test_Vector_Euclid(self):
        """Test scalar products and related functions for the class Vector."""

        # Square and norm
        v1 = PS.Vector([0., 1., -2.])
        self.assertEqual(v1.square(), 5.)
        v2 = PS.Vector([3., 4., 0.])
        self.assertEqual(abs(v2), 5.)
        v2n = 0.2 * v2
        v2.normalize()
        self.assertEqual(v2, v2n)
        self.assertEqual(PS.Vector.dot(v1, v2), 0.8)
        v3 = PS.Vector([random.random() for _ in range(3)])
        w = PS.Vector([random.random() for _ in range(3)])
        v3p = v3.project_onto(w)
        v3t = v3.component_orthogonal_to(w)
        self.assertAlmostEqual(v3, v3p + v3t)
        self.assertAlmostEqual(PS.Vector.dot(v3t, w), 0.)
        self.assertAlmostEqual(abs(v3p+w), abs(v3p)+abs(w))

    def test_Vector_Minkowski(self):
        """Test the relativistic norm."""

        # Square and norm
        v1 = PS.LorentzVector([1, 0, 0, +1])
        v2 = PS.LorentzVector([1, 0, 0, -1])
        self.assertEqual(v1.square(), 0)
        self.assertEqual(v2.square(), 0)
        self.assertEqual(v1.dot(v2), 2)

#===============================================================================
# Test the phase-space generators
#===============================================================================
class PhaseSpaceGeneratorsTest(unittest.TestCase):
    """ Test various phase-space generators."""

    def setUp(self):
        pass

    def test_flat_invertible_phase_space(self):
        """ Tests the flat invertible phase-space."""
        
        E_cm  = 5000.0
    
        # Try to run the above for a 2->8.
        my_PS_generator = PS.FlatInvertiblePhasespace(
            [0.]*2, [100. + 10.*i for i in range(8)],beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0))
        # Try to run the above for a 2->1.    
        #    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [5000.0])
        
        random_variables = [random.random() for _ in range(my_PS_generator.nDimPhaseSpace())]
    
        momenta, wgt = my_PS_generator.generateKinematics(E_cm, random_variables)
       
        #print "\n ========================="
        #print " ||    PS generation    ||"
        #print " ========================="    
        #print "\nRandom variables :\n",random_variables
        #print "\n%s\n"%my_PS_generator.nice_momenta_string(
        #                    momenta, recompute_mass=True, n_initial=my_PS_generator.n_initial)
        #print "Phase-space weight : %.16e\n"%wgt,
    
        variables_reconstructed, wgt_reconstructed = \
                                             my_PS_generator.invertKinematics(E_cm, momenta)
    
        #print "\n ========================="
        #print " || Kinematic inversion ||"
        #print " ========================="
        #print "\nReconstructed random variables :\n",variables_reconstructed
        differences = [abs(variables_reconstructed[i]-random_variables[i]) 
                                        for i in range(len(variables_reconstructed))]

        self.assertLess(max(differences[i]/random_variables[i] for i in range(len(differences))), 1.0e-10)
        self.assertLess(abs(wgt-wgt_reconstructed)/abs(wgt), 1.0e-10)
        
        #print "Reconstructed weight = %.16e"%wgt_reconstructed
        #if differences:
        #    print "\nMax. relative diff. in reconstructed variables = %.3e"%\
        #        max(differences[i]/random_variables[i] for i in range(len(differences)))
        #print "Rel. diff. in PS weight = %.3e\n"%((wgt_reconstructed-wgt)/wgt)

#===============================================================================
# Test the phase-space mappers
#===============================================================================

class CollinearVariablesTest(unittest.TestCase):
    """Test class for variables describing internal collinear structure."""

    my_mapping = PS.ElementaryMappingCollinearFinal()
    n_children = 3

    def setUp(self):
        pass

    def test_variables(self):
        """Test determination of collinear variables and reverse mapping."""

        # Generate n_children+1 random massless vectors
        # (the light-cone direction is also random)
        my_PS_point = dict()
        for i in range(self.n_children+1):
            my_PS_point[i] = PS.LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            )
            my_PS_point[i][0] = math.sqrt(-my_PS_point[i].square())
        # Compute collinear variables
        variables = dict()
        self.my_mapping.get_collinear_variables(
            my_PS_point, self.n_children, range(self.n_children),
            variables
        )
        # Compute total momentum
        total_momentum = PS.LorentzVector(4, 0.)
        for i in range(self.n_children):
            total_momentum += my_PS_point[i]
        # Compute new phase space point
        new_PS_point = dict()
        new_PS_point[self.n_children] = my_PS_point[self.n_children]
        self.my_mapping.set_collinear_variables(
            new_PS_point, self.n_children, range(self.n_children),
            total_momentum, variables
        )
        for i in range(self.n_children):
            for j in range(4):
                self.assertAlmostEqual(my_PS_point[i][j], new_PS_point[i][j])
