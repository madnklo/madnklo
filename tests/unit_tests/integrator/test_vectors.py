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
"""Unit test library for vectors."""

import tests.unit_tests as unittest
import random
from madgraph.integrator.vectors import Vector, LorentzVector

#=========================================================================================
# Test Vectors
#=========================================================================================

class VectorsTest(unittest.TestCase):
    """Test class for Vector and LorentzVector."""

    def test_Vector_basic(self):
        """Test the basic operations for vectors."""

        # Component-wise addition and subtraction
        v1 = Vector([1., 0., 2., 3.])
        v2 = Vector(4*[0.5, ])
        self.assertEqual(v1 + v2, Vector([1.5, 0.5, 2.5, 3.5]))
        v3 = Vector([1., 2., 0.5, -1.])
        self.assertEqual(v1 - v3, Vector([0., -2., 1.5, 4.]))
        v3 += Vector([1., 0., 0., 0.])
        self.assertEqual(v3, Vector([2., 2., 0.5, -1.]))
        v3 -= Vector([0.5, 1.5, 0., 0.])
        self.assertEqual(v3, Vector([1.5, 0.5, 0.5, -1.]))

        # Multiplication and division by scalars
        self.assertEqual(v2 * 2., Vector(4*[1., ]))
        self.assertEqual(v1 / 4., Vector([0.25, 0., 0.5, 0.75]))
        v3 *= 3
        self.assertEqual(v3, Vector([4.5, 1.5, 1.5, -3.]))
        v3 /= 1.5
        self.assertEqual(v3, Vector([3., 1., 1., -2.]))
        self.assertEqual(3 * v1, v1 * 3)

        # Negative
        self.assertEqual(-v3, Vector([-3., -1., -1., 2.]))

    def test_Vector_Euclid(self):
        """Test scalar products and related functions for the class Vector."""

        # Square and norm
        v1 = Vector([0., 1., -2.])
        self.assertEqual(v1.square(), 5.)
        v2 = Vector([3., 4., 0.])
        self.assertEqual(abs(v2), 5.)
        v2n = 0.2 * v2
        v2.normalize()
        self.assertEqual(v2, v2n)
        self.assertEqual(Vector.dot(v1, v2), 0.8)
        v3 = Vector([random.random() for _ in range(3)])
        w = Vector([random.random() for _ in range(3)])
        v3p = v3.project_onto(w)
        v3t = v3.component_orthogonal_to(w)
        self.assertEqual(v3, v3p + v3t)
        self.assertAlmostEqual(Vector.dot(v3t, w), 0.)
        self.assertAlmostEqual(abs(v3p+w), abs(v3p)+abs(w))

    def test_Vector_Minkowski(self):
        """Test the relativistic norm."""

        # Square and norm
        v1 = LorentzVector([1, 0, 0, +1])
        v2 = LorentzVector([1, 0, 0, -1])
        self.assertEqual(v1.square(), 0)
        self.assertEqual(v2.square(), 0)
        self.assertEqual(v1.dot(v2), 2)

        # Test the rotoboost
        #   Definition
        v1 = LorentzVector([0, ] + [random.random() for _ in range(3)])
        v2 = LorentzVector([0, ] + [random.random() for _ in range(3)])
        m = random.random()
        v1.set_square(m**2)
        v2.set_square(m**2)
        self.assertEqual(v2.rotoboost(v2, v1), v1)
        #   Inverse
        v3 = LorentzVector([0, ] + [random.random() for _ in range(3)])
        v3.set_square(m**2)
        v4 = LorentzVector([random.random() for _ in range(4)])
        v4_old = v4.get_copy()
        v4.rotoboost(v1, v3)
        v4.rotoboost(v3, v1)
        self.assertEqual(v4, v4_old)
