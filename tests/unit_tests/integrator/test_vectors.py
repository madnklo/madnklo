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

    n_random_tests = 10
    seed = 42

    def test_Vector_basic(self):
        """Test the basic operations for vectors."""

        # Component-wise addition and subtraction
        v1 = Vector([1., 0., 2., 3.])
        v2 = Vector(4*[0.5, ])
        v12 = Vector([1.5, 0.5, 2.5, 3.5])
        v1pv2 = v1 + v2
        self.assertEqual(v1pv2, v12)
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

        random.seed(self.seed)
        for _ in range(self.n_random_tests):

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

        random.seed(self.seed)
        for _ in range(self.n_random_tests):

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
            # WARNING: this is delicate because of numerics
            self.assertEqual(v4, v4_old)

            # Test the boost_from_to in its two implementations
            #    Generate two random vectors with equal squares (potentially massless)
            v5 = LorentzVector([0, ] + [random.random() for _ in range(3)])
            v6 = LorentzVector([0, ] + [random.random() for _ in range(3)])
            massive = random.choice([True, False])
            if massive: m = random.random()
            else: m = 0
            v5.set_square(m**2)
            v6.set_square(m**2)
            #    Test that both versions send v5 in v6
            v5_boosted_1 = LorentzVector(v5)
            v5_boosted_2 = LorentzVector(v5)
            boostvector_5_6 = LorentzVector.boost_vector_from_to(v5, v6)
            self.assertEqual(v5_boosted_1.boost_from_to(v5, v6), v6)
            self.assertEqual(v5_boosted_2.boost(boostvector_5_6), v6)
            #    Test that both versions recover v5 after inverse
            boostvector_6_5 = LorentzVector.boost_vector_from_to(v6, v5)
            self.assertEqual(v5_boosted_1.boost_from_to(v6, v5), v5)
            self.assertEqual(v5_boosted_2.boost(boostvector_6_5), v5)
            #    Test that the action of the two on a different vector is the same
            v7 = LorentzVector([0, ] + [random.random() for _ in range(3)])
            v7.set_square(0 if random.choice([True, False]) else random.random())
            v7_boosted_1 = LorentzVector(v7)
            v7_boosted_1.boost_from_to(v5, v6)
            v7_boosted_2 = LorentzVector(v7)
            v7_boosted_2.boost(boostvector_5_6)
            self.assertEqual(v7_boosted_1, v7_boosted_2)
