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
"""Unit test library for mappings."""

import madgraph.integrator.walkers as walkers
import madgraph.core.subtraction as subtraction
import madgraph.core.base_objects as base_objects
import madgraph.various.misc as misc
from madgraph.integrator.vectors import Vector, LorentzVector
from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList

import math
import random
import os

import tests.unit_tests as unittest
import tests.input_files.simple_qcd as simple_qcd

pjoin = os.path.join

#=========================================================================================
# Shorthands for initial and final state
#=========================================================================================

INITIAL = base_objects.Leg.INITIAL
FINAL   = base_objects.Leg.FINAL

assert subtraction.SubtractionLeg.INITIAL == INITIAL
assert subtraction.SubtractionLeg.FINAL   == FINAL

#=========================================================================================
# AssertDictAlmostEqual
#=========================================================================================

def assertDictAlmostEqual(test, dict1, dict2):

    test.assertEqual(sorted(dict1.keys()), sorted(dict2.keys()))
    for key in dict1.keys():
        if isinstance(dict1[key], Vector):
            test.assertEqual(dict1[key].__class__, dict2[key].__class__)
            for i in range(max(len(dict1[key]), len(dict2[key]))):
                test.assertAlmostEqual(dict1[key][i], dict2[key][i])
        else:
            test.assertAlmostEqual(dict1[key], dict2[key])

#=========================================================================================
# Generation of random momenta
#=========================================================================================

def random_momentum(square=None):
    """Generate a random momentum, either massive or massless."""

    foo = LorentzVector([0., ] + [(2 * random.random() - 1) for _ in range(3)])
    if square is None:
        foo.set_square(random.random())
    else:
        foo.set_square(square)
    return foo

#=========================================================================================
# Test the phase-space walkers
#=========================================================================================

class WalkersTest(unittest.TestCase):
    """Test class for walkers."""

    # Parameters
    #=====================================================================================

    # Verbosity (silent = 0, max = 4)
    verbosity = 0
    # Random seed to make tests deterministic
    seed = 42
    # Number of PS points the invertibility test is run for (more = stronger, slower test)
    n_test_invertible = 3
    # Number of PS points the approach_limit test is run for (more = stronger, slower test)
    n_test_approach = 5
    # Number of PS points for the low_level_approach_limit test (more = stronger, slower test)
    n_test_low_level_approach = 100
    # Values of the parameter in approach_limit (more, smaller = stronger, slower test)
    parameter_values = [0.1 ** i for i in range(5)]

    # Setup
    #=====================================================================================

    # IRSubtraction module
    irs = subtraction.IRSubtraction(
        simple_qcd.model, coupling_types=('QCD', ), n_unresolved=None )
    # Not aligning the initial state with z avoids numerical issues
    initial_along_z = False
    # Separator for output
    stars = "*" * 90

    # Functions
    #=====================================================================================

    @classmethod
    def generate_PS_point(cls, process):
        """Generate a phase-space point to test the walker."""

        model = process.get('model')
        # Generate random vectors
        my_PS_point = LorentzVectorDict()
        legs_FS = tuple(
            subtraction.SubtractionLeg(leg)
            for leg in process['legs']
            if leg['state'] == FINAL )
        for leg in legs_FS:
            my_PS_point[leg.n] = random_momentum(
                simple_qcd.masses[model.get_particle(leg.pdg)['mass'] ] )
        total_momentum = LorentzVector()
        for leg in legs_FS: total_momentum += my_PS_point[leg.n]
        legs_IS = tuple(
            subtraction.SubtractionLeg(leg)
            for leg in process['legs']
            if leg['state'] == INITIAL )
        if len(legs_IS) == 1:
            my_PS_point[legs_IS[0].n] = total_momentum
        elif len(legs_IS) == 2:
            if cls.initial_along_z:
                bv = total_momentum.boostVector()
                for key in my_PS_point.keys():
                    my_PS_point[key].boost(-bv)
                total_momentum.boost(-bv)
                E = total_momentum[0]
                my_PS_point[1] = LorentzVector([E/2., 0., 0., +E/2.])
                my_PS_point[2] = LorentzVector([E/2., 0., 0., -E/2.])
            else:
                E = abs(total_momentum.square()) ** 0.5
                rest_momentum = LorentzVector([E, 0., 0., 0.])
                my_PS_point[1] = LorentzVector([E/2., 0., 0., +E/2.])
                my_PS_point[2] = LorentzVector([E/2., 0., 0., -E/2.])
                my_PS_point[1].rotoboost(rest_momentum, total_momentum)
                my_PS_point[2].rotoboost(rest_momentum, total_momentum)
        else: raise BaseException
        return my_PS_point

    def _test_invertible(
        self, walker, process,
        max_unresolved_in_elementary, max_unresolved_in_combination):
        """Check that the walker and its inverse yield the same result.

        :param walker: Mapping walker to be tested
        :type walker: walkers.VirtualWalker

        :param process: The physical process the walker will be tested for
        :type process: base_objects.Process

        :param max_unresolved_in_elementary: Maximum number of unresolved particles
        within the same elementary operator
        :type max_unresolved_in_elementary: positive integer

        :param max_unresolved_in_combination: Maximum number of unresolved particles
        within a combination of elementary operators
        :type max_unresolved_in_combination: positive integer
        """

        random.seed(self.seed)
        if self.verbosity > 2:
            print "\n" + self.stars * (self.verbosity - 2)
        if self.verbosity > 0:
            tmp_str = "test_invertible for " + walker.__class__.__name__
            tmp_str += " with " + process.nice_string()
            print tmp_str
        if self.verbosity > 2:
            print self.stars * (self.verbosity - 2) + "\n"
        my_operators = self.irs.get_all_elementary_operators(
            process, max_unresolved_in_elementary)
        my_combinations = self.irs.get_all_combinations(
            my_operators, max_unresolved_in_combination)
        my_counterterms = [
            self.irs.get_counterterm(combination, process)
            for combination in my_combinations ]

        # For each counterterm
        for i in range(len(my_counterterms)):
            if self.verbosity > 3: print "\n" + self.stars * (self.verbosity - 3)
            if self.verbosity > 1: print "Considering counterterm", my_counterterms[i]
            if self.verbosity > 3: print self.stars * (self.verbosity - 3) + "\n"
            for j in range(self.n_test_invertible):
                if self.verbosity > 2:
                    print "Phase space point #", j+1
                my_PS_point = self.generate_PS_point(process)
                # Compute collinear variables
                res_dict1 = walker.walk_to_lower_multiplicity(
                    my_PS_point, my_counterterms[i],
                    compute_kinematic_variables=True, compute_jacobian=True )
                currs1, ME1, kin = (
                    res_dict1['currents'], res_dict1['matrix_element'],
                    res_dict1['kinematic_variables'] )
                if self.verbosity > 3:
                    print "Walking down"
                    for curr in currs1:
                        print curr['higher_PS_point']
                    print ME1[1]
                res_dict2 = walker.walk_to_higher_multiplicity(
                    ME1[1], my_counterterms[i], kin,
                    compute_jacobian=True )
                currs2, ME2 = res_dict2['currents'], res_dict2['matrix_element']
                if self.verbosity > 3:
                    print "Walking up"
                    print ME2[1]
                    for curr in currs2:
                        print curr['higher_PS_point']
                # Check currents
                self.assertEqual(len(currs1), len(currs2))
                for i_curr in range(len(currs1)):
                    c1, c2 = currs1[i_curr], currs2[i_curr]
                    self.assertEqual(c1['stroll_currents'], c2['stroll_currents'])
                    self.assertDictEqual(c1['higher_PS_point'], c2['higher_PS_point'])
                    self.assertDictEqual(c1['lower_PS_point'], c2['lower_PS_point'])
                    assertDictAlmostEqual(self, c1['stroll_vars'], c2['stroll_vars'])
                # Check MEs
                self.assertEqual(ME1[0], ME2[0])
                self.assertDictEqual(ME1[1], ME2[1])

    def _test_approach_limit(
        self, walker, process,
        max_unresolved_in_elementary, max_unresolved_in_combination):
        """Check that the walker is capable of approaching limits.

        :param walker: Mapping walker to be tested
        :type walker: walkers.VirtualWalker

        :param process: The physical process the walker will be tested for
        :type process: base_objects.Process

        :param max_unresolved_in_elementary: Maximum number of unresolved particles
        within the same elementary operator
        :type max_unresolved_in_elementary: positive integer

        :param max_unresolved_in_combination: Maximum number of unresolved particles
        within a combination of elementary operators
        :type max_unresolved_in_combination: positive integer
        """

        if self.verbosity > 2:
            print "\n" + self.stars * (self.verbosity - 2)
        if self.verbosity > 0:
            tmp_str = "test_approach_limit for " + walker.__class__.__name__
            tmp_str += " with " + process.nice_string()
            print tmp_str
        if self.verbosity > 2:
            print self.stars * (self.verbosity - 2) + "\n"
        random.seed(self.seed)
        # Generate all counterterms for this process, and separate the non-singular one
        my_operators = self.irs.get_all_elementary_operators(
            process, max_unresolved_in_elementary)
        my_combinations = self.irs.get_all_combinations(
            my_operators, max_unresolved_in_combination)
        my_counterterms = [
            self.irs.get_counterterm(combination, process)
            for combination in my_combinations ]
        # Get all legs in the FS and the model to check masses after approach_limit
        legs_FS = tuple(
            subtraction.SubtractionLeg(leg)
            for leg in process['legs']
            if leg['state'] == FINAL )
        model = process.get('model')
        # For each counterterm
        for ct in my_counterterms:
            if not ct.is_singular():
                continue
            if self.verbosity > 3: print "\n" + self.stars * (self.verbosity - 3)
            if self.verbosity > 1: print "Considering counterterm", ct
            if self.verbosity > 3: print self.stars * (self.verbosity - 3) + "\n"
            ss = ct.reconstruct_complete_singular_structure()
            for j in range(self.n_test_invertible):
                if self.verbosity > 2:
                    print "Phase space point #", j+1
                # Generate random vectors
                my_PS_point = self.generate_PS_point(process)
                if self.verbosity > 3:
                    print "Starting phase space point:\n", my_PS_point, "\n"
                squares = {key: my_PS_point[key].square() for key in my_PS_point.keys()}
                # Compute collinear variables
                for alpha in self.parameter_values:
                    new_PS_point = walker.approach_limit(my_PS_point, ss, alpha, process)
                    if self.verbosity > 4:
                        print "New PS point for", alpha, ":\n", new_PS_point
                    for leg in legs_FS:
                        if model.get_particle(leg.pdg)['mass'].lower() == 'zero':
                            self.assertLess(
                                abs(new_PS_point[leg.n].square()),
                                math.sqrt(new_PS_point[leg.n].eps) )
                        else:
                            self.assertAlmostEqual(
                                new_PS_point[leg.n].square(),
                                squares[leg.n] )

    def _test_low_level_approach_limit(self, process, low_level_limit):

        model = process.get('model')
        legs = process.get('legs')
        for j in range(self.n_test_low_level_approach):
            if self.verbosity > 0:
                print "Phase space point #", j + 1
            my_PS_point = self.generate_PS_point(process)
            clean_momenta_dict = subtraction.create_momenta_dict(process)
            new_PS_point = walkers.low_level_approach_limit(
                my_PS_point, low_level_limit, 10 ** (-8*random.random()), clean_momenta_dict,
                verbose=True )
            # Sanity checks on masses and energy positivity
            for leg in legs:
                pdg = leg['id']
                n = leg['number']
                if model.get_particle(pdg)['mass'].lower() == 'zero':
                    self.assertLess(
                        abs(new_PS_point[n].square()),
                        new_PS_point[n].eps)
                else:
                    self.assertAlmostEqual(
                        new_PS_point[n].square(),
                        my_PS_point[n].square())
                self.assertTrue(
                    new_PS_point[n][0] > 0 or abs(new_PS_point[n][0]) < new_PS_point[n].eps)

    # Processes
    #=====================================================================================

    # H > u u~ d d~
    H_to_uuxddx_legs = base_objects.LegList([
        base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
        base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
        base_objects.Leg({'number': 3, 'id': -1, 'state': FINAL}),
        base_objects.Leg({'number': 4, 'id':  2, 'state': FINAL}),
        base_objects.Leg({'number': 5, 'id': -2, 'state': FINAL}),
    ])
    H_to_uuxddx = base_objects.Process({
        'legs': H_to_uuxddx_legs,
        'model': simple_qcd.model,
        'n_loops': 0
    })

    # H > u u~ d d~ g
    H_to_uuxddxg_legs = base_objects.LegList([
        base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
        base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
        base_objects.Leg({'number': 3, 'id': -1, 'state': FINAL}),
        base_objects.Leg({'number': 4, 'id':  2, 'state': FINAL}),
        base_objects.Leg({'number': 5, 'id': -2, 'state': FINAL}),
        base_objects.Leg({'number': 6, 'id': 21, 'state': FINAL}),
    ])
    H_to_uuxddxg = base_objects.Process({
        'legs': H_to_uuxddxg_legs,
        'model': simple_qcd.model,
        'n_loops': 0
    })

    # H > u u~ d d~ H
    H_to_uuxddxH_legs = base_objects.LegList([
        base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
        base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
        base_objects.Leg({'number': 3, 'id': -1, 'state': FINAL}),
        base_objects.Leg({'number': 4, 'id':  2, 'state': FINAL}),
        base_objects.Leg({'number': 5, 'id': -2, 'state': FINAL}),
        base_objects.Leg({'number': 6, 'id': 25, 'state': FINAL}),
    ])
    H_to_uuxddxH = base_objects.Process({
        'legs': H_to_uuxddxH_legs,
        'model': simple_qcd.model,
        'n_loops': 0
    })

    # H > b b~ u u~ g g g
    H_to_bbxuuxggg_legs = base_objects.LegList([
        base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
        base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
        base_objects.Leg({'number': 3, 'id': -1, 'state': FINAL}),
        base_objects.Leg({'number': 4, 'id':  2, 'state': FINAL}),
        base_objects.Leg({'number': 5, 'id': -2, 'state': FINAL}),
        base_objects.Leg({'number': 6, 'id': 21, 'state': FINAL}),
        base_objects.Leg({'number': 7, 'id': 21, 'state': FINAL}),
        base_objects.Leg({'number': 8, 'id': 21, 'state': FINAL}),
    ])
    H_to_bbxuuxggg = base_objects.Process({
        'legs': H_to_bbxuuxggg_legs,
        'model': simple_qcd.model,
        'n_loops': 0
    })

    # H > q q~ g g H
    H_to_qqxggH_legs = base_objects.LegList([
        base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
        base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
        base_objects.Leg({'number': 3, 'id': -1, 'state': FINAL}),
        base_objects.Leg({'number': 4, 'id': 21, 'state': FINAL}),
        base_objects.Leg({'number': 5, 'id': 21, 'state': FINAL}),
        base_objects.Leg({'number': 6, 'id': 25, 'state': FINAL}),
    ])
    H_to_qqxggH = base_objects.Process({
        'legs': H_to_qqxggH_legs,
        'model': simple_qcd.model,
        'n_loops': 0
    })

    # q q~ > g g H
    # HACK: one gluon more to avoid issue with final-only soft mapping
    qqx_to_ggH_legs = base_objects.LegList([
        base_objects.Leg({'number': 1, 'id':  1, 'state': INITIAL}),
        base_objects.Leg({'number': 2, 'id': -1, 'state': INITIAL}),
        base_objects.Leg({'number': 3, 'id': 21, 'state': FINAL}),
        base_objects.Leg({'number': 4, 'id': 21, 'state': FINAL}),
        base_objects.Leg({'number': 5, 'id': 21, 'state': FINAL}),
        base_objects.Leg({'number': 6, 'id': 25, 'state': FINAL}),
    ])
    qqx_to_ggH = base_objects.Process({
        'legs': qqx_to_ggH_legs,
        'model': simple_qcd.model,
        'n_loops': 0
    })

    # Test low-level approach limit
    #=====================================================================================

    def test_low_level_approach_limit(self):

        K = subtraction.SingularStructure
        C = subtraction.CollStructure
        S = subtraction.SoftStructure
        def L(n, state=FINAL):
            return subtraction.SubtractionLeg(n, 0, state)
        limit1 = [
            (walkers.mappings.FinalGroupingMapping, K(C(L(3), L(4)), L(2), L(5), L(6)), 0.,),
            (walkers.mappings.FinalGroupingMapping, K(C(L(5), L(9)), L(2), L(6)), 1.,)
        ]
        self._test_low_level_approach_limit(self.H_to_bbxuuxggg, limit1)

    # Test NLO walkers
    #=====================================================================================

    def test_FinalRescalingOneWalker_invertible(self):

        walker = walkers.FinalRescalingOneWalker()
        self._test_invertible(walker, self.H_to_uuxddx, 1, 1)
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)

    def test_FinalLorentzOneWalker_invertible(self):

        walker = walkers.FinalLorentzOneWalker()
        self._test_invertible(walker, self.H_to_uuxddx, 1, 1)
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)

    def test_FinalRescalingNLOWalker_invertible(self):

        walker = walkers.FinalRescalingNLOWalker()
        self._test_invertible(walker, self.H_to_uuxddx, 1, 1)
        self._test_invertible(walker, self.H_to_uuxddxg, 1, 1)
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)
        self._test_invertible(walker, self.H_to_qqxggH, 1, 1)

    def test_FinalRescalingNLOWalker_approach_limit(self):

        walker = walkers.FinalRescalingNLOWalker()
        self._test_approach_limit(walker, self.H_to_uuxddx, 1, 1)
        self._test_approach_limit(walker, self.H_to_uuxddxg, 1, 1)
        self._test_approach_limit(walker, self.H_to_uuxddxH, 1, 1)
        self._test_approach_limit(walker, self.H_to_qqxggH, 1, 1)

    def test_FinalLorentzNLOWalker_invertible(self):

        walker = walkers.FinalLorentzNLOWalker()
        self._test_invertible(walker, self.H_to_uuxddxg, 1, 1)
        self._test_invertible(walker, self.H_to_qqxggH, 1, 1)

    def test_FinalLorentzNLOWalker_approach_limit(self):

        walker = walkers.FinalLorentzNLOWalker()
        self._test_approach_limit(walker, self.H_to_uuxddxH, 1, 1)
        self._test_approach_limit(walker, self.H_to_qqxggH, 1, 1)

    def test_LorentzNLOWalker_invertible(self):

        walker = walkers.LorentzNLOWalker()
        self._test_approach_limit(walker, self.H_to_uuxddxg, 1, 1)
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)
        self._test_invertible(walker, self.H_to_qqxggH, 1, 1)
        self._test_invertible(walker, self.qqx_to_ggH, 1, 1)

    def test_LorentzNLOWalker_approach_limit(self):

        walker = walkers.LorentzNLOWalker()
        self._test_approach_limit(walker, self.H_to_uuxddxg, 1, 1)
        # The following fails because the soft mapping does not handle massive particles
        # self._test_approach_limit(walker, self.H_to_uuxddxH, 1, 1)
        # self._test_approach_limit(walker, self.H_to_qqxggH, 1, 1)
        # self._test_approach_limit(walker, self.qqx_to_ggH, 1, 1)

    def test_SoftBeamsRecoilNLOWalker_invertible(self):

        walker = walkers.SoftBeamsRecoilNLOWalker()
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)
        # TODO understand why this fails
        # self._test_invertible(walker, self.H_to_qqxggH, 1, 1)
        self._test_invertible(walker, self.qqx_to_ggH, 1, 1)

    def test_SoftBeamsRecoilNLOWalker_approach_limit(self):

        walker = walkers.SoftBeamsRecoilNLOWalker()
        self._test_approach_limit(walker, self.H_to_qqxggH, 1, 1)
        self._test_approach_limit(walker, self.qqx_to_ggH, 1, 1)

    # Test disjoint walkers
    #=====================================================================================

    def test_FinalLorentzDisjointWalker_invertible(self):

        walker = walkers.FinalLorentzDisjointWalker()
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 2)

    def test_FinalLorentzDisjointWalker_approach_limit(self):

        walker = walkers.FinalLorentzDisjointWalker()
        self._test_approach_limit(walker, self.H_to_uuxddxH, 1, 2)

    def test_FinalGroupingDisjointWalker_invertible(self):

        walker = walkers.FinalGroupingDisjointWalker()
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 2)

    def test_FinalGroupingDisjointWalker_approach_limit(self):

        walker = walkers.FinalGroupingDisjointWalker()
        self._test_approach_limit(walker, self.H_to_uuxddxH, 1, 2)
