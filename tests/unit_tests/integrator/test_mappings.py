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

import madgraph.integrator.mappings as mappings
import madgraph.core.subtraction as subtraction
import madgraph.core.base_objects as base_objects
import madgraph.integrator.phase_space_generators as PS
from madgraph.integrator.vectors import Vector, LorentzVector
from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList

import copy
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
# Test final-collinear variables
#=========================================================================================

class FinalCollinearVariablesTest(unittest.TestCase):
    """Test class for variables describing internal collinear structure."""

    n_children = 3
    massive = True

    def test_final_collinear_variables_away_from_limit(self):
        """Test determination of collinear variables and reverse mapping,
        for completely generic input values.
        """

        # Make test deterministic by setting seed
        random.seed(42)
        # Check many times
        for _ in range(100):
            # Generate n_children random Lorentz vectors
            my_PS_point = LorentzVectorDict()
            parent = self.n_children
            children = tuple(range(self.n_children))
            for i in children:
                if (self.massive): my_PS_point[i] = random_momentum()
                else: my_PS_point[i] = random_momentum(0)
            # Generate two random light-cone directions
            na = random_momentum(0)
            nb = random_momentum(0)
            # Compute collinear variables
            variables = dict()
            mappings.FinalCollinearVariables.get(
                my_PS_point, children, na, nb, variables )
            # Compute new phase space point
            new_PS_point = LorentzVectorDict()
            new_PS_point[parent] = sum(my_PS_point[child] for child in children)
            mappings.FinalCollinearVariables.set(
                new_PS_point, parent, children, na, nb, variables )
            # Check the two phase-space points are equal
            self.assertDictEqual(my_PS_point, new_PS_point)

    def test_final_collinear_variables_close_to_limit(self):
        """Test determination of collinear variables and reverse mapping
        in a typical collinear situation.
        """

        # Make test deterministic by setting seed
        random.seed(42)
        # Check many times
        for _ in range(20):
            parent = self.n_children
            children = tuple(range(self.n_children))
            # Generate children squares
            if self.massive:
                squares = {i: random.random() for i in children}
            else:
                squares = {i: 0. for i in children}
            # Generate children random starting directions
            directions = {
                i: Vector([random.random() for _ in range(3)])
                for i in children }
            coll_direction = Vector(3*[0., ])
            for n in directions.values():
                coll_direction += n
            coll_lorentz = LorentzVector([0., ] + list(coll_direction)).set_square(0)
            # Generate values for the parameter that describes approach to limit
            pars = (math.pow(0.25, i) for i in range(14))
            # This is one way like any other to approach the limit
            for par in pars:
                new_directions = {
                    i: (par*n + (1-par)*coll_direction)
                    for (i, n) in directions.items() }
                my_PS_point = {
                    i: LorentzVector([0., ] + list(n)).set_square(squares[i])
                    for (i, n) in new_directions.items() }
                na, nb = mappings.FinalCollinearVariables.collinear_and_reference(
                    coll_lorentz )
                # Compute collinear variables
                variables = dict()
                mappings.FinalCollinearVariables.get(
                    my_PS_point, children, na, nb, variables)
                # Compute new phase space point
                new_PS_point = LorentzVectorDict()
                new_PS_point[parent] = sum(my_PS_point[child] for child in children)
                mappings.FinalCollinearVariables.set(
                    new_PS_point, parent, children, na, nb, variables)
                # Check the two phase-space points are equal
                self.assertDictEqual(my_PS_point, new_PS_point)

#=========================================================================================
# Test initial-collinear variables
#=========================================================================================

class InitialCollinearVariablesTest(unittest.TestCase):
    """Test class for variables describing internal initial-collinear structure."""

    n_children = 3
    massive = False

    def test_initial_collinear_variables_away_from_limit(self):
        """Test initial-collinear variables getting and setting,
        for completely generic input values.
        """

        # Make test deterministic by setting seed
        random.seed(42)
        # Check many times
        for _ in range(100):
            # Generate n_children random Lorentz vectors
            my_PS_point = LorentzVectorDict()
            parent = self.n_children
            is_child = 0
            fs_children = tuple(range(1, self.n_children))
            my_PS_point[is_child] = random_momentum(0)
            # my_PS_point[is_child] = LorentzVector([1.,0.,0.,1])
            for i in fs_children:
                if self.massive: my_PS_point[i] = random_momentum()
                else: my_PS_point[i] = random_momentum(0)
            # Generate two random light-cone directions
            na, nb = mappings.InitialCollinearVariables.collinear_and_reference(
                my_PS_point[is_child])
            # Compute collinear variables
            variables = dict()
            mappings.InitialCollinearVariables.get(
                my_PS_point, fs_children, is_child, na, nb, variables )
            # Compute new phase space point
            new_PS_point = LorentzVectorDict()
            new_PS_point[is_child] = my_PS_point[is_child]
            mappings.InitialCollinearVariables.set(
                new_PS_point, is_child, fs_children, na, nb, variables )
            # Check the two phase-space points are equal
            self.assertDictEqual(my_PS_point, new_PS_point)

    def test_initial_collinear_variables_close_to_limit(self):
        """Test initial-collinear variables getting and setting,
        in a typical collinear situation.
        """

        # Make test deterministic by setting seed
        random.seed(42)
        # Check many times
        for _ in range(20):
            is_child = 0
            fs_children = tuple(range(1, self.n_children))
            # Generate children squares
            if self.massive:
                squares = {i: random.random() for i in fs_children}
            else:
                squares = {i: 0. for i in fs_children}
            # Generate children random starting directions
            directions = {
                i: Vector([random.random() for _ in range(3)])
                for i in fs_children }
            coll_direction = Vector([random.random() for _ in range(3)])
            # coll_direction = Vector([0.,0.,1.])
            coll_lorentz = LorentzVector([0., ] + list(coll_direction)).set_square(0)
            na, nb = mappings.InitialCollinearVariables.collinear_and_reference(coll_lorentz)
            # Generate values for the parameter that describes approach to limit
            pars = (math.pow(0.25, i) for i in range(14))
            # This is one way like any other to approach the limit
            for par in pars:
                new_directions = {
                    i: (par*n + (1-par)*coll_direction)
                    for (i, n) in directions.items() }
                my_PS_point = {
                    i: LorentzVector([0., ] + list(n)).set_square(squares[i])
                    for (i, n) in new_directions.items() }
                my_PS_point[is_child] = coll_lorentz
                # Compute collinear variables
                variables = dict()
                mappings.InitialCollinearVariables.get(
                    my_PS_point, fs_children, is_child, na, nb, variables )
                # Compute new phase space point
                new_PS_point = LorentzVectorDict()
                new_PS_point[is_child] = my_PS_point[is_child]
                mappings.InitialCollinearVariables.set(
                    new_PS_point, is_child, fs_children,
                    na, nb, variables )
                # Check the two phase-space points are equal
                self.assertDictEqual(my_PS_point, new_PS_point)

#=========================================================================================
# Test soft variables
#=========================================================================================

class SoftVariablesTest(unittest.TestCase):
    """Test class for variables describing internal collinear structure."""

    # This is trivial if one stores the full momentum
    # Nevertheless checking and supporting more complicated scenarios

    my_mapping = mappings.ElementaryMappingSoft()
    n_soft = 3

    def test_soft_variables_away_from_limit(self):
        """Test determination of soft variables and reverse mapping,
        for completely generic input values.
        """

        # Make test deterministic by setting seed
        random.seed(42)
        # Check many times
        for _ in range(100):
            softs = tuple(range(self.n_soft))
            # Generate n_soft random massless vectors
            my_PS_point = LorentzVectorDict()
            for i in softs:
                my_PS_point[i] = random_momentum(0)
            # Compute soft variables
            variables = dict()
            mappings.SoftVariables.get(my_PS_point, softs, variables)
            # Compute new phase space point
            new_PS_point = LorentzVectorDict()
            mappings.SoftVariables.set(new_PS_point, softs, variables)
            # Check the two phase-space points are equal
            self.assertDictEqual(my_PS_point, new_PS_point)

    def test_soft_variables_close_to_limit(self):
        """Test determination of soft variables and reverse mapping
        in a typical soft situation.
        """

        # Make test deterministic by setting seed
        random.seed(42)
        # Check many times
        for _ in range(20):
            softs = tuple(range(self.n_soft))
            # Generate n_soft random massless vectors
            my_PS_point = LorentzVectorDict()
            for i in softs:
                my_PS_point[i] = random_momentum(0)
            # Generate values for the parameter that describes approach to limit
            pars = (math.pow(0.25, i) for i in range(8))
            # This is one way like any other to approach the limit
            for par in pars:
                old_PS_point = {
                    i: par*my_PS_point[i]
                    for i in my_PS_point.keys() }
                # Compute collinear variables
                variables = dict()
                mappings.SoftVariables.get(old_PS_point, softs, variables)
                # Compute new phase space point
                new_PS_point = LorentzVectorDict()
                mappings.SoftVariables.set(new_PS_point, softs, variables)
                self.assertDictEqual(old_PS_point, new_PS_point)

#=========================================================================================
# MappingsTest
#=========================================================================================

class MappingsTest(unittest.TestCase):
    """Collection of functions to test mappings."""

    # Leaving the IS not aligned with z avoids numerical issues
    v1 = LorentzVector([4., 0., 0., math.sqrt(3.)])
    v2 = LorentzVector([4., 1., 1., 1.])

    seed = 42
    n_tests = 20
    verbose = False

    @staticmethod
    def randomize(pars):

        # Randomly generate the structure of the final state
        n_coll_sets = random.randint(
            pars.get('min_coll_sets', 0), pars.get('max_coll_sets', 0) )
        n_soft_sets = random.randint(
            pars.get('min_soft_sets', 0), pars.get('max_soft_sets', 0) )
        n_recoilers = random.randint(pars['min_recoilers'], pars['max_recoilers'])
        n_unchanged = random.randint(0, pars['max_unchanged'])
        n_coll_per_set = [
            random.randint(
                pars['min_unresolved_per_set']+1, pars['max_unresolved_per_set']+1)
            for _ in range(n_coll_sets) ]
        n_soft_per_set = [
            random.randint(
                pars['min_unresolved_per_set'], pars['max_unresolved_per_set'])
            for _ in range(n_soft_sets) ]
        initial_sets = pars.get('initial_sets', 0)
        n_collinears_1 = 1
        n_collinears_2 = 1
        if initial_sets == 1:
            if bool(random.getrandbits(1)):
                n_collinears_1 = random.randint(
                    pars['min_unresolved_per_set']+1, pars['max_unresolved_per_set']+1)
            else:
                n_collinears_2 = random.randint(
                    pars['min_unresolved_per_set']+1, pars['max_unresolved_per_set']+1)
        elif initial_sets > 1:
            n_collinears_1 = random.randint(
                pars['min_unresolved_per_set']+1, pars['max_unresolved_per_set']+1)
            n_collinears_2 = random.randint(
                pars['min_unresolved_per_set']+1, pars['max_unresolved_per_set']+1)
        # Compute the total number of particles in the final state
        n_tot = sum((n_collinears_1 + n_collinears_2,
                     sum(n_coll_per_set) + sum(n_soft_per_set),
                     n_recoilers + n_unchanged))
        pars['n_tot'] = n_tot
        # Initialise trivial entries {i: i} of the momenta dictionary
        pars['momenta_dict'] = subtraction.bidict()
        for i in range(1, n_tot+1):
            pars['momenta_dict'][i] = frozenset((i, ))
        # Randomly generate nonconsecutive numbers of the final state
        shuffled_numbers = range(3, n_tot+1)
        random.shuffle(shuffled_numbers)
        pars['structure'] = []
        pars['parents'] = []
        n_set = 1
        # Build the structure of the set collinear to 1 (if any)
        if n_collinears_1 > 1:
            numbers_collinears_1 = shuffled_numbers[:n_collinears_1-1]
            shuffled_numbers = shuffled_numbers[n_collinears_1-1:]
            legs_in_this_set = [
                subtraction.SubtractionLeg(1, 21, INITIAL) ] + [
                subtraction.SubtractionLeg(i, 21, FINAL)
                for i in numbers_collinears_1 ]
            legs_in_this_set = subtraction.SubtractionLegSet(legs_in_this_set)
            pars['structure'] += [ subtraction.CollStructure(legs=legs_in_this_set), ]
            pars['momenta_dict'][n_tot + n_set] = frozenset([1, ] + numbers_collinears_1)
            pars['parents'].append(n_tot + n_set)
            n_set += 1
        # Build the structure of the set collinear to 2 (if any)
        if n_collinears_2 > 1:
            numbers_collinears_2 = shuffled_numbers[:n_collinears_2-1]
            shuffled_numbers = shuffled_numbers[n_collinears_2-1:]
            legs_in_this_set = [
                subtraction.SubtractionLeg(2, 21, INITIAL) ] + [
                subtraction.SubtractionLeg(i, 21, FINAL)
                for i in numbers_collinears_2 ]
            legs_in_this_set = subtraction.SubtractionLegSet(legs_in_this_set)
            pars['structure'] += [ subtraction.CollStructure(legs=legs_in_this_set), ]
            pars['momenta_dict'][n_tot + n_set] = frozenset([2, ] + numbers_collinears_2)
            pars['parents'].append(n_tot + n_set)
            n_set += 1
        # Build the structure of final-collinear sets
        for n_in_this_set in n_coll_per_set:
            numbers_in_this_set = shuffled_numbers[:n_in_this_set]
            shuffled_numbers = shuffled_numbers[n_in_this_set:]
            legs_in_this_set = subtraction.SubtractionLegSet(
                subtraction.SubtractionLeg(i, 21, FINAL)
                for i in numbers_in_this_set)
            pars['structure'] += [subtraction.CollStructure(legs=legs_in_this_set), ]
            if n_in_this_set != 1:
                pars['momenta_dict'][n_tot + n_set] = frozenset(numbers_in_this_set)
                pars['parents'].append(n_tot + n_set)
            else:
                pars['parents'].append(numbers_in_this_set[0])
            n_set += 1
        # Build the structure of soft sets
        for n_in_this_set in n_soft_per_set:
            numbers_in_this_set = shuffled_numbers[:n_in_this_set]
            shuffled_numbers = shuffled_numbers[n_in_this_set:]
            legs_in_this_set = subtraction.SubtractionLegSet(
                subtraction.SubtractionLeg(i, 21, FINAL)
                for i in numbers_in_this_set )
            pars['structure'] += [subtraction.SoftStructure(legs=legs_in_this_set), ]
        # Include recoilers and unchanged
        pars['structure'] += [
            subtraction.SubtractionLeg(i, 21, FINAL)
            for i in shuffled_numbers[:n_recoilers]]
        pars['structure'] = subtraction.SingularStructure(*pars['structure'])
        pars['unchanged'] = shuffled_numbers[n_recoilers:]
        # print "The random structure for this test is %s, with %s unchanged" % (
        #     str(pars['structure']), str(pars['unchanged']) )
        assert len(pars['unchanged']) == n_unchanged

    @staticmethod
    def generate_PS(pars):
        """Generate a valid random phase-space point for the mapping test."""

        masses = []
        for i in range(3, pars['n_tot']+1):
            # Recoilers
            if i in pars['structure'].legs:
                if pars['supports_massive_recoilers']: masses.append(random.random())
                else: masses.append(0.)
            # Soft particles or unchanged
            else:
                masses.append(random.random())
        E_cm = sum(masses) / random.random()
        gen = PS.FlatInvertiblePhasespace(
            (0., 0.), masses, beam_Es=(E_cm/2., E_cm/2.), beam_types=(0, 0) )
        randoms = [random.random() for _ in range(gen.nDimPhaseSpace())]
        my_PS_point, _ = gen.generateKinematics(E_cm, randoms)
        return my_PS_point.to_dict()

    def _test_invertible(self, pars):
        """Test mapping and inverse."""

        # Make test deterministic by setting seed
        random.seed(self.seed)
        # Check many times
        for _ in range(self.n_tests):
            # Generate a random setup
            MappingsTest.randomize(pars)
            if self.verbose:
                print "Structure:", pars['structure']
            my_PS_point = MappingsTest.generate_PS(pars)
            squared_masses = dict()
            for parent in pars['parents']:
                if pars['masses']:
                    p = sum(my_PS_point[child] for child in pars['momenta_dict'][parent])
                    squared_masses['m2' + str(parent)] = random.random()*p.square()
                else:
                    squared_masses['m2' + str(parent)] = 0.
            # Rotate it to avoid zero components
            for key in my_PS_point.keys():
                my_PS_point[key].rotoboost(MappingsTest.v1, MappingsTest.v2)
            if self.verbose:
                print "Starting PS point:\n", my_PS_point
            # Compute collinear variables
            variables = dict()
            low_PS_point, low_jac = pars['mapping'].map_to_lower_multiplicity(
                my_PS_point, pars['structure'], pars['momenta_dict'], squared_masses,
                variables, True )
            if self.verbose:
                print "Mapped PS point:\n", low_PS_point
                print "with jacobian:", low_jac
            high_PS_point, high_jac = pars['mapping'].map_to_higher_multiplicity(
                low_PS_point, pars['structure'], pars['momenta_dict'],
                variables, True )
            if self.verbose:
                print "Unmapped PS point:\n", high_PS_point
                print "with jacobian:", high_jac
            assertDictAlmostEqual(self, my_PS_point, high_PS_point)
            self.assertAlmostEqual(low_jac, high_jac)

    # Test masses mappings
    #=====================================================================================

    def test_FinalZeroMassesMapping_invertible(self):
        """Test if FinalZeroMassesMapping is invertible."""

        pars = {
            'mapping': mappings.FinalZeroMassesMapping(),
            'min_coll_sets': 2, 'max_coll_sets': 5,
            'min_recoilers': 0, 'max_recoilers': 0,
            'max_unchanged': 0, 'masses': None,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 0,
            'supports_massive_recoilers': False}
        self._test_invertible(pars)

    def test_FinalMassesMapping_invertible(self):
        """Test if FinalMassesMapping is invertible."""

        pars = {
            'mapping': mappings.FinalMassesMapping(),
            'min_coll_sets': 2, 'max_coll_sets': 5,
            'min_recoilers': 0, 'max_recoilers': 0,
            'max_unchanged': 0, 'masses': True,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 0,
            'supports_massive_recoilers': False}
        self._test_invertible(pars)

    # Test final-collinear mappings
    #=====================================================================================

    def test_FinalRescalingOneMapping_invertible(self):
        """Test if FinalRescalingOneMapping is invertible."""

        pars = {
            'mapping': mappings.FinalRescalingOneMapping(),
            'min_coll_sets': 1, 'max_coll_sets': 1,
            'min_recoilers': 1, 'max_recoilers': 3,
            'max_unchanged': 3, 'masses': None,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': False}
        self._test_invertible(pars)

    def test_FinalLorentzOneMapping_invertible(self):
        """Test if FinalLorentzOneMapping is invertible."""

        pars = {
            'mapping': mappings.FinalLorentzOneMapping(),
            'min_coll_sets': 1, 'max_coll_sets': 1,
            'min_recoilers': 1, 'max_recoilers': 5,
            'max_unchanged': 3, 'masses': None,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': True}
        self._test_invertible(pars)


    def test_FinalGroupingMapping_invertible(self):
        """Test if FinalGroupingMapping is invertible."""

        pars = {
            'mapping': mappings.FinalGroupingMapping(),
            'max_unchanged': 3, 'masses': True,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 3,
            'supports_massive_recoilers': True}
        pars.update({
            'min_coll_sets': 1, 'max_coll_sets': 3,
            'min_recoilers': 1, 'max_recoilers': 4,})
        self._test_invertible(pars)
        pars.update({
            'min_coll_sets': 2, 'max_coll_sets': 3,
            'min_recoilers': 0, 'max_recoilers': 3,})
        self._test_invertible(pars)

    def test_FinalLorentzMapping_invertible(self):
        """Test if FinalLorentzMapping is invertible."""

        pars = {
            'mapping': mappings.FinalLorentzMapping(),
            'max_unchanged': 3, 'masses': True,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 3,
            'supports_massive_recoilers': True}
        pars.update({
            'min_coll_sets': 1, 'max_coll_sets': 3,
            'min_recoilers': 1, 'max_recoilers': 4,})
        self._test_invertible(pars)
        pars.update({
            'min_coll_sets': 2, 'max_coll_sets': 3,
            'min_recoilers': 0, 'max_recoilers': 3,})
        self._test_invertible(pars)

    # Test initial-collinear mappings
    #=====================================================================================

    def test_InitialLorentzOneMapping_invertible(self):
        """Test if InitialLorentzOneMapping is invertible."""

        pars = {
            'mapping': mappings.InitialLorentzOneMapping(),
            'initial_sets': 1,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'min_recoilers': 1, 'max_recoilers': 5,
            'max_unchanged': 3, 'masses': None,
            'supports_massive_recoilers': True}
        self._test_invertible(pars)

    # Test soft mappings
    #=====================================================================================

    def test_SoftVsFinalMapping_invertible(self):
        """Test if SoftVsFinalMapping is invertible."""

        pars = {
            'mapping': mappings.SoftVsFinalMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 3,
            'min_recoilers': 2, 'max_recoilers': 5,
            'max_unchanged': 3, 'masses': None,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': False}
        self._test_invertible(pars)

#=========================================================================================
# Test the phase-space walkers
#=========================================================================================

class WalkersTest(unittest.TestCase):
    """Test class for walkers."""

    # Parameters
    #=====================================================================================

    # Verbosity (silent = 0, max = 4)
    verbosity = 1
    # Random seed to make tests deterministic
    seed = 42
    # Number of PS points the invertibility test is run for (more = stronger, slower test)
    n_test_invertible = 3
    # Number of PS points the approach_limit test is run for (more = stronger, slower test)
    n_test_approach = 5
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
                E = abs(total_momentum)
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
        :type walker: mappings.VirtualWalker

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
                (currs1, ME1, jac1, kin)  = (
                    res_dict1['currents'], res_dict1['matrix_element'],
                    res_dict1['jacobian'], res_dict1['kinematic_variables'] )
                if self.verbosity > 3:
                    print "Walking down"
                    for curr in currs1:
                        print curr[1]
                    print ME1[1]
                res_dict2 = walker.walk_to_higher_multiplicity(
                    ME1[1], my_counterterms[i], kin,
                    compute_jacobian=True )
                (currs2, ME2, jac2)  = (
                    res_dict2['currents'], res_dict2['matrix_element'],
                    res_dict2['jacobian'] )
                if self.verbosity > 3:
                    print "Walking up"
                    print ME2[1]
                    for curr in currs2:
                        print curr[1]
                    print "Jacobians:", jac1, jac2
                # Check currents
                self.assertEqual(len(currs1), len(currs2))
                for i_curr in range(len(currs1)):
                    self.assertEqual(currs1[i_curr][0], currs2[i_curr][0])
                    self.assertDictEqual(currs1[i_curr][1], currs2[i_curr][1])
                # Check MEs
                self.assertEqual(ME1[0], ME2[0])
                self.assertDictEqual(ME1[1], ME2[1])
                # Check mapping variables
                self.assertAlmostEqual(jac1, jac2)

    def _test_approach_limit(
        self, walker, process,
        max_unresolved_in_elementary, max_unresolved_in_combination):
        """Check that the walker is capable of approaching limits.

        :param walker: Mapping walker to be tested
        :type walker: mappings.VirtualWalker

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
                                math.sqrt(new_PS_point[leg.n].eps()) )
                        else:
                            self.assertAlmostEqual(
                                new_PS_point[leg.n].square(),
                                squares[leg.n] )

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

    # Test NLO walkers
    #=====================================================================================

    def test_FinalRescalingOneWalker_invertible(self):

        walker = mappings.FinalRescalingOneWalker()
        self._test_invertible(walker, self.H_to_uuxddx, 1, 1)
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)

    def test_FinalLorentzOneWalker_invertible(self):

        walker = mappings.FinalLorentzOneWalker()
        self._test_invertible(walker, self.H_to_uuxddx, 1, 1)
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)

    def test_FinalRescalingNLOWalker_invertible(self):

        walker = mappings.FinalRescalingNLOWalker()
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)
        self._test_invertible(walker, self.H_to_qqxggH, 1, 1)

    def test_FinalRescalingNLOWalker_approach_limit(self):

        walker = mappings.FinalRescalingNLOWalker()
        self._test_approach_limit(walker, self.H_to_uuxddxH, 1, 1)
        self._test_approach_limit(walker, self.H_to_qqxggH, 1, 1)

    def test_FinalLorentzNLOWalker_invertible(self):

        walker = mappings.FinalLorentzNLOWalker()
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)
        self._test_invertible(walker, self.H_to_qqxggH, 1, 1)

    def test_FinalLorentzNLOWalker_approach_limit(self):

        walker = mappings.FinalLorentzNLOWalker()
        self._test_approach_limit(walker, self.H_to_uuxddxH, 1, 1)
        self._test_approach_limit(walker, self.H_to_qqxggH, 1, 1)

    def test_LorentzNLOWalker_invertible(self):

        walker = mappings.LorentzNLOWalker()
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 1)
        self._test_invertible(walker, self.H_to_qqxggH, 1, 1)
        self._test_invertible(walker, self.qqx_to_ggH, 1, 1)

    def test_LorentzNLOWalker_approach_limit(self):

        walker = mappings.LorentzNLOWalker()
        self._test_approach_limit(walker, self.H_to_qqxggH, 1, 1)
        self._test_approach_limit(walker, self.qqx_to_ggH, 1, 1)

    # Test disjoint walkers
    #=====================================================================================

    def test_FinalLorentzDisjointWalker_invertible(self):

        walker = mappings.FinalLorentzDisjointWalker()
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 2)

    def test_FinalLorentzDisjointWalker_approach_limit(self):

        walker = mappings.FinalLorentzDisjointWalker()
        self._test_approach_limit(walker, self.H_to_uuxddxH, 1, 2)

    def test_FinalGroupingDisjointWalker_invertible(self):

        walker = mappings.FinalGroupingDisjointWalker()
        self._test_invertible(walker, self.H_to_uuxddxH, 1, 2)

    def test_FinalGroupingDisjointWalker_approach_limit(self):

        walker = mappings.FinalGroupingDisjointWalker()
        self._test_approach_limit(walker, self.H_to_uuxddxH, 1, 2)
