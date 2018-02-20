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

    test.assertEqual(dict1.keys(), dict2.keys())
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
            # Compute total momentum
            total_momentum = LorentzVector()
            for i in children:
                total_momentum += my_PS_point[i]
            # Compute new phase space point
            new_PS_point = LorentzVectorDict()
            mappings.FinalCollinearVariables.set(
                new_PS_point, children, total_momentum, na, nb, variables )
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
                # Compute total momentum
                total_momentum = LorentzVector()
                for i in children:
                    total_momentum += my_PS_point[i]
                # Compute new phase space point
                new_PS_point = LorentzVectorDict()
                mappings.FinalCollinearVariables.set(
                    new_PS_point, children, total_momentum, na, nb, variables )
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
            is_child = 0
            fs_children = tuple(range(1, self.n_children))
            my_PS_point[is_child] = random_momentum(0)
            # my_PS_point[is_child] = LorentzVector([1.,0.,0.,1])
            for i in fs_children:
                if self.massive: my_PS_point[i] = random_momentum()
                else: my_PS_point[i] = random_momentum(0)
            # Generate two random light-cone directions
            na, nb = mappings.InitialCollinearVariables.collinear_and_reference(my_PS_point[is_child])
            # Compute collinear variables
            variables = dict()
            mappings.InitialCollinearVariables.get(
                my_PS_point, fs_children, is_child, na, nb, variables )
            # Compute new phase space point
            new_PS_point = LorentzVectorDict()
            mappings.InitialCollinearVariables.set(
                new_PS_point, fs_children, is_child, my_PS_point[is_child], na, nb, variables )
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
                mappings.InitialCollinearVariables.set(
                    new_PS_point, fs_children, is_child,
                    my_PS_point[is_child], na, nb, variables )
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
# MappingTest
#=========================================================================================

class MappingTest(object):
    """Collection of functions to test mappings."""

    # Leaving the IS not aligned with z avoids numerical issues
    v1 = LorentzVector([4., 0., 0., math.sqrt(3.)])
    v2 = LorentzVector([4., 1., 1., 1.])

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
            random.randint(2, pars['max_unresolved_per_set']+1)
            for _ in range(n_coll_sets) ]
        n_soft_per_set = [
            random.randint(1, pars['max_unresolved_per_set'])
            for _ in range(n_soft_sets) ]
        initial_sets = pars.get('initial_sets', 0)
        n_collinears_1 = 1
        n_collinears_2 = 1
        if initial_sets == 1:
            if bool(random.getrandbits(1)):
                n_collinears_1 = random.randint(2, pars['max_unresolved_per_set']+1)
            else:
                n_collinears_2 = random.randint(2, pars['max_unresolved_per_set']+1)
        elif initial_sets > 1:
            n_collinears_1 = random.randint(2, pars['max_unresolved_per_set']+1)
            n_collinears_2 = random.randint(2, pars['max_unresolved_per_set']+1)
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
                pars['momenta_dict'][n_tot + n_set + 1] = frozenset(numbers_in_this_set)
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

    @staticmethod
    def test_invertible(pars, test):
        """Test mapping and inverse."""

        # Make test deterministic by setting seed
        random.seed(42)
        # Check many times
        for _ in range(20):
            # Generate a random setup
            MappingTest.randomize(pars)
            my_PS_point = MappingTest.generate_PS(pars)
            # Rotate it to avoid zero components
            for key in my_PS_point.keys():
                my_PS_point[key].rotoboost(MappingTest.v1, MappingTest.v2)
            # I know what I'm doing
            old_PS_point = copy.deepcopy(my_PS_point)
            # Compute collinear variables
            variables = dict()
            lres = pars['mapping'].map_to_lower_multiplicity(
                my_PS_point, pars['structure'], pars['momenta_dict'], variables )
            hres = pars['mapping'].map_to_higher_multiplicity(
                my_PS_point, pars['structure'], pars['momenta_dict'], variables )
            # print pars['structure'], pars['unchanged']
            # print my_PS_point
            # print old_PS_point
            assertDictAlmostEqual(test, my_PS_point, old_PS_point)
            assertDictAlmostEqual(test, lres, hres)

#=========================================================================================
# Test final-collinear mappings
#=========================================================================================

class FinalRescalingMappingOneTest(unittest.TestCase):
    """Test class for FinalRescalingMappingOne."""

    # Test settings
    pars = {
        'mapping': mappings.FinalRescalingMappingOne(),
        'min_coll_sets': 1, 'max_coll_sets': 1,
        'min_recoilers': 1, 'max_recoilers': 3,
        'max_unchanged': 3,
        'max_unresolved_per_set': 4,
        'supports_massive_recoilers': False }

    def test_FinalRescalingMappingOne_invertible(self):

        MappingTest.test_invertible(self.pars, self)

class FinalLorentzMappingOneTest(unittest.TestCase):
    """Test class for FinalLorentzMappingOne."""

    # Test settings
    pars = {
        'mapping': mappings.FinalLorentzMappingOne(),
        'min_coll_sets': 1, 'max_coll_sets': 1,
        'min_recoilers': 1, 'max_recoilers': 5,
        'max_unchanged': 3,
        'max_unresolved_per_set': 4,
        'supports_massive_recoilers': True }

    def test_FinalLorentzMappingOne_invertible(self):

        MappingTest.test_invertible(self.pars, self)

#=========================================================================================
# Test initial-collinear mappings
#=========================================================================================

class InitialLorentzMappingOneTest(unittest.TestCase):
    """Test class for InitialLorentzMappingOne."""

    # Test settings
    pars = {
        'mapping': mappings.InitialLorentzMappingOne(),
        'initial_sets': 1,
        'max_unresolved_per_set': 4,
        'min_recoilers': 1, 'max_recoilers': 5,
        'max_unchanged': 3,
        'supports_massive_recoilers': True}

    def test_InitialLorentzMappingOne_invertible(self):

        MappingTest.test_invertible(self.pars, self)

#=========================================================================================
# Test soft mappings
#=========================================================================================

class SomogyietalSoftTest(unittest.TestCase):
    """Test class for MappingSomogyietalSoft."""

    # Test settings
    pars = {
        'mapping': mappings.MappingSomogyietalSoft(),
        'min_soft_sets': 1, 'max_soft_sets': 3,
        'min_recoilers': 2, 'max_recoilers': 5,
        'max_unchanged': 3,
        'max_unresolved_per_set': 4,
        'supports_massive_recoilers': False }

    def test_SomogyietalSoft_invertible(self):

        MappingTest.test_invertible(self.pars, self)

#=========================================================================================
# Test the phase-space walkers
#=========================================================================================

# Functions
#=========================================================================================

class WalkerTest(object):
    """Test class for walkers."""

    verbose = True
    irs = subtraction.IRSubtraction(
        simple_qcd.model, coupling_types=('QCD', ), n_unresolved=1 )
    parameter_values = [math.pow(0.1, i) for i in range(5)]
    max_unresolved = None
    initial_along_z = False # Leaving the IS not aligned with z avoids numerical issues
    # v1 = LorentzVector([4, 0, 0, 3])
    # v2 = LorentzVector([4, 1, 1, 1])

    @classmethod
    def generate_PS_point(cls, process):
        """Generate a phase-space point to test the walker."""

        model = process.get('model')
        # is_masses = tuple(simple_qcd.masses[model.get_particle(pdg)['mass']]
        #                   for pdg in process.get_initial_ids())
        # fs_masses = tuple(simple_qcd.masses[model.get_particle(pdg)['mass']]
        #                   for pdg in process.get_final_ids())
        # E_cm = sum(fs_masses) / random.random()
        # gen = PS.FlatInvertiblePhasespace(
        #     is_masses, fs_masses, beam_Es=(E_cm / 2., E_cm / 2.), beam_types=(0, 0) )
        # randoms = [random.random() for _ in range(gen.nDimPhaseSpace())]
        # my_PS_point, _ = gen.generateKinematics(E_cm, randoms)
        # for key in my_PS_point.keys():
        #     my_PS_point[key].rotoboost(cls.v1, cls.v2)
        # return my_PS_point.to_dict()

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
        else: raise
        return my_PS_point

    @classmethod
    def test_invertible(cls, test, walker, process):
        """Test walk and inverse."""

        if cls.verbose:
            print "Testing process", process.nice_string()
        my_operators = cls.irs.get_all_elementary_operators(process)
        my_combinations = cls.irs.get_all_combinations(my_operators, cls.max_unresolved)
        my_counterterms = [
            cls.irs.get_counterterm(combination, process)
            for combination in my_combinations ]

        # For each counterterm
        for i in range(len(my_counterterms)):
            if cls.verbose:
                print "Considering counterterm", my_counterterms[i]
            my_PS_point = cls.generate_PS_point(process)
            # Compute collinear variables
            res_dict1 = walker.walk_to_lower_multiplicity(
                my_PS_point, my_counterterms[i], True )
            (currs1, ME1, mv1, kin)  = (
                res_dict1['currents'], res_dict1['matrix_element'],
                res_dict1['mapping_variables'], res_dict1['kinematic_variables'] )
            if cls.verbose:
                print "Walking down"
                for curr in currs1:
                    print curr[1]
                print ME1[1]
            res_dict2 = walker.walk_to_higher_multiplicity(
                ME1[1], my_counterterms[i], kin )
            (currs2, ME2, mv2)  = (
                res_dict2['currents'], res_dict2['matrix_element'],
                res_dict2['mapping_variables'] )
            if cls.verbose:
                print "Walking up"
                print ME2[1]
                for curr in currs2:
                    print curr[1]
            # Check currents
            test.assertEqual(len(currs1), len(currs2))
            for i_curr in range(len(currs1)):
                test.assertEqual(currs1[i_curr][0], currs2[i_curr][0])
                test.assertDictEqual(currs1[i_curr][1], currs2[i_curr][1])
            # Check MEs
            test.assertEqual(ME1[0], ME2[0])
            test.assertDictEqual(ME1[1], ME2[1])
            # Check mapping variables
            assertDictAlmostEqual(test, mv1, mv2)

    @classmethod
    def test_approach_limit(cls, test, walker, process):
        """Test limit approach."""

        counterterms = cls.irs.get_all_counterterms(process)[0][1:]

        legs_FS = tuple(
            subtraction.SubtractionLeg(leg)
            for leg in process['legs']
            if leg['state'] == FINAL )

        # For each counterterm
        for ct in counterterms:
            if cls.verbose:
                print "\n" + "*" * 100
                print "Considering counterterm", ct
                print "*" * 100 + "\n"
            # Generate random vectors
            my_PS_point = cls.generate_PS_point(process)
            if cls.verbose:
                print "Starting phase space point:\n", my_PS_point, "\n"
            squares = {key: my_PS_point[key].square() for key in my_PS_point.keys()}
            # Compute collinear variables
            for alpha in cls.parameter_values:
                new_PS_point = my_PS_point.get_copy()
                walker.approach_limit(new_PS_point, ct, alpha, counterterms)
                if cls.verbose:
                    print "New PS point for", alpha, ":\n", new_PS_point
                for leg in legs_FS:
                    if process.get('model').get_particle(leg.pdg)['mass'].lower() == 'zero':
                        test.assertLess(
                            abs(new_PS_point[leg.n].square()),
                            math.sqrt(new_PS_point[leg.n].eps()) )
                    else:
                        test.assertAlmostEqual(new_PS_point[leg.n].square(), squares[leg.n])

# Processes
#=========================================================================================

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

# HACK: one gluon more to avoid issue with soft mapping

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

# Tests
#=========================================================================================

class FlatCollinearWalkerTest(unittest.TestCase):
    """Test class for FlatCollinearWalker."""

    walker = mappings.FlatCollinearWalker()

    def test_FlatCollinearWalker_invertible(self):

        WalkerTest.max_unresolved = 1
        random.seed(42)
        WalkerTest.test_invertible(self, self.walker, H_to_uuxddxH)

class FinalNLOWalkerTest(unittest.TestCase):
    """Test class for FinalNLOWalker."""

    walker = mappings.FinalNLOWalker()

    def test_FinalNLOWalker_invertible(self):

        WalkerTest.max_unresolved = 1
        WalkerTest.test_invertible(self, self.walker, H_to_uuxddxH)
        WalkerTest.test_invertible(self, self.walker, H_to_qqxggH)

    def test_FinalNLOWalker_approach_limit(self):

        WalkerTest.test_approach_limit(self, self.walker, H_to_qqxggH)

    # def test_sc_approach_limit(self):
    #
    #     # Set up a soft-collinear counterterm
    #     sc_ll = base_objects.LegList([
    #         base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
    #         base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
    #         base_objects.Leg({'number': 5, 'id': -1, 'state': FINAL}), ])
    #     sc_rp = base_objects.Process({
    #         'legs': sc_ll, 'model': simple_qcd.model, 'n_loops': 0 })
    #     sc_ss = subtraction.CollStructure(
    #         legs=subtraction.SubtractionLegSet((
    #             subtraction.SubtractionLeg(3, -1, FINAL), )),
    #         substructures=[
    #             subtraction.SoftStructure(
    #                 subtraction.SubtractionLeg(4, 21, FINAL)) ] )
    #     sc_md = subtraction.bidict({i: frozenset((i, )) for i in range(1, 5)})
    #     sc_md[5] = frozenset((3, 4, ))
    #     sc_ct = subtraction.Counterterm(
    #         process=sc_rp,
    #         nodes=[
    #             subtraction.CountertermNode(
    #                 current=subtraction.Current({
    #                     'singular_structure': sc_ss }) ) ],
    #         prefactor=1,
    #         momenta_dict=sc_md )
    #
    #     # Start from a random phase space point
    #     # The Higgs is going to have a random mass, but it doesn't matter
    #     PS_point = LorentzVectorDict()
    #     PS_point[1] = LorentzVector()
    #     for i in range(2, 5):
    #         PS_point[i] = random_momentum(0)
    #         PS_point[1] += PS_point[i]
    #
    #     hike_down = self.walker.walk_to_lower_multiplicity(PS_point, sc_ct, True)
    #     lower_PS_point = hike_down['resulting_PS_point']
    #     starting_variables = hike_down['kinematic_variables']
    #
    #     print sc_ss
    #     print sc_ct
    #
    #     ratios = []
    #     flucts = []
    #     for par in range(10):
    #         x = math.pow(0.1, par)
    #         hike_up = self.walker.approach_limit(
    #             lower_PS_point, sc_ct, starting_variables, x )
    #         ratios_vec = LorentzVector([
    #             hike_up['resulting_PS_point'][3][i] / hike_up['resulting_PS_point'][4][i]
    #             for i in range(4) ])
    #         ratios.append(ratios_vec.view(type=Vector).square())
    #         norm_ratios = ratios_vec.view(type=Vector)
    #         norm_ratios.normalize()
    #         flucts.append(abs(norm_ratios - Vector(4 * [0.5, ])))
    #
    #     print ratios
    #
    #     # Skipping the first few values, not close enough to the limit
    #     # Numerical effects can make the test fail for the deep IR region
    #     for i in range(2, len(ratios)-2):
    #         self.assertLess(
    #             abs(ratios[i + 2] / ratios[i + 1] - 10.),
    #             abs(ratios[i + 1] / ratios[i] - 10.) )
    #     for i in range(2, len(flucts)-1):
    #         self.assertLess(flucts[i+1], flucts[i])

class NLOWalkerTest(unittest.TestCase):
    """Test class for NLOWalker."""

    walker = mappings.NLOWalker()

    def test_NLOWalker_invertible(self):

        WalkerTest.max_unresolved = 1
        # Make test deterministic by setting seed
        random.seed(42)
        # Check many times
        for _ in range(5):
            WalkerTest.test_invertible(self, self.walker, H_to_uuxddxH)
            WalkerTest.test_invertible(self, self.walker, H_to_qqxggH)
            WalkerTest.test_invertible(self, self.walker, qqx_to_ggH)

    def test_NLOWalker_approach_limit(self):

        # Make test deterministic by setting seed
        random.seed(42)
        # Check many times
        for _ in range(5):
            WalkerTest.test_approach_limit(self, self.walker, H_to_qqxggH)
            WalkerTest.test_approach_limit(self, self.walker, qqx_to_ggH)
