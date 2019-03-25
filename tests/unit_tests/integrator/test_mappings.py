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
import madgraph.various.misc as misc
from madgraph.integrator.vectors import Vector, LorentzVector
from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList

import math
import random
import os

import tests.unit_tests as unittest

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
    verbose = 0
    n_skip = 0

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
            pars['structure'] += [ subtraction.CollStructure(legs=legs_in_this_set), ]
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
        leg_ns = [leg.n for leg in pars['structure'].legs]
        for i in range(3, pars['n_tot']+1):
            # Recoilers
            if i in leg_ns:
                if pars.get('supports_massive_recoilers', True):
                    masses.append(random.random())
                else: masses.append(0.)
            # Soft particles or unchanged
            else:
                if pars.get('supports_massive_unresolved', True):
                    masses.append(random.random())
                else: masses.append(0.)
        M = sum(masses)
        if M > 0.:
            E_cm = sum(masses) / random.random()
        else:
            E_cm = random.random()
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
        for n_test in range(pars.get('n_tests', self.n_tests)):
            # Generate a random setup
            MappingsTest.randomize(pars)
            if self.verbose:
                print "Test" , n_test+1, ": Structure" , str(pars['structure'])
            my_PS_point = MappingsTest.generate_PS(pars)
            squared_masses = dict()
            for parent in pars['parents']:
                if pars['masses']:
                    p = sum(my_PS_point[child] for child in pars['momenta_dict'][parent])
                    squared_masses['m2' + str(parent)] = random.random()*p.square()
                else:
                    squared_masses['m2' + str(parent)] = 0.
            if n_test < self.n_skip:
                print "Test", n_test+1, "skipped"
                continue
            # Rotate it to avoid zero components
            for key in my_PS_point.keys():
                my_PS_point[key].rotoboost(MappingsTest.v1, MappingsTest.v2)
            if self.verbose > 1:
                misc.sprint("Starting PS point:\n" + str(my_PS_point))
            # Compute collinear variables
            variables = dict()
            low_PS_point, low_vars = pars['mapping'].map_to_lower_multiplicity(
                my_PS_point, pars['structure'], pars['momenta_dict'], squared_masses,
                variables, pars.get('jacobian', True) )
            if self.verbose > 1:
                misc.sprint("Mapped PS point:\n" + str(low_PS_point))
                misc.sprint("with variables: " + str(low_vars))
            high_PS_point, high_vars = pars['mapping'].map_to_higher_multiplicity(
                low_PS_point, pars['structure'], pars['momenta_dict'],
                variables, pars.get('jacobian', True) )
            if self.verbose > 1:
                misc.sprint("Unmapped PS point:\n" + str(high_PS_point))
                misc.sprint("with variables: " + str(high_vars))
            assertDictAlmostEqual(self, my_PS_point, high_PS_point)
            assertDictAlmostEqual(self, low_vars, high_vars)

    def _test_equal_lower(self, pars):
        """Test whether the map_to_lower_multiplicity functions
         of two mappings are the same.
         """

        # Make test deterministic by setting seed
        random.seed(self.seed)
        # Check many times
        for _ in range(pars.get('n_tests', self.n_tests)):
            # Generate a random setup
            MappingsTest.randomize(pars)
            if self.verbose:
                misc.sprint("Structure: " + str(pars['structure']))
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
            if self.verbose > 1:
                misc.sprint("Starting PS point:\n" + str(my_PS_point))
            # Map to lower multiplicity using mapping1 and mapping2
            low_PS_point1, low_vars1 = pars['mapping1'].map_to_lower_multiplicity(
                my_PS_point, pars['structure'], pars['momenta_dict'], squared_masses,
                None, pars.get('jacobian', True) )
            low_PS_point2, low_vars2 = pars['mapping2'].map_to_lower_multiplicity(
                my_PS_point, pars['structure'], pars['momenta_dict'], squared_masses,
                None, pars.get('jacobian', True) )
            if self.verbose > 1:
                print pars['mapping1'].__class__.__name__ + " yields"
                print str(low_PS_point1)
                print "variables:", low_vars1
                print pars['mapping2'].__class__.__name__ + " yields"
                print str(low_PS_point2)
                print "variables:", low_vars2
            assertDictAlmostEqual(self, low_PS_point1, low_PS_point2)
            assertDictAlmostEqual(self, low_vars1, low_vars2)

    def _test_associative(self, pars):
        """Test whether a mapping is associative."""

        # Make test deterministic by setting seed
        random.seed(self.seed)
        # Check many times
        for _ in range(pars.get('n_tests', self.n_tests)):
            # Generate a random setup
            MappingsTest.randomize(pars)
            structure = pars['structure']
            momenta_dict = pars['momenta_dict']
            if self.verbose:
                print "Structure:", structure
            # Randomly divide the mapping in two steps
            half_1_substructures = []
            half_2_substructures = []
            half_1_legs = [leg for leg in structure.legs]
            half_2_legs = [leg for leg in structure.legs]
            for substructure in structure.substructures:
                if bool(random.getrandbits(1)):
                    half_1_substructures.append(substructure)
                    if substructure.name() == 'C':
                        children = frozenset([leg.n for leg in substructure.legs])
                        parent = momenta_dict.inv[children]
                        state = subtraction.SubtractionLeg.FINAL
                        if substructure.legs.has_initial_state_leg():
                            state = subtraction.SubtractionLeg.INITIAL
                        parent_leg = subtraction.SubtractionLeg(parent, 0, state)
                        half_2_legs.append(parent_leg)
                else:
                    half_2_substructures.append(substructure)
                    half_1_legs += substructure.legs
            half_1 = subtraction.SingularStructure(
                substructures=half_1_substructures, legs=half_1_legs )
            half_2 = subtraction.SingularStructure(
                substructures=half_2_substructures, legs=half_2_legs )
            if self.verbose:
                print "Half 1:", half_1
                print "Half 2:", half_2
            my_PS_point = MappingsTest.generate_PS(pars)
            squared_masses = dict()
            for leg in structure.legs:
                squared_masses['m2' + str(leg.n)] = my_PS_point[leg.n].square()
            for parent in pars['parents']:
                if pars['masses']:
                    p = sum(my_PS_point[child] for child in momenta_dict[parent])
                    squared_masses['m2' + str(parent)] = random.random()*p.square()
                else:
                    squared_masses['m2' + str(parent)] = 0.
            # Rotate it to avoid zero components
            for key in my_PS_point.keys():
                my_PS_point[key].rotoboost(MappingsTest.v1, MappingsTest.v2)
            if self.verbose > 1:
                print "Starting PS point:\n" + str(my_PS_point)
            # Perform the mapping in one go
            low_PS_point, low_vars = pars['mapping'].map_to_lower_multiplicity(
                my_PS_point, structure, momenta_dict, squared_masses,
                compute_jacobian=False )
            if self.verbose > 1:
                print "1 PS point:\n" + str(low_PS_point)
            squared_masses_1 = {key: val for key, val in squared_masses.items()}
            for leg in half_1_legs:
                squared_masses_1['m2' + str(leg.n)] = my_PS_point[leg.n].square()
            # Perform the mapping in two steps
            half_PS_point, half_vars = pars['mapping'].map_to_lower_multiplicity(
                my_PS_point, half_1, momenta_dict, squared_masses_1,
                compute_jacobian=False )
            if self.verbose > 1:
                print "1*0.5 PS point:\n" + str(half_PS_point)
            squared_masses_2 = {key: val for key, val in squared_masses.items()}
            for leg in half_2_legs:
                squared_masses_2['m2' + str(leg.n)] = half_PS_point[leg.n].square()
            full_PS_point, full_vars = pars['mapping'].map_to_lower_multiplicity(
                half_PS_point, half_2, momenta_dict, squared_masses_2,
                compute_jacobian=False )
            if self.verbose > 1:
                print "2*0.5 PS point:\n" + str(full_PS_point)
            assertDictAlmostEqual(self, low_PS_point, full_PS_point)

    def _test_commutative(self, pars):
        """Test whether a mapping is commutative."""

        # Make test deterministic by setting seed
        random.seed(self.seed)
        # Check many times
        for _ in range(pars.get('n_tests', self.n_tests)):
            # Generate a random setup
            MappingsTest.randomize(pars)
            structure = pars['structure']
            momenta_dict = pars['momenta_dict']
            if self.verbose:
                print "Structure:", structure
            # Randomly divide the mapping in two steps
            half_A_substructures = []
            half_B_substructures = []
            recoilers = [leg for leg in structure.legs]
            half_A_legs_1 = list(recoilers)
            half_A_legs_2 = list(recoilers)
            half_B_legs_1 = list(recoilers)
            half_B_legs_2 = list(recoilers)
            for substructure in structure.substructures:
                if bool(random.getrandbits(1)):
                    half_A_substructures.append(substructure)
                    half_B_legs_1 += substructure.legs
                    if substructure.name() == 'C':
                        children = frozenset([leg.n for leg in substructure.legs])
                        parent = momenta_dict.inv[children]
                        state = subtraction.SubtractionLeg.FINAL
                        if substructure.legs.has_initial_state_leg():
                            state = subtraction.SubtractionLeg.INITIAL
                        parent_leg = subtraction.SubtractionLeg(parent, 0, state)
                        half_B_legs_2.append(parent_leg)
                else:
                    half_B_substructures.append(substructure)
                    half_A_legs_1 += substructure.legs
                    if substructure.name() == 'C':
                        children = frozenset([leg.n for leg in substructure.legs])
                        parent = momenta_dict.inv[children]
                        state = subtraction.SubtractionLeg.FINAL
                        if substructure.legs.has_initial_state_leg():
                            state = subtraction.SubtractionLeg.INITIAL
                        parent_leg = subtraction.SubtractionLeg(parent, 0, state)
                        half_A_legs_2.append(parent_leg)
            half_A_1 = subtraction.SingularStructure(
                substructures=half_A_substructures, legs=half_A_legs_1 )
            half_A_2 = subtraction.SingularStructure(
                substructures=half_A_substructures, legs=half_A_legs_2 )
            half_B_1 = subtraction.SingularStructure(
                substructures=half_B_substructures, legs=half_B_legs_1 )
            half_B_2 = subtraction.SingularStructure(
                substructures=half_B_substructures, legs=half_B_legs_2 )
            if self.verbose:
                print "Path 1:"
                print half_A_1
                print half_B_2
                print "Path 2:"
                print half_B_1
                print half_A_2
            my_PS_point = MappingsTest.generate_PS(pars)
            squared_masses = dict()
            for leg in structure.legs:
                squared_masses['m2' + str(leg.n)] = my_PS_point[leg.n].square()
            for parent in pars['parents']:
                if pars['masses']:
                    p = sum(my_PS_point[child] for child in momenta_dict[parent])
                    squared_masses['m2' + str(parent)] = random.random()*p.square()
                else:
                    squared_masses['m2' + str(parent)] = 0.
            # Rotate it to avoid zero components
            for key in my_PS_point.keys():
                my_PS_point[key].rotoboost(MappingsTest.v1, MappingsTest.v2)
            if self.verbose > 1:
                print "Starting PS point:\n" + str(my_PS_point)
            # Perform the mapping A then B
            squared_masses_A_1 = {key: val for key, val in squared_masses.items()}
            for leg in half_A_legs_1:
                squared_masses_A_1['m2' + str(leg.n)] = my_PS_point[leg.n].square()
            PS_point_A_1, _ = pars['mapping'].map_to_lower_multiplicity(
                my_PS_point, half_A_1, momenta_dict, squared_masses_A_1,
                compute_jacobian=False )
            if self.verbose > 1:
                print "A 1 PS point:\n" + str(PS_point_A_1)
            squared_masses_B_2 = {key: val for key, val in squared_masses.items()}
            for leg in half_B_legs_2:
                squared_masses_B_2['m2' + str(leg.n)] = PS_point_A_1[leg.n].square()
            PS_point_B_2, _ = pars['mapping'].map_to_lower_multiplicity(
                PS_point_A_1, half_B_2, momenta_dict, squared_masses_B_2,
                compute_jacobian=False )
            if self.verbose > 1:
                print "B 2 PS point:\n" + str(PS_point_B_2)
            # Perform the mapping B then A
            squared_masses_B_1 = {key: val for key, val in squared_masses.items()}
            for leg in half_B_legs_1:
                squared_masses_B_1['m2' + str(leg.n)] = my_PS_point[leg.n].square()
            PS_point_B_1, _ = pars['mapping'].map_to_lower_multiplicity(
                my_PS_point, half_B_1, momenta_dict, squared_masses_B_1,
                compute_jacobian=False )
            if self.verbose > 1:
                print "B 1 PS point:\n" + str(PS_point_B_1)
            squared_masses_A_2 = {key: val for key, val in squared_masses.items()}
            for leg in half_A_legs_2:
                squared_masses_A_2['m2' + str(leg.n)] = PS_point_B_1[leg.n].square()
            PS_point_A_2, _ = pars['mapping'].map_to_lower_multiplicity(
                PS_point_B_1, half_A_2, momenta_dict, squared_masses_A_2,
                compute_jacobian=False )
            if self.verbose > 1:
                print "B 2 PS point:\n" + str(PS_point_B_2)
            assertDictAlmostEqual(self, PS_point_A_2, PS_point_B_2)

    # Test masses mappings
    #=====================================================================================

    def test_FinalZeroMassesMapping_invertible(self):
        """Test if FinalZeroMassesMapping is invertible."""

        pars = {
            'mapping': mappings.FinalZeroMassesMapping(),
            'min_coll_sets': 2, 'max_coll_sets': 5,
            'min_recoilers': 0, 'max_recoilers': 0,
            'max_unchanged': 0, 'masses': False,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 0, }
        self._test_invertible(pars)

    def test_FinalMassesMapping_invertible(self):
        """Test if FinalMassesMapping is invertible."""

        pars = {
            'mapping': mappings.FinalMassesMapping(),
            'min_coll_sets': 2, 'max_coll_sets': 5,
            'min_recoilers': 0, 'max_recoilers': 0,
            'max_unchanged': 0, 'masses': True,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 0, }
        self._test_invertible(pars)

    def test_FinalMasses_reduces_to_FinalZeroMasses(self):

        pars = {
            'mapping1': mappings.FinalMassesMapping(),
            'mapping2': mappings.FinalZeroMassesMapping(),
            'min_coll_sets': 2, 'max_coll_sets': 5,
            'min_recoilers': 0, 'max_recoilers': 0,
            'max_unchanged': 0, 'masses': False,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 0, }
        self._test_equal_lower(pars)

    # Test final-collinear mappings
    #=====================================================================================

    def test_FinalRescalingOneMapping_invertible(self):
        """Test if FinalRescalingOneMapping is invertible."""

        pars = {
            'mapping': mappings.FinalRescalingOneMapping(),
            'min_coll_sets': 1, 'max_coll_sets': 1,
            'min_recoilers': 1, 'max_recoilers': 3,
            'max_unchanged': 3, 'masses': False,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': False, }
        self._test_invertible(pars)

    def test_FinalLorentzOneMapping_invertible(self):
        """Test if FinalLorentzOneMapping is invertible."""

        pars = {
            'mapping': mappings.FinalLorentzOneMapping(),
            'min_coll_sets': 1, 'max_coll_sets': 1,
            'min_recoilers': 1, 'max_recoilers': 5,
            'max_unchanged': 3, 'masses': False,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4, }
        self._test_invertible(pars)

    def test_FinalLorentzOne_reduces_to_FinalRescalingOne(self):

        pars = {
            'mapping1': mappings.FinalRescalingOneMapping(),
            'mapping2': mappings.FinalLorentzOneMapping(),
            'min_coll_sets': 1, 'max_coll_sets': 1,
            'min_recoilers': 1, 'max_recoilers': 1,
            'max_unchanged': 0, 'masses': False,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 1,
            'supports_massive_recoilers': False }
        self._test_equal_lower(pars)

    def test_FinalGroupingMapping_invertible(self):
        """Test if FinalGroupingMapping is invertible."""

        pars = {
            'mapping': mappings.FinalGroupingMapping(),
            'max_unchanged': 3, 'masses': True,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 3, }
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
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 3, }
        pars.update({
            'min_coll_sets': 1, 'max_coll_sets': 3,
            'min_recoilers': 1, 'max_recoilers': 4,})
        self._test_invertible(pars)
        pars.update({
            'min_coll_sets': 2, 'max_coll_sets': 3,
            'min_recoilers': 0, 'max_recoilers': 3,})
        self._test_invertible(pars)

    def test_FinalLorentz_reduces_to_FinalLorentzOne(self):
        """Test that FinalLorentzMapping reduces to FinalLorentzOneMapping
        for a single collinear set.
        """

        pars = {
            'mapping1': mappings.FinalLorentzMapping(),
            'mapping2': mappings.FinalLorentzOneMapping(),
            'min_coll_sets': 1, 'max_coll_sets': 1,
            'min_recoilers': 1, 'max_recoilers': 4,
            'max_unchanged': 1, 'masses': False,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4, }
        self._test_equal_lower(pars)

    def test_FinalGrouping_reduces_to_FinalRescalingOne(self):
        """Test that FinalGroupingMapping reduces to FinalRescalingOneMapping
        for a single collinear set and massless recoilers.
        """

        pars = {
            'mapping1': mappings.FinalGroupingMapping(),
            'mapping2': mappings.FinalRescalingOneMapping(),
            'min_coll_sets': 1, 'max_coll_sets': 1,
            'min_recoilers': 1, 'max_recoilers': 4,
            'max_unchanged': 3, 'masses': False,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': False }
        self._test_equal_lower(pars)

    def test_FinalLorentz_equal_FinalGrouping(self):
        """Test that FinalLorentzMapping and FinalGroupingMapping
        are the same for zero or a single massive/massless recoiler.
        """

        pars = {
            'mapping1': mappings.FinalLorentzMapping(),
            'mapping2': mappings.FinalGroupingMapping(),
            'min_coll_sets': 1, 'max_coll_sets': 5,
            'min_recoilers': 1, 'max_recoilers': 1,
            'max_unchanged': 1, 'masses': False,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'n_tests': 5, }
        self._test_equal_lower(pars)
        pars['supports_massive_recoilers'] = False
        self._test_equal_lower(pars)
        pars.update({'min_coll_sets': 2, 'min_recoilers': 0, 'max_recoilers': 0})
        self._test_equal_lower(pars)

    def test_FinalGroupingMapping_associative(self):
        """Test if FinalGroupingMapping is associative."""

        pars = {
            'mapping': mappings.FinalGroupingMapping(),
            'max_unchanged': 3, 'masses': True,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 3, }
        pars.update({
            'min_coll_sets': 2, 'max_coll_sets': 4,
            'min_recoilers': 1, 'max_recoilers': 4,})
        self._test_associative(pars)
        pars.update({
            'min_coll_sets': 2, 'max_coll_sets': 5,
            'min_recoilers': 0, 'max_recoilers': 3,})
        self._test_associative(pars)

    def test_FinalGroupingMapping_commutative(self):
        """Test if FinalGroupingMapping is associative."""

        pars = {
            'mapping': mappings.FinalGroupingMapping(),
            'max_unchanged': 3, 'masses': True,
            'min_unresolved_per_set': 0, 'max_unresolved_per_set': 2, }
        pars.update({
            'min_coll_sets': 3, 'max_coll_sets': 4,
            'min_recoilers': 1, 'max_recoilers': 4,})
        self._test_commutative(pars)
        pars.update({
            'min_coll_sets': 3, 'max_coll_sets': 5,
            'min_recoilers': 0, 'max_recoilers': 3,})
        self._test_commutative(pars)

    # Test initial-collinear mappings
    #=====================================================================================

    def test_InitialLorentzOneMapping_invertible(self):
        """Test if InitialLorentzOneMapping is invertible."""

        pars = {
            'mapping': mappings.InitialLorentzOneMapping(),
            'initial_sets': 1,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'min_recoilers': 1, 'max_recoilers': 5,
            'max_unchanged': 3, 'masses': False, }
        self._test_invertible(pars)

    # Test soft mappings
    #=====================================================================================

    def test_FinalAssociativeSoftMappingZero_invertible(self):
        """Test if FinalAssociativeSoftMappingZero is invertible."""

        # WARNING
        # This mapping may produce a negative jacobian,
        # in which case there are two solutions for alpha in the inverse mapping
        # and the wrong one is picked.
        # This causes this test to fail.
        pars = {
            'mapping': mappings.FinalAssociativeSoftMappingZero(),
            'min_soft_sets': 1, 'max_soft_sets': 5,
            'min_recoilers': 2, 'max_recoilers': 5,
            'max_unchanged': 0,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': False, }
        self._test_invertible(pars)

    def test_FinalAssociativeSoftMappingZero_associative(self):
        """Test if FinalAssociativeSoftMappingZero is associative."""

        pars = {
            'mapping': mappings.FinalAssociativeSoftMappingZero(),
            'min_soft_sets': 1, 'max_soft_sets': 6,
            'min_recoilers': 2, 'max_recoilers': 4,
            'max_unchanged': 1,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 3,
            'supports_massive_recoilers': False, 'supports_massive_unresolved': False, }
        self._test_associative(pars)

    def test_FinalAssociativeSoftMappingZero_commutative(self):
        """Test if FinalAssociativeSoftMappingZero is commutative."""

        pars = {
            'mapping': mappings.FinalAssociativeSoftMappingZero(),
            'min_soft_sets': 1, 'max_soft_sets': 5,
            'min_recoilers': 2, 'max_recoilers': 4,
            'max_unchanged': 1,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 3,
            'supports_massive_recoilers': False, 'supports_massive_unresolved': False, }
        self._test_commutative(pars)

    def test_FinalAssociativeSoftMapping_invertible(self):
        """Test if FinalAssociativeSoftMapping is invertible."""

        # WARNING
        # This mapping may produce a negative jacobian,
        # in which case there are two solutions for alpha in the inverse mapping
        # and the wrong one is picked.
        # This causes this test to fail.
        pars = {
            'mapping': mappings.FinalAssociativeSoftMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 5,
            'min_recoilers': 2, 'max_recoilers': 5,
            'max_unchanged': 3,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4, }
        self._test_invertible(pars)

    def test_FinalAssociativeSoftMapping_associative(self):
        """Test if FinalAssociativeSoftMapping is associative."""

        pars = {
            'mapping': mappings.FinalAssociativeSoftMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 6,
            'min_recoilers': 2, 'max_recoilers': 4,
            'max_unchanged': 2,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 3, }
        self._test_associative(pars)

    def test_FinalAssociativeSoftMapping_commutative(self):
        """Test if FinalAssociativeSoftMapping is commutative."""

        pars = {
            'mapping': mappings.FinalAssociativeSoftMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 5,
            'min_recoilers': 2, 'max_recoilers': 4,
            'max_unchanged': 2,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 3, }
        self._test_commutative(pars)

    def test_SoftVsFinalPureRescalingMapping_invertible(self):
        """Test if SoftVsFinalPureRescalingMapping is invertible."""

        pars = {
            'mapping': mappings.SoftVsFinalPureRescalingMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 3,
            'min_recoilers': 2, 'max_recoilers': 5,
            'max_unchanged': 3,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': False, }
        self._test_invertible(pars)

    def test_SoftVsFinal_BoostThenRescale_invertible(self):
        """Test if SoftVsFinalBoostThenRescaleMapping is invertible."""

        pars = {
            'mapping': mappings.SoftVsFinalBoostThenRescaleMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 3,
            'min_recoilers': 2, 'max_recoilers': 5,
            'max_unchanged': 3, 'masses': None,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': False, 'jacobian': False, }
        self._test_invertible(pars)

    def test_SoftVsFinal_RescaleThenBoost_reduces_to_PureRescaling(self):

        pars = {
            'mapping1': mappings.SoftVsFinalPureRescalingMapping(),
            'mapping2': mappings.SoftVsFinalRescaleThenBoostMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 3,
            'min_recoilers': 2, 'max_recoilers': 5,
            'max_unchanged': 3, 'masses': False,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': False, 'jacobian': False, }
        self._test_equal_lower(pars)

    def test_SoftVsFinal_BoostThenRescale_reduces_to_PureRescaling(self):

        pars = {
            'mapping1': mappings.SoftVsFinalPureRescalingMapping(),
            'mapping2': mappings.SoftVsFinalBoostThenRescaleMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 3,
            'min_recoilers': 2, 'max_recoilers': 5,
            'max_unchanged': 3, 'masses': False,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': False, 'jacobian': False, }
        self._test_equal_lower(pars)

    def test_SoftVsFinal_BoostThenRescale_equal_RescaleThenBoost(self):

        pars = {
            'mapping1': mappings.SoftVsFinalRescaleThenBoostMapping(),
            'mapping2': mappings.SoftVsFinalBoostThenRescaleMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 3,
            'min_recoilers': 2, 'max_recoilers': 5,
            'max_unchanged': 3, 'masses': False,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4,
            'supports_massive_recoilers': True, 'jacobian': False, }
        self._test_equal_lower(pars)

    def test_SoftVsInitialMapping_invertible(self):
        """Test if SoftVsInitialMapping is invertible."""

        pars = {
            'mapping': mappings.SoftVsInitialMapping(),
            'min_soft_sets': 1, 'max_soft_sets': 5,
            'min_recoilers': 2, 'max_recoilers': 5,
            'max_unchanged': 3,
            'min_unresolved_per_set': 1, 'max_unresolved_per_set': 4, }
        self._test_invertible(pars)

    # TODO Confirm that indeed the SoftVsInitialMapping is neither associative or commutative
    # def test_SoftVsInitialMapping_associative(self):
    #     """Test if SoftVsInitialMapping is associative."""
    #
    #     pars = {
    #         'mapping': mappings.SoftVsInitialMapping(),
    #         'min_soft_sets': 1, 'max_soft_sets': 6,
    #         'min_recoilers': 2, 'max_recoilers': 4,
    #         'max_unchanged': 2,
    #         'min_unresolved_per_set': 1, 'max_unresolved_per_set': 3, }
    #     self._test_associative(pars)
    #
    # def test_SoftVsInitialMapping_commutative(self):
    #     """Test if SoftVsInitialMapping is commutative."""
    #
    #     pars = {
    #         'mapping': mappings.SoftVsInitialMapping(),
    #         'min_soft_sets': 1, 'max_soft_sets': 5,
    #         'min_recoilers': 2, 'max_recoilers': 4,
    #         'max_unchanged': 2,
    #         'min_unresolved_per_set': 1, 'max_unresolved_per_set': 3, }
    #     self._test_commutative(pars)