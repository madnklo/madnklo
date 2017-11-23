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
# Generation of random momenta
#=========================================================================================

def random_momentum(massive=False):
    """Generate a random momentum, either massive or massless."""

    foo = LorentzVector([0., ] + [(2 * random.random() - 1) for _ in range(3)])
    if massive:
        foo.set_square(random.random())
    else:
        foo.set_square(0)
    return foo

#=========================================================================================
# Test collinear variables
#=========================================================================================

class CollinearVariablesTest(unittest.TestCase):
    """Test class for variables describing internal collinear structure."""

    n_children = 3
    massive = True

    def test_collinear_variables_away_from_limit(self):
        """Test determination of collinear variables and reverse mapping,
        for completely generic input values.
        """

        # Generate n_children random massless vectors
        my_PS_point = LorentzVectorDict()
        children = tuple(range(self.n_children))
        for i in children:
            my_PS_point[i] = random_momentum(self.massive)
        # Generate two random light-cone directions
        na = random_momentum(False)
        nb = random_momentum(False)
        # Compute collinear variables
        variables = dict()
        mappings.CollinearVariables.get(
            my_PS_point, children, na, nb, variables )
        # Compute total momentum
        total_momentum = LorentzVector()
        for i in children:
            total_momentum += my_PS_point[i]
        # Compute new phase space point
        new_PS_point = LorentzVectorDict()
        mappings.CollinearVariables.set(
            new_PS_point, children, total_momentum, na, nb, variables )
        # Check the two phase-space points are equal
        self.assertDictEqual(my_PS_point, new_PS_point)

    def test_collinear_variables_close_to_limit(self):
        """Test determination of collinear variables and reverse mapping
        in a typical collinear situation.
        """

        children = tuple(range(self.n_children))
        # Generate n_children random starting directions
        directions = {
            i: Vector([random.random() for _ in range(3)])
            for i in children }
        coll_direction = Vector(3*[0., ])
        for n in directions.values():
            coll_direction += n
        # Generate values for the parameter that describes approach to limit
        pars = (math.pow(0.25, i) for i in range(14))
        # This is one way like any other to approach the limit
        for par in pars:
            new_directions = {
                i: (par*n + (1-par)*coll_direction)
                for (i, n) in directions.items() }
            my_PS_point = {
                i: LorentzVector([0., ] + list(n)).set_square(0)
                for (i, n) in new_directions.items() }
            # na = LorentzVector([1, 0, 0, +1])
            # nb = LorentzVector([1, 0, 0, -1])
            na = LorentzVector([0., ] + list( coll_direction)).set_square(0)
            nb = LorentzVector([0., ] + list(-coll_direction)).set_square(0)
            # Compute collinear variables
            variables = dict()
            mappings.CollinearVariables.get(my_PS_point, children, na, nb, variables)
            # Compute total momentum
            total_momentum = LorentzVector()
            for i in children:
                total_momentum += my_PS_point[i]
            # Compute new phase space point
            new_PS_point = LorentzVectorDict()
            mappings.CollinearVariables.set(
                new_PS_point, children, total_momentum, na, nb, variables )
            # Check the two phase-space points are equal
            self.assertDictEqual(my_PS_point, new_PS_point)

#=========================================================================================
# Test collinear mappings
#=========================================================================================

class FFCollinearMappingTest(object):
    """Test class for final-final collinear mappings."""

    @staticmethod
    def randomize(pars):

        # Randomly generate the structure of the final state
        n_collinear_sets = random.randint(1, pars['max_collinear_sets'])
        real_min_recoilers = pars['min_recoilers']
        if n_collinear_sets == 1 and pars['min_recoilers'] == 0:
            real_min_recoilers = 1
        n_recoilers = random.randint(real_min_recoilers, pars['max_recoilers'])
        n_unchanged = random.randint(0, pars['max_unchanged'])
        n_collinears_per_set = [
            random.randint(2, pars['max_collinears_per_set'])
            for _ in range(n_collinear_sets) ]
        # Compute the total number of particles in the final state
        n_tot = sum(n_collinears_per_set) + n_recoilers + n_unchanged
        # Initialise trivial entries {i: i} of the momenta dictionary
        pars['momenta_dict'] = subtraction.bidict()
        for i in range(1, n_tot+3):
            pars['momenta_dict'][i] = frozenset((i, ))
        # Randomly generate nonconsecutive numbers of the final state
        shuffled_numbers = range(2, n_tot+2)
        random.shuffle(shuffled_numbers)
        pars['structure'] = []
        n_set = 1
        for n_in_this_set in n_collinears_per_set:
            numbers_in_this_set = shuffled_numbers[:n_in_this_set]
            shuffled_numbers = shuffled_numbers[n_in_this_set:]
            legs_in_this_set = subtraction.SubtractionLegSet(
                subtraction.SubtractionLeg(i, 21, FINAL)
                for i in numbers_in_this_set)
            pars['structure'] += [subtraction.CollStructure(legs=legs_in_this_set), ]
            if n_in_this_set != 1:
                pars['momenta_dict'][n_tot + n_set + 1] = frozenset(numbers_in_this_set)
            n_set += 1
        pars['structure'] += [
            subtraction.SubtractionLeg(i, 21, FINAL)
            for i in shuffled_numbers[:n_recoilers]]
        pars['structure'] = subtraction.SingularStructure(*pars['structure'])
        pars['unchanged'] = shuffled_numbers[n_recoilers:]
        assert len(pars['unchanged']) == n_unchanged
        # print "The random structure for this test is", str(pars['structure'])
        # print "with particles number", pars['unchanged'], "unchanged"

    @staticmethod
    def test_invertible(pars, test):
        """Test mapping and inverse."""

        # Generate a random phase space point
        my_PS_point = LorentzVectorDict()
        for i in pars['unchanged']:
            my_PS_point[i] = random_momentum(True)
        for leg in pars['structure'].legs:
            my_PS_point[leg.n] = random_momentum(
                pars['supports_massive_recoilers'] )
        for set in pars['structure'].substructures:
            leg_ns = [leg.n for leg in set.get_all_legs()]
            # If only one particle is going collinear,
            # interpret its momentum as a collective one and make it massive
            if len(leg_ns) == 1:
                my_PS_point[leg_ns[0]] = random_momentum(True)
            else:
                for i in leg_ns:
                    my_PS_point[i] = random_momentum(False)
        # I know what I'm doing
        old_PS_point = copy.deepcopy(my_PS_point)
        # Compute collinear variables
        variables = dict()
        pars['mapping'].map_to_lower_multiplicity(
            my_PS_point, pars['structure'], pars['momenta_dict'], variables )
        pars['mapping'].map_to_higher_multiplicity(
            my_PS_point, pars['structure'], pars['momenta_dict'], variables )
        test.assertDictEqual(my_PS_point, old_PS_point)

class FFRescalingMappingOneTest(unittest.TestCase):
    """Test class for FFRescalingMappingOne."""

    # Test settings
    pars = {
        'mapping': mappings.FFRescalingMappingOne(),
        'max_collinear_sets': 1,
        'max_collinears_per_set': 4,
        'min_recoilers': 1,
        'max_recoilers': 3,
        'max_unchanged': 3,
        'supports_massive_parents': True,
        'supports_massive_recoilers': False }

    def setUp(self):

        FFCollinearMappingTest.randomize(self.pars)

    def test_FFRescalingMappingOne_invertible(self):

        FFCollinearMappingTest.test_invertible(self.pars, self)

class FFLorentzMappingOneTest(unittest.TestCase):
    """Test class for FFRescalingMappingOne."""

    # Test settings
    pars = {
        'mapping': mappings.FFLorentzMappingOne(),
        'max_collinear_sets': 1,
        'max_collinears_per_set': 4,
        'min_recoilers': 1,
        'max_recoilers': 5,
        'max_unchanged': 3,
        'supports_massive_parents': True,
        'supports_massive_recoilers': True }

    def setUp(self):

        FFCollinearMappingTest.randomize(self.pars)

    def test_FFLorentzMappingOne_invertible(self):

        FFCollinearMappingTest.test_invertible(self.pars, self)

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

        softs = tuple(range(self.n_soft))
        # Generate n_soft random massless vectors
        my_PS_point = LorentzVectorDict()
        for i in softs:
            my_PS_point[i] = random_momentum(False)
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

        softs = tuple(range(self.n_soft))
        # Generate n_soft random massless vectors
        my_PS_point = LorentzVectorDict()
        for i in softs:
            my_PS_point[i] = random_momentum(False)
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
# Soft mappings
#=========================================================================================

class SomogyietalSoftTest(unittest.TestCase):
    """Test class for MappingSomogyietalSoft."""

    # This specific mapping only for massless,
    # but 'massive' will be useful in the future so might as well keep it
    mapping = mappings.MappingSomogyietalSoft()
    n_soft = (2, 3, )
    n_recoilers = 3
    massive = False
    structure = subtraction.SingularStructure()
    momenta_dict = subtraction.bidict()

    def setUp(self):

        n_tot = sum(self.n_soft) + self.n_recoilers
        for i in range(n_tot):
            self.momenta_dict[i] = frozenset((i,))
        self.structure = []
        n_soft_so_far = 0
        n_bunch = 0
        for n_in_this_bunch in self.n_soft:
            this_bunch_numbers = tuple(range(
                    n_soft_so_far, n_soft_so_far + n_in_this_bunch
                ))
            this_bunch_legs = subtraction.SubtractionLegSet(
                subtraction.SubtractionLeg(
                    i, 21, FINAL
                ) for i in this_bunch_numbers
            )
            self.structure += [subtraction.SoftStructure(legs=this_bunch_legs), ]
            n_soft_so_far += n_in_this_bunch
            n_bunch += 1
        self.structure += [
                subtraction.SubtractionLeg(
                    i, 21, FINAL
                )
                for i in range(n_soft_so_far, n_tot)
            ]
        self.structure = subtraction.SingularStructure(*self.structure)

    def test_soft_map_invertible(self):
        """Test mapping and inverse."""

        # Generate n_soft random massless vectors plus
        # n_recoilers (massive = True, False) random vectors
        my_PS_point = LorentzVectorDict()
        n_soft_so_far = 0
        for n_in_this_bunch in self.n_soft:
            for i in range(
                n_soft_so_far, n_soft_so_far + n_in_this_bunch
            ):
                my_PS_point[i] = random_momentum(False)
                my_PS_point[i].set_square(0)
            n_soft_so_far += n_in_this_bunch
        for i in range(
            n_soft_so_far, n_soft_so_far + self.n_recoilers
        ):
            my_PS_point[i] = random_momentum(self.massive)
        # I know what I'm doing
        old_PS_point = copy.deepcopy(my_PS_point)
        # Compute collinear variables
        variables = dict()
        self.mapping.map_to_lower_multiplicity(
            my_PS_point, self.structure, self.momenta_dict, variables )
        self.mapping.map_to_higher_multiplicity(
            my_PS_point, self.structure, self.momenta_dict, variables )
        for i in my_PS_point.keys():
            self.assertEqual(my_PS_point[i], old_PS_point[i])

#=========================================================================================
# Test the phase-space walkers
#=========================================================================================

# Define processes for phase-space walker tests
#=========================================================================================

# Subtraction instance

my_subtraction = subtraction.IRSubtraction(
    simple_qcd.model, coupling_types=('QCD', ), n_unresolved=1 )

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

def walk_invertible(test, walker, process, max_unresolved=None, verbose=False):
    """Test walk and inverse."""

    my_operators = my_subtraction.get_all_elementary_operators(process)
    my_combinations = my_subtraction.get_all_combinations(my_operators, max_unresolved)
    my_counterterms = [
        my_subtraction.get_counterterm(combination, process)
        for combination in my_combinations ]

    # For each counterterm
    for i in range(len(my_counterterms)):
        if verbose:
            print "Considering counterterm", my_counterterms[i]
        # Generate random vectors
        my_PS_point = LorentzVectorDict()
        legs_FS = (
            subtraction.SubtractionLeg(leg)
            for leg in process['legs']
            if leg['state'] == FINAL )
        for leg in legs_FS:
            my_PS_point[leg.n] = random_momentum(
                process.get('model').get_particle(leg.pdg)['mass'].lower() != 'zero' )
        # Compute collinear variables
        res_dict1 = walker.walk_to_lower_multiplicity(
            my_PS_point, my_counterterms[i], True )
        (currs1, ME1, jac1, kin)  = (
            res_dict1['currents'], res_dict1['matrix_element'],
            res_dict1['jacobian'], res_dict1['kinematic_variables'] )
        if verbose:
            print "Walking down"
            for curr in currs1:
                print curr[1]
            print ME1[1]
        res_dict2 = walker.walk_to_higher_multiplicity(
            ME1[1], my_counterterms[i], kin )
        (currs2, ME2, jac2)  = (
            res_dict2['currents'], res_dict2['matrix_element'],
            res_dict2['jacobian'] )
        if verbose:
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
        # Check jacobians
        test.assertAlmostEqual(jac1*jac2, 1.)

def approach_limit(test, walker, process, verbose=False):
    """Test walk and inverse."""

    elementary_operators = my_subtraction.get_all_elementary_operators(process)
    elementary_structures = [
        subtraction.SingularOperatorList([op, ]).simplify()
        for op in elementary_operators ]
    elementary_counterterms = [
        my_subtraction.get_counterterm(combination, process)
        for combination in elementary_structures ]

    # For each counterterm
    for ct in elementary_counterterms:
        if verbose:
            print "\n", "*"*100
            print "Considering counterterm", ct
            print "*" * 100, "\n"
        # Generate random vectors
        my_PS_point = LorentzVectorDict()
        legs_FS = tuple(
            subtraction.SubtractionLeg(leg)
            for leg in process['legs']
            if leg['state'] == FINAL )
        for leg in legs_FS:
            my_PS_point[leg.n] = random_momentum(
                process.get('model').get_particle(leg.pdg)['mass'].lower() != 'zero' )
        if verbose:
            print "Starting phase space point:\n", my_PS_point, "\n"
        squares = {key: my_PS_point[key].square() for key in my_PS_point.keys()}
        # Compute collinear variables
        res_dict = walker.walk_to_lower_multiplicity(my_PS_point, ct, True)
        base_PS_point = res_dict['matrix_element'][1]
        base_kinematic_variables = res_dict['kinematic_variables']
        if verbose:
            print "Baseline PS point:"
            print base_PS_point, "\n"
            print "Starting kinematic variables:"
            print base_kinematic_variables, "\n"
        for alpha in [math.pow(0.1, i) for i in range(5)]:
            res_dict = walker.approach_limit(
                base_PS_point, ct, base_kinematic_variables, alpha)
            new_PS_point = res_dict['currents'][0][1]
            if verbose:
                print "New PS point for", alpha, ":\n", new_PS_point
            for leg in legs_FS:
                if process.get('model').get_particle(leg.pdg)['mass'].lower() == 'zero':
                    test.assertLess(
                        abs(new_PS_point[leg.n].square()),
                        math.sqrt(new_PS_point[leg.n].eps()) )
                else:
                    test.assertAlmostEqual(new_PS_point[leg.n].square(), squares[leg.n])

class FlatCollinearWalkerTest(unittest.TestCase):
    """Test class for FlatCollinearWalker."""

    walker = mappings.FlatCollinearWalker()

    def test_FlatCollinearWalker_invertible(self):

        walk_invertible(self, self.walker, H_to_uuxddxH, 2)

class FFNLOWalkerTest(unittest.TestCase):
    """Test class for FlatCollinearWalker."""

    walker = mappings.FFNLOWalker()

    def test_FFNLOWalker_invertible(self):

        walk_invertible(self, self.walker, H_to_uuxddxH, 1)
        walk_invertible(self, self.walker, H_to_qqxggH, 1)

    def test_FFNLOWalker_approach_limit(self):

        approach_limit(self, self.walker, H_to_qqxggH, True)

    def test_sc_approach_limit(self):

        # Set up a soft-collinear counterterm
        sc_ll = base_objects.LegList([
            base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
            base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
            base_objects.Leg({'number': 5, 'id': -1, 'state': FINAL}), ])
        sc_rp = base_objects.Process({
            'legs': sc_ll, 'model': simple_qcd.model, 'n_loops': 0 })
        sc_ss = subtraction.CollStructure(
            legs=subtraction.SubtractionLegSet((
                subtraction.SubtractionLeg(3, -1, FINAL), )),
            substructures=[
                subtraction.SoftStructure(
                    subtraction.SubtractionLeg(4, 21, FINAL)) ] )
        sc_md = subtraction.bidict({i: frozenset((i, )) for i in range(1, 5)})
        sc_md[5] = frozenset((3, 4, ))
        sc_ct = subtraction.Counterterm(
            process=sc_rp,
            nodes=[
                subtraction.CountertermNode(
                    current=subtraction.Current({
                        'singular_structure': sc_ss }) ) ],
            prefactor=1,
            momenta_dict=sc_md )

        # Start from a random phase space point
        # The Higgs is going to have a random mass, but it doesn't matter
        PS_point = LorentzVectorDict()
        PS_point[1] = LorentzVector()
        for i in range(2, 5):
            PS_point[i] = random_momentum(False)
            PS_point[1] += PS_point[i]

        hike_down = self.walker.walk_to_lower_multiplicity(PS_point, sc_ct, True)
        lower_PS_point = hike_down['resulting_PS_point']
        starting_variables = hike_down['kinematic_variables']

        ratios = []
        flucts = []
        for par in range(10):
            x = math.pow(0.1, par)
            hike_up = self.walker.approach_limit(
                lower_PS_point, sc_ct, starting_variables, x )
            ratios_vec = LorentzVector([
                hike_up['resulting_PS_point'][3][i] / hike_up['resulting_PS_point'][4][i]
                for i in range(4) ])
            ratios.append(ratios_vec.view(type=Vector).square())
            norm_ratios = ratios_vec.view(type=Vector)
            norm_ratios.normalize()
            flucts.append(abs(norm_ratios - Vector(4 * [0.5, ])))

        # Skipping the first few values, not close enough to the limit
        # Numerical effects can make the test fail for the deep IR region
        for i in range(2, len(ratios)-2):
            self.assertLess(
                abs(ratios[i + 2] / ratios[i + 1] - 10.),
                abs(ratios[i + 1] / ratios[i] - 10.) )
        for i in range(2, len(flucts)-1):
            self.assertLess(flucts[i+1], flucts[i])
