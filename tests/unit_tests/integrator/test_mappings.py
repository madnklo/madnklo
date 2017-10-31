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
# Final-final collinear mappings
#=========================================================================================

class CollinearVariablesTest(unittest.TestCase):
    """Test class for variables describing internal collinear structure."""

    my_mapping = mappings.FinalFinalCollinearMapping()
    n_children = 3

    def test_collinear_variables_away_from_limit(self):
        """Test determination of collinear variables and reverse mapping,
        for completely generic input values.
        """

        # Generate n_children+1 random massless vectors
        # (the light-cone direction is also random)
        my_PS_point = LorentzVectorDict()
        for i in range(self.n_children+1):
            my_PS_point[i] = LorentzVector(
                [0., ] + [random.random() for _ in range(3)] )
            my_PS_point[i].set_square(0)
        # Compute collinear variables
        variables = dict()
        self.my_mapping.get_collinear_variables(
            my_PS_point, self.n_children, range(self.n_children),
            variables )
        # Compute total momentum
        total_momentum = LorentzVector()
        for i in range(self.n_children):
            total_momentum += my_PS_point[i]
        # Compute new phase space point
        new_PS_point = LorentzVectorDict()
        new_PS_point[self.n_children] = my_PS_point[self.n_children]
        self.my_mapping.set_collinear_variables(
            new_PS_point, self.n_children, range(self.n_children),
            total_momentum, variables )
        # Check the two phase-space points are equal
        for i in range(self.n_children):
            self.assertEqual(my_PS_point[i], new_PS_point[i])

    def test_collinear_variables_close_to_limit(self):
        """Test determination of collinear variables and reverse mapping
        in a typical collinear situation.
        """

        # Generate n_children random starting directions
        directions = {
            i: Vector([random.random() for _ in range(3)])
            for i in range(self.n_children)
        }
        coll_direction = Vector(3*[0., ])
        for n in directions.values():
            coll_direction += n
        # Generate values for the parameter that describes approach to limit
        pars = (math.pow(0.25, i) for i in range(6))
        # This is one way like any other to approach the limit
        for par in pars:
            new_directions = {
                i: (par*n + (1-par)*coll_direction)
                for (i, n) in directions.items()
            }
            my_PS_point = {
                i: LorentzVector([0., ] + list(n)).set_square(0)
                for (i, n) in new_directions.items()
            }
            my_PS_point[self.n_children] = LorentzVector(
                [0., ] + list(coll_direction)
            ).set_square(0)
            # Compute collinear variables
            variables = dict()
            self.my_mapping.get_collinear_variables(
                my_PS_point, self.n_children, range(self.n_children),
                variables
            )
            # Compute total momentum
            total_momentum = LorentzVector(4*[0., ])
            for i in range(self.n_children):
                total_momentum += my_PS_point[i]
            # Compute new phase space point
            new_PS_point = LorentzVectorDict()
            new_PS_point[self.n_children] = my_PS_point[self.n_children]
            self.my_mapping.set_collinear_variables(
                new_PS_point, self.n_children, range(self.n_children),
                total_momentum, variables
            )
            # Check the two phase-space points are equal
            for i in range(self.n_children):
                self.assertEqual(my_PS_point[i], new_PS_point[i])

class FFRescalingMappingOneTest(unittest.TestCase):
    """Test class for FFRescalingMappingOneTest."""

    mapping = mappings.FFLorentzMappingOne()
    n_collinear = (2, )
    n_recoilers = 3
    massive = False
    structure = subtraction.SingularStructure()
    momenta_dict = subtraction.bidict()

    def setUp(self):

        n_tot = sum(self.n_collinear) + self.n_recoilers
        for i in range(n_tot):
            self.momenta_dict[i] = frozenset((i,))
        self.structure = []
        n_collinear_so_far = 0
        n_bunch = 0
        for n_in_this_bunch in self.n_collinear:
            this_bunch_numbers = tuple(range(
                    n_collinear_so_far, n_collinear_so_far + n_in_this_bunch
                ))
            this_bunch_legs = subtraction.SubtractionLegSet(
                subtraction.SubtractionLeg(
                    i, 21, FINAL
                ) for i in this_bunch_numbers
            )
            self.structure += [subtraction.CollStructure(legs=this_bunch_legs), ]
            if n_in_this_bunch != 1:
                self.momenta_dict[n_tot + n_bunch] = frozenset(this_bunch_numbers)
            n_collinear_so_far += n_in_this_bunch
            n_bunch += 1
        self.structure += [
                subtraction.SubtractionLeg(
                    i, 21, FINAL
                )
                for i in range(n_collinear_so_far, n_tot)
            ]
        self.structure = subtraction.SingularStructure(*self.structure)

    def test_collinear_map_invertible(self):
        """Test mapping and inverse."""

        # Generate n_collinear random massive vectors plus
        # n_recoilers (massive = True, False) random vectors
        my_PS_point = LorentzVectorDict()
        n_collinear_so_far = 0
        for n_in_this_bunch in self.n_collinear:
            for i in range(
                n_collinear_so_far, n_collinear_so_far + n_in_this_bunch
            ):
                my_PS_point[i] = LorentzVector(
                    [0., ] + [random.random() for _ in range(3)]
                )
                my_PS_point[i].set_square(0)
            if n_in_this_bunch == 1:
                my_PS_point[n_collinear_so_far].set_square(random.random())
            n_collinear_so_far += n_in_this_bunch
        for i in range(
            n_collinear_so_far, n_collinear_so_far + self.n_recoilers
        ):
            my_PS_point[i] = LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            )
            if self.massive:
                my_PS_point[i].set_square(random.random())
            else:
                my_PS_point[i].set_square(0)
        # I know what I'm doing
        old_PS_point = copy.deepcopy(my_PS_point)
        # Compute collinear variables
        variables = dict()
        print my_PS_point
        self.mapping.map_to_lower_multiplicity(
            my_PS_point, self.structure, self.momenta_dict, variables
        )
        print my_PS_point
        self.mapping.map_to_higher_multiplicity(
            my_PS_point, self.structure, self.momenta_dict, variables
        )
        print my_PS_point
        for i in my_PS_point.keys():
            self.assertEqual(my_PS_point[i], old_PS_point[i])

#=========================================================================================
# Soft mappings
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

        # Generate n_soft random massless vectors
        my_PS_point = LorentzVectorDict()
        for i in range(self.n_soft):
            my_PS_point[i] = LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            )
            my_PS_point[i].set_square(0)
        # Compute soft variables
        variables = dict()
        self.my_mapping.get_soft_variables(
            my_PS_point, range(self.n_soft), variables
        )
        # Compute new phase space point
        new_PS_point = LorentzVectorDict()
        self.my_mapping.set_soft_variables(
            new_PS_point, range(self.n_soft), variables
        )
        # Check the two phase-space points are equal
        for i in range(self.n_soft):
            self.assertEqual(my_PS_point[i], new_PS_point[i])


    def test_soft_variables_close_to_limit(self):
        """Test determination of soft variables and reverse mapping
        in a typical soft situation.
        """

        # Generate n_soft random massless vectors
        my_PS_point = LorentzVectorDict()
        for i in range(self.n_soft):
            my_PS_point[i] = LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            )
            my_PS_point[i].set_square(0)
        # Generate values for the parameter that describes approach to limit
        pars = (math.pow(0.25, i) for i in range(8))
        # This is one way like any other to approach the limit
        for par in pars:
            old_PS_point = {
                i: par*my_PS_point[i]
                for i in my_PS_point.keys()
            }
            # Compute collinear variables
            variables = dict()
            self.my_mapping.get_soft_variables(
                old_PS_point, range(self.n_soft), variables
            )
            # Compute new phase space point
            new_PS_point = LorentzVectorDict()
            self.my_mapping.set_soft_variables(
                new_PS_point, range(self.n_soft), variables
            )
            for i in range(self.n_soft):
                self.assertEqual(old_PS_point[i], new_PS_point[i])

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
                my_PS_point[i] = LorentzVector(
                    [0., ] + [random.random() for _ in range(3)]
                )
                my_PS_point[i].set_square(0)
            n_soft_so_far += n_in_this_bunch
        for i in range(
            n_soft_so_far, n_soft_so_far + self.n_recoilers
        ):
            my_PS_point[i] = LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            )
            if self.massive:
                my_PS_point[i].set_square(random.random())
            else:
                my_PS_point[i].set_square(0)
        # I know what I'm doing
        old_PS_point = copy.deepcopy(my_PS_point)
        # Compute collinear variables
        variables = dict()
        self.mapping.map_to_lower_multiplicity(
            my_PS_point, self.structure, self.momenta_dict, variables
        )
        self.mapping.map_to_higher_multiplicity(
            my_PS_point, self.structure, self.momenta_dict, variables
        )
        for i in my_PS_point.keys():
            self.assertEqual(my_PS_point[i], old_PS_point[i])

#=========================================================================================
# Test the phase-space walkers
#=========================================================================================

# Define processes for phase-space walker tests
#=========================================================================================

# Subtraction instance

my_subtraction = subtraction.IRSubtraction(
    simple_qcd.model, coupling_types=('QCD'), n_unresolved=1 )

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

def walk_invertible(test, walker, process, max_unresolved=None):
    """Test walk and inverse."""

    my_operators = my_subtraction.get_all_elementary_operators(process)
    my_combinations = my_subtraction.get_all_combinations(my_operators, max_unresolved)
    my_counterterms = [
        my_subtraction.get_counterterm(combination, process)
        for combination in my_combinations ]

    # For each counterterm
    for i in range(len(my_counterterms)):
        print "Testing counterterm", my_counterterms[i]
        # Generate random vectors
        my_PS_point = LorentzVectorDict()
        legs_FS = (
            subtraction.SubtractionLeg(leg)
            for leg in process['legs']
            if leg['state'] == FINAL
        )
        leg_numbers = (leg.n for leg in legs_FS)
        for j in leg_numbers:
            my_PS_point[j] = LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            ).set_square(0)
        # Compute collinear variables
        res_dict1 = walker.walk_to_lower_multiplicity(
            my_PS_point, my_counterterms[i], True
        )
        (currs1, ME1, jac1, kin)  = (
            res_dict1['currents'], res_dict1['matrix_element'],
            res_dict1['jacobian'], res_dict1['kinematic_variables'] )

        res_dict2 = walker.walk_to_higher_multiplicity(
            ME1[1], my_counterterms[i], kin
        )
        (currs2, ME2, jac2)  = (
            res_dict2['currents'], res_dict2['matrix_element'],
            res_dict2['jacobian']
        )

        # Check currents
        test.assertEqual(len(currs1), len(currs2))
        for i_curr in range(len(currs1)):
            test.assertEqual(currs1[i_curr][0], currs2[i_curr][0])
            test.assertEqual(
                currs1[i_curr][1].keys(), currs2[i_curr][1].keys()
            )
            for i_part in currs1[i_curr][1].keys():
                test.assertEqual(
                    currs1[i_curr][1][i_part],
                    currs2[i_curr][1][i_part]
                )
        # Check MEs
        test.assertEqual(ME1[0], ME2[0])
        test.assertEqual(ME1[1].keys(), ME2[1].keys())
        for i_part in ME1[1].keys():
            test.assertEqual(ME1[1][i_part], ME2[1][i_part])
        # Check jacobians
        test.assertAlmostEqual(jac1*jac2, 1.)

class FlatCollinearWalkerTest(unittest.TestCase):
    """Test class for FlatCollinearWalker."""

    walker = mappings.FlatCollinearWalker()

    def test_FlatCollinearWalker_invertible(self):

        walk_invertible(self, self.walker, H_to_uuxddxH, 2)

class SimpleNLOWalkerTest(unittest.TestCase):
    """Test class for FlatCollinearWalker."""

    walker = mappings.SimpleNLOWalker()

    def test_SimpleNLOWalker_invertible(self):

        # Want to make the test a bit more serious,
        # so use combinations of NLO elementary operators up to 2 unresolved particles

        walk_invertible(self, self.walker, H_to_uuxddxH, 2)
        walk_invertible(self, self.walker, H_to_qqxggH, 2)

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
            PS_point[i] = LorentzVector([0., ] + [random.random() for _ in range(3)])
            PS_point[i].set_square(0)
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
