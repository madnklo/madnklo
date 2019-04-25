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

"""Unit test library for the routines of the core library related to writing
color information for diagrams.
"""

import random

import madgraph.core.base_objects as base_objects
import madgraph.core.subtraction as sub

import tests.unit_tests as unittest
import tests.input_files.simple_qcd as simple_qcd

#=========================================================================================
# Shorthands for initial and final state
#=========================================================================================

INITIAL = base_objects.Leg.INITIAL
FINAL   = base_objects.Leg.FINAL

assert sub.SubtractionLeg.INITIAL == INITIAL
assert sub.SubtractionLeg.FINAL   == FINAL

#=========================================================================================
# Test ordered tuples
#=========================================================================================

class SubsetTest(unittest.TestCase):
    """Test class for the function is_subset."""

    def test_subset(self):
        """Test the function is_subset."""

        self.assertTrue(sub.is_subset(
            (3,),
            (3,)
        ))
        self.assertTrue(sub.is_subset(
            (3, 5, 91),
            (3, 3, 5, 7, 91, 123)
        ))
        self.assertFalse(sub.is_subset(
            (3, 3),
            (3, 5, 91)
        ))
        self.assertFalse(sub.is_subset(
            (3, 7),
            (1, 5, 42)
        ))
        self.assertTrue(sub.is_subset(
            "aabccddfkk",
            "aabbbbccddeeffffhhhhkkkk"
        ))
        self.assertFalse(sub.is_subset(
            "aabccddfkk",
            "aabbbbcddeeffffhhhhkkkk"
        ))

    def test_disjoint(self):
        """Test the function is_disjoint."""

        # Batch test
        for _ in range(100):
            na = random.randint(0, 20)
            nb = random.randint(0, 20)
            a = sorted([random.randint(1, 100) for _ in range(na)])
            b = sorted([random.randint(1, 100) for _ in range(nb)])
            self.assertEqual(sub.is_disjoint(a, b), set(a).isdisjoint(set(b)))

        # Pathological cases
        a = ()
        b = (1, 2)
        self.assertTrue(sub.is_disjoint(a, b))
        self.assertTrue(sub.is_disjoint(a, a))

    def test_union(self):
        """Test the function union."""

        a = "abel"
        b = "bio"
        c = "klmmmnnnn"
        self.assertEqual(
            sub.union(a),
            tuple(iter(a))
        )
        self.assertEqual(
            sub.union(a, b),
            tuple(iter("abbeilo"))
        )
        self.assertEqual(
            sub.union(a, b, b, c),
            tuple(iter("abbbeiikllmmmnnnnoo"))
        )

    def test_difference(self):
        """Test the function difference."""

        # Test without multiplicities
        for _ in range(100):
            na = random.randint(0, 50)
            nb = random.randint(0, 30)
            aset = set([random.randint(1, 100) for _ in range(na)])
            bset = set([random.randint(1, 100) for _ in range(nb)])
            a = sorted(aset)
            b = sorted(bset)
            self.assertEqual(
                set(sub.difference(a, b)),
                aset.difference(bset)
            )

        # Test with multiplicities
        a = "abbello"
        b = "bello"
        self.assertEqual(
            sub.difference(a, b),
            tuple(iter("ab"))
        )

#=========================================================================================
# Test SubtractionLegSet
#=========================================================================================

SLeg = sub.SubtractionLeg
SLegSet = sub.SubtractionLegSet

class SubtractionLegSetTest(unittest.TestCase):
    """Test class for SubtractionLegSet."""

    def test_split_by_pdg_abs(self):
        """Test the method split_by_pdg_abs."""

        subtraction_leg_set = SLegSet((
            SLeg(0, +1, INITIAL),
            SLeg(1, -1, INITIAL),
            SLeg(2, 25, FINAL),
            SLeg(3, +5, FINAL),
            SLeg(4, +1, FINAL),
            SLeg(5, +5, FINAL),
            SLeg(6, -1, FINAL),
            SLeg(7, +2, FINAL),
            SLeg(8, -5, FINAL),
        ))
        legs_by_pdg, pdg_counts, rest = subtraction_leg_set.split_by_pdg_abs([1,2,3,4,5], FINAL)
        target_legs_by_pdg = {
            1: [SLegSet([subtraction_leg_set[4], ]), SLegSet([subtraction_leg_set[6], ])],
            2: [SLegSet(), SLegSet([subtraction_leg_set[7], ])],
            3: [SLegSet(), SLegSet()],
            4: [SLegSet(), SLegSet()],
            5: [SLegSet([subtraction_leg_set[8], ]),
                SLegSet([subtraction_leg_set[3], subtraction_leg_set[5], ])]
        }
        target_pdg_counts = {(1, 1): [1], (0, 1): [2], (0, 0): [3, 4], (1, 2): [5]}
        target_rest = SLegSet(subtraction_leg_set[0:3])
        self.assertDictEqual(target_legs_by_pdg, legs_by_pdg)
        self.assertDictEqual(target_pdg_counts, pdg_counts)
        self.assertEqual(target_rest, rest)

    def test_map_leg_numbers(self):
        """Test the method map_leg_numbers."""

        # Basic test
        set_1a = SLegSet([SLeg(2, 21, FINAL), SLeg(6, -1, FINAL), SLeg(7, +2, FINAL), ])
        set_1b = SLegSet([SLeg(1, 21, FINAL), SLeg(5, -1, FINAL), SLeg(6, +2, FINAL), ])
        map_1b = set_1a.map_leg_numbers(set_1b)
        target_map_1 = {2: 1, 6: 5, 7: 6}
        self.assertDictEqual(target_map_1, map_1b)
        set_1c = SLegSet([SLeg(1, 21, FINAL), SLeg(5, -1, FINAL), ])
        map_1c = set_1a.map_leg_numbers(set_1c)
        self.assertIsNone(map_1c)
        set_1d = SLegSet([SLeg(1, 21, FINAL), SLeg(5, -1, FINAL), SLeg(6, -2, FINAL), ])
        map_1d = set_1a.map_leg_numbers(set_1d)
        self.assertIsNone(map_1d)
        map_1d = set_1a.map_leg_numbers(set_1d, [[2, ], ])
        self.assertDictEqual(target_map_1, map_1d)
        set_1e = SLegSet([SLeg(1, 21, FINAL), SLeg(5, -1, FINAL), SLeg(6, +3, FINAL), ])
        map_1e = set_1a.map_leg_numbers(set_1e)
        self.assertIsNone(map_1e)
        map_1e = set_1a.map_leg_numbers(set_1e, [[2, 3, ], ])
        self.assertDictEqual(target_map_1, map_1e)
        map_1e = set_1a.map_leg_numbers(set_1e, [[1, 2, 3, ], ])
        self.assertDictEqual(target_map_1, map_1e)

        # Example test: d dx u
        set_2a = SLegSet([SLeg(2, +1, FINAL), SLeg(6, -1, FINAL), SLeg(7, +2, FINAL), ])
        set_2b = SLegSet([SLeg(2, +4, FINAL), SLeg(6, -2, FINAL), SLeg(7, +2, FINAL), ])
        map_2b = set_2a.map_leg_numbers(set_2b, [range(1, 6), ])
        target_maps_2 = [
            {2: 7, 6: 6, 7: 2},
            {2: 7, 6: 2, 7: 6},
        ]
        self.assertTrue(map_2b in target_maps_2)
        set_2c = SLegSet(set_2b + (SLeg(8, +2, FINAL), ))
        map_2c = set_2a.map_leg_numbers(set_2c, [range(1, 6), ])
        self.assertIsNone(map_2c)

        # Example test: d dx s sx
        set_3a = SLegSet([
            SLeg(12, +1, FINAL), SLeg(11, -1, FINAL),
            SLeg(17, -3, FINAL), SLeg(19, +3, FINAL),
        ])
        set_3b = SLegSet([
            SLeg(21, +4, FINAL), SLeg(24, +3, FINAL),
            SLeg(22, -4, FINAL), SLeg(23, -3, FINAL),
        ])
        map_3b = set_3a.map_leg_numbers(set_3b, [range(1, 6), ])
        target_maps_3 = [
            {12: 21, 11: 22, 17: 24, 19: 23},
            {12: 21, 11: 22, 17: 23, 19: 24},
            {12: 22, 11: 21, 17: 24, 19: 23},
            {12: 22, 11: 21, 17: 23, 19: 24},
            {12: 23, 11: 24, 17: 21, 19: 22},
            {12: 23, 11: 24, 17: 22, 19: 21},
            {12: 24, 11: 23, 17: 21, 19: 22},
            {12: 24, 11: 23, 17: 22, 19: 21},
        ]
        self.assertTrue(map_3b in target_maps_3)
        set_3c = SLegSet([
            SLeg(21, +4, FINAL), SLeg(24, +4, FINAL),
            SLeg(22, -4, FINAL), SLeg(23, -4, FINAL),
        ])
        map_3c = set_3a.map_leg_numbers(set_3c, [range(1, 6), ])
        self.assertIsNone(map_3c)
        set_3d = SLegSet([SLeg(21, +4, FINAL), SLeg(24, +3, FINAL), ])
        map_3d = set_3a.map_leg_numbers(set_3d, [range(1, 6), ])
        self.assertIsNone(map_3d)

        # Example test: INITIAL d FINAL d
        set_4a = SLegSet([SLeg(1, +1, INITIAL), SLeg(2, +1, FINAL), ])
        set_4b = SLegSet([SLeg(4, +4, INITIAL), SLeg(1, +4, FINAL), ])
        map_4b = set_4a.map_leg_numbers(set_4b, [range(1, 6), ])
        target_map_4 = {1: 4, 2: 1}
        self.assertDictEqual(map_4b, target_map_4)
        set_4c = SLegSet([SLeg(4, +4, FINAL), SLeg(1, +4, FINAL), ])
        map_4c = set_4a.map_leg_numbers(set_4c, [range(1, 6), ])
        self.assertIsNone(map_4c)

#=========================================================================================
# Test SingularStructure
#=========================================================================================

class SingularStructureTest(unittest.TestCase):
    """Test class for SingularStructure."""

    def test_singular_structure_init(self):
        """Test initialization of SingularStructure."""

        subtraction_leg_set = sub.SubtractionLegSet((
            sub.SubtractionLeg(1,  1, INITIAL),
            sub.SubtractionLeg(3, 25, FINAL)
        ))
        self.assertEqual(
            sub.SingularStructure(legs=subtraction_leg_set),
            sub.SingularStructure(*subtraction_leg_set)
        )

    def test_act_on(self):
        """Test action of operators."""

        leg1 = base_objects.Leg({'number': 1, 'id':  1, 'state': INITIAL})
        leg4 = base_objects.Leg({'number': 4, 'id': 21, 'state': FINAL})
        leg5 = base_objects.Leg({'number': 5, 'id': 21, 'state': FINAL})

        C14  = sub.CollStructure(leg1, leg4)
        C45  = sub.CollStructure(leg4, leg5)
        C145 = sub.CollStructure(leg1, leg4, leg5)
        S4   = sub.SoftStructure(leg4)
        S45  = sub.SoftStructure(leg4, leg5)

        CC_list = sub.SingularStructure(C14, C145)
        CC_simple = CC_list.nest()
        CC_benchmark = sub.SingularStructure(
            sub.CollStructure(
                sub.CollStructure(leg1, leg4),
                leg5 ) )
        self.assertEqual(CC_simple,CC_benchmark)

        SC_list = sub.SingularStructure(S4, C145)
        SC_simple = SC_list.nest()
        self.assertEqual(str(SC_simple), '(C(S(4),1,5),)')

        SS_list = sub.SingularStructure(S4, S45)
        SS_simple = SS_list.nest()
        self.assertTrue(SS_simple.is_void)

        C45S45_list = sub.SingularStructure(C45, S45)
        C45S45_simple = C45S45_list.nest()
        C45S45_benchmark = sub.SingularStructure(
            sub.SoftStructure(
                sub.CollStructure(leg4, leg5) ) )
        self.assertEqual(C45S45_simple,C45S45_benchmark)

        S45C45_list = sub.SingularStructure(S45, C45)
        S45C45_simple = S45C45_list.nest()
        self.assertTrue(S45C45_simple.is_void)

    def test_decompose(self):
        """Test decomposition of singular structures."""

        a_structure = sub.SingularStructure(sub.CollStructure(
            sub.CollStructure(
                sub.SubtractionLeg(1,  1, INITIAL),
                sub.SubtractionLeg(4, 21, FINAL)
            ),
            sub.SoftStructure(
                sub.SubtractionLeg(14,  2, FINAL),
                sub.SubtractionLeg(11, -2, FINAL),
            ),
            sub.SubtractionLeg(5, 21, FINAL)
        ))
        sub_1 = sub.CollStructure(
            sub.SubtractionLeg(1,  1, INITIAL),
            sub.SubtractionLeg(4, 21, FINAL) )
        sub_2 = sub.SoftStructure(
            sub.SubtractionLeg(14, 2, FINAL),
            sub.SubtractionLeg(11, -2, FINAL) )
        sub_3 = sub.CollStructure(
            sub.SubtractionLeg(1,   1, INITIAL),
            sub.SubtractionLeg(4,  21, FINAL),
            sub.SubtractionLeg(5,  21, FINAL),
            sub.SubtractionLeg(14,  2, FINAL),
            sub.SubtractionLeg(11, -2, FINAL) )
        dec = a_structure.decompose()
        ben = [sub_1, sub_2, sub_3]
        self.assertEqual(dec, ben)

    def test_count_unresolved(self):
        """Test counting of unresolved particles."""

        leg1 = base_objects.Leg({'number': 1, 'id':  1, 'state': INITIAL})
        leg4 = base_objects.Leg({'number': 4, 'id': 21, 'state': FINAL})
        leg5 = base_objects.Leg({'number': 5, 'id': 21, 'state': FINAL})

        C14  = sub.CollStructure(leg1, leg4)
        C145 = sub.CollStructure(leg1, leg4, leg5)
        S4   = sub.SoftStructure(leg4)
        S5   = sub.SoftStructure(leg5)
        S45  = sub.SoftStructure(leg4, leg5)

        list1 = sub.SingularStructure(C145).nest()
        self.assertEqual(list1.count_unresolved(), 2)

        list2 = sub.SingularStructure(S45).nest()
        self.assertEqual(list2.count_unresolved(), 2)

        list3 = sub.SingularStructure(S45, C145).nest()
        self.assertEqual(list3.count_unresolved(), 2)

        list4 = sub.SingularStructure(S4, S5).nest()
        self.assertEqual(list4.count_unresolved(), 2)

        list5 = sub.SingularStructure(S4, S5, C145).nest()
        self.assertEqual(list5.count_unresolved(), 2)

        list6 = sub.SingularStructure(C14).nest()
        self.assertEqual(list6.count_unresolved(), 1)

        list7 = sub.SingularStructure(S4, C14).nest()
        self.assertEqual(list7.count_unresolved(), 1)

    def test_from_string(self):
        """Test reconstruction of singular structures from a string."""

        mylegs = base_objects.LegList([
            base_objects.Leg({'number': 1, 'id':  1, 'state': INITIAL}),
            base_objects.Leg({'number': 2, 'id': -1, 'state': INITIAL}),
            base_objects.Leg({'number': 3, 'id':  1, 'state': FINAL}),
            base_objects.Leg({'number': 4, 'id': -1, 'state': FINAL}),
            base_objects.Leg({'number': 5, 'id': 21, 'state': FINAL}),
            base_objects.Leg({'number': 6, 'id': 21, 'state': FINAL}),
        ])

        myprocess = base_objects.Process({
            'legs': mylegs, 'model': simple_qcd.model })

        string0 = "C(3,6)"
        ss0 = sub.CollStructure(legs=[mylegs[2], mylegs[5]])
        ss0_reco = sub.SingularStructure.from_string(string0, myprocess)
        self.assertEqual(ss0, ss0_reco)

        string1 = "C(S(5),3,6)"
        ss1 = sub.CollStructure(
            substructures=[sub.SoftStructure(legs=(mylegs[4],))],
            legs=[mylegs[2], mylegs[5]] )
        ss1_reco = sub.SingularStructure.from_string(string1, myprocess)
        self.assertEqual(ss1, ss1_reco)

        string2 = "C(1,3,6)"
        ss2_reco = sub.SingularStructure.from_string(string2, myprocess)
        self.assertEqual(string2, str(ss2_reco))

        string3 = "(C(3,6),C(4,5),)"
        ss3_reco = sub.SingularStructure.from_string(string3, myprocess)
        self.assertEqual(string3, str(ss3_reco))
        string3b = "(C(3,6),C(4,5))"
        ss3b_reco = sub.SingularStructure.from_string(string3b, myprocess)
        self.assertEqual(str(ss3b_reco), str(ss3_reco))

        string4 = "(C((3,6),C(4,5),)"
        ss4_reco = sub.SingularStructure.from_string(string4, myprocess)
        self.assertIsNone(ss4_reco)
        string5 = "C(1,x,6)"
        ss5_reco = sub.SingularStructure.from_string(string5, myprocess)
        self.assertIsNone(ss5_reco)
        string6 = "(XY(1,6),S(5))"
        ss6_reco = sub.SingularStructure.from_string(string6, myprocess)
        self.assertIsNone(ss6_reco)
        string7 = "C(3,6,S(5)),)"
        ss7_reco = sub.SingularStructure.from_string(string7, myprocess)
        self.assertIsNone(ss7_reco)

    def test_map_leg_number(self):
        """Test matching of singular structures."""

        mylegs = SLegSet([
            SLeg(0,  1, INITIAL),
            SLeg(1, -1, INITIAL),
            SLeg(2,  2, FINAL),
            SLeg(3, -2, FINAL),
            SLeg(4, 21, FINAL),
            SLeg(5, 21, FINAL),
            SLeg(6,  2, FINAL),
            SLeg(7,  3, FINAL),
            SLeg(8, -3, FINAL),
            SLeg(9,  1, FINAL),
        ])

        # Match C(d1,g2)
        ss1 = sub.CollStructure(legs=[SLeg(1, 1, FINAL), SLeg(2, 21, FINAL)])
        ss1a = sub.CollStructure(legs=[mylegs[2], mylegs[5]])
        map1a = ss1.map_leg_numbers(ss1a, equivalent_pdg_sets=[range(1,6)])
        self.assertDictEqual({1: 2, 2: 5}, map1a)
        ss1b = sub.SoftStructure(legs=[mylegs[2], mylegs[5]])
        map1b = ss1.map_leg_numbers(ss1b, equivalent_pdg_sets=[range(1,6)])
        self.assertIsNone(map1b)
        ss1c = sub.CollStructure(legs=[mylegs[4], mylegs[3]])
        map1c = ss1.map_leg_numbers(ss1c, equivalent_pdg_sets=[range(1,6)])
        self.assertDictEqual({1: 3, 2: 4}, map1c)
        ss1d = sub.SoftStructure(legs=[mylegs[4], mylegs[5]])
        map1d = ss1.map_leg_numbers(ss1d, equivalent_pdg_sets=[range(1,6)])
        self.assertIsNone(map1d)

        # Match S(C(u1, u~2))
        ss2 = sub.SoftStructure(substructures=[
            sub.CollStructure(legs=[SLeg(1, +1, FINAL), SLeg(2, -1, FINAL)])
        ])
        ss2a = sub.SoftStructure(substructures=[
            sub.CollStructure(legs=[mylegs[3], mylegs[6]])
        ])
        map2a = ss2.map_leg_numbers(ss2a, equivalent_pdg_sets=[range(1,6)])
        maps2a = [{1: 3, 2: 6}, {1: 6, 2: 3}]
        self.assertTrue(map2a in maps2a)
        ss2b = sub.SoftStructure(
            substructures=[sub.CollStructure(legs=[mylegs[3], mylegs[6]])],
            legs=[mylegs[4]]
        )
        map2b = ss2.map_leg_numbers(ss2b, equivalent_pdg_sets=[range(1,6)])
        self.assertIsNone(map2b)
        ss2c = sub.CollStructure(substructures=[
            sub.CollStructure(legs=[mylegs[3], mylegs[6]])
        ])
        map2c = ss2.map_leg_numbers(ss2c, equivalent_pdg_sets=[range(1,6)])
        self.assertIsNone(map2c)
        ss2d = sub.SoftStructure(substructures=[
            sub.CollStructure(legs=[mylegs[7], mylegs[8]])
        ])
        map2d = ss2.map_leg_numbers(ss2d, equivalent_pdg_sets=[range(1,6)])
        maps2d = [{1: 7, 2: 8}, {1: 8, 2: 7}]
        self.assertTrue(map2d in maps2d)
        ss2e = sub.SoftStructure(substructures=[
            sub.CollStructure(legs=[mylegs[2], mylegs[6]])
        ])
        map2e = ss2.map_leg_numbers(ss2e, equivalent_pdg_sets=[range(1,6)])
        self.assertIsNone(map2e)

        # Match (C(u1, u~2, g3), C(d4, d~5))
        ss3 = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[SLeg(1, +1, FINAL), SLeg(2, -1, FINAL), SLeg(3, 21, FINAL)]),
            sub.CollStructure(legs=[SLeg(4, +2, FINAL), SLeg(5, -2, FINAL)]),
        ])
        ss3a = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[mylegs[8], mylegs[7]]),
            sub.CollStructure(legs=[mylegs[4], mylegs[3], mylegs[6]]),
        ])
        map3a = ss3.map_leg_numbers(ss3a, equivalent_pdg_sets=[range(1,6)])
        maps3a = [
            {1: 3, 2: 6, 3: 4, 4: 7, 5: 8},
            {1: 3, 2: 6, 3: 4, 4: 8, 5: 7},
            {1: 6, 2: 3, 3: 4, 4: 7, 5: 8},
            {1: 6, 2: 3, 3: 4, 4: 8, 5: 7},
        ]
        self.assertTrue(map3a in maps3a)
        ss3b = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[mylegs[8], mylegs[5], mylegs[7]]),
            sub.CollStructure(legs=[mylegs[2], mylegs[3]]),
        ])
        map3b = ss3.map_leg_numbers(ss3b, equivalent_pdg_sets=[range(1,6)])
        maps3b = [
            {1: 7, 2: 8, 3: 5, 4: 2, 5: 3},
            {1: 7, 2: 8, 3: 5, 4: 3, 5: 2},
            {1: 8, 2: 7, 3: 5, 4: 2, 5: 3},
            {1: 8, 2: 7, 3: 5, 4: 3, 5: 2},
        ]
        self.assertTrue(map3b in maps3b)
        ss3c = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[mylegs[8], mylegs[5], mylegs[7]]),
            sub.SoftStructure(legs=[mylegs[2], mylegs[3]]),
        ])
        map3c = ss3.map_leg_numbers(ss3c, equivalent_pdg_sets=[range(1,6)])
        self.assertIsNone(map3c)
        ss3d = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[mylegs[8], mylegs[5], mylegs[7]]),
            sub.CollStructure(legs=[mylegs[2], mylegs[4], mylegs[3]]),
        ])
        map3d = ss3.map_leg_numbers(ss3d, equivalent_pdg_sets=[range(1,6)])
        self.assertIsNone(map3d)
        ss3e = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[mylegs[8], mylegs[5], mylegs[7]]),
            sub.CollStructure(legs=[mylegs[2], mylegs[3]]),
            sub.CollStructure(legs=[mylegs[4], mylegs[6]]),
        ])
        map3e = ss3.map_leg_numbers(ss3e, equivalent_pdg_sets=[range(1,6)])
        self.assertIsNone(map3e)

#=========================================================================================
# Test Counterterm (and CountertermNode)
#=========================================================================================

class CountertermTest(unittest.TestCase):

    def random_children(
        self, parent_number, momenta_dict, prob_parent, max_children):
        """Branch a leaf in a momenta dictionary."""

        n_children = random.randint(2, max_children)
        children = set()
        for i_child in range(n_children):
            max_n_guess = len(momenta_dict.keys()) + 1
            while True:
                child_number = random.randint(1, max_n_guess)
                if child_number in momenta_dict.keys():
                    max_n_guess += 1
                    continue
                children.add(child_number)
                momenta_dict[child_number] = frozenset([child_number, ])
                break
        momenta_dict[parent_number] = frozenset(children)
        for child_number in children:
            if random.random() < prob_parent/len(children):
                self.random_children(
                    child_number, momenta_dict, prob_parent**2, max_children)

    def test_get_ancestor(self):
        """Test the reconstruction of an ancestor within a momentum dictionary."""

        random.seed(42)
        for i in range(100):
            # Generate starting tree
            n_start = random.randint(1, 4)
            momenta_dict = sub.bidict({
                i: frozenset([i, ]) for i in range(1, n_start+1) })
            for key in momenta_dict.keys():
                self.random_children(key, momenta_dict, 0.2, 4)
            # Branch one leaf
            leaves_before = [key for (key, val) in momenta_dict.items() if len(val) == 1]
            ancestor = leaves_before[random.randrange(len(leaves_before))]
            self.random_children(ancestor, momenta_dict, 0.8, 4)
            leaves_after  = [key for (key, val) in momenta_dict.items() if len(val) == 1]
            ancestor_leaves = set(leaves_after).difference(set(leaves_before))
            found_ancestor = sub.Counterterm.get_ancestor(frozenset(ancestor_leaves), momenta_dict)
            self.assertEqual(ancestor, found_ancestor)

    def test_split_loops_flat(self):
        """Test the assignment of loop numbers within a flat counterterm."""

        # Set up a flat counterterm
        # Momentum dictionary: particles from 1 to 3
        md = sub.bidict({i: frozenset((i, )) for i in range(1, 5)})
        # SingularStructures:  S(4), S(5)
        s4 = sub.SoftStructure(sub.SubtractionLeg(4, 21, FINAL))
        s5 = sub.SoftStructure(sub.SubtractionLeg(5, 21, FINAL))
        # Reduced process: 1 > 2 3 (H > d d~)
        ll = base_objects.LegList([
            base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
            base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
            base_objects.Leg({'number': 3, 'id': -1, 'state': FINAL}), ])
        rp = sub.Current({'legs': ll, 'model': simple_qcd.model})
        node4 = sub.CountertermNode(current=sub.Current({'singular_structure': s4}) )
        node5 = sub.CountertermNode(current=sub.Current({'singular_structure': s5}) )
        ct = sub.Counterterm(
            process=rp, nodes=[node4, node5],
            prefactor=1, momenta_dict=md )

        # Sprinkle 2 loops over this counterterm
        ct_all = ct.distribute_loops(2)

        # Build the benchmark by hand
        # Node #1
        ct_0 = ct.get_copy('n_loops')
        ct_0.current.set('n_loops', 0)
        ct_1 = ct.get_copy('n_loops')
        ct_1.current.set('n_loops', 1)
        ct_2 = ct.get_copy('n_loops')
        ct_2.current.set('n_loops', 2)
        # Eliminate determined after #1
        ct_2.nodes[0].current.set('n_loops', 0)
        ct_2.nodes[1].current.set('n_loops', 0)
        # Node #2
        ct_00 = ct_0.get_copy('n_loops')
        ct_00.nodes[0].current.set('n_loops', 0)
        ct_01 = ct_0.get_copy('n_loops')
        ct_01.nodes[0].current.set('n_loops', 1)
        ct_02 = ct_0.get_copy('n_loops')
        ct_02.nodes[0].current.set('n_loops', 2)
        ct_10 = ct_1.get_copy('n_loops')
        ct_10.nodes[0].current.set('n_loops', 0)
        ct_11 = ct_1.get_copy('n_loops')
        ct_11.nodes[0].current.set('n_loops', 1)
        # All are determined after #2
        ct_00.nodes[1].current.set('n_loops', 2)
        ct_01.nodes[1].current.set('n_loops', 1)
        ct_02.nodes[1].current.set('n_loops', 0)
        ct_10.nodes[1].current.set('n_loops', 1)
        ct_11.nodes[1].current.set('n_loops', 0)
        # Build list of all loop combinations
        ct_bench = (ct_2, ct_00, ct_01, ct_02, ct_10, ct_11)

        # Convert to strings, discard order
        ct_bench_str = frozenset(str(cti) for cti in ct_bench)
        ct_all_str = frozenset(str(cti) for cti in ct_all)

        # Check
        self.assertEqual(ct_all_str, ct_bench_str)

    def test_split_loops_nested(self):
        """Test the assignment of loop numbers within a nested counterterm."""

        # Set up a nested counterterm
        # Momentum dictionary: particles from 1 to 5
        md = sub.bidict({i: frozenset((i, )) for i in range(1, 6)})
        # Momentum dictionary: 6 is parent of (3, 4) and 7 of (5, 6)
        md[6]  = frozenset((3, 4, ))
        md[7]  = frozenset((5, 6, ))
        # SingularStructures:  C(3, 4), C(5, 6)
        c34 = sub.CollStructure(
            sub.SubtractionLeg(3, -1, FINAL),
            sub.SubtractionLeg(4, 21, FINAL) )
        c56 = sub.CollStructure(
            sub.SubtractionLeg(5, 21, FINAL),
            sub.SubtractionLeg(6, -1, FINAL) )
        # Reduced process: 1 > 2 10 (H > d d~)
        ll = base_objects.LegList([
            base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
            base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
            base_objects.Leg({'number': 7, 'id': -1, 'state': FINAL}), ])
        rp = sub.Current({'legs': ll, 'model': simple_qcd.model})
        node56 = sub.CountertermNode(
            current=sub.Current({'singular_structure': c56}),
            nodes=[
                sub.CountertermNode(
                    current=sub.Current({'singular_structure': c34}) ), ] )
        ct = sub.Counterterm(
            process=rp, nodes=[node56, ],
            prefactor=1, momenta_dict=md )

        # Sprinkle 2 loops over this counterterm
        ct_all = ct.distribute_loops(2)

        # Build the benchmark by hand
        # Node #1
        ct_0 = ct.get_copy('n_loops')
        ct_0.current.set('n_loops', 0)
        ct_1 = ct.get_copy('n_loops')
        ct_1.current.set('n_loops', 1)
        ct_2 = ct.get_copy('n_loops')
        ct_2.current.set('n_loops', 2)
        # Eliminate determined after #1
        ct_2.nodes[0].current.set('n_loops', 0)
        ct_2.nodes[0].nodes[0].current.set('n_loops', 0)
        # Node #2
        ct_00 = ct_0.get_copy('n_loops')
        ct_00.nodes[0].current.set('n_loops', 0)
        ct_01 = ct_0.get_copy('n_loops')
        ct_01.nodes[0].current.set('n_loops', 1)
        ct_02 = ct_0.get_copy('n_loops')
        ct_02.nodes[0].current.set('n_loops', 2)
        ct_10 = ct_1.get_copy('n_loops')
        ct_10.nodes[0].current.set('n_loops', 0)
        ct_11 = ct_1.get_copy('n_loops')
        ct_11.nodes[0].current.set('n_loops', 1)
        # All are determined after #2
        ct_00.nodes[0].nodes[0].current.set('n_loops', 2)
        ct_01.nodes[0].nodes[0].current.set('n_loops', 1)
        ct_02.nodes[0].nodes[0].current.set('n_loops', 0)
        ct_10.nodes[0].nodes[0].current.set('n_loops', 1)
        ct_11.nodes[0].nodes[0].current.set('n_loops', 0)
        # Build list of all loop combinations
        ct_bench = (ct_2, ct_00, ct_01, ct_02, ct_10, ct_11)

        # Convert to strings, discard order
        ct_bench_str = frozenset(str(cti) for cti in ct_bench)
        ct_all_str = frozenset(str(cti) for cti in ct_all)

        # Check
        self.assertEqual(ct_all_str, ct_bench_str)

    def test_split_loops(self):
        """Test the assignment of loop numbers within a current."""

        # Set up a counterterm that is complicated enough
        # Momentum dictionary: particles from 1 to 7
        md = sub.bidict({i: frozenset((i, )) for i in range(1, 8)})
        # Momentum dictionary: 8 is parent of (4, 5), 9 of (6, 8), 10 of (7, 9)
        md[8]  = frozenset((4, 5, ))
        md[9]  = frozenset((6, 8, ))
        md[10] = frozenset((7, 9, ))
        # SingularStructures:  S(3), C(4, 5), C(6, 8), C(7, 9),
        s3  = sub.SoftStructure(sub.SubtractionLeg(3, 21, FINAL))
        c45 = sub.CollStructure(
            sub.SubtractionLeg(5, 21, FINAL),
            sub.SubtractionLeg(4, -1, FINAL) )
        c68 = sub.CollStructure(
            sub.SubtractionLeg(6, 21, FINAL),
            sub.SubtractionLeg(8, -1, FINAL) )
        c79 = sub.CollStructure(
            sub.SubtractionLeg(7, 21, FINAL),
            sub.SubtractionLeg(9, -1, FINAL) )
        # Reduced process: 1 > 2 10 (H > d d~)
        ll = base_objects.LegList([
            base_objects.Leg({'number': 1, 'id': 25, 'state': INITIAL}),
            base_objects.Leg({'number': 2, 'id':  1, 'state': FINAL}),
            base_objects.Leg({'number': 10, 'id': -1, 'state': FINAL}), ])
        rp = sub.Current({
            'legs': ll, 'model': simple_qcd.model})
        node79 = sub.CountertermNode(
            current=sub.Current({'singular_structure': c79}),
            nodes=[
                sub.CountertermNode(
                    current=sub.Current({'singular_structure': c68}),
                    nodes=[
                        sub.CountertermNode(
                            current=sub.Current({'singular_structure': c45}) ) ] ) ] )
        node3 = sub.CountertermNode(current=sub.Current({'singular_structure': s3}) )
        ct = sub.Counterterm(
            process=rp, nodes=[node3, node79],
            prefactor=1, momenta_dict=md )

        # Sprinkle 2 loops over this counterterm
        ct_all = ct.distribute_loops(2)

        # Painfully build the benchmark by hand
        # Node #1
        ct_0 = ct.get_copy('n_loops')
        ct_0.current.set('n_loops', 0)
        ct_1 = ct.get_copy('n_loops')
        ct_1.current.set('n_loops', 1)
        ct_2 = ct.get_copy('n_loops')
        ct_2.current.set('n_loops', 2)
        # Eliminate determined after #1
        ct_2.nodes[0].current.set('n_loops', 0)
        ct_2.nodes[1].current.set('n_loops', 0)
        ct_2.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_2.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        # Node #2
        ct_00 = ct_0.get_copy('n_loops')
        ct_00.nodes[0].current.set('n_loops', 0)
        ct_01 = ct_0.get_copy('n_loops')
        ct_01.nodes[0].current.set('n_loops', 1)
        ct_02 = ct_0.get_copy('n_loops')
        ct_02.nodes[0].current.set('n_loops', 2)
        ct_10 = ct_1.get_copy('n_loops')
        ct_10.nodes[0].current.set('n_loops', 0)
        ct_11 = ct_1.get_copy('n_loops')
        ct_11.nodes[0].current.set('n_loops', 1)
        # Eliminate determined after #2
        ct_02.nodes[1].current.set('n_loops', 0)
        ct_02.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_02.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        ct_11.nodes[1].current.set('n_loops', 0)
        ct_11.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_11.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        # Node #3
        ct_000 = ct_00.get_copy('n_loops')
        ct_000.nodes[1].current.set('n_loops', 0)
        ct_001 = ct_00.get_copy('n_loops')
        ct_001.nodes[1].current.set('n_loops', 1)
        ct_002 = ct_00.get_copy('n_loops')
        ct_002.nodes[1].current.set('n_loops', 2)
        ct_010 = ct_01.get_copy('n_loops')
        ct_010.nodes[1].current.set('n_loops', 0)
        ct_011 = ct_01.get_copy('n_loops')
        ct_011.nodes[1].current.set('n_loops', 1)
        ct_100 = ct_10.get_copy('n_loops')
        ct_100.nodes[1].current.set('n_loops', 0)
        ct_101 = ct_10.get_copy('n_loops')
        ct_101.nodes[1].current.set('n_loops', 1)
        # Eliminate determined after #3
        ct_002.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_002.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        ct_011.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_011.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        ct_101.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_101.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        # Node #4
        ct_0000 = ct_000.get_copy('n_loops')
        ct_0000.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_0001 = ct_000.get_copy('n_loops')
        ct_0001.nodes[1].nodes[0].current.set('n_loops', 1)
        ct_0002 = ct_000.get_copy('n_loops')
        ct_0002.nodes[1].nodes[0].current.set('n_loops', 2)
        ct_0010 = ct_001.get_copy('n_loops')
        ct_0010.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_0011 = ct_001.get_copy('n_loops')
        ct_0011.nodes[1].nodes[0].current.set('n_loops', 1)
        ct_0100 = ct_010.get_copy('n_loops')
        ct_0100.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_0101 = ct_010.get_copy('n_loops')
        ct_0101.nodes[1].nodes[0].current.set('n_loops', 1)
        ct_1000 = ct_100.get_copy('n_loops')
        ct_1000.nodes[1].nodes[0].current.set('n_loops', 0)
        ct_1001 = ct_100.get_copy('n_loops')
        ct_1001.nodes[1].nodes[0].current.set('n_loops', 1)
        # All are determined after #4
        ct_0000.nodes[1].nodes[0].nodes[0].current.set('n_loops', 2)
        ct_0001.nodes[1].nodes[0].nodes[0].current.set('n_loops', 1)
        ct_0002.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        ct_0010.nodes[1].nodes[0].nodes[0].current.set('n_loops', 1)
        ct_0011.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        ct_0100.nodes[1].nodes[0].nodes[0].current.set('n_loops', 1)
        ct_0101.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        ct_1000.nodes[1].nodes[0].nodes[0].current.set('n_loops', 1)
        ct_1001.nodes[1].nodes[0].nodes[0].current.set('n_loops', 0)
        # Build list of all loop combinations
        ct_bench = (
            ct_2, ct_02, ct_11, ct_002, ct_011, ct_101, ct_0000, ct_0001, ct_0002,
            ct_0010, ct_0011, ct_0100, ct_0101, ct_1000, ct_1001 )

        # Convert to strings, discard order
        ct_bench_str = frozenset(str(cti) for cti in ct_bench)
        ct_all_str = frozenset(str(cti) for cti in ct_all)

        # Finally check
        self.assertEqual(ct_all_str, ct_bench_str)

#=========================================================================================
# Test IRSubtraction
#=========================================================================================

class IRSubstractionTest(unittest.TestCase):

    subtraction = sub.IRSubtraction(
        simple_qcd.model, coupling_types=('QCD', ), )

    def test_parent_PDGs(self):
        """Test determination of parent PDGs."""

        children1 = sub.SubtractionLegSet((
            sub.SubtractionLeg(1,  1, FINAL),
            sub.SubtractionLeg(2, 21, FINAL),
            sub.SubtractionLeg(3, -1, FINAL),
            sub.SubtractionLeg(4, -2, FINAL),
            sub.SubtractionLeg(5,  1, FINAL),
            sub.SubtractionLeg(8, -1, FINAL),
            sub.SubtractionLeg(9,  2, FINAL),
        ))
        self.assertEqual(self.subtraction.parent_PDGs_from_legs(children1), [21])

        children2 = sub.SubtractionLegSet((
            sub.SubtractionLeg(1,  1, INITIAL),
            sub.SubtractionLeg(2, 21, FINAL),
            sub.SubtractionLeg(3,  1, FINAL),
            sub.SubtractionLeg(4, -2, FINAL),
        ))
        self.assertEqual(self.subtraction.parent_PDGs_from_legs(children2), [2])

        children3 = sub.SubtractionLegSet((
            sub.SubtractionLeg(1,  1, INITIAL),
            sub.SubtractionLeg(2, 21, FINAL),
            sub.SubtractionLeg(3, -1, FINAL),
        ))
        self.assertEqual(self.subtraction.parent_PDGs_from_legs(children3), [])

#=========================================================================================
# Test counterterm generation
#=========================================================================================

class NLOSubtractionTest(unittest.TestCase):
    """Test class for the subtraction module."""

    mylegs = base_objects.LegList()
    myprocess = base_objects.Process()
    mysubtraction = None

    def setUp(self):

        # Setting up a process and its subtraction

        self.mylegs = base_objects.LegList([
            base_objects.Leg({'number': 1, 'id':  1, 'state': INITIAL}),
            base_objects.Leg({'number': 2, 'id': -1, 'state': INITIAL}),
            base_objects.Leg({'number': 3, 'id': 25, 'state': FINAL}),
            base_objects.Leg({'number': 4, 'id': 21, 'state': FINAL}),
            base_objects.Leg({'number': 5, 'id': 21, 'state': FINAL})
        ])

        self.myprocess = base_objects.Process({
            'legs': self.mylegs,
            'model': simple_qcd.model } )

        self.mysubtraction = sub.IRSubtraction(
            simple_qcd.model, coupling_types=('QCD', ) )

    def test_generation_of_elementary_structures(self):
        """Test generation of all elementary operators for selected process."""

        elem_structures_target = [
            sub.SoftStructure(self.mylegs[3]),
            sub.SoftStructure(self.mylegs[4]),
            sub.CollStructure(self.mylegs[3], self.mylegs[4]),
            sub.CollStructure(self.mylegs[0], self.mylegs[3]),
            sub.CollStructure(self.mylegs[0], self.mylegs[4]),
            sub.CollStructure(self.mylegs[1], self.mylegs[3]),
            sub.CollStructure(self.mylegs[1], self.mylegs[4]),
        ]

        elem_structures = self.mysubtraction.get_all_elementary_structures(
            self.myprocess, 1)

        self.assertEqual(
            set(str(op) for op in elem_structures),
            set(str(op) for op in elem_structures_target)
        )

    def test_structure_combinations(self):
        """Test the generation of all operator combinations for one selected process."""

        target_NLO_combos = [
            # 0
            '()',
            # 1
            '(S(4),)', '(S(5),)',
            '(C(4,5),)', '(C(1,4),)', '(C(1,5),)', '(C(2,4),)', '(C(2,5),)',
            # 2
            '(C(S(4),5),)', '(C(S(4),1),)', '(C(S(4),2),)',
            '(C(S(5),4),)', '(C(S(5),1),)', '(C(S(5),2),)'
        ]

        elem_structures = self.mysubtraction.get_all_elementary_structures(
            self.myprocess, 1)
        combos = self.mysubtraction.get_all_combinations(elem_structures, 1)

        self.assertEqual(
            set(target_NLO_combos),
            set(str(combo) for combo in combos)
        )

        all_currents = []

        for combo in combos:
            ct = self.mysubtraction.get_counterterm(combo, self.myprocess)
            for ct_n_loops in ct.distribute_loops(1):
                all_currents += ct_n_loops.get_all_currents()

        currents_to_store = []
        for current in all_currents:
            current.discard_leg_numbers()
            if current not in currents_to_store:
                currents_to_store += [current, ]

        self.assertEqual(len(currents_to_store), 14)

class NNLOSubtractionTest(unittest.TestCase):
    """Test class for the subtraction module at NNLO."""

    mylegs = base_objects.LegList()
    myprocess = base_objects.Process()
    mysubtraction = None

    def setUp(self):

        # Setting up a process and its subtraction

        self.mylegs = base_objects.LegList([
            base_objects.Leg({'number': 1, 'id':  25, 'state': INITIAL}),
            base_objects.Leg({'number': 2, 'id':  21, 'state': FINAL}),
            base_objects.Leg({'number': 3, 'id':  21, 'state': FINAL}),
            base_objects.Leg({'number': 4, 'id':  21, 'state': FINAL}),
            base_objects.Leg({'number': 5, 'id':   1, 'state': FINAL}),
            base_objects.Leg({'number': 6, 'id':  -1, 'state': FINAL}),
        ])

        self.myprocess = base_objects.Process({
            'legs': self.mylegs,
            'model': simple_qcd.model } )

        self.mysubtraction = sub.IRSubtraction(
            simple_qcd.model, coupling_types=('QCD', ), )

    def test_generation_of_elementary_structures_NN(self):
        """Test generation of all elementary operators for selected process."""

        elem_structures_target = [
            sub.SoftStructure(self.mylegs[1]),
            sub.SoftStructure(self.mylegs[2]),
            sub.SoftStructure(self.mylegs[3]),
            sub.SoftStructure(self.mylegs[1], self.mylegs[2]),
            sub.SoftStructure(self.mylegs[1], self.mylegs[3]),
            sub.SoftStructure(self.mylegs[2], self.mylegs[3]),
            sub.SoftStructure(self.mylegs[4], self.mylegs[5]),
            sub.CollStructure(self.mylegs[1], self.mylegs[2], self.mylegs[3]),
            sub.CollStructure(self.mylegs[1], self.mylegs[2], self.mylegs[4]),
            sub.CollStructure(self.mylegs[1], self.mylegs[2], self.mylegs[5]),
            sub.CollStructure(self.mylegs[1], self.mylegs[3], self.mylegs[4]),
            sub.CollStructure(self.mylegs[1], self.mylegs[3], self.mylegs[5]),
            sub.CollStructure(self.mylegs[1], self.mylegs[4], self.mylegs[5]),
            sub.CollStructure(self.mylegs[2], self.mylegs[3], self.mylegs[4]),
            sub.CollStructure(self.mylegs[2], self.mylegs[3], self.mylegs[5]),
            sub.CollStructure(self.mylegs[2], self.mylegs[4], self.mylegs[5]),
            sub.CollStructure(self.mylegs[3], self.mylegs[4], self.mylegs[5]),
            sub.CollStructure(self.mylegs[1], self.mylegs[2]),
            sub.CollStructure(self.mylegs[1], self.mylegs[3]),
            sub.CollStructure(self.mylegs[1], self.mylegs[4]),
            sub.CollStructure(self.mylegs[1], self.mylegs[5]),
            sub.CollStructure(self.mylegs[2], self.mylegs[3]),
            sub.CollStructure(self.mylegs[2], self.mylegs[4]),
            sub.CollStructure(self.mylegs[2], self.mylegs[5]),
            sub.CollStructure(self.mylegs[3], self.mylegs[4]),
            sub.CollStructure(self.mylegs[3], self.mylegs[5]),
            sub.CollStructure(self.mylegs[4], self.mylegs[5]),
        ]

        elem_structures = self.mysubtraction.get_all_elementary_structures(
            self.myprocess, 2)

        self.assertEqual(
            set(str(op) for op in elem_structures),
            set(str(op) for op in elem_structures_target)
        )

    def test_structure_combinations_NN(self):
        """Test the generation of all operator combinations for one selected process."""

        elem_structures = self.mysubtraction.get_all_elementary_structures(
            self.myprocess, 2)
        combos = self.mysubtraction.get_all_combinations(elem_structures, 2)

        for combo in combos:
            ct = self.mysubtraction.get_counterterm(combo, self.myprocess)
            # Remove the outermost (which is for beams) layer from counterterms,
            # except from the empty counterterm
            if str(combo) == "()":
                self.assertEqual(str(combo), str(ct))
            else:
                if str(combo) != str(ct)[1:-2]:
                    print ct.nice_string()
                self.assertEqual(str(combo), str(ct)[1:-2])

class HiggsN3LOSubtractionTest(unittest.TestCase):
    """Test the generation of counterterms for Higgs production at N3LO."""

    ggglegs = base_objects.LegList()
    ggg = base_objects.Process()
    mysubtraction = None

    def setUp(self):

        # g g > g g g H

        self.ggglegs = base_objects.LegList([
            base_objects.Leg({'number': 1, 'id': 21, 'state': INITIAL}),
            base_objects.Leg({'number': 2, 'id': 21, 'state': INITIAL}),
            base_objects.Leg({'number': 3, 'id': 21, 'state': FINAL}),
            base_objects.Leg({'number': 4, 'id': 21, 'state': FINAL}),
            base_objects.Leg({'number': 5, 'id': 21, 'state': FINAL}),
            base_objects.Leg({'number': 6, 'id': 25, 'state': FINAL}),
        ])

        self.ggg = base_objects.Process({
            'legs': self.ggglegs,
            'model': simple_qcd.model
        })

        self.mysubtraction = sub.IRSubtraction(
            simple_qcd.model, coupling_types=('QCD', ), )

    def test_generation_of_elementary_structures(self):
        """Test generation of all elementary operators for Higgs RRR."""

        # Shorthand
        ggglegs = self.ggglegs

        elem_structures_target = [
            # single-unresolved
            sub.SoftStructure(ggglegs[2]),
            sub.SoftStructure(ggglegs[3]),
            sub.SoftStructure(ggglegs[4]),
            sub.CollStructure(ggglegs[2], ggglegs[3]),
            sub.CollStructure(ggglegs[2], ggglegs[4]),
            sub.CollStructure(ggglegs[3], ggglegs[4]),
            sub.CollStructure(ggglegs[0], ggglegs[2]),
            sub.CollStructure(ggglegs[0], ggglegs[3]),
            sub.CollStructure(ggglegs[0], ggglegs[4]),
            sub.CollStructure(ggglegs[1], ggglegs[2]),
            sub.CollStructure(ggglegs[1], ggglegs[3]),
            sub.CollStructure(ggglegs[1], ggglegs[4]),
            # double-unresolved
            sub.SoftStructure(ggglegs[2], ggglegs[3]),
            sub.SoftStructure(ggglegs[2], ggglegs[4]),
            sub.SoftStructure(ggglegs[3], ggglegs[4]),
            sub.CollStructure(ggglegs[2], ggglegs[3], ggglegs[4]),
            sub.CollStructure(ggglegs[0], ggglegs[2], ggglegs[3]),
            sub.CollStructure(ggglegs[0], ggglegs[2], ggglegs[4]),
            sub.CollStructure(ggglegs[0], ggglegs[3], ggglegs[4]),
            sub.CollStructure(ggglegs[1], ggglegs[2], ggglegs[3]),
            sub.CollStructure(ggglegs[1], ggglegs[2], ggglegs[4]),
            sub.CollStructure(ggglegs[1], ggglegs[3], ggglegs[4]),
            # triple-unresolved
            sub.SoftStructure(ggglegs[2], ggglegs[3], ggglegs[4]),
            sub.CollStructure(ggglegs[0], ggglegs[2], ggglegs[3], ggglegs[4]),
            sub.CollStructure(ggglegs[1], ggglegs[2], ggglegs[3], ggglegs[4]),
        ]

        elem_structures = self.mysubtraction.get_all_elementary_structures(self.ggg, 3)

        self.assertEqual(
            set(str(op) for op in elem_structures),
            set(str(op) for op in elem_structures_target)
        )

    def test_structure_combinations(self):
        """Test the generation of all operator combinations for g g > g g g H."""

        elem_structures = self.mysubtraction.get_all_elementary_structures(self.ggg, 3)
        combos = self.mysubtraction.get_all_combinations(elem_structures, 3)
        # The number of counterterms may turn out to be incorrect,
        # but test that the generation goes through
        self.assertEqual(1398, len(combos))
