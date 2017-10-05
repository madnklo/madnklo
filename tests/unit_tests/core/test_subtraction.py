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
# Test SingularStructure and SingularOperator
#=========================================================================================

class SingularStructureOperatorTest(unittest.TestCase):
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

        C14  = sub.CollOperator(leg1, leg4)
        C145 = sub.CollOperator(leg1, leg4, leg5)
        S4   = sub.SoftOperator(leg4)
        S45  = sub.SoftOperator(leg4, leg5)

        CC_list = sub.SingularOperatorList([C14, C145])
        CC_simple = CC_list.simplify()
        CC_benchmark = sub.SingularStructure(sub.CollStructure(
            sub.CollStructure(
                sub.SubtractionLeg(1,  1, INITIAL),
                sub.SubtractionLeg(4, 21, FINAL)
            ),
            sub.SubtractionLeg(5, 21, FINAL)
        ))
        self.assertEqual(CC_simple,CC_benchmark)

        SC_list = sub.SingularOperatorList([S4, C145])
        SC_simple = SC_list.simplify()
        self.assertEqual(str(SC_simple), '(C(S(4),1,5),)')

        SS_list = sub.SingularOperatorList([S4, S45])
        SS_simple = SS_list.simplify()
        self.assertEqual(SS_simple.is_void, True)

    def test_count_unresolved(self):
        """Test counting of unresolved particles."""

        leg1 = base_objects.Leg({'number': 1, 'id':  1, 'state': INITIAL})
        leg4 = base_objects.Leg({'number': 4, 'id': 21, 'state': FINAL})
        leg5 = base_objects.Leg({'number': 5, 'id': 21, 'state': FINAL})

        C14  = sub.CollOperator(leg1, leg4)
        C145 = sub.CollOperator(leg1, leg4, leg5)
        S4   = sub.SoftOperator(leg4)
        S5   = sub.SoftOperator(leg5)
        S45  = sub.SoftOperator(leg4, leg5)

        list1 = sub.SingularOperatorList([C145, ]).simplify()
        self.assertEqual(list1.count_unresolved(), 2)

        list2 = sub.SingularOperatorList([S45, ]).simplify()
        self.assertEqual(list2.count_unresolved(), 2)

        list3 = sub.SingularOperatorList([S45, C145, ]).simplify()
        self.assertEqual(list3.count_unresolved(), 2)

        list4 = sub.SingularOperatorList([S4, S5, ]).simplify()
        self.assertEqual(list4.count_unresolved(), 2)

        list5 = sub.SingularOperatorList([S4, S5, C145, ]).simplify()
        self.assertEqual(list5.count_unresolved(), 2)

        list6 = sub.SingularOperatorList([C14, ]).simplify()
        self.assertEqual(list6.count_unresolved(), 1)

        list7 = sub.SingularOperatorList([S4, C14, ]).simplify()
        self.assertEqual(list7.count_unresolved(), 1)

#=========================================================================================
# Test Counterterm (and CountertermNode)
#=========================================================================================

class CountertermTest(unittest.TestCase):

    # TODO Move test_split_loops and test_split_orders here
    pass

#=========================================================================================
# Test IRSubtraction
#=========================================================================================

class IRSubstractionTest(unittest.TestCase):

    subtraction = sub.IRSubtraction(
        simple_qcd.model, coupling_types=('QCD', ), n_unresolved=2 )

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

    def test_split_loops(self):
        """Test the assignment of loop numbers within a current."""

        leg1 = base_objects.Leg({'number': 1, 'id':  1, 'state': INITIAL})
        leg2 = base_objects.Leg({'number': 2, 'id': -1, 'state': INITIAL})
        leg3 = base_objects.Leg({'number': 3, 'id': 21, 'state': FINAL})
        leg4 = base_objects.Leg({'number': 4, 'id': 21, 'state': FINAL})
        leg5 = base_objects.Leg({'number': 5, 'id': 21, 'state': FINAL})
        leg6 = base_objects.Leg({'number': 6, 'id': 22, 'state': FINAL})

        leg_list = base_objects.LegList([leg1, leg2, leg3, leg4, leg5, leg6])

        myprocess = base_objects.Process({
            'legs': leg_list,
            'model': simple_qcd.model
        })

        C14  = sub.CollOperator(leg1, leg4)
        C145 = sub.CollOperator(leg1, leg4, leg5)
        S3   = sub.SoftOperator(leg3)
        S5   = sub.SoftOperator(leg5)
        S45  = sub.SoftOperator(leg4, leg5)

        ct1 = self.subtraction.get_counterterm(
            sub.SingularOperatorList([C145, S3]).simplify(),
            myprocess
        )

        ct1_0 = ct1.get_copy()
        ct1_0.current.set('n_loops', 0)
        ct1_1 = ct1.get_copy()
        ct1_1.current.set('n_loops', 1)
        ct1_2 = ct1.get_copy()
        ct1_2.current.set('n_loops', 2)

        ct1_00 = ct1_0.get_copy()
        ct1_00.nodes[0].current.set('n_loops', 0)
        ct1_01 = ct1_0.get_copy()
        ct1_01.nodes[0].current.set('n_loops', 1)
        ct1_02 = ct1_0.get_copy()
        ct1_02.nodes[0].current.set('n_loops', 2)
        ct1_10 = ct1_1.get_copy()
        ct1_10.nodes[0].current.set('n_loops', 0)
        ct1_11 = ct1_1.get_copy()
        ct1_11.nodes[0].current.set('n_loops', 1)
        ct1_20 = ct1_2.get_copy()
        ct1_20.nodes[0].current.set('n_loops', 0)

        ct1_002 = ct1_00.get_copy()
        ct1_002.nodes[1].current.set('n_loops', 2)
        ct1_011 = ct1_01.get_copy()
        ct1_011.nodes[1].current.set('n_loops', 1)
        ct1_020 = ct1_02.get_copy()
        ct1_020.nodes[1].current.set('n_loops', 0)
        ct1_101 = ct1_10.get_copy()
        ct1_101.nodes[1].current.set('n_loops', 1)
        ct1_110 = ct1_11.get_copy()
        ct1_110.nodes[1].current.set('n_loops', 0)
        ct1_200 = ct1_20.get_copy()
        ct1_200.nodes[1].current.set('n_loops', 0)
        ct1_bench = (ct1_200, ct1_101, ct1_110, ct1_002, ct1_011, ct1_020)
        ct1_bench_str = frozenset(str(ct) for ct in ct1_bench)
        ct1_all = ct1.split_loops(2)
        ct1_all_str = frozenset(str(ct) for ct in ct1_all)
        self.assertEqual(ct1_all_str, ct1_bench_str)

    def test_split_orders(self):
        """Test the assignment of coupling orders within a current."""

        pass

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
            simple_qcd.model, coupling_types=('QCD', ), n_unresolved=1 )

    def test_generation_of_elementary_operators(self):
        """Test generation of all elementary operators for selected process."""

        elem_operators_target = [
            sub.SoftOperator(self.mylegs[3]),
            sub.SoftOperator(self.mylegs[4]),
            sub.CollOperator(self.mylegs[3], self.mylegs[4]),
            sub.CollOperator(self.mylegs[0], self.mylegs[3]),
            sub.CollOperator(self.mylegs[0], self.mylegs[4]),
            sub.CollOperator(self.mylegs[1], self.mylegs[3]),
            sub.CollOperator(self.mylegs[1], self.mylegs[4]),
        ]

        elem_operators = self.mysubtraction.get_all_elementary_operators(self.myprocess)

        self.assertEqual(
            set(str(op) for op in elem_operators),
            set(str(op) for op in elem_operators_target)
        )

    def test_operator_combinations(self):
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

        elem_operators = self.mysubtraction.get_all_elementary_operators(self.myprocess)
        combos = self.mysubtraction.get_all_combinations(elem_operators)

        self.assertEqual(
            set(target_NLO_combos),
            set(str(combo) for combo in combos)
        )

        all_currents = []

        for combo in combos:
            ct = self.mysubtraction.get_counterterm(combo, self.myprocess)
            for ct_n_loops in ct.split_loops(1):
                all_currents += ct_n_loops.get_all_currents()

        currents_to_store = []
        for current in all_currents:
            current.discard_leg_numbers()
            if current not in currents_to_store:
                currents_to_store += [current, ]

        self.assertEqual(len(currents_to_store), 14)

class HiggsN3LOSubtractionTest(unittest.TestCase):
    """Test the generation of counterterms for Higgs production at N3LO."""

    ggglegs = base_objects.LegList()
    ggg = base_objects.Process()
    mysubtraction = None

    def setUp(self):

        # g g > g g g H

        self.ggglegs = base_objects.LegList([
            base_objects.Leg(
                    {'number': 1, 'id': 21, 'state': base_objects.Leg.INITIAL}),
            base_objects.Leg(
                    {'number': 2, 'id': 21, 'state': base_objects.Leg.INITIAL}),
            base_objects.Leg(
                    {'number': 3, 'id': 21, 'state': base_objects.Leg.FINAL}),
            base_objects.Leg(
                    {'number': 4, 'id': 21, 'state': base_objects.Leg.FINAL}),
            base_objects.Leg(
                    {'number': 5, 'id': 21, 'state': base_objects.Leg.FINAL}),
            base_objects.Leg(
                    {'number': 6, 'id': 25, 'state': base_objects.Leg.FINAL}),
        ])

        self.ggg = base_objects.Process({
            'legs': self.ggglegs,
            'model': simple_qcd.model
        })

        self.mysubtraction = sub.IRSubtraction(
            simple_qcd.model, coupling_types=('QCD', ), n_unresolved=3 )

    def test_generation_of_elementary_operators(self):
        """Test generation of all elementary operators for Higgs RRR."""

        # Shorthand
        ggglegs = self.ggglegs

        elem_operators_target = [
            # single-unresolved
            sub.SoftOperator(ggglegs[2]),
            sub.SoftOperator(ggglegs[3]),
            sub.SoftOperator(ggglegs[4]),
            sub.CollOperator(ggglegs[2], ggglegs[3]),
            sub.CollOperator(ggglegs[2], ggglegs[4]),
            sub.CollOperator(ggglegs[3], ggglegs[4]),
            sub.CollOperator(ggglegs[0], ggglegs[2]),
            sub.CollOperator(ggglegs[0], ggglegs[3]),
            sub.CollOperator(ggglegs[0], ggglegs[4]),
            sub.CollOperator(ggglegs[1], ggglegs[2]),
            sub.CollOperator(ggglegs[1], ggglegs[3]),
            sub.CollOperator(ggglegs[1], ggglegs[4]),
            # double-unresolved
            sub.SoftOperator(ggglegs[2], ggglegs[3]),
            sub.SoftOperator(ggglegs[2], ggglegs[4]),
            sub.SoftOperator(ggglegs[3], ggglegs[4]),
            sub.CollOperator(ggglegs[2], ggglegs[3], ggglegs[4]),
            sub.CollOperator(ggglegs[0], ggglegs[2], ggglegs[3]),
            sub.CollOperator(ggglegs[0], ggglegs[2], ggglegs[4]),
            sub.CollOperator(ggglegs[0], ggglegs[3], ggglegs[4]),
            sub.CollOperator(ggglegs[1], ggglegs[2], ggglegs[3]),
            sub.CollOperator(ggglegs[1], ggglegs[2], ggglegs[4]),
            sub.CollOperator(ggglegs[1], ggglegs[3], ggglegs[4]),
            # triple-unresolved
            sub.SoftOperator(ggglegs[2], ggglegs[3], ggglegs[4]),
            sub.CollOperator(ggglegs[0], ggglegs[2], ggglegs[3], ggglegs[4]),
            sub.CollOperator(ggglegs[1], ggglegs[2], ggglegs[3], ggglegs[4]),
        ]

        elem_operators = self.mysubtraction.get_all_elementary_operators(self.ggg)

        self.assertEqual(
            set(str(op) for op in elem_operators),
            set(str(op) for op in elem_operators_target)
        )

    def test_operator_combinations(self):
        """Test the generation of all operator combinations for g g > g g g H."""

        elem_operators = self.mysubtraction.get_all_elementary_operators(self.ggg)
        combos = self.mysubtraction.get_all_combinations(elem_operators)
        # The number of counterterms may turn out to be incorrect,
        # but test that the generation goes through
        self.assertEqual(len(combos), 4752)
