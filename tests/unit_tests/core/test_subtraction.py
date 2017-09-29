################################################################################
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
################################################################################

"""Unit test library for the routines of the core library related to writing
color information for diagrams.
"""

import random

import madgraph.core.base_objects as base_objects
import madgraph.core.subtraction as sub

import tests.unit_tests as unittest
import tests.input_files.simple_qcd as simple_qcd

#===============================================================================
# Test ordered tuples
#===============================================================================

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


#===============================================================================
# Test counterterm generation
#===============================================================================

class NLOSubtractionTest(unittest.TestCase):
    """Test class for the subtraction module."""

    mylegs = base_objects.LegList()
    myprocess = base_objects.Process()
    mysubtraction = None

    def setUp(self):

        # Setting up a process and its subtraction

        self.mylegs = base_objects.LegList([
            base_objects.Leg(
                    {'number': 1, 'id': 1, 'state': base_objects.Leg.INITIAL}),
            base_objects.Leg(
                    {'number': 2, 'id': -1, 'state': base_objects.Leg.INITIAL}),
            base_objects.Leg(
                    {'number': 3, 'id': 25, 'state': base_objects.Leg.FINAL}),
            base_objects.Leg(
                    {'number': 4, 'id': 21, 'state': base_objects.Leg.FINAL}),
            base_objects.Leg(
                    {'number': 5, 'id': 21, 'state': base_objects.Leg.FINAL})
        ])

        self.myprocess = base_objects.Process({
            'legs': self.mylegs,
            'model': simple_qcd.model
        })

        self.mysubtraction = sub.IRSubtraction(
                simple_qcd.model,
                orders = {'QCD': 2}
        )

    def test_singular_structure_init(self):
        """Test initialization of SingularStructure."""

        subtraction_leg_set = sub.SubtractionLegSet(
                self.mylegs[0], self.mylegs[2]
        )
        self.assertEqual(
                sub.SingularStructure(subtraction_leg_set),
                sub.SingularStructure(leg for leg in subtraction_leg_set)
        )

    def test_parent_PDGs(self):
        """Test determination of parent PDGs."""

        children1 = sub.SubtractionLegSet(
            sub.SubtractionLeg(1,  1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(2, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(3, -1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(4, -2, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(5,  1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(8, -1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(9,  2, sub.SubtractionLeg.FINAL),
        )
        self.assertEqual(self.mysubtraction.parent_PDGs_from_legs(children1), [21])

        children2 = sub.SubtractionLegSet(
            sub.SubtractionLeg(1,  1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(2, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(3,  1, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(4, -2, sub.SubtractionLeg.FINAL),
        )
        self.assertEqual(self.mysubtraction.parent_PDGs_from_legs(children2), [2])

        children3 = sub.SubtractionLegSet(
            sub.SubtractionLeg(1, 1, sub.SubtractionLeg.INITIAL),
            sub.SubtractionLeg(2, 21, sub.SubtractionLeg.FINAL),
            sub.SubtractionLeg(3, -1, sub.SubtractionLeg.FINAL),
        )
        self.assertEqual(self.mysubtraction.parent_PDGs_from_legs(children3), [])

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

        if self.mysubtraction.orders == {'QCD': 1}:
            self.assertEqual(
                    set(str(op) for op in elem_operators),
                    set(str(op) for op in elem_operators_target)
            )

    def test_act_on(self):
        """Test action of operators."""

        elem_operators = self.mysubtraction.get_all_elementary_operators(
            self.myprocess)

        C14 = sub.CollOperator(self.mylegs[0], self.mylegs[3])
        C145 = sub.CollOperator(self.mylegs[0], self.mylegs[3], self.mylegs[4])
        S4 = sub.SoftOperator(self.mylegs[3])
        S45 = sub.SoftOperator(self.mylegs[3], self.mylegs[4])

        CC_list = sub.SingularOperatorList([C14, C145])
        CC_simple = CC_list.simplify()
        CC_benchmark = sub.SingularStructure(sub.CollStructure(
                sub.CollStructure(
                        sub.SubtractionLeg(1,  1, sub.SubtractionLeg.INITIAL),
                        sub.SubtractionLeg(4, 21, sub.SubtractionLeg.FINAL)
                ),
                sub.SubtractionLeg(5, 21, sub.SubtractionLeg.FINAL)
        ))
        self.assertEqual(str(CC_simple), '(C(C(1,4),5),)')
        # self.assertEqual(CC_simple,CC_benchmark)

        SC_list = sub.SingularOperatorList([S4, C145])
        SC_simple = SC_list.simplify()
        self.assertEqual(str(SC_simple), '(C(S(4),1,5),)')

        SS_list = sub.SingularOperatorList([S4, S45])
        SS_simple = SS_list.simplify()
        self.assertEqual(SS_simple.is_void(), True)

    def test_count_unresolved(self):
        """Test counting of unresolved particles."""

        elem_operators = self.mysubtraction.get_all_elementary_operators(self.myprocess)

        C14  = sub.CollOperator(self.mylegs[0], self.mylegs[3])
        C145 = sub.CollOperator(self.mylegs[0], self.mylegs[3], self.mylegs[4])
        S4   = sub.SoftOperator(self.mylegs[3])
        S5   = sub.SoftOperator(self.mylegs[4])
        S45  = sub.SoftOperator(self.mylegs[3], self.mylegs[4])

        list1 = sub.SingularOperatorList([C145,]).simplify()
        self.assertEqual(
                self.mysubtraction.count_unresolved(list1),
                2
        )
        list2 = sub.SingularOperatorList([S45,]).simplify()
        self.assertEqual(
                self.mysubtraction.count_unresolved(list2),
                2
        )
        list3 = sub.SingularOperatorList([S45, C145,]).simplify()
        self.assertEqual(
                self.mysubtraction.count_unresolved(list3),
                2
        )
        list4 = sub.SingularOperatorList([S4, S5,]).simplify()
        self.assertEqual(
                self.mysubtraction.count_unresolved(list4),
                2
        )
        list5 = sub.SingularOperatorList([S4, S5, C145,]).simplify()
        self.assertEqual(
                self.mysubtraction.count_unresolved(list5),
                2
        )
        list6 = sub.SingularOperatorList([C14,]).simplify()
        self.assertEqual(
                self.mysubtraction.count_unresolved(list6),
                1
        )
        list7 = sub.SingularOperatorList([S4, C14,]).simplify()
        self.assertEqual(
                self.mysubtraction.count_unresolved(list7),
                1
        )

    def test_operator_combinations(self):
        """Test the generation of all elementary operators
        for one selected process.
        """

        target_combos = [
            # 0
            '()',
            # 1
            '(S(4),)', '(S(5),)',
            '(C(4,5),)', '(C(1,4),)', '(C(1,5),)', '(C(2,4),)', '(C(2,5),)',
            # 2
            '(S(4),S(5),)',
            '(C(S(4),5),)', '(C(S(4),1),)', '(C(S(4),2),)',
            '(C(1,5),S(4),)', '(C(2,5),S(4),)',
            '(C(S(5),4),)', '(C(S(5),1),)', '(C(S(5),2),)',
            '(C(1,4),S(5),)', '(C(2,4),S(5),)',
            '(C(1,4),C(2,5),)', '(C(1,5),C(2,4),)',
            # 3
            '(C(S(5),1),S(4),)', '(C(S(5),2),S(4),)',
            '(C(S(4),1),S(5),)', '(C(S(4),2),S(5),)',
            '(C(2,5),C(S(4),1),)', '(C(1,5),C(S(4),2),)',
            '(C(2,4),C(S(5),1),)', '(C(1,4),C(S(5),2),)',
            # 4
            '(C(S(4),1),C(S(5),2),)', '(C(S(4),2),C(S(5),1),)',
        ]

        elem_operators = self.mysubtraction.get_all_elementary_operators(self.myprocess)

        combos = self.mysubtraction.get_all_combinations(elem_operators)
        if self.mysubtraction.orders == {'QCD': 1}:
            self.assertEqual(
                    set(target_combos),
                    set(str(combo) for combo in combos)
            )

        target_filtered_NLO_combos = [
            # 0
            '()',
            # 1
            '(S(4),)', '(S(5),)',
            '(C(4,5),)', '(C(1,4),)', '(C(1,5),)', '(C(2,4),)', '(C(2,5),)',
            # 2
            '(C(S(4),5),)', '(C(S(4),1),)', '(C(S(4),2),)',
            '(C(S(5),4),)', '(C(S(5),1),)', '(C(S(5),2),)'
        ]

        filtered_NLO_combos = self.mysubtraction.filter_combinations(combos)
        if self.mysubtraction.orders == {'QCD': 1}:
            self.assertEqual(
                    set(target_filtered_NLO_combos),
                    set(str(combo) for combo in filtered_NLO_combos)
            )

        # sub_legs = [sub.SubtractionLeg(leg) for leg in self.mylegs]
        #
        # target_NLO_counterterms = [
        #     # S(4)
        #     sub.Current({
        #         'parent_subtraction_leg': sub.SoftStructure(sub_legs[4]),
        #         'singular_structure': sub.SoftStructure(sub_legs[3])
        #         }),
        #     # S(5)
        #     sub.Current({
        #         'parent_subtraction_leg': sub.SoftStructure(sub_legs[4]),
        #         'singular_structure': sub.SoftStructure(sub_legs[4])
        #         }),
        #     # TODO Keep going....
        # ]

        all_currents = []

        for combo in filtered_NLO_combos:
            print '-'*80
            print combo
            print "Prefactor:", combo.get_subtraction_prefactor()
            ct = self.mysubtraction.get_counterterm(combo, self.myprocess)
            print ct
            print ""
            for ct_n_loops in self.mysubtraction.split_loops(ct, 1):
                print ct_n_loops
                all_currents += ct_n_loops.get_all_currents()
        print len(filtered_NLO_combos)

        currents_to_store = []
        for current in all_currents:
            current.discard_leg_numbers()
            if current not in currents_to_store:
                currents_to_store += [current, ]

        for current in currents_to_store:
            print str(current)

        # self.assertEqual(set(elementary_NLO_currents), set(target_elementary_NLO_currents))
        #
        #
        # target_mapped_elementary_NLO_currents = {
        #     # TODO
        #     # ProcessKey(DefiningCurrent) : (DefiningCurrent, [list of mapped currents])
        #                                          }
        #
        # mapped_elementary_NLO_currents = self.mysubtraction.group_currents(elementary_NLO_currents)
        #
        # self.assertEqual(target_mapped_elementary_NLO_currents, mapped_elementary_NLO_currents)
