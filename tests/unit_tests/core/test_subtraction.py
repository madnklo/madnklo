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
            sub.SubtractionLeg(1,   1, INITIAL),
            sub.SubtractionLeg(4,  21, FINAL),
            sub.SubtractionLeg(5,  21, FINAL),
            sub.SubtractionLeg(14,  2, FINAL),
            sub.SubtractionLeg(11, -2, FINAL) )
        sub_2 = sub.CollStructure(
            sub.SubtractionLeg(1,  1, INITIAL),
            sub.SubtractionLeg(4, 21, FINAL) )
        sub_3 = sub.SoftStructure(
            sub.SubtractionLeg(14, 2, FINAL),
            sub.SubtractionLeg(11, -2, FINAL) )
        dec = a_structure.decompose()
        ben = {sub_1, sub_2, sub_3}
        dec_str = set(el.__str__(True, True, True) for el in dec)
        ben_str = set(el.__str__(True, True, True) for el in ben)
        self.assertEqual(dec_str, ben_str)

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
        ct_all = ct.split_loops(2)

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
        ct_all = ct.split_loops(2)

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
        ct_all = ct.split_loops(2)

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
        self.assertEqual(len(combos), 1026)
