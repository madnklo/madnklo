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
"""Unit test library for the various phase-space generation/handling features."""

import madgraph.integrator.phase_space_generators as PS
import madgraph.core.subtraction as subtraction
import madgraph.core.base_objects as base_objects
import madgraph.core.color_algebra as color

import copy
import math
import random

import tests.unit_tests as unittest

#===============================================================================
# Test Vectors
#===============================================================================

class VectorsTest(unittest.TestCase):
    """Test class for BaseVector, Vector and LorentzVector."""

    def setUp(self):
        pass

    def test_Vector_basic(self):
        """Test the basic operations for vectors."""

        # Component-wise addition and subtraction
        v1 = PS.Vector([1., 0., 2., 3.])
        v2 = PS.Vector(4*[0.5, ])
        self.assertEqual(v1 + v2, PS.Vector([1.5, 0.5, 2.5, 3.5]))
        v3 = PS.Vector([1., 2., 0.5, -1.])
        self.assertEqual(v1 - v3, PS.Vector([0., -2., 1.5, 4.]))
        v3 += PS.Vector([1., 0., 0., 0.])
        self.assertEqual(v3, PS.Vector([2., 2., 0.5, -1.]))
        v3 -= PS.Vector([0.5, 1.5, 0., 0.])
        self.assertEqual(v3, PS.Vector([1.5, 0.5, 0.5, -1.]))

        # Multiplication and division by scalars
        self.assertEqual(v2 * 2., PS.Vector(4*[1., ]))
        self.assertEqual(v1 / 4., PS.Vector([0.25, 0., 0.5, 0.75]))
        v3 *= 3
        self.assertEqual(v3, PS.Vector([4.5, 1.5, 1.5, -3.]))
        v3 /= 1.5
        self.assertEqual(v3, PS.Vector([3., 1., 1., -2.]))
        self.assertEqual(3 * v1, v1 * 3)

        # Negative
        self.assertEqual(-v3, PS.Vector([-3., -1., -1., 2.]))

    def test_Vector_Euclid(self):
        """Test scalar products and related functions for the class Vector."""

        # Square and norm
        v1 = PS.Vector([0., 1., -2.])
        self.assertEqual(v1.square(), 5.)
        v2 = PS.Vector([3., 4., 0.])
        self.assertEqual(abs(v2), 5.)
        v2n = 0.2 * v2
        v2.normalize()
        self.assertEqual(v2, v2n)
        self.assertEqual(PS.Vector.dot(v1, v2), 0.8)
        v3 = PS.Vector([random.random() for _ in range(3)])
        w = PS.Vector([random.random() for _ in range(3)])
        v3p = v3.project_onto(w)
        v3t = v3.component_orthogonal_to(w)
        self.assertEqual(v3, v3p + v3t)
        self.assertAlmostEqual(PS.Vector.dot(v3t, w), 0.)
        self.assertAlmostEqual(abs(v3p+w), abs(v3p)+abs(w))

    def test_Vector_Minkowski(self):
        """Test the relativistic norm."""

        # Square and norm
        v1 = PS.LorentzVector([1, 0, 0, +1])
        v2 = PS.LorentzVector([1, 0, 0, -1])
        self.assertEqual(v1.square(), 0)
        self.assertEqual(v2.square(), 0)
        self.assertEqual(v1.dot(v2), 2)

#===============================================================================
# Test the phase-space generators
#===============================================================================
class PhaseSpaceGeneratorsTest(unittest.TestCase):
    """ Test various phase-space generators."""

    def setUp(self):
        pass

    def test_flat_invertible_phase_space(self):
        """ Tests the flat invertible phase-space."""
        
        E_cm  = 5000.0
    
        # Try to run the above for a 2->8.
        my_PS_generator = PS.FlatInvertiblePhasespace(
            [0.]*2, [100. + 10.*i for i in range(8)],beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0))
        # Try to run the above for a 2->1.    
        #    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [5000.0])
        
        random_variables = [random.random() for _ in range(my_PS_generator.nDimPhaseSpace())]
    
        momenta, wgt = my_PS_generator.generateKinematics(E_cm, random_variables)
       
        #print "\n ========================="
        #print " ||    PS generation    ||"
        #print " ========================="    
        #print "\nRandom variables :\n",random_variables
        #print "\n%s\n"%my_PS_generator.nice_momenta_string(
        #                    momenta, recompute_mass=True, n_initial=my_PS_generator.n_initial)
        #print "Phase-space weight : %.16e\n"%wgt,
    
        variables_reconstructed, wgt_reconstructed = \
                                             my_PS_generator.invertKinematics(E_cm, momenta)
    
        #print "\n ========================="
        #print " || Kinematic inversion ||"
        #print " ========================="
        #print "\nReconstructed random variables :\n",variables_reconstructed
        differences = [abs(variables_reconstructed[i]-random_variables[i]) 
                                        for i in range(len(variables_reconstructed))]

        self.assertLess(max(differences[i]/random_variables[i] for i in range(len(differences))), 1.0e-10)
        self.assertLess(abs(wgt-wgt_reconstructed)/abs(wgt), 1.0e-10)
        
        #print "Reconstructed weight = %.16e"%wgt_reconstructed
        #if differences:
        #    print "\nMax. relative diff. in reconstructed variables = %.3e"%\
        #        max(differences[i]/random_variables[i] for i in range(len(differences)))
        #print "Rel. diff. in PS weight = %.3e\n"%((wgt_reconstructed-wgt)/wgt)

#===============================================================================
# Test the phase-space mappers
#===============================================================================

class CollinearVariablesTest(unittest.TestCase):
    """Test class for variables describing internal collinear structure."""

    my_mapping = PS.ElementaryMappingCollinearFinal()
    n_children = 3

    def setUp(self):
        pass

    def test_variables_away_from_limit(self):
        """Test determination of collinear variables and reverse mapping,
        for completely generic input values.
        """

        # Generate n_children+1 random massless vectors
        # (the light-cone direction is also random)
        my_PS_point = dict()
        for i in range(self.n_children+1):
            my_PS_point[i] = PS.LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            )
            my_PS_point[i][0] = math.sqrt(-my_PS_point[i].square())
        # Compute collinear variables
        variables = dict()
        self.my_mapping.get_collinear_variables(
            my_PS_point, self.n_children, range(self.n_children),
            variables
        )
        # Compute total momentum
        total_momentum = PS.LorentzVector(4*[0., ])
        for i in range(self.n_children):
            total_momentum += my_PS_point[i]
        # Compute new phase space point
        new_PS_point = dict()
        new_PS_point[self.n_children] = my_PS_point[self.n_children]
        self.my_mapping.set_collinear_variables(
            new_PS_point, self.n_children, range(self.n_children),
            total_momentum, variables
        )
        for i in range(self.n_children):
            for j in range(4):
                self.assertAlmostEqual(my_PS_point[i][j], new_PS_point[i][j])


    def test_variables_close_to_limit(self):
        """Test determination of collinear variables and reverse mapping
        in a typical collinear situation.
        """

        # Generate n_children random starting directions
        directions = {
            i: PS.Vector([random.random() for _ in range(3)])
            for i in range(self.n_children)
        }
        coll_direction = PS.Vector(3*[0., ])
        for n in directions.values():
            coll_direction += n
        # Generate values for the parameter that describes approach to limit
        pars = (math.pow(0.25, i) for i in range(10))
        # This is one way like any other to approach the limit
        for par in pars:
            new_directions = {
                i: (par*n + (1-par)*coll_direction)
                for (i, n) in directions.items()
            }
            my_PS_point = {
                i: PS.LorentzVector([math.sqrt(n.square()),] + list(n))
                for (i, n) in new_directions.items()
            }
            my_PS_point[self.n_children] = PS.LorentzVector(
                [math.sqrt(n.square()), ] + list(n)
            )
            # print "Phase space point", my_PS_point
            # Compute collinear variables
            variables = dict()
            self.my_mapping.get_collinear_variables(
                my_PS_point, self.n_children, range(self.n_children),
                variables
            )
            # Compute total momentum
            total_momentum = PS.LorentzVector(4*[0., ])
            for i in range(self.n_children):
                total_momentum += my_PS_point[i]
            # Compute new phase space point
            new_PS_point = dict()
            new_PS_point[self.n_children] = my_PS_point[self.n_children]
            self.my_mapping.set_collinear_variables(
                new_PS_point, self.n_children, range(self.n_children),
                total_momentum, variables
            )
            for i in range(self.n_children):
                for j in range(4):
                    self.assertAlmostEqual(my_PS_point[i][j], new_PS_point[i][j])

class CataniSeymourFFOneTest(unittest.TestCase):
    """Test class for MappingCataniSeymourFFOne."""

    mapping = PS.MappingCataniSeymourFFOne()
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
            this_bunch_legs = (
                subtraction.SubtractionLeg(
                    i, 21, subtraction.SubtractionLeg.FINAL
                ) for i in this_bunch_numbers
            )
            self.structure += [subtraction.CollStructure(this_bunch_legs), ]
            if n_in_this_bunch != 1:
                self.momenta_dict[n_tot + n_bunch] = frozenset(this_bunch_numbers)
            n_collinear_so_far += n_in_this_bunch
            n_bunch += 1
        self.structure += [
                subtraction.SubtractionLeg(
                    i, 21, subtraction.SubtractionLeg.FINAL
                )
                for i in range(n_collinear_so_far, n_tot)
            ]
        self.structure = subtraction.SingularStructure(self.structure)

    def test_map_invertible(self):
        """Test mapping and inverse."""

        # Generate n_collinear random massive vectors plus
        # n_recoilers (massive = True, False) random vectors
        my_PS_point = dict()
        n_collinear_so_far = 0
        for n_in_this_bunch in self.n_collinear:
            for i in range(
                n_collinear_so_far, n_collinear_so_far + n_in_this_bunch
            ):
                my_PS_point[i] = PS.LorentzVector(
                    [0., ] + [random.random() for _ in range(3)]
                )
                my_PS_point[i][0] = math.sqrt(-my_PS_point[i].square())
            if n_in_this_bunch == 1:
                my_PS_point[n_collinear_so_far][0] = math.sqrt(
                    random.random() - my_PS_point[n_collinear_so_far].square()
                )
            n_collinear_so_far += n_in_this_bunch
        for i in range(
            n_collinear_so_far, n_collinear_so_far + self.n_recoilers
        ):
            my_PS_point[i] = PS.LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            )
            if self.massive:
                my_PS_point[i][0] = math.sqrt(
                    random.random() - my_PS_point[i].square()
                )
            else:
                my_PS_point[i][0] = math.sqrt(-my_PS_point[i].square())
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
            for j in range(4):
                self.assertAlmostEqual(my_PS_point[i][j], old_PS_point[i][j])

class FlatCollinearWalkerTest(unittest.TestCase):
    """Test class for FlatCollinearWalker."""

    walker = PS.FlatCollinearWalker()

    mypartlist = base_objects.ParticleList()
    myinterlist = base_objects.InteractionList()
    mymodel = base_objects.Model()
    mylegs = base_objects.LegList()
    process = base_objects.Process()
    subtraction = None

    def setUp(self):
        # Setting up a dumb model

        # A gluon
        self.mypartlist.append(base_objects.Particle({'name': 'g',
                                                      'antiname': 'g',
                                                      'spin': 3,
                                                      'color': 8,
                                                      'mass': 'zero',
                                                      'width': 'zero',
                                                      'texname': 'g',
                                                      'antitexname': 'g',
                                                      'line': 'curly',
                                                      'charge': 0.,
                                                      'pdg_code': 21,
                                                      'propagating': True,
                                                      'is_part': True,
                                                      'self_antipart': True}))

        # A quark U and its antiparticle
        self.mypartlist.append(base_objects.Particle({'name': 'u',
                                                      'antiname': 'u~',
                                                      'spin': 2,
                                                      'color': 3,
                                                      'mass': 'zero',
                                                      'width': 'zero',
                                                      'texname': 'u',
                                                      'antitexname': '\bar u',
                                                      'line': 'straight',
                                                      'charge': 2. / 3.,
                                                      'pdg_code': 2,
                                                      'propagating': True,
                                                      'is_part': True,
                                                      'self_antipart': False}))
        antiu = copy.copy(self.mypartlist[1])
        antiu.set('is_part', False)

        # A quark D and its antiparticle
        self.mypartlist.append(base_objects.Particle({'name': 'd',
                                                      'antiname': 'd~',
                                                      'spin': 2,
                                                      'color': 3,
                                                      'mass': 'zero',
                                                      'width': 'zero',
                                                      'texname': 'u',
                                                      'antitexname': '\bar u',
                                                      'line': 'straight',
                                                      'charge': -1. / 3.,
                                                      'pdg_code': 1,
                                                      'propagating': True,
                                                      'is_part': True,
                                                      'self_antipart': False}))
        antid = copy.copy(self.mypartlist[2])
        antid.set('is_part', False)

        # A photon
        self.mypartlist.append(base_objects.Particle({'name': 'a',
                                                      'antiname': 'a',
                                                      'spin': 3,
                                                      'color': 1,
                                                      'mass': 'zero',
                                                      'width': 'zero',
                                                      'texname': '\gamma',
                                                      'antitexname': '\gamma',
                                                      'line': 'wavy',
                                                      'charge': 0.,
                                                      'pdg_code': 22,
                                                      'propagating': True,
                                                      'is_part': True,
                                                      'self_antipart': True}))

        # A Higgs
        self.mypartlist.append(base_objects.Particle({'name': 'h',
                                                      'antiname': 'h',
                                                      'spin': 1,
                                                      'color': 1,
                                                      'mass': 'mh',
                                                      'width': 'wh',
                                                      'texname': 'h',
                                                      'antitexname': 'h',
                                                      'line': 'dashed',
                                                      'charge': 0.,
                                                      'pdg_code': 25,
                                                      'propagating': True,
                                                      'is_part': True,
                                                      'self_antipart': True}))

        # 3 gluon vertiex
        self.myinterlist.append(base_objects.Interaction({
            'id': 1,
            'particles': base_objects.ParticleList(
                [self.mypartlist[0]] * 3
            ),
            'color': [color.ColorString([color.f(0, 1, 2)])],
            'lorentz': ['L1'],
            'couplings': {(0, 0): 'G'},
            'orders': {'QCD': 1}}))

        # 4 gluon vertex
        self.myinterlist.append(base_objects.Interaction({
            'id': 2,
            'particles': base_objects.ParticleList(
                [self.mypartlist[0]] * 4
            ),
            'color': [color.ColorString([color.f(-1, 0, 2),
                                         color.f(-1, 1, 3)]),
                      color.ColorString([color.f(-1, 0, 3),
                                         color.f(-1, 1, 2)]),
                      color.ColorString([color.f(-1, 0, 1),
                                         color.f(-1, 2, 3)])],
            'lorentz': ['L(p1,p2,p3)', 'L(p2,p3,p1)', 'L3'],
            'couplings': {(0, 0): 'G^2',
                          (1, 1): 'G^2',
                          (2, 2): 'G^2'},
            'orders': {'QCD': 2}}))

        # Gluon couplings to up and down quarks
        self.myinterlist.append(base_objects.Interaction({
            'id': 3,
            'particles': base_objects.ParticleList(
                [self.mypartlist[1],
                 antiu,
                 self.mypartlist[0]]),
            'color': [color.ColorString([color.T(2, 0, 1)])],
            'lorentz': ['L1'],
            'couplings': {(0, 0): 'GQQ'},
            'orders': {'QCD': 1}}))

        self.myinterlist.append(base_objects.Interaction({
            'id': 4,
            'particles': base_objects.ParticleList(
                [self.mypartlist[2],
                 antid,
                 self.mypartlist[0]]),
            'color': [color.ColorString([color.T(2, 0, 1)])],
            'lorentz': ['L1'],
            'couplings': {(0, 0): 'GQQ'},
            'orders': {'QCD': 1}}))

        # Photon coupling to up
        self.myinterlist.append(base_objects.Interaction({
            'id': 5,
            'particles': base_objects.ParticleList(
                [self.mypartlist[1],
                 antiu,
                 self.mypartlist[3]]),
            'color': [color.ColorString([color.T(0, 1)])],
            'lorentz': ['L1'],
            'couplings': {(0, 0): 'GQED'},
            'orders': {'QED': 1}}))

        self.mymodel.set('particles', self.mypartlist)
        self.mymodel.set('interactions', self.myinterlist)
        self.mymodel.set('name', "sm4test")

        # Setting up a process and its subtraction

        self.mylegs = base_objects.LegList([
            base_objects.Leg(
                {'number': 1, 'id': 25, 'state': base_objects.Leg.INITIAL}),
            base_objects.Leg(
                {'number': 2, 'id':  1, 'state': base_objects.Leg.FINAL}),
            base_objects.Leg(
                {'number': 3, 'id': -1, 'state': base_objects.Leg.FINAL}),
            base_objects.Leg(
                {'number': 4, 'id':  2, 'state': base_objects.Leg.FINAL}),
            base_objects.Leg(
                {'number': 5, 'id': -2, 'state': base_objects.Leg.FINAL}),
            base_objects.Leg(
                {'number': 6, 'id': 25, 'state': base_objects.Leg.FINAL}),
        ])

        self.process = base_objects.Process({
            'legs': self.mylegs,
            'model': self.mymodel,
            'n_loops': 0
        })

        self.subtraction = subtraction.IRSubtraction(
            self.mymodel,
            orders={'QCD': 1}
        )

        operators = self.subtraction.get_all_elementary_operators(self.process)
        self.combinations = self.subtraction.get_all_combinations(operators)
        self.counterterms = [
            self.subtraction.get_counterterm(combination, self.process)
            for combination in self.combinations
        ]
        for i in self.counterterms:
            print str(i)

    def test_walk_invertible(self):
        """Test walk and inverse."""

        # For each counterterm
        for i in range(len(self.counterterms)):
            # Generate random vectors
            my_PS_point = dict()
            legs_FS = (
                subtraction.SubtractionLeg(leg)
                for leg in self.process['legs']
                if leg['state'] == base_objects.Leg.FINAL
            )
            leg_numbers = (leg.n for leg in legs_FS)
            for j in leg_numbers:
                my_PS_point[j] = PS.LorentzVector(
                    [0., ] + [random.random() for _ in range(3)]
                )
                my_PS_point[j][0] = math.sqrt(-my_PS_point[j].square())
            # Compute collinear variables
            res_dict1 = self.walker.walk_to_lower_multiplicity(
                my_PS_point, self.counterterms[i], True
            )
            (currs1, ME1, jac1, kin)  = (
                res_dict1['currents'], res_dict1['matrix_element'],
                res_dict1['jacobian'], res_dict1['kinematic_variables'] )
            
            res_dict2 = self.walker.walk_to_higher_multiplicity(
                ME1[1], self.counterterms[i], kin
            )
            (currs2, ME2, jac2)  = (
                res_dict2['currents'], res_dict2['matrix_element'],
                res_dict2['jacobian'] )
            
            # Check currents
            self.assertEqual(len(currs1), len(currs2))
            for i_curr in range(len(currs1)):
                self.assertEqual(
                    currs1[i_curr][0], currs2[i_curr][0]
                )
                self.assertEqual(
                    currs1[i_curr][1].keys(), currs2[i_curr][1].keys()
                )
                for i_part in currs1[i_curr][1].keys():
                    for i_comp in range(4):
                        self.assertAlmostEqual(
                            currs1[i_curr][1][i_part][i_comp],
                            currs2[i_curr][1][i_part][i_comp]
                        )
            # Check MEs
            self.assertEqual(
                ME1[0], ME2[0]
            )
            self.assertEqual(
                ME1[1].keys(), ME2[1].keys()
            )
            for i_part in ME1[1].keys():
                for i_comp in range(4):
                    self.assertAlmostEqual(
                        ME1[1][i_part][i_comp],
                        ME2[1][i_part][i_comp]
                    )
            # Check jacobians
            self.assertAlmostEqual(jac1*jac2, 1.)
            print ME1[1]
