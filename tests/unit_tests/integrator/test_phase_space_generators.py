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
import madgraph.various.misc as misc
import models.import_ufo as import_ufo
from madgraph import MG5DIR

import copy
import math
import random
import os

import tests.unit_tests as unittest
import tests.input_files.simple_qcd as simple_qcd

pjoin = os.path.join

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

        # Test the rotoboost
        #   Definition
        v1 = PS.LorentzVector([0, ] + [random.random() for _ in range(3)])
        v2 = PS.LorentzVector([0, ] + [random.random() for _ in range(3)])
        m = random.random()
        v1.set_square(m**2)
        v2.set_square(m**2)
        self.assertEqual(v2.rotoboost(v2, v1), v1)
        #   Inverse
        v3 = PS.LorentzVector([0, ] + [random.random() for _ in range(3)])
        v3.set_square(m**2)
        v4 = PS.LorentzVector([random.random() for _ in range(4)])
        v4_old = v4.get_copy()
        v4.rotoboost(v1, v3)
        v4.rotoboost(v3, v1)
        self.assertEqual(v4, v4_old)

#===============================================================================
# Test the phase-space generators
#===============================================================================

class PhaseSpaceGeneratorsTest(unittest.TestCase):
    """ Test various phase-space generators."""

    def setUp(self):
        """Instantiate a model,
        which will be useful for the non-flat phase-space generator test.
        """
        
        model_with_params_set = import_ufo.import_model(
                pjoin(MG5DIR,'models','sm'), prefix=True,
                complex_mass_scheme = False )
        model_with_params_set.pass_particles_name_in_mg_default()
        model_with_params_set.set_parameters_and_couplings(
                param_card = pjoin(MG5DIR,'models','sm','restrict_default.dat'),
                complex_mass_scheme=False)
        self.model = model_with_params_set

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
        #print "\n%s\n"%momenta.__str__(n_initial=my_PS_generator.n_initial)
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

    def test_multi_channel_phase_space(self):
        """ Test the multichannel phase-space that is aligned along specific s- and t-channels."""
        
        # A specific sets of s- and t-channels for this test:

        ####################################################################
        # a) A simple unique massless photon s-channel from e+ e- > d d~ / z
        ####################################################################
        
        massless_photon_schannel_specifier = (
            # s-channels first:
            base_objects.VertexList([
                base_objects.Vertex({
                    'id': 15,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': -1,
                            'number': 4,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': 1,
                            'number': 3,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': 22,
                            'number': -1,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        })
                    ])
                }),
            ]),
            # t-channels then:
            base_objects.VertexList([
                base_objects.Vertex({
                    'id': 34,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 11,
                            'number': 1,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': 22,
                            'number': -1,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': -11,
                            'number': -2,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        })
                    ])
                }),
            ])
        )                          
        
        ####################################################################
        # a) A simple unique massive Z-boson s-channel from e+ e- > d d~ / a
        ####################################################################
        
        massive_zboson_schannel_specifier = (
            # s-channels first:
            base_objects.VertexList([
                base_objects.Vertex({
                    'id': 22,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': -1,
                            'number': 4,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': 1,
                            'number': 3,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': 23,
                            'number': -1,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                }),
            ]),
            # t-channels then:
            base_objects.VertexList([
                base_objects.Vertex({
                    'id': 40,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 11,
                            'number': 1,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': 23,
                            'number': -1,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': -11,
                            'number': -2,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                }),
            ]),
        )
      
        ###############################################################################
        # c) A complicated fully decayed VBF topology: 
        #    from: generate u c > h > u c e+ e- mu+ mu- $$ c u / a s d s~ d~ QCD=0 --LO
        ###############################################################################
        vbf_topology_s_and_t_channel_specifier = (
            # s-channels first:
            base_objects.VertexList([
                base_objects.Vertex({
                    'id': 41,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 13,
                            'number': 8,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': -13,
                            'number': 7,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': 23,
                            'number': -1,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                }),
                base_objects.Vertex({
                    'id': 40,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 11,
                            'number': 6,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': -11,
                            'number': 5,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': 23,
                            'number': -2,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                }),
                base_objects.Vertex({
                    'id': 13,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 23,
                            'number': -1,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': 23,
                            'number': -2,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': 25,
                            'number': -3,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                })
            ]),
            # t-channels then:
            base_objects.VertexList([
                base_objects.Vertex({
                    'id': 63,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': -2,
                            'number': 1,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': 2,
                            'number': 3,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': 23,
                            'number': -4,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                }),
                base_objects.Vertex({
                    'id': 13,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 23,
                            'number': -4,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None  
                        }),
                        base_objects.Leg({
                            'id': 25,
                            'number': -3,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': 23,
                            'number': -5,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                }),
                base_objects.Vertex({
                    'id': 64,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 23,
                            'number': -5,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None  
                        }),
                        base_objects.Leg({
                            'id': 4,
                            'number': 4,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': -4,
                            'number': -6,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                }),
            ]),
        )

        # The following topologies have been defined above
        #    massless_photon_schannel_specifier
        #    massive_zboson_schannel_specifier
        #    vbf_topology_s_and_t_channel_specifier
        
        # Example of a nice way to printout what these topologies are:        
        def print_topology(topology_to_print):
            return 's-channels:\n%s\nand t-channels:\n%s'%\
                    ( ', '.join('%s > %d'%(
                        ' '.join('%d'%leg['number'] for leg in vertex['legs'][:-1]),
                        vertex['legs'][-1]['number']) for vertex in topology_to_print[0]), 
                      ', '.join('%s > %d'%(
                        ' '.join('%d'%leg['number'] for leg in vertex['legs'][:-1]),
                        vertex['legs'][-1]['number']) for vertex in topology_to_print[1])
                    )

        E_cm  = 5000.0
        
        # Example on how to retrieve a numerical value for a model parameter
        MH = self.model.get('parameter_dict')['mdl_MH']
                
        # Now try it on the two e+ e- s-channel topologies. First the massless one.
        misc.sprint('='*100)
        misc.sprint('Now considering the following topology:\n'+print_topology(massless_photon_schannel_specifier))
        my_epem_PS_generator_massless = PS.MultiChannelPhasespace([0.]*2, [0.]*2, 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            model=self.model, topology=massless_photon_schannel_specifier)
        
        # Random variables set to None meaning that we ask the PS generator to generated them uniformly
        # themselves. Otherwise, it should be an instance of madgraph.integrator.integrands.DimensionList
        random_variables = None
        
        PS_point, wgt, xb_1, xb_2 =  my_epem_PS_generator_massless.get_PS_point(random_variables)
        misc.sprint('Generated the following PS point:')
        misc.sprint(str(PS_point))
        misc.sprint('jacobian = %.16e'%wgt)
        misc.sprint('xb_1 = %.16e'%xb_1)
        misc.sprint('xb_2 = %.16e'%xb_2)
        #
        # TODO: perform some (comparison) test of the PS point generated
        #

        # Now try it on the two e+ e- s-channel topologies. Then the massive one.
        misc.sprint('='*100)
        misc.sprint('Now considering the following topology:\n'+print_topology(massive_zboson_schannel_specifier))
        my_epem_PS_generator_massive = PS.MultiChannelPhasespace([0.]*2, [0.]*2, 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            model=self.model, topology=massive_zboson_schannel_specifier)

        PS_point, wgt, xb_1, xb_2 =  my_epem_PS_generator_massive.get_PS_point(random_variables)
        misc.sprint('Generated the following PS point:')
        misc.sprint(str(PS_point))
        misc.sprint('jacobian = %.16e'%wgt)
        misc.sprint('xb_1 = %.16e'%xb_1)
        misc.sprint('xb_2 = %.16e'%xb_2)
        #
        # TODO: perform some (comparison) test of the PS point generated
        #
        
        # Now try this pp multichannel PS generator on the hadronic VBF topology
        misc.sprint('='*100)
        misc.sprint('Now considering the following topology:\n'+print_topology(vbf_topology_s_and_t_channel_specifier))
        my_pp_PS_generator_VBF = PS.MultiChannelPhasespace([0.]*2, [0.]*6,
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(1,1), 
              model=self.model, topology=vbf_topology_s_and_t_channel_specifier)
        
        PS_point, wgt, xb_1, xb_2 =  my_pp_PS_generator_VBF.get_PS_point(random_variables)
        misc.sprint('Generated the following PS point:')
        misc.sprint(str(PS_point))
        misc.sprint('jacobian = %.16e'%wgt)
        misc.sprint('xb_1 = %.16e'%xb_1)
        misc.sprint('xb_2 = %.16e'%xb_2)
        #
        # TODO: perform some (comparison) test of the PS point generated
        #

#===============================================================================
# Final-final collinear mappings
#===============================================================================

class CollinearVariablesTest(unittest.TestCase):
    """Test class for variables describing internal collinear structure."""

    my_mapping = PS.ElementaryMappingCollinearFinal()
    n_children = 3

    def test_collinear_variables_away_from_limit(self):
        """Test determination of collinear variables and reverse mapping,
        for completely generic input values.
        """

        # Generate n_children+1 random massless vectors
        # (the light-cone direction is also random)
        my_PS_point = PS.LorentzVectorDict()
        for i in range(self.n_children+1):
            my_PS_point[i] = PS.LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            )
            my_PS_point[i].set_square(0)
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
        new_PS_point = PS.LorentzVectorDict()
        new_PS_point[self.n_children] = my_PS_point[self.n_children]
        self.my_mapping.set_collinear_variables(
            new_PS_point, self.n_children, range(self.n_children),
            total_momentum, variables
        )
        # Check the two phase-space points are equal
        for i in range(self.n_children):
            self.assertEqual(my_PS_point[i], new_PS_point[i])

    def test_collinear_variables_close_to_limit(self):
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
        pars = (math.pow(0.25, i) for i in range(8))
        # This is one way like any other to approach the limit
        for par in pars:
            new_directions = {
                i: (par*n + (1-par)*coll_direction)
                for (i, n) in directions.items()
            }
            my_PS_point = {
                i: PS.LorentzVector([0., ] + list(n)).set_square(0)
                for (i, n) in new_directions.items()
            }
            my_PS_point[self.n_children] = PS.LorentzVector(
                [0., ] + list(coll_direction)
            ).set_square(0)
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
            new_PS_point = PS.LorentzVectorDict()
            new_PS_point[self.n_children] = my_PS_point[self.n_children]
            self.my_mapping.set_collinear_variables(
                new_PS_point, self.n_children, range(self.n_children),
                total_momentum, variables
            )
            # Check the two phase-space points are equal
            for i in range(self.n_children):
                self.assertEqual(my_PS_point[i], new_PS_point[i])

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

    def test_collinear_map_invertible(self):
        """Test mapping and inverse."""

        # Generate n_collinear random massive vectors plus
        # n_recoilers (massive = True, False) random vectors
        my_PS_point = PS.LorentzVectorDict()
        n_collinear_so_far = 0
        for n_in_this_bunch in self.n_collinear:
            for i in range(
                n_collinear_so_far, n_collinear_so_far + n_in_this_bunch
            ):
                my_PS_point[i] = PS.LorentzVector(
                    [0., ] + [random.random() for _ in range(3)]
                )
                my_PS_point[i].set_square(0)
            if n_in_this_bunch == 1:
                my_PS_point[n_collinear_so_far].set_square(random.random())
            n_collinear_so_far += n_in_this_bunch
        for i in range(
            n_collinear_so_far, n_collinear_so_far + self.n_recoilers
        ):
            my_PS_point[i] = PS.LorentzVector(
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

#===============================================================================
# Soft mappings
#===============================================================================

class SoftVariablesTest(unittest.TestCase):
    """Test class for variables describing internal collinear structure."""

    # This is trivial if one stores the full momentum
    # Nevertheless checking and supporting more complicated scenarios

    my_mapping = PS.ElementaryMappingSoft()
    n_soft = 3

    def test_soft_variables_away_from_limit(self):
        """Test determination of soft variables and reverse mapping,
        for completely generic input values.
        """

        # Generate n_soft random massless vectors
        my_PS_point = PS.LorentzVectorDict()
        for i in range(self.n_soft):
            my_PS_point[i] = PS.LorentzVector(
                [0., ] + [random.random() for _ in range(3)]
            )
            my_PS_point[i].set_square(0)
        # Compute soft variables
        variables = dict()
        self.my_mapping.get_soft_variables(
            my_PS_point, range(self.n_soft), variables
        )
        # Compute new phase space point
        new_PS_point = PS.LorentzVectorDict()
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
        my_PS_point = PS.LorentzVectorDict()
        for i in range(self.n_soft):
            my_PS_point[i] = PS.LorentzVector(
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
            new_PS_point = PS.LorentzVectorDict()
            self.my_mapping.set_soft_variables(
                new_PS_point, range(self.n_soft), variables
            )
            for i in range(self.n_soft):
                self.assertEqual(old_PS_point[i], new_PS_point[i])

class SomogyietalSoftTest(unittest.TestCase):
    """Test class for MappingSomogyietalSoft."""

    # This specific mapping only for massless,
    # but 'massive' will be useful in the future so might as well keep it
    mapping = PS.MappingSomogyietalSoft()
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
            this_bunch_legs = (
                subtraction.SubtractionLeg(
                    i, 21, subtraction.SubtractionLeg.FINAL
                ) for i in this_bunch_numbers
            )
            self.structure += [subtraction.SoftStructure(this_bunch_legs), ]
            n_soft_so_far += n_in_this_bunch
            n_bunch += 1
        self.structure += [
                subtraction.SubtractionLeg(
                    i, 21, subtraction.SubtractionLeg.FINAL
                )
                for i in range(n_soft_so_far, n_tot)
            ]
        self.structure = subtraction.SingularStructure(self.structure)

    def test_soft_map_invertible(self):
        """Test mapping and inverse."""

        # Generate n_soft random massless vectors plus
        # n_recoilers (massive = True, False) random vectors
        my_PS_point = PS.LorentzVectorDict()
        n_soft_so_far = 0
        for n_in_this_bunch in self.n_soft:
            for i in range(
                n_soft_so_far, n_soft_so_far + n_in_this_bunch
            ):
                my_PS_point[i] = PS.LorentzVector(
                    [0., ] + [random.random() for _ in range(3)]
                )
                my_PS_point[i].set_square(0)
            n_soft_so_far += n_in_this_bunch
        for i in range(
            n_soft_so_far, n_soft_so_far + self.n_recoilers
        ):
            my_PS_point[i] = PS.LorentzVector(
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

#===============================================================================
# Test the phase-space walkers
#===============================================================================

# Define processes for phase-space walker tests
#===============================================================================

# Subtraction instance

my_subtraction = subtraction.IRSubtraction(
    simple_qcd.model,
    orders={'QCD': 1}
)

# H > u u~ d d~ H

H_to_uuxddxH_legs = base_objects.LegList([
    base_objects.Leg({'number': 1, 'id': 25, 'state': base_objects.Leg.INITIAL}),
    base_objects.Leg({'number': 2, 'id':  1, 'state': base_objects.Leg.FINAL}),
    base_objects.Leg({'number': 3, 'id': -1, 'state': base_objects.Leg.FINAL}),
    base_objects.Leg({'number': 4, 'id':  2, 'state': base_objects.Leg.FINAL}),
    base_objects.Leg({'number': 5, 'id': -2, 'state': base_objects.Leg.FINAL}),
    base_objects.Leg({'number': 6, 'id': 25, 'state': base_objects.Leg.FINAL}),
])

H_to_uuxddxH = base_objects.Process({
    'legs': H_to_uuxddxH_legs,
    'model': simple_qcd.model,
    'n_loops': 0
})

# H > q q~ g g H

H_to_qqxggH_legs = base_objects.LegList([
    base_objects.Leg({'number': 1, 'id': 25, 'state': base_objects.Leg.INITIAL}),
    base_objects.Leg({'number': 2, 'id':  1, 'state': base_objects.Leg.FINAL}),
    base_objects.Leg({'number': 3, 'id': -1, 'state': base_objects.Leg.FINAL}),
    base_objects.Leg({'number': 4, 'id': 21, 'state': base_objects.Leg.FINAL}),
    base_objects.Leg({'number': 5, 'id': 21, 'state': base_objects.Leg.FINAL}),
    base_objects.Leg({'number': 6, 'id': 25, 'state': base_objects.Leg.FINAL}),
])

H_to_qqxggH = base_objects.Process({
    'legs': H_to_qqxggH_legs,
    'model': simple_qcd.model,
    'n_loops': 0
})

def walk_invertible(test, walker, process, orders_for_filtering=None):
    """Test walk and inverse."""

    original_orders = my_subtraction.orders
    my_operators = my_subtraction.get_all_elementary_operators(process)
    my_combinations = my_subtraction.get_all_combinations(my_operators)
    if orders_for_filtering is not None:
        my_subtraction.orders = orders_for_filtering
    my_combinations = my_subtraction.filter_combinations(my_combinations)
    my_subtraction.orders = original_orders
    my_counterterms = [
        my_subtraction.get_counterterm(combination, process)
        for combination in my_combinations
    ]

    # For each counterterm
    for i in range(len(my_counterterms)):
        print "Testing counterterm", my_counterterms[i]
        # Generate random vectors
        my_PS_point = PS.LorentzVectorDict()
        legs_FS = (
            subtraction.SubtractionLeg(leg)
            for leg in process['legs']
            if leg['state'] == base_objects.Leg.FINAL
        )
        leg_numbers = (leg.n for leg in legs_FS)
        for j in leg_numbers:
            my_PS_point[j] = PS.LorentzVector(
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

    walker = PS.FlatCollinearWalker()

    def test_FlatCollinearWalker_invertible(self):

        walk_invertible(self, self.walker, H_to_uuxddxH, {'QCD': 2})

class SimpleNLOWalkerTest(unittest.TestCase):
    """Test class for FlatCollinearWalker."""

    walker = PS.SimpleNLOWalker()

    def test_SimpleNLOWalker_invertible(self):

        # Want to make the test a bit more serious,
        # so use combinations of NLO elementary operators up to 2 unresolved particles

        walk_invertible(self, self.walker, H_to_uuxddxH, {'QCD': 2})
        walk_invertible(self, self.walker, H_to_qqxggH, {'QCD': 2})
