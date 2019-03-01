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
"""Unit test library for the various phase-space generation/handling features."""

import madgraph.integrator.phase_space_generators as PS
import madgraph.core.subtraction as subtraction
import madgraph.core.base_objects as base_objects
import madgraph.integrator.vectors as vectors
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

#=========================================================================================
# Shorthands for initial and final state
#=========================================================================================

INITIAL = base_objects.Leg.INITIAL
FINAL   = base_objects.Leg.FINAL

assert subtraction.SubtractionLeg.INITIAL == INITIAL
assert subtraction.SubtractionLeg.FINAL   == FINAL

#=========================================================================================
# Test Vectors
#=========================================================================================

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
        self.assertAlmostEqual(v1.square(), 0, 10)
        self.assertAlmostEqual(v2.square(), 0, 10)
        self.assertAlmostEqual(v1.dot(v2), 2, 10)

        # Test the rotoboost
        #   Definition
        v1 = PS.LorentzVector([0, ] + [random.random() for _ in range(3)])
        v2 = PS.LorentzVector([0, ] + [random.random() for _ in range(3)])
        m = random.random()
        v1.set_square(m**2)
        v2.set_square(m**2)
        for x, y in zip(v2.rotoboost(v2, v1), v1):
            self.assertAlmostEqual(x,y,10)
        #   Inverse
        v3 = PS.LorentzVector([0, ] + [random.random() for _ in range(3)])
        v3.set_square(m**2)
        v4 = PS.LorentzVector([random.random() for _ in range(4)])
        v4_old = v4.get_copy()
        v4.rotoboost(v1, v3)
        v4.rotoboost(v3, v1)
        for x, y in zip(v4, v4_old):
            self.assertAlmostEqual(x,y,10)

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

        # A specific sets of s- and t-channels for this test:

        ####################################################################
        # a) A simple unique massless photon s-channel from u u~ > d d~ / z w+ QCD=0
        ####################################################################
        
        self.massless_photon_schannel_specifier = (
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
                            'id': 2,
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
                            'id': -2,
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
        # a) A simple unique massive Z-boson s-channel from u u~ > d d~ / a w+ QCD = 0
        ####################################################################
        
        self.massive_zboson_schannel_specifier = (
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
                            'id': 2,
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
                            'id': -2,
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
        
        ####################################################################
        # b) A simple unique massive Z-boson t-channel from generate u d > u d /a w+ w- QCD=0 --LO
        ####################################################################
        
        self.massive_zboson_tchannel_specifier = (
            # s-channels first:
            base_objects.VertexList([]),
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
                            'number': -1,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                }),
                base_objects.Vertex({
                    'id': 22,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 23,
                            'number': -1,
                            'state': False,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': 1,
                            'number': 4,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': -1,
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
        
        ####################################################################
        # b) A simple unique massless photon t-channel from u d > u d / z
        ####################################################################
        
        self.massless_photon_tchannel_specifier = (
            # s-channels first:
            base_objects.VertexList([]),
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
                            'id': 22,
                            'number': -1,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                    ])
                }),
                base_objects.Vertex({
                    'id': 22,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 22,
                            'number': -1,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': 1,
                            'number': 4,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': -1,
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
        self.vbf_topology_s_and_t_channel_specifier = (
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


        ###############################################################################
        # d) A complicated fully decayed VBF topology: 
        #    from: generate e- e+ > h > e+ e- mu+ mu- ta+ ta- $$ e+ e- \ a QCD=0 --diagram_filter --LO
        ###############################################################################
        # where diagram filter removes the first three diagrams
        # import model sm-dario
        self.vbf_topology_s_and_t_channel_specifier2 = (
            # s-channels first:
            base_objects.VertexList([
                base_objects.Vertex({
                    'id': 42,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 15,
                            'number': 8,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None  
                        }),
                        base_objects.Leg({
                            'id': -15,
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
                    'id': 41,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 13,
                            'number': 6,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': -13,
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
                    'id': 40,
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                             'id': -11,
                             'number': 1,
                             'state': False,
                             'from_group': True,
                             'loop_line': False,
                             'onshell': None   
                        }),
                        base_objects.Leg({
                            'id': 11,
                            'number': 4,
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
                    'id': 40,
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
                            'id': -11,
                            'number': 3,
                            'state': True,
                            'from_group': True,
                            'loop_line': False,
                            'onshell': None
                        }),
                        base_objects.Leg({
                            'id': 11,
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
        
        
    def test_flat_invertible_phase_space(self):
        """ Tests the flat invertible phase-space."""
        
        E_cm  = 5000.0
    
        # Try to run the above for a 2->8.
        my_PS_generator = PS.FlatInvertiblePhasespace(
            [0.]*2, [100. + 10.*i for i in range(8)],beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0))
        # Try to run the above for a 2->1.    
        #    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [5000.0])
        
        random_variables = [random.random() for _ in range(my_PS_generator.nDimPhaseSpace())]
    
#        import time
#        start = time.time()
        n_loops = 1
        for _ in range(n_loops):
            momenta, wgt = my_PS_generator.generateKinematics(E_cm, random_variables)
#        end = time.time()
#        misc.sprint('Time per call',(end-start)/float(n_loops))
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
            
    def test_matrix_element_integration(self):
        """ Test the multicannel phase-space integration over simple s- and t- channel matrix elements."""
        
        pass
        
    def test_path_generator(self):
        
        import madgraph.interface.madgraph_interface as madgraph_interface
        import madgraph.iolibs.drawing_eps as draw
        import madgraph.core.base_objects as base_objects
        import madgraph.core.drawing as draw_lib
        
        my_topology = self.vbf_topology_s_and_t_channel_specifier
        
        max_leg_nr = 3
        for channel in my_topology:
            for vertex in channel:
                for leg in vertex.get('legs'):
                    leg_nr = leg.get('number')
                    if leg_nr > max_leg_nr:
                        max_leg_nr = leg_nr
        
        nr_final = max_leg_nr - 2
        E_cm = 5000
        my_PS_generator = my_epem_PS_generator = PS.SingleChannelPhasespace([0.]*2, [0.]*nr_final, beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            model=self.model, topology=my_topology)
        
        my_random_path = my_PS_generator.get_random_path()
        ##misc.sprint('\nThis is a random path ', my_random_path)
        """
        #misc.sprint('\n',my_topology)
        my_command = madgraph_interface.MadGraphCmd()
        #my_command.draw(line='diagrams')
        #[{'vertices': my_topology, 'orders': {'WEIGHTED': 8, 'QED': 4, 'QCD': 0}}]
        #misc.sprint(my_topology)
        options = {'horizontal': False, 'add_gap': 0, 'external': 0, 'max_size': 1.5, 'contract_non_propagating': True}
        options = draw_lib.DrawOption(options)
        my_diag = base_objects.VertexList([vertex for vertex in my_topology[0]]+[vertex for vertex in my_topology[1]])
        #for vertex in my_topology[1]:
        #my_diag = my_diag.append(my_topology[1][0])
        #misc.sprint(my_diag)
        for vertex in my_diag:
            leg1 = vertex.get('legs')[0].copy()
            leg2 = vertex.get('legs')[1].copy()
            leg3 = vertex.get('legs')[2].copy()
            if leg1['number']>leg2['number']:
                vertex['legs'][0] = leg2
                vertex['legs'][1] = leg1
            nr = leg3['number']
            leg1 = vertex.get('legs')[0].copy()
            leg2 = vertex.get('legs')[1].copy()
            leg3['number'] = leg1['number']
            for vert in my_diag:
                for leg in vert['legs']:
                    if leg['number'] == nr:
                        leg['number'] = leg3['number']
        #my_diag[-1]['legs'][-1]['number'] == 2
        diag = base_objects.Diagram({'vertices': my_diag, 'orders': {'WEIGHTED': 12, 'QED': 6, 'QCD': 0}})
        misc.sprint(diag)
        
        #diag = base_objects.Diagram({'vertices': base_objects.VertexList()})
        diags = base_objects.DiagramList([diag])
        filename='test_path_generator.eps'
        plot = draw.MultiEpsDiagramDrawer(diagramlist=diags,filename=filename,model=self.model,amplitude=None,legend=None,diagram_type=None)
        plot.draw(opt=options)
        #my_command
        my_command.exec_cmd('open %s' % filename)
        """    
    
    def test_single_channel_wgt_reconstruction(self):
        """ Test single channel weight reconstruction. """  
        
        E_cm = 5000
        nr_final = 6
        
        SCPS = PS.SingleChannelPhasespace([0.]*2, [0.]*nr_final, beam_Es =(E_cm/2.,E_cm/2.), beam_types=(1,1),
            model=self.model, topology=self.vbf_topology_s_and_t_channel_specifier)

        path = SCPS.get_random_path()
        random_variables = SCPS.dimensions.random_sample()  
        PS_point, wgt, xb_1, xb_2 = SCPS.get_PS_point(random_variables,path=path)
        reconstructed_variables, reconstructed_wgt = SCPS.get_PS_point(PS_point,path=path)
        
        ##misc.sprint('\n Random variables       : ',random_variables,'\n Reconstructed variables: ',reconstructed_variables)
        ##misc.sprint('\n Direct weight       : %.6e'%wgt + '\n Reconstructed weight: %.6e'%reconstructed_wgt)
        ##misc.sprint('\n Ratio direct/reconstructed : %.6e'%(wgt/reconstructed_wgt))
        # The inversion is for now coded up only for the invariants. The angles are for now not reconstructed and are None.
        differences = [abs(reconstructed_variables[i]-random_variables[i]) 
                                        for i in range(len(reconstructed_variables)) if reconstructed_variables[i] is not None]
        self.assertLess(max(differences[i]/max(random_variables[i],1.0e-10) for i in range(len(differences))), 1.0e-10)
        self.assertLess(abs(wgt-reconstructed_wgt)/max(abs(wgt),1.0e-10), 1.0e-10)

    
    def test_phase_space_volume(self):
        """ Test the Singlechannel phase-space that is aligned along specific s- and t-channels."""
        
        import madgraph.integrator.vegas3_integrator as vegas3
        import madgraph.integrator.integrands as integrands

        class IntegrandForTest(integrands.VirtualIntegrand):
            """An integrand for this phase-space volume test."""
            def __init__(self, phase_space_generator):
                super(IntegrandForTest, self).__init__(phase_space_generator.get_dimensions())
                self.phase_space_generator = phase_space_generator
                self.counter               = 0
                #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
                #    self.my_random_path = self.phase_space_generator.generate_random_path()
            
            def __call__(self, continuous_inputs, discrete_inputs, **opts):
                #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
                #    PS_point, wgt, x1, x2 = self.phase_space_generator.get_PS_point(continuous_inputs,self.my_random_path)
                #else:
                self.counter += 1
                PS_point, wgt, x1, x2 = self.phase_space_generator.get_PS_point(continuous_inputs)
                return wgt
        
        def analytical_phase_space_vol(E_cm,n):
            return math.pow((math.pi/2.0),n-1)*(math.pow((E_cm**2),n-2)/(math.factorial(n-1)*math.factorial(n-2)))
        
        verbose = False 

        E_cm = 5000
        nr_final = 6
        
        #Analytical
        if verbose: misc.sprint('\n%d-body phase space analytical: %.4e'%(nr_final,analytical_phase_space_vol(E_cm,nr_final)))
        analytical_PS_volume = analytical_phase_space_vol(E_cm,nr_final)

        #SCPS
        my_epem_PS_generator = PS.SingleChannelPhasespace([0.]*2, [0.]*nr_final, beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            model=self.model, topology=self.vbf_topology_s_and_t_channel_specifier)
        my_integrand = IntegrandForTest(my_epem_PS_generator)
        my_integrator = vegas3.Vegas3Integrator(my_integrand, n_points_survey=400, n_points_refine=400, accuracy_target=None)
        # Finally integrate
        if verbose: misc.sprint('\nSCPS %d-body phase space '%nr_final + 'SCPS: Final result: %.4e +/- %.2e'%my_integrator.integrate())
        SCPS_volume = my_integrator.integrate()
        #misc.sprint(abs(SCPS_volume[0]-analytical_PS_volume)/abs(SCPS_volume[1]), 5.0)
        self.assertTrue(abs(SCPS_volume[0]-analytical_PS_volume)/abs(SCPS_volume[1])<5.0)

        #FLATPS
        my_PS_generator = PS.FlatInvertiblePhasespace([0.]*2, [0.]*nr_final, beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0))
        my_integrand = IntegrandForTest(my_PS_generator)
        my_integrator = vegas3.Vegas3Integrator(my_integrand, n_points_survey=100, n_points_refine=100, accuracy_target=None)
        # Finally integrate
        if verbose: misc.sprint('\n FLATPS %d-body phase space '%nr_final + 'FLATPS: Final result: %.4e +/- %.2e'%my_integrator.integrate())
        FLATPS_volume = my_integrator.integrate()
        self.assertTrue(abs(FLATPS_volume[0]-analytical_PS_volume)/abs(FLATPS_volume[1])<5.0)

        return 0

    @staticmethod
    def compare_PS_point(a,b,threshold=1.0e-10):
        for i_vec, (a_vec, b_vec) in enumerate(zip(a,b)):
            if max(abs(a_vec_el-b_vec_el)/max(abs(a_vec_el), 1.0e-10) for a_vec_el, b_vec_el in zip(a_vec, b_vec)) > threshold:
                misc.sprint('Lorentz vector #%d differ: %s vs %s'%(i_vec+1, str(a_vec), str(b_vec)))
                return False
        return True

    def test_single_channel_phase_space(self):
        """ Test the single channel phase-space that is aligned along specific s- and t-channels."""
        
        verbose = False 

        # The following topologies have been defined above
        #    massless_photon_schannel_specifier
        #    massive_zboson_schannel_specifier
        #    vbf_topology_s_and_t_channel_specifier
        
        # Example of a nice way to printout what these topologies are:        
        def print_topology(topology_to_print):
            if  topology_to_print[0] == None:
                return 'no s-channels \nand t-channels:\n%s'%\
                    ( ', '.join('%s > %d'%(
                        ' '.join('%d'%leg['number'] for leg in vertex['legs'][:-1]),
                        vertex['legs'][-1]['number']) for vertex in topology_to_print[1])
                    )
            else:
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
        
        """
        # first the massless photon t-channel
        misc.sprint('='*100)
        misc.sprint('Now considering the following topology:\n'+print_topology(self.massless_photon_tchannel_specifier))
        my_epem_PS_generator_massive = PS.SingleChannelPhasespace([0.]*2, [0.]*2, 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(1,1),
            model=self.model, topology=self.massless_photon_tchannel_specifier)
        
        # Random variables set to None meaning that we ask the PS generator to generated them uniformly
        # themselves. Otherwise, it should be an instance of madgraph.integrator.integrands.DimensionList
        random_variables = None
        
        out_data = open('Myoutput_massless.dat','w')
        for n_point in range(10**4):
            PS_point, wgt, xb_1, xb_2 =  my_epem_PS_generator_massive.get_PS_point(random_variables)
            out_data.write('%.16f, '%(PS_point[1]-PS_point[3]).square()+'%.16f\n'%(xb_1*xb_2))
        out_data.close()
        
        misc.sprint('Generated the following PS point:')
        misc.sprint(str(PS_point))
        misc.sprint('jacobian = %.16e'%wgt)
        misc.sprint('xb_1 = %.16e'%xb_1)
        misc.sprint('xb_2 = %.16e'%xb_2)
        #
        # TODO: perform some (comparison) test of the PS point generated
        #
        """
         
        # Now try this pp Singlechannel PS generator on the hadronic VBF topology
        if verbose: misc.sprint('='*100)
        if verbose: misc.sprint('Now considering the following topology:\n'+print_topology(self.vbf_topology_s_and_t_channel_specifier2))
        my_pp_PS_generator_VBF = PS.SingleChannelPhasespace([0.]*2, [0.]*6,
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0), 
              model=self.model, topology=self.vbf_topology_s_and_t_channel_specifier2)

        #my_pp_PS_generator_VBF.path = my_pp_PS_generator_VBF.get_random_path()
        my_pp_PS_generator_VBF.path = [[1, 0, 2], [(0, 1)]]
        if verbose: misc.sprint('parametrisation path=',my_pp_PS_generator_VBF.path)

        
        # Random variables set to None meaning that we ask the PS generator to generated them uniformly
        random_variables = [ 0.24102803, 0.57037533, 0.48274477, 0.44783153, 0.44533712, 0.69884539, 
                             0.85848863, 0.33310392, 0.0262083,  0.89275099, 0.88319001, 0.8931308,
                             0.91795734, 0.29064606]
        # random_variables = my_pp_PS_generator_VBF.dimensions.random_sample()
        if verbose: misc.sprint('Input random variables:',random_variables)
        PS_point, wgt, xb_1, xb_2 =  my_pp_PS_generator_VBF.get_PS_point(random_variables)
        if verbose: misc.sprint('Generated the following PS point:')
        if verbose: misc.sprint(str(PS_point))
        if verbose: misc.sprint('jacobian = %.16e'%wgt)
        if verbose: misc.sprint('xb_1 = %.16e'%xb_1[0])
        if verbose: misc.sprint('xb_2 = %.16e'%xb_2[0]) 
        self.assertAlmostEqual(wgt, 4.9941144626860431e+18, 10)
        self.assertAlmostEqual(xb_1[0], 1.0, 10)
        self.assertAlmostEqual(xb_2[0], 1.0, 10)        
        self.assertTrue(self.compare_PS_point(PS_point,
            vectors.LorentzVectorList([
                vectors.LorentzVector([2.5000000000000000e+03,   0.0000000000000000e+00,    0.0000000000000000e+00,    2.5000000000000000e+03]),
                vectors.LorentzVector([2.5000000000000000e+03,   0.0000000000000000e+00,    0.0000000000000000e+00,   -2.5000000000000000e+03]),
                vectors.LorentzVector([1.3777142427916745e+03,   4.3262266420292940e+02,    1.2993259345667730e+03,   -1.5061965665602452e+02]),
                vectors.LorentzVector([2.1710049623024338e+03,  -1.5089883488272969e+03,   -1.3423949792402009e+03,   -7.9636199630316389e+02]),
                vectors.LorentzVector([5.6619767333718691e+01,   4.2460556304040125e+01,    2.1934323050718518e+01,    3.0360907160073381e+01]),
                vectors.LorentzVector([8.4753536614633481e+02,   5.9434471364340163e+02,    1.6833840509615015e+01,    6.0397614195151914e+02]),
                vectors.LorentzVector([1.7222182422139423e+02,   1.5636728285898403e+02,   -2.2821767218714310e+01,    6.8474787553216402e+01]),
                vectors.LorentzVector([3.7490383720444157e+02,   2.8319313181793984e+02,    2.7122648331808783e+01,    2.4416981629437822e+02]),
            ])
        ))
        
        # the massive Z-boson t-channel
        if verbose: misc.sprint('='*100)
        if verbose: misc.sprint('Now considering the following topology:\n'+print_topology(self.massive_zboson_tchannel_specifier))
        my_epem_PS_generator_massive = PS.SingleChannelPhasespace([0.]*2, [0.]*2, 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(1,1),
            model=self.model, topology=self.massive_zboson_tchannel_specifier)

        my_epem_PS_generator_massive.path = [[], []]
        if verbose: misc.sprint('parametrisation path=',my_epem_PS_generator_massive.path)        
        
        # Random variables set to None meaning that we ask the PS generator to generated them uniformly
        # themselves. Otherwise, it should be an instance of madgraph.integrator.integrands.DimensionList
        #random_variables = my_epem_PS_generator_massive.dimensions.random_sample()
        random_variables = [0.73239213, 0.48891956, 0.30386151, 0.96030843]   
        if verbose: misc.sprint('Input random variables:',random_variables) 
        
        PS_point, wgt, xb_1, xb_2 =  my_epem_PS_generator_massive.get_PS_point(random_variables)
        
        if verbose: misc.sprint('Generated the following PS point:')
        if verbose: misc.sprint(str(PS_point))
        if verbose: misc.sprint('jacobian = %.16e'%wgt)
        if verbose: misc.sprint('xb_1 = %.16e'%xb_1[0])
        if verbose: misc.sprint('xb_2 = %.16e'%xb_2[0])
        self.assertAlmostEqual(wgt, 1.1239946709787594e+00, 10)
        self.assertAlmostEqual(xb_1[0], 8.2572891098195511e-01, 10)
        self.assertAlmostEqual(xb_2[0], 5.9210665139700069e-01, 10)        
        self.assertTrue(self.compare_PS_point(PS_point,
            vectors.LorentzVectorList([
                vectors.LorentzVector([1.7480696146807511e+03,    0.0000000000000000e+00,    0.0000000000000000e+00,    1.7480696146807511e+03]),
                vectors.LorentzVector([1.7480696146807511e+03,    0.0000000000000000e+00,    0.0000000000000000e+00,   -1.7480696146807511e+03]),
                vectors.LorentzVector([1.7480696146807511e+03,    1.5582119406049219e+03,   -3.9686365896323360e+02,    6.8572746927672858e+02]),
                vectors.LorentzVector([1.7480696146807511e+03,   -1.5582119406049219e+03,    3.9686365896323360e+02,   -6.8572746927672858e+02]),
            ])
        ))
                
        # Now try it on the two e+ e- s-channel topologies. First the massless one.
        if verbose: misc.sprint('='*100)
        if verbose: misc.sprint('Now considering the following topology:\n'+print_topology(self.massless_photon_tchannel_specifier))
        my_epem_PS_generator_massless = PS.SingleChannelPhasespace([0.]*2, [0.]*2, 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(1,1),
            model=self.model, topology=self.massless_photon_tchannel_specifier)
        
        my_epem_PS_generator_massless.path = [[], []]
        if verbose: misc.sprint('parametrisation path=',my_epem_PS_generator_massless.path)        
        
        # Random variables set to None meaning that we ask the PS generator to generated them uniformly
        # themselves. Otherwise, it should be an instance of madgraph.integrator.integrands.DimensionList
        #random_variables = my_epem_PS_generator_massless.dimensions.random_sample()
        random_variables = [0.29901866, 0.60439549, 0.08707295, 0.43305739]   
        if verbose: misc.sprint('Input random variables:',random_variables) 

        PS_point, wgt, xb_1, xb_2 =  my_epem_PS_generator_massless.get_PS_point(random_variables)
       
        if verbose: misc.sprint('Generated the following PS point:')
        if verbose: misc.sprint(str(PS_point))
        if verbose: misc.sprint('jacobian = %.16e'%wgt)
        if verbose: misc.sprint('xb_1 = %.16e'%xb_1[0])
        if verbose: misc.sprint('xb_2 = %.16e'%xb_2[0])
       
        self.assertAlmostEqual(wgt, 9.4479392122324497e-07, 10)
        self.assertAlmostEqual(xb_1[0], 7.0260342125405872e-01, 10)
        self.assertAlmostEqual(xb_2[0], 8.6022283345192141e-01, 10)        
        self.assertTrue(self.compare_PS_point(PS_point,
            vectors.LorentzVectorList([
                vectors.LorentzVector([1.9435719465461336e+03,    0.0000000000000000e+00,    0.0000000000000000e+00,    1.9435719465461336e+03]),
                vectors.LorentzVector([1.9435719465461336e+03,    0.0000000000000000e+00,    0.0000000000000000e+00,   -1.9435719465461336e+03]),
                vectors.LorentzVector([1.9435719465461339e+03,   -9.4651597058939463e-01,    4.2338396204149914e-01,    1.9435716699557713e+03]),
                vectors.LorentzVector([1.9435719465461339e+03,    9.4651597058939463e-01,   -4.2338396204149914e-01,   -1.9435716699557713e+03]),
            ])
        ))

        # Now try it on the two e+ e- s-channel topologies. Then the massive one.
        if verbose: misc.sprint('='*100)
        if verbose: misc.sprint('Now considering the following topology:\n'+print_topology(self.massive_zboson_schannel_specifier))
        my_epem_PS_generator_massive = PS.SingleChannelPhasespace([0.]*2, [0.]*2, 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(1,1),
            model=self.model, topology=self.massive_zboson_schannel_specifier)
       
        my_epem_PS_generator_massive.path = [[], []]
        if verbose: misc.sprint('parametrisation path=',my_epem_PS_generator_massive.path)        
        
        # Random variables set to None meaning that we ask the PS generator to generated them uniformly
        # themselves. Otherwise, it should be an instance of madgraph.integrator.integrands.DimensionList
        #random_variables = my_epem_PS_generator_massive.dimensions.random_sample()
        random_variables = [0.12194472, 0.55495504, 0.74388797, 0.01441303]   
        if verbose: misc.sprint('Input random variables:',random_variables)

        PS_point, wgt, xb_1, xb_2 =  my_epem_PS_generator_massive.get_PS_point(random_variables)

        if verbose: misc.sprint('Generated the following PS point:')
        if verbose: misc.sprint(str(PS_point))
        if verbose: misc.sprint('jacobian = %.16e'%wgt)
        if verbose: misc.sprint('xb_1 = %.16e'%xb_1[0])
        if verbose: misc.sprint('xb_2 = %.16e'%xb_2[0])

        self.assertAlmostEqual(wgt, 3.6086874339847757e-04, 10)
        self.assertAlmostEqual(xb_1[0], 8.8709098163260083e-04, 10)
        self.assertAlmostEqual(xb_2[0], 3.7681867561422167e-01, 10)        
        self.assertTrue(self.compare_PS_point(PS_point,
            vectors.LorentzVectorList([
                vectors.LorentzVector([4.5707798079766739e+01,    0.0000000000000000e+00,    0.0000000000000000e+00,    4.5707798079766739e+01]),
                vectors.LorentzVector([4.5707798079766739e+01,    0.0000000000000000e+00,    0.0000000000000000e+00,   -4.5707798079766739e+01]),
                vectors.LorentzVector([4.5707798079766739e+01,   -3.9737978862046255e+01,   -3.6085309654448117e+00,   -2.2295164173688413e+01]),
                vectors.LorentzVector([4.5707798079766739e+01,    3.9737978862046255e+01,    3.6085309654448117e+00,    2.2295164173688413e+01]),
            ])
        ))

    def test_multi_channel_phase_space(self):  
       
        verbose = False 

        E_cm  = 500
        
        topologies = [self.massive_zboson_schannel_specifier,self.massive_zboson_tchannel_specifier]  
        
        if verbose: misc.sprint('='*100)
        if verbose: misc.sprint('Consider topologies: massive z boson, s-&t-channel')

        my_epem_PS_generator_massive = PS.MultiChannelPhasespace([0.]*2, [0.]*2, 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            model=self.model, topologies=topologies)
         
        # Random variables set to None meaning that we ask the PS generator to generated them uniformly
        # themselves. Otherwise, it should be an instance of madgraph.integrator.integrands.DimensionList
        #random_variables = my_epem_PS_generator_massive.dimensions.random_sample()
        random_variables = [0.19995697, 0.70414753]   
        if verbose: misc.sprint('Input random variables:',random_variables)

        alphas = [0.25,0.75]
        PS_point, wgt, xb_1, xb_2 =  my_epem_PS_generator_massive.get_PS_point(random_variables,adaptive_wgts=alphas,channel_nr=0)
        
        if verbose: misc.sprint('Generated the following PS point for channel #1')
        if verbose: misc.sprint(str(PS_point))
        if verbose: misc.sprint('jacobian = %.16e'%wgt)
        if verbose: misc.sprint('xb_1 = %.16e'%xb_1[0])
        if verbose: misc.sprint('xb_2 = %.16e'%xb_2[0])
        self.assertAlmostEqual(wgt, 1.5707963267948963e+00, 10)
        self.assertAlmostEqual(xb_1[0], 1.0000000000000000e+00, 10)
        self.assertAlmostEqual(xb_2[0], 1.0000000000000000e+00, 10)        
        self.assertTrue(self.compare_PS_point(PS_point,
            vectors.LorentzVectorList([
                vectors.LorentzVector([2.5000000000000000e+02,    0.0000000000000000e+00,    0.0000000000000000e+00,    2.5000000000000000e+02]),
                vectors.LorentzVector([2.5000000000000000e+02,    0.0000000000000000e+00,    0.0000000000000000e+00,   -2.5000000000000000e+02]),
                vectors.LorentzVector([2.5000000000000000e+02,    5.6821540761378650e+01,    1.9174164269299396e+02,    1.5002151499999999e+02]),
                vectors.LorentzVector([2.5000000000000000e+02,   -5.6821540761378650e+01,   -1.9174164269299396e+02,   -1.5002151499999999e+02]),
            ])
        ))

        PS_point, wgt, xb_1, xb_2 =  my_epem_PS_generator_massive.get_PS_point(random_variables,adaptive_wgts=alphas,channel_nr=1)

        if verbose: misc.sprint('Generated the following PS point for channel #2')
        if verbose: misc.sprint(str(PS_point))
        if verbose: misc.sprint('jacobian = %.16e'%wgt)
        if verbose: misc.sprint('xb_1 = %.16e'%xb_1[0])
        if verbose: misc.sprint('xb_2 = %.16e'%xb_2[0])
        self.assertAlmostEqual(wgt, 1.5707963267948963e+00, 10)
        self.assertAlmostEqual(xb_1[0], 1.0000000000000000e+00, 10)
        self.assertAlmostEqual(xb_2[0], 1.0000000000000000e+00, 10)        
        self.assertTrue(self.compare_PS_point(PS_point,
            vectors.LorentzVectorList([
                vectors.LorentzVector([2.5000000000000000e+02,    0.0000000000000000e+00,    0.0000000000000000e+00,    2.5000000000000000e+02]),
                vectors.LorentzVector([2.5000000000000000e+02,    0.0000000000000000e+00,    0.0000000000000000e+00,   -2.5000000000000000e+02]),
                vectors.LorentzVector([2.5000000000000000e+02,   -5.6821540761378650e+01,   -1.9174164269299396e+02,    1.5002151499999999e+02]),
                vectors.LorentzVector([2.5000000000000000e+02,    5.6821540761378650e+01,    1.9174164269299396e+02,   -1.5002151499999999e+02]),
            ])
        ))

        return 0
