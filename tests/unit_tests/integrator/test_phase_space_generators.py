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
"""Unit test library for phase-space generation."""

import madgraph.integrator.phase_space_generators as PS
import madgraph.core.base_objects as base_objects
import madgraph.various.misc as misc
import models.import_ufo as import_ufo
from madgraph import MG5DIR

import random
import os

import tests.unit_tests as unittest

pjoin = os.path.join

#=========================================================================================
# Test the phase-space generators
#=========================================================================================

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
            [0.]*2, [100. + 10.*i for i in range(8)],
            beam_Es =(E_cm/2., E_cm/2.), beam_types=(0, 0) )
        # Try to run the above for a 2->1.    
        #    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [5000.0])
        
        random_variables = [random.random() for _ in range(my_PS_generator.nDimPhaseSpace())]

#        import time
#        start = time.time()
#        n_loops = 1
#        for _ in range(n_loops):
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
