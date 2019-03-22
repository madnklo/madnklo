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
# Shorthands for initial and final state
#=========================================================================================

INITIAL = base_objects.Leg.INITIAL
FINAL   = base_objects.Leg.FINAL

assert subtraction.SubtractionLeg.INITIAL == INITIAL
assert subtraction.SubtractionLeg.FINAL   == FINAL

#===============================================================================
# Test the phase-space generators
#=========================================================================================

class PhaseSpaceGeneratorsTest(unittest.TestCase):
    """ Test various phase-space generators."""

    verbosity = 0

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
        if self.verbosity > 1:
            print "\n ========================="
            print " ||    PS generation    ||"
            print " ========================="
            print "\nRandom variables :\n",random_variables
            print "\n%s\n"%momenta.__str__(n_initial=my_PS_generator.n_initial)
            print "Phase-space weight : %.16e\n" % wgt
    
        variables_reconstructed, wgt_reconstructed = \
            my_PS_generator.invertKinematics(E_cm, momenta)
    
        if self.verbosity > 1:
            print "\n ========================="
            print " || Kinematic inversion ||"
            print " ========================="
            print "\nReconstructed random variables :\n", variables_reconstructed
            print "\nReconstructed weight : %.16e\n" % wgt_reconstructed
        differences = [abs(variables_reconstructed[i]-random_variables[i])
                       for i in range(len(variables_reconstructed))]

        self.assertLess(max(differences[i]/random_variables[i] for i in range(len(differences))), 1.0e-10)
        self.assertLess(abs(wgt-wgt_reconstructed)/abs(wgt), 1.0e-10)
        
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
#        misc.sprint('='*100)
#        misc.sprint('Now considering the following topology:\n'+print_topology(massless_photon_schannel_specifier))
        my_epem_PS_generator_massless = PS.MultiChannelPhasespace([0.]*2, [0.]*2, 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            model=self.model, topology=massless_photon_schannel_specifier)
        
        # Random variables set to None meaning that we ask the PS generator to generated them uniformly
        # themselves. Otherwise, it should be an instance of madgraph.integrator.integrands.DimensionList
        random_variables = None
        
        
        PS_point, wgt, x1s, x2s =  my_epem_PS_generator_massless.get_PS_point(random_variables)

 #       misc.sprint('Generated the following PS point:')
 #       misc.sprint(str(PS_point))
 #       misc.sprint('jacobian = %.16e'%wgt)
 #       misc.sprint('xb_1 = %.16e'%x1s[0])
 #       misc.sprint('xb_2 = %.16e'%x2s[0])
        #
        # TODO: perform some (comparison) test of the PS point generated
        #

        # Now try it on the two e+ e- s-channel topologies. Then the massive one.
#       misc.sprint('='*100)
#        misc.sprint('Now considering the following topology:\n'+print_topology(massive_zboson_schannel_specifier))
        my_epem_PS_generator_massive = PS.MultiChannelPhasespace([0.]*2, [0.]*2, 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            model=self.model, topology=massive_zboson_schannel_specifier)

        PS_point, wgt, x1s, x2s =  my_epem_PS_generator_massive.get_PS_point(random_variables)
#        misc.sprint('Generated the following PS point:')
#        misc.sprint(str(PS_point))
#        misc.sprint('jacobian = %.16e'%wgt)
#        misc.sprint('xb_1 = %.16e'%x1s[0])
#        misc.sprint('xb_2 = %.16e'%x2s[0])
        #
        # TODO: perform some (comparison) test of the PS point generated
        #
        
        # Now try this pp multichannel PS generator on the hadronic VBF topology
#        misc.sprint('='*100)
#        misc.sprint('Now considering the following topology:\n'+print_topology(vbf_topology_s_and_t_channel_specifier))
        my_pp_PS_generator_VBF = PS.MultiChannelPhasespace([0.]*2, [0.]*6,
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(1,1), 
              model=self.model, topology=vbf_topology_s_and_t_channel_specifier)
        
        PS_point, wgt, x1s, x2s =  my_pp_PS_generator_VBF.get_PS_point(random_variables)
#        misc.sprint('Generated the following PS point:')
#        misc.sprint(str(PS_point))
#        misc.sprint('jacobian = %.16e'%wgt)
#        misc.sprint('xb_1 = %.16e'%x1s[0])
#        misc.sprint('xb_2 = %.16e'%x2s[0])
        #
        # TODO: perform some (comparison) test of the PS point generated
        #
