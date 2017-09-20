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

import copy
import fractions
import os

import madgraph.core.base_objects as base_objects
import madgraph.loop.loop_base_objects as loop_base_objects
import madgraph.core.subtraction as sub
import madgraph.core.color_algebra as color
import madgraph.various.misc as misc
import madgraph.core.contributions as contributions
import madgraph.interface.master_interface as cmd
import madgraph.various.misc as misc
import madgraph.iolibs.export_ME7 as export_ME7
import madgraph.iolibs.save_load_object as save_load_object
import madgraph.integrator.phase_space_generators as phase_space_generators
import models.import_ufo as import_ufo
from madgraph import MG4DIR, MG5DIR

import tests.unit_tests as unittest

pjoin = os.path.join

root_path = os.path.dirname(os.path.realpath( __file__ ))

class ME7ContributionTest(unittest.TestCase):
    """Test class for functionalities related to contributions."""
    
    mypartlist = base_objects.ParticleList()
    myinterlist = base_objects.InteractionList()
    mymodel = loop_base_objects.LoopModel()
    mylegs = base_objects.LegList()
    myprocess = base_objects.Process()
    # To be set during setUp
    LO_contributions  = None
    NLO_contributions = None
 
    @misc.mute_logger()
    def setUp(self):

        # Setting up a dumb model

        # A gluon
        self.mypartlist.append(base_objects.Particle({'name':'g',
                      'antiname':'g',
                      'spin':3,
                      'color':8,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'g',
                      'antitexname':'g',
                      'line':'curly',
                      'charge':0.,
                      'pdg_code':21,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':True}))

        # A quark U and its antiparticle
        self.mypartlist.append(base_objects.Particle({'name':'u',
                      'antiname':'u~',
                      'spin':2,
                      'color':3,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'u',
                      'antitexname':'\bar u',
                      'line':'straight',
                      'charge':2. / 3.,
                      'pdg_code':2,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':False}))
        antiu = copy.copy(self.mypartlist[1])
        antiu.set('is_part', False)

        # A quark D and its antiparticle
        self.mypartlist.append(base_objects.Particle({'name':'d',
                      'antiname':'d~',
                      'spin':2,
                      'color':3,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'u',
                      'antitexname':'\bar u',
                      'line':'straight',
                      'charge':-1. / 3.,
                      'pdg_code':1,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':False}))
        antid = copy.copy(self.mypartlist[2])
        antid.set('is_part', False)

        # A photon
        self.mypartlist.append(base_objects.Particle({'name':'a',
                      'antiname':'a',
                      'spin':3,
                      'color':1,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'\gamma',
                      'antitexname':'\gamma',
                      'line':'wavy',
                      'charge':0.,
                      'pdg_code':22,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':True}))

        # A Higgs
        self.mypartlist.append(base_objects.Particle({'name':'h',
                      'antiname':'h',
                      'spin':1,
                      'color':1,
                      'mass':'mh',
                      'width':'wh',
                      'texname':'h',
                      'antitexname':'h',
                      'line':'dashed',
                      'charge':0.,
                      'pdg_code':25,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':True}))

        # 3 gluon vertiex
        self.myinterlist.append(base_objects.Interaction({
                      'id': 1,
                      'particles': base_objects.ParticleList(
                              [self.mypartlist[0]] * 3
                      ),
                      'color': [color.ColorString([color.f(0, 1, 2)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'G'},
                      'orders':{'QCD':1}}))

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
                      'lorentz':['L(p1,p2,p3)', 'L(p2,p3,p1)', 'L3'],
                      'couplings':{(0, 0):'G^2',
                                   (1, 1):'G^2',
                                   (2, 2):'G^2'},
                      'orders':{'QCD':2}}))

        # Gluon couplings to up and down quarks
        self.myinterlist.append(base_objects.Interaction({
                      'id': 3,
                      'particles': base_objects.ParticleList(
                                            [self.mypartlist[1],
                                             antiu,
                                             self.mypartlist[0]]),
                      'color': [color.ColorString([color.T(2, 0, 1)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQQ'},
                      'orders':{'QCD':1}}))

        self.myinterlist.append(base_objects.Interaction({
                      'id': 4,
                      'particles': base_objects.ParticleList(
                                            [self.mypartlist[2],
                                             antid,
                                             self.mypartlist[0]]),
                      'color': [color.ColorString([color.T(2, 0, 1)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQQ'},
                      'orders':{'QCD':1}}))

        # Photon coupling to up
        self.myinterlist.append(base_objects.Interaction({
                      'id': 5,
                      'particles': base_objects.ParticleList(
                                            [self.mypartlist[1],
                                             antiu,
                                             self.mypartlist[3]]),
                      'color': [color.ColorString([color.T(0, 1)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQED'},
                      'orders':{'QED':1}}))

        self.mymodel.set('particles', self.mypartlist)
        self.mymodel.set('interactions', self.myinterlist)
        self.mymodel.set('name', "sm4test")

        # Gluon couplings to a Higgs
        self.myinterlist.append(base_objects.Interaction({
                      'id': 6,
                      'particles': base_objects.ParticleList(
                                            [self.mypartlist[4],
                                             self.mypartlist[0],
                                             self.mypartlist[0]]),
                      'color': [color.ColorString([color.Tr(1, 2)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GGH'},
            'orders':{'QCD':2, 'QED':1}}))

        # Setting up the process p p > h j j and its subtraction
        self.mylegs = base_objects.MultiLegList([
            base_objects.MultiLeg(
                    {'ids': [1,2,-1,-2,21], 'state': base_objects.Leg.INITIAL}),
            base_objects.MultiLeg(
                    {'ids': [1,2,-1,-2,21], 'state': base_objects.Leg.INITIAL}),
            base_objects.MultiLeg(
                    {'ids': [25], 'state': base_objects.Leg.FINAL}),
            base_objects.MultiLeg(
                    {'ids': [21,1,-1], 'state': base_objects.Leg.FINAL})
        ])

        self.myprocdef = base_objects.ProcessDefinition({
            'legs': self.mylegs,
            'model': self.mymodel
        })

        # Generate only Born LO contributions
        with misc.TMP_directory(debug=False) as tmp_path:
            
            # Generate the output for this.
            self.madgraph_cmd = cmd.MasterCmd(main='MadGraph')
            self.madgraph_cmd._curr_model = self.mymodel
            self.madgraph_cmd._export_dir = pjoin(tmp_path,'ME7ContributionTest_LO')

            # Generate contributions
            generation_options = {'ME7_definition': True, 
                                  'diagram_filter': False, 
                                  'LO': True, 
                                  'NNLO': [], 
                                  'optimize': False, 
                                  'NLO': [], 
                                  'loop_induced': [],
                                  'ignore_contributions' : []}

            self.madgraph_cmd.add_contributions(self.myprocdef, generation_options)
            LO_contributions = self.madgraph_cmd._curr_contribs
            LO_contributions.apply_method_to_all_contribs(
                    'generate_amplitudes', log='Generate diagrams for')

            self.exporter = export_ME7.ME7Exporter(self.madgraph_cmd, False, group_subprocesses=True)
            self.exporter.pass_information_from_cmd(self.madgraph_cmd)
            self.exporter.copy_template(self.madgraph_cmd._curr_model)
            self.exporter.export(True, args=[])
            # We don't want to finalize and output the model which is not complete.
            # The export above was enough to have fully functional contributions to test
#            self.exporter.finalize(['nojpeg'], self.madgraph_cmd.history)
            self.LO_contributions = self.madgraph_cmd._curr_contribs         


        # Generate all NLO contributions 
        with misc.TMP_directory(debug=False) as tmp_path:

            # Generate the output for this.
            self.madgraph_cmd = cmd.MasterCmd(main='MadGraph')
            self.madgraph_cmd._curr_model = self.mymodel
            self.madgraph_cmd._export_dir = pjoin(tmp_path,'ME7ContributionTest_LO')

            # Generate contributions
            generation_options = {'ME7_definition': True, 
                                  'diagram_filter': False, 
                                  'LO': True, 
                                  'NNLO': [], 
                                  'optimize': False, 
                                  'NLO': ['QCD'], 
                                  'loop_induced': [],
                                  'ignore_contributions' : ['V']}

            self.madgraph_cmd.add_contributions(self.myprocdef, generation_options)  

            self.madgraph_cmd._curr_contribs.apply_method_to_all_contribs(
                    'generate_amplitudes', log='Generate diagrams for')

            self.exporter = export_ME7.ME7Exporter(self.madgraph_cmd, False, group_subprocesses=True)
            self.exporter.pass_information_from_cmd(self.madgraph_cmd)
            self.exporter.copy_template(self.madgraph_cmd._curr_model)
            self.exporter.export(True, args=[])
            # The export above was enough to have fully functional contributions to test            
#            self.exporter.finalize(['nojpeg'], self.madgraph_cmd.history)
            self.NLO_contributions = self.madgraph_cmd._curr_contribs 

    def test_subtraction_current_generation_and_access(self):
        """Test the generation of counterterms in single real-emission type of contributions."""
        
        real_emission_contrib = self.NLO_contributions.get_contributions_of_type(contributions.Contribution_R)[0]
        
        accessors_dict = contributions.MEAccessorDict()
        
        all_currents = real_emission_contrib.get_all_necessary_subtraction_currents(accessors_dict)
        # Print all currents
        misc.sprint('\n'.join([str(current) for current in all_currents]))
        
        with misc.TMP_directory(debug=False) as tmp_path:
            
            # For now test the handling of the first current only.
            all_currents = all_currents[:1]
            
            real_emission_contrib.add_current_accessors(accessors_dict, tmp_path, all_currents)
            # Print all accessor keys
            #misc.sprint('\n'.join([str(key) for key in accessors_dict]))
            misc.sprint("A total of %d accessor keys have been generated."%len(accessors_dict))
            
            # Try and call one current

            # We must initialize the model parameters for this first
            model_with_params_set = import_ufo.import_model(
                pjoin(MG5DIR,'models','loop_sm'), prefix=True,
                complex_mass_scheme = False )
            model_with_params_set.pass_particles_name_in_mg_default()
            model_with_params_set.set_parameters_and_couplings(
                param_card = pjoin(MG5DIR,'models','loop_sm','restrict_default.dat'), 
                scale=0.118, 
                complex_mass_scheme=False)
            accessors_dict.synchronize(model = model_with_params_set)
    
            one_current = all_currents[0]
            # Use a random PS point
            a_PS_point = dict((id, phase_space_generators.LorentzVector([10.0,2.0,3.0,4.0])) for id in range(10))
            misc.sprint('Call of the current %s:'%str(one_current))
            specific_eval, all_evals = accessors_dict(
                one_current,
                a_PS_point,
                hel_config=None
            )
            misc.sprint('specific_eval is:\n',specific_eval)
            misc.sprint('all_evals are:\n',all_evals)
