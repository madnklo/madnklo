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
import unittest
import logging
import random
import multiprocessing

logger = logging.getLogger('test_PS_volume')

import madgraph.integrator.functions as functions
import madgraph.integrator.integrands as integrands
import madgraph.integrator.phase_space_generators as psgen
import madgraph.integrator.vegas3_integrator as vegas
import madgraph.integrator.mappings as mappings
import madgraph.core.subtraction as sub
import madgraph.various.cluster as cluster
import madgraph.various.misc as misc

#=========================================================================================
# Global settings
#=========================================================================================

# Parallelization on/off
this_cluster = cluster.MultiCore(nb_core=multiprocessing.cpu_count())
# this_cluster = cluster.MultiCore(nb_core=1)

# Default seeds for the tests, more seeds = stronger check, but slower
default_seeds = (1, 2 ,3, )

#=========================================================================================
# Direct PS integrand
#=========================================================================================

class DirectPS(functions.VirtualFunction):

    def __init__(self, final_masses, cut_function=None):

        self.generator = psgen.FlatInvertiblePhasespace(
            [0., ] * 2, final_masses, (0.5, 0.5), beam_types=(0, 0) )
        self.cut_function = cut_function

    def __call__(self, continuous_inputs, discrete_inputs, **opts):

        ps, wgt = self.generator.generateKinematics(1., continuous_inputs)
        if self.cut_function:
            ps_dict = ps.to_dict()
            if not self.cut_function(ps_dict):
                return {'weight': 0.}
        return {'weight': wgt}

#=========================================================================================
# Mapped PS integrands
#=========================================================================================

class MappedPS(functions.VirtualFunction):

    @staticmethod
    def leg(i):

        return sub.SubtractionLeg(i, 0, sub.SubtractionLeg.FINAL)

    @staticmethod
    def structure(i):

        return sub.CollStructure(legs=[MappedPS.leg(i), ])

    def __init__(self, mapping, intermediate_masses, final_masses):

        self.mapping = mapping
        assert len(intermediate_masses) == len(final_masses)
        self.generator = psgen.FlatInvertiblePhasespace(
            [0., ] * 2, intermediate_masses, (0.5, 0.5), beam_types=(0, 0) )
        self.momenta_dict = sub.bidict({1: frozenset((1, )), 2: frozenset((2, ))})
        legs = []
        substructures = []
        self.final_masses = dict()
        self.intermediate_masses = dict()
        for i in range(len(intermediate_masses)):
            self.momenta_dict[i+3] = frozenset((i+3, ))
            if intermediate_masses[i] == final_masses[i]:
                legs.append(MappedPS.leg(i+3))
            else:
                substructures.append(MappedPS.structure(i+3))
                self.final_masses['s' + str(i+3)] = final_masses[i]**2
                self.intermediate_masses['s' + str(i+3)] = intermediate_masses[i] ** 2
        self.singular_structure = sub.SingularStructure(
            legs=legs, substructures=substructures )
        print self.singular_structure

class LowerMappedPS(MappedPS):

    def __call__(self, continuous_inputs, discrete_inputs, **opts):

        ps, wgt = self.generator.generateKinematics(1., continuous_inputs)
        ps_dict = ps.to_dict()
        try:
            res = self.mapping.map_to_lower_multiplicity(
                ps_dict, self.singular_structure, self.momenta_dict,
                None, True, self.final_masses )
            return {'weight': wgt / res['jacobian']}
        except mappings.FailedMapping:
            return {'weight': 0.}

class HigherMappedPS(MappedPS):

    def __call__(self, continuous_inputs, discrete_inputs, **opts):

        ps, wgt = self.generator.generateKinematics(1., continuous_inputs)
        ps_dict = ps.to_dict()
        try:
            res = self.mapping.map_to_higher_multiplicity(
                ps_dict, self.singular_structure, self.momenta_dict,
                self.final_masses, True )
            return {'weight': wgt * res['jacobian']}
        except mappings.FailedMapping:
            return {'weight': 0.}

#=========================================================================================
# TestPSVolume
#=========================================================================================

class PSVolumeTest(object):
    """Collection of functions to test phase space volumes."""

    verbosity = 1
    survey_n_iterations = 5
    refine_n_iterations = 5
    survey_n_points = 2000
    refine_n_points = 10000
    max_sigmas = 5

    @staticmethod
    def randomize(pars):

        n_changing = random.randint(
            pars.get('min_changing', 1), pars.get('max_changing', 1) )
        n_unchanged = random.randint(
            pars.get('min_unchanged', 1), pars.get('max_unchanged', 1) )
        if pars.get('higher_massive', True):
            masses_higher_changing = [random.random() for _ in range(n_changing)]
        else:
            masses_higher_changing = [0., ] * n_changing
        if pars.get('lower_massive', False):
            masses_lower_changing = [random.random() for _ in range(n_changing)]
        else:
            masses_lower_changing = [0., ] * n_changing
        if pars.get('recoilers_massive', False):
            masses_recoilers = [random.random() for _ in range(n_unchanged)]
        else:
            masses_recoilers = [0., ] * n_unchanged
        E_higher = (sum(masses_higher_changing) + sum(masses_recoilers)) / random.random()
        E_lower  = (sum(masses_lower_changing ) + sum(masses_recoilers)) / random.random()
        E = max([E_higher, E_lower])
        masses_higher_changing = [m / E for m in masses_higher_changing]
        masses_lower_changing = [m / E for m in masses_lower_changing]
        masses_recoilers = [m / E for m in masses_recoilers]
        if pars['mapto'] == "lower":
            pars['intermediate_masses'] = masses_higher_changing + masses_recoilers
            pars['final_masses'] = masses_lower_changing + masses_recoilers
        elif pars['mapto'] == "higher":
            pars['intermediate_masses'] = masses_lower_changing + masses_recoilers
            pars['final_masses'] = masses_higher_changing + masses_recoilers
        else:
            raise ValueError

    @staticmethod
    def needs_benchmark_integration(final_masses):

        for mass in final_masses:
            if mass != 0.:
                return True
        return False

    @staticmethod
    def setup(pars):

        cut = None
        if pars['mapto'] == "lower":
            mapped_PS = LowerMappedPS(
                pars['mapping'], pars['intermediate_masses'], pars['final_masses'] )
            if pars.get('check_cut', False):
                cut = lambda ps: pars['mapping'].can_map_to_higher_multiplicity(
                    ps, mapped_PS.singular_structure, mapped_PS.momenta_dict,
                    mapped_PS.intermediate_masses )
        elif pars['mapto'] == "higher":
            mapped_PS = HigherMappedPS(
                pars['mapping'], pars['intermediate_masses'], pars['final_masses'] )
            if pars.get('check_cut', False):
                cut = lambda ps: pars['mapping'].can_map_to_lower_multiplicity(
                    ps, mapped_PS.singular_structure, mapped_PS.momenta_dict,
                    None, masses=mapped_PS.intermediate_masses )
        else:
            raise ValueError
        if (PSVolumeTest.needs_benchmark_integration(pars['final_masses'])
            or pars.get('check_cut', False)):
            pars['benchmark_integrand'] = integrands.VirtualIntegrand()
            direct_PS = DirectPS(pars['final_masses'], cut)
            pars['benchmark_integrand'].set_dimensions(
                direct_PS.generator.get_dimensions() )
            pars['benchmark_integrand'].function_list = functions.FunctionList()
            pars['benchmark_integrand'].function_list.append(direct_PS)
            pars['benchmark_integrator'] = vegas.Vegas3Integrator(
                pars['benchmark_integrand'],
                cluster=this_cluster,
                verbosity=PSVolumeTest.verbosity-1,
                survey_n_iterations=PSVolumeTest.survey_n_iterations,
                survey_n_points=PSVolumeTest.survey_n_points,
                refine_n_iterations=PSVolumeTest.refine_n_iterations,
                refine_n_points=PSVolumeTest.refine_n_points )
        pars['integrand'] = integrands.VirtualIntegrand()
        pars['integrand'].set_dimensions(mapped_PS.generator.get_dimensions())
        pars['integrand'].function_list = functions.FunctionList()
        pars['integrand'].function_list.append(mapped_PS)
        # Setup the integrator
        pars['integrator'] = vegas.Vegas3Integrator(
            pars['integrand'],
            verbosity=PSVolumeTest.verbosity-1,
            cluster=this_cluster,
            survey_n_iterations=PSVolumeTest.survey_n_iterations,
            survey_n_points=PSVolumeTest.survey_n_points,
            refine_n_iterations=PSVolumeTest.refine_n_iterations,
            refine_n_points=PSVolumeTest.refine_n_points )

    @staticmethod
    def test_PS_volume(pars, test):
        """Compute the phase-space volume."""

        for seed in pars.get('seeds', default_seeds):

            print "-" * 200
            if PSVolumeTest.verbosity > 0:
                print "Seed", seed
            random.seed(seed)
            PSVolumeTest.randomize(pars)
            if PSVolumeTest.verbosity > 0:
                print pars['intermediate_masses'], ">", pars['final_masses']
            PSVolumeTest.setup(pars)
            # Compute the benchmark phase space
            if pars.has_key('benchmark_integrator'):
                bmk, bmkerr = pars['benchmark_integrator'].integrate()
                print "Benchmark: %.10e +/- %.2e" % (bmk, bmkerr)
            else:
                bmk = psgen.FlatInvertiblePhasespace.get_flatWeights(
                    1., len(pars['final_masses']))
                bmkerr = 0.
                print "Benchmark: ", bmk
            # Compute the mapped phase space
            res, reserr = pars['integrator'].integrate()
            print 'Result: %.10e +/- %.2e' % (res, reserr)
            err = (bmkerr**2 + reserr**2) ** 0.5
            print 'Difference: %.2e std deviations' % ((res-bmk)/err)
            test.assertLess(
                abs(res-bmk), PSVolumeTest.max_sigmas*err, "Jacobian does not check out")

#=========================================================================================
# Test invariant mass mappings
#=========================================================================================

class FinalZeroMassesMappingPSVolumeTest(unittest.TestCase):
    """Test class for the jacobian of FinalZeroMassesMapping."""

    # Test settings
    pars_lower = {
        'mapping': mappings.FinalZeroMassesMapping(),
        'mapto': 'lower', 'check_cut': False,
        'min_changing': 2, 'max_changing': 5,
        'min_unchanged': 0, 'max_unchanged': 0,
        'higher_massive': True, 'lower_massive': False,
        'recoilers_massive': False, }
    pars_higher = {
        'mapping': mappings.FinalZeroMassesMapping(),
        'mapto': 'higher', 'check_cut': False,
        'min_changing': 2, 'max_changing': 5,
        'min_unchanged': 0, 'max_unchanged': 0,
        'higher_massive': True, 'lower_massive': False,
        'recoilers_massive': False, }

    def test_FinalZeroMassesMapping_PSvolume_lower(self):

        PSVolumeTest.test_PS_volume(self.pars_lower, self)

    def test_FinalZeroMassesMapping_PSvolume_higher(self):

        PSVolumeTest.test_PS_volume(self.pars_higher, self)

class FinalMassesMappingPSVolumeTest(unittest.TestCase):
    """Test class for the jacobian of FinalMassesMapping."""

    # Test settings
    pars_lower = {
        'mapping': mappings.FinalMassesMapping(),
        'mapto': 'lower', 'check_cut': False,
        'min_changing': 2, 'max_changing': 5,
        'min_unchanged': 0, 'max_unchanged': 0,
        'higher_massive': True, 'lower_massive': True,
        'recoilers_massive': False, }
    pars_higher = {
        'mapping': mappings.FinalMassesMapping(),
        'mapto': 'higher', 'check_cut': False,
        'min_changing': 2, 'max_changing': 5,
        'min_unchanged': 0, 'max_unchanged': 0,
        'higher_massive': True, 'lower_massive': True,
        'recoilers_massive': False, }

    def test_FinalMassesMapping_PSvolume_lower(self):

        PSVolumeTest.test_PS_volume(self.pars_lower, self)

    def test_FinalMassesMapping_PSvolume_higher(self):

        PSVolumeTest.test_PS_volume(self.pars_higher, self)

class FinalRescalingOneMappingPSVolumeTest(unittest.TestCase):
    """Test class for the jacobian of FinalRescalingOneMapping."""

    # Test settings
    pars_lower = {
        'mapping': mappings.FinalRescalingOneMapping(),
        'mapto': 'lower', 'check_cut': False,
        'min_changing': 1, 'max_changing': 1,
        'min_unchanged': 1, 'max_unchanged': 4,
        'higher_massive': True, 'lower_massive': False,
        'recoilers_massive': False, }
    pars_higher = {
        'mapping': mappings.FinalRescalingOneMapping(),
        'mapto': 'higher', 'check_cut': False,
        'min_changing': 1, 'max_changing': 1,
        'min_unchanged': 1, 'max_unchanged': 4,
        'higher_massive': True, 'lower_massive': False,
        'recoilers_massive': False, }

    def test_FinalRescalingOneMapping_PSvolume_lower(self):

        PSVolumeTest.test_PS_volume(self.pars_lower, self)

    def test_FinalRescalingOneMapping_PSvolume_higher(self):

        PSVolumeTest.test_PS_volume(self.pars_higher, self)

class FinalLorentzOneMappingPSVolumeTest(unittest.TestCase):
    """Test class for the jacobian of FinalLorentzOneMapping."""

    # Test settings
    pars_lower = {
        'mapping': mappings.FinalLorentzOneMapping(),
        'mapto': 'lower', 'check_cut': True,
        'min_changing': 1, 'max_changing': 1,
        'min_unchanged': 1, 'max_unchanged': 4,
        'higher_massive': True, 'lower_massive': False,
        'recoilers_massive': True, }
    pars_higher = {
        'mapping': mappings.FinalLorentzOneMapping(),
        'mapto': 'higher', 'check_cut': False,
        'min_changing': 1, 'max_changing': 1,
        'min_unchanged': 1, 'max_unchanged': 4,
        'higher_massive': True, 'lower_massive': False,
        'recoilers_massive': True, }

    def test_FinalLorentzOneMapping_PSvolume_lower(self):

        PSVolumeTest.test_PS_volume(self.pars_lower, self)

    def test_FinalLorentzOneMapping_PSvolume_higher(self):

        PSVolumeTest.test_PS_volume(self.pars_higher, self)

class FinalGroupingMappingPSVolumeTest(unittest.TestCase):
    """Test class for the jacobian of FinalGroupingMapping."""

    # Test settings
    # For min_unchanged=0 this is already tested with FinalMassesMapping
    # Skipping to avoid the case of min_changing=1, min_unchanged=0 which fails
    pars_lower = {
        'mapping': mappings.FinalGroupingMapping(),
        'mapto': 'lower', 'check_cut': False,
        'min_changing': 1, 'max_changing': 4,
        'min_unchanged': 1, 'max_unchanged': 4,
        'higher_massive': True, 'lower_massive': True,
        'recoilers_massive': True, }
    pars_higher = {
        'mapping': mappings.FinalGroupingMapping(),
        'mapto': 'higher', 'check_cut': False,
        'min_changing': 1, 'max_changing': 4,
        'min_unchanged': 1, 'max_unchanged': 4,
        'higher_massive': True, 'lower_massive': True,
        'recoilers_massive': True, }

    def test_FinalGroupingMapping_PSvolume_lower(self):

        PSVolumeTest.test_PS_volume(self.pars_lower, self)

    def test_FinalGroupingMapping_PSvolume_higher(self):

        PSVolumeTest.test_PS_volume(self.pars_higher, self)

class FinalLorentzMappingPSVolumeTest(unittest.TestCase):
    """Test class for the jacobian of FinalLorentzMapping."""

    # Test settings
    pars_lower = {
        'mapping': mappings.FinalLorentzMapping(),
        'mapto': 'lower', 'check_cut': True,
        'min_changing': 1, 'max_changing': 3,
        'min_unchanged': 1, 'max_unchanged': 4,
        'higher_massive': True, 'lower_massive': True,
        'recoilers_massive': True, }
    pars_higher = {
        'mapping': mappings.FinalLorentzMapping(),
        'mapto': 'higher', 'check_cut': True,
        'min_changing': 1, 'max_changing': 2,
        'min_unchanged': 1, 'max_unchanged': 2,
        'higher_massive': True, 'lower_massive': True,
        'recoilers_massive': True, }

    def test_FinalLorentzMapping_PSvolume_lower(self):

        PSVolumeTest.test_PS_volume(self.pars_lower, self)

    def test_FinalLorentzMapping_PSvolume_higher(self):

        PSVolumeTest.test_PS_volume(self.pars_higher, self)
