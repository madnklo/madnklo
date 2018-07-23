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

"""Unit test library for the NNLO subtraction currents"""

import copy
import fractions
import tests.input_files.simple_qcd as simple_qcd

#import madgraph.core.subtraction as sub
import madgraph.iolibs.template_files.subtraction.colorful.NNLO.local_currents as thecurrents
import tests.unit_tests as unittest


class NNLOCurrentsTest(unittest.TestCase):
    """Test class for the NNLO currents"""

    def setUp(self):
        self.ggg = thecurrents.QCD_final_collinear_0_ggg(simple_qcd.model)
        self.z1 = 0.20438592914045145
        self.z2 = 0.03438168578192958
        self.z3 = 0.3312881257808562
        self.s12 = 0.06778492406087142
        self.s13 = 0.6531483241035966
        self.s23 = 0.10987224288268904
        self.s123 = self.s12 + self.s13 + self.s23
        self.t123 = (2 * (self.z1*self.s23 - self.z2*self.s13)/(self.z1 + self.z2)
                    + (self.z1 - self.z2)*self.s12/(self.z1 + self.z2)
                     )

    def test_g_2_ggg_coef_gmunu_permuted(self):
        """Numerical test of the permuted coefficient of g_mu_nu"""
        computed = self.ggg.evaluate_permuted_gmunu(self.z1,self.z2,self.z3,self.s12,self.s13,self.s23)
        expected = -41201.58049492838
        self.assertAlmostEqual(expected, computed)

    def test_g_2_ggg_coef_gmunu(self):
        """Numerical test of the unpermuted coefficient of g_mu_nu"""
        computed = self.ggg.evaluate_gmunu(self.z1,self.z2,self.z3,self.s12,self.s13,self.s23,self.s123,self.t123)
        expected = -4488.262231052243
        self.assertAlmostEqual(expected, computed)

    def test_g_2_ggg_coef_k1k1(self):
        """Numerical test of the coefficient of g_mu_nu"""
        computed = self.ggg.evaluatecoefk1k1(self.z1,self.z2,self.z3,self.s12,self.s13,self.s23)
        expected = 6750.75554682
        self.assertAlmostEqual(expected, computed)

    def test_g_2_ggg_coef_k2k2(self):
        """Numerical test of the coefficient of g_mu_nu"""
        computed = self.ggg.evaluatecoefk2k2(self.z1, self.z2, self.z3, self.s12, self.s13, self.s23)
        expected = -2831.39768039
        self.assertAlmostEqual(expected, computed)

    def test_g_2_ggg_coef_k3k3(self):
        """Numerical test of the coefficient of g_mu_nu"""
        computed = self.ggg.evaluatecoefk3k3(self.z1, self.z2, self.z3, self.s12, self.s13, self.s23)
        expected = 5336.16427029
        self.assertAlmostEqual(expected, computed)








if __name__ == '__main__':
    unittest.main()
