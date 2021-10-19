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
"""Common factors and cuts for subtraction schemes."""

from scipy import special


#damping_factors = [alpha_soft, beta_FF, beta_FI, beta_IF, beta_II]  ~ beta_jr

alpha = 0.0
beta_FF = 0.0
beta_FI = 0.0
beta_IF = 0.0
beta_II = 0.0

damping_factors = [alpha, beta_FF, beta_FI, beta_IF, beta_II]



Eulergamma = 0.57721566490153286061

def polygamma(self):

    return special.polygamma(1, 2 +  self)

def A1(self):

    return Eulergamma + special.psi(1 + self)

def A2(self):

    return Eulergamma + special.psi(2 + self)


def no_factor(**opts):
    return 1.0

def no_cut(**opts):
    return False
