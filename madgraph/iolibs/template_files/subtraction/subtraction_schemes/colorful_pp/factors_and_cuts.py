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
# TODO Ideally these functions should be part of colorful and colorful_pp

y_0_prime   = 0.5

def cut_initial_coll(**opts):

    pA    = opts['pA']
    pR    = opts['pR']
    Q     = opts['Q']
    y_0p  = (2.*pA.dot(pR))/Q.square()
    # Include the counterterm only up to y_0_prime
    return y_0p > y_0_prime

def no_factor(**opts):
    return 1.0

def no_cut(**opts):
    return False