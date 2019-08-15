'''Explicit expressions for the integrals of the modified partial-fractionned soft scheme at NLO
Computations performed in ND's
research/Active/Subtractions/Cbars/NLO_cbars/
    integrated_hard_coll_DS.nb
    integrated_modified_DS.nb
TODO put this in some common place
'''
import logging
import math
from madgraph import MG5DIR
from madgraph.various import misc
from madgraph.core.base_objects import EpsilonExpansion


logger=logging.getLogger("SubtractionCurrents")

import madgraph.various.math_tools.mpl as MPL

# =============================================
#          Harcoded Expressions
# =============================================

def Soft_FF_distributed_modified_eikonal(y0):
    '''Integral of the modified partial fractionned eikonal divided by the overall normalization of
    $(\mu^2/\tilde{s})^(\epsilon)$ (included in the `evaluate_integrated_current` function)

    :param y0: y0*stilde is the upper value of the invariant mass of the collinear pair,
    where stilde is the invariant mass of the reduced emitting dipole
    :type y0: float
    :return: numerical value of the finite part of the integral
    :rtype: EpsilonExpansion
    '''
    return -2*(2. - math.pi**2/4. + MPL.G([0,-1],y0) - math.log(y0) - \
    y0*math.log(y0) - math.log(y0)**2/2. + math.log(1 + y0) + \
    y0*math.log(1 + y0))

