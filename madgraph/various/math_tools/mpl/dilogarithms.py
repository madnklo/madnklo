"""Implementation of the MPLs that can be expressed in terms of simple dilogarithms.
The evaluation is based on the mpmath implementation"""

from mpmath import polylog
from math import log
import logging

logger=logging.getLogger(__name__)

def li2(x):
    """Float-returning dilogarithm"""
    return float(polylog(2,x))

def mpl_dilog(entries,x):
    """We use the formula between general real weight-two MPL and logs and dilogs"""
    a = entries[0]
    b = entries[1]
    if a==b:
        if a==0.:
            return 0.5*log(x)**2
        else:
            return 0.5*log(1-x/a)**2
    elif a == 0.:
        return -li2(x/b)
    elif b == 0.:
        return log(x)*log(1-x/a)+li2(x/a)
    else:
        return -li2(b/(b-a)) + li2((b - x)/(b-a)) + log((x-a)/(b-a))*log(1 - x/b)


def mpl_dilog_conditions(entries,x):
    """Check the reality condition on the entries
    We do not check for singularities as these already raise errors in the actual functions
    """
    if len(entries) != 2:
        return False
    a = entries[0]
    b = entries[1]

    if a == b: # G(a,a,x) == 1/2 log^2(1-x/a) or log(x)^2
        # if x == b:
        #     msg = "Trying to evaluate log(0)^2"
        #     logger.error(msg)
        #     raise ValueError(msg)
        # el
        if b == 0.:
            if x>=0:
                return True
        else:
            if x/b <= 1.:
                return True
    elif a == 0: # G(0,b,x) == - Li2(x/b) (b==0 already caught before)
        if x/b <= 1.:
            return True
    elif b == 0: #G(a,0,x) = int log(t)/(t-a) has a log divergence at x=a
        if x>=0:
            if (a<=0 or a>=x):
                return True
            # elif a==x:
            #     msg = "Trying to evaluate log(0)"
            #     logger.error(msg)
            #     raise ValueError(msg)
    else: # General case: int log(1-t/b)/(t-a) dt
        if x/b <= 1. and (a/x <= 0. or a/x >= 1.):
            return True

    logger.warning("Trying to evaluate a weight 2 polylog outside of the reality condition. Something is most likely wrong")
    return False





