"""Implementation of the MPLs that can be expressed in terms of simple dilogarithms.
The evaluation is based on the mpmath implementation"""

from mpmath import polylog
from math import log

def li2(x):
    """Float-returning dilogarithm"""
    return float(polylog(2,x))

def mpl_dilog(entries,x):
    """We use the formula between general real weight-two MPL and logs and dilogs"""
    a = entries[0]
    b = entries[1]
    if a==b:
        return 0.5*log(x)**2
    elif a == 0.:
        return -li2(x/b)
    elif b == 0.:
        return log(x)*log(1-x/a)+li2(x/a)
    else:
        return -li2(b/(b-a)) + li2((b - x)/(b-a)) + log((x-a)/(b-a))*log(1 - x/b)


def mpl_dilog_conditions(entries,x):
    """Check the reality condition on the entries"""
    if len(entries) != 2:
        return False
    condition = True
    if b!=0:
        condition = condition and x/b < 1.
    else:
        condition = condition and x > 0.
    if a != 0:
        condition = condition and a/x > 0.
    else:
        condition = condition and x > 0.
