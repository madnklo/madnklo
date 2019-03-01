"""Implementation of weight one MPLs in terms of python math.log"""
from math import log

def mpl_log_conditions(entries,x):
    """Enforce that logs are weight one and evaluated for positive entries"""
    return (len(entries)==1) and (entries[0]==0. or x/entries < 1.)

def mpl_log(entries,x):
    """Return G([a],x) in terms of log"""
    if entries[0] == 0.:
        return log(x)
    else:
        return log(1-x/entries[0])