"""Implementation of weight one MPLs in terms of python math.log"""
from math import log
import logging
logger=logging.getLogger(__name__)

def mpl_log_conditions(entries,x):
    """Enforce that logs are weight one and evaluated for positive entries"""
    if not len(entries)==1:
        return False
    elif entries[0] == 0.:
        if x<=0:
            logger.warning("Trying to evaluate a log with negative argument. This is probably bad")
            return False
        else:
            return True
    elif (x/entries[0]>1):
        logger.warning("Trying to evaluate a log with negative argument. This is probably bad")
        return False
    else:
        return True

def mpl_log(entries,x):
    """Return G([a],x) in terms of log"""
    if entries[0] == 0.:
        return log(x)
    else:
        return log(1-x/entries[0])