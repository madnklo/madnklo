"""Main dispatcher module that centralizes all polylogarithm evaluations and calls the relevant specialized tools"""

import logging
from numbers import Number
import zero_weight_mpl
import python_log
import dilogarithms
import ginac_mpl

logger=logging.getLogger(__name__)

class MPL_implementation(object):
    """Representation of one type of implementation of some subset of MPL that verifies conditions
    checked in compatibility conditions.
    """
    def __init__(self,implementation,conditions):
        self.implementation = implementation
        self.compatibility_conditions=conditions

    def __call__(self, entries, x):
        """Check the compatibility conditions. If possible return the value of the MPL, if not return None
        For now we enforce a hard check that all entries are real (numeric types)
        """
        logger.debug("Calling the MPL implementation %s"%self.implementation.__name__)
        all_reals = isinstance(x,Number) and all([isinstance(a,Number) for a in entries])
        if not all_reals:
            raise NotImplementedError("The entries ({},{}) are incompatible with all MPL"
                                      "implementations available".format(entries, x))
        if self.compatibility_conditions(entries,x):
            return self.implementation(entries,x)
        else:
            return None

class MPL_implementation_list(list):
    """List of MPL implementations"""
    def add_implementation(self,implementation,conditions):
        self.append(MPL_implementation(implementation,conditions))


mpl_implementations = MPL_implementation_list()
mpl_implementations.add_implementation(zero_weight_mpl.zero_weight_mpl, zero_weight_mpl.zero_weight_mpl_conditions)
mpl_implementations.add_implementation(python_log.mpl_log, python_log.mpl_log_conditions)
mpl_implementations.add_implementation(dilogarithms.mpl_dilog,dilogarithms.mpl_dilog_conditions)
mpl_implementations.add_implementation(ginac_mpl.GinacMPL.G, ginac_mpl.GinacMPL.G_conditions)





