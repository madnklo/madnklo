"""This module implement a singleton class that evaluates MPLs on entries based on a list of implementations
that are tried sequentially for compatibility"""

import logging
logging.getLogger(__name__)


class MPL_dispatcher(object):
    """Dispatcher class for MPLs"""
    def __init__(self,implementations):
        self.implementations=implementations
    def __call__(self,entries,x):
        """Sequentially try the different implementations of MPLs to evaluate them"""
        for implementation in self.implementations:
            result = implementation(entries,x)
            if result is not None:
                return result

        raise NotImplementedError("The entries ({},{}) are incompatible with all MPL"
                                  "implementations available".format(entries, x))
