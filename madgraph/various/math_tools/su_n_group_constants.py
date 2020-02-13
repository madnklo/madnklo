"""Casimir operators of SU(N) groups
This is somewhat redundant with the madgraph.core.color_algebra module but here we just want numerical
constants to avoid redundant calculations
"""
# Let's avoid fraction issues
from __future__ import division

class SU_N_group:
    def __init__(self,nc=3,tr=0.5):
        self.NC = nc
        self.CF = (nc**2-1)/(2*nc)
        self.CA = nc
        self.TR = tr

    def casimir(self, representation_dimension):
        """Return the Casimir invariant of a representation with dimension representation_dimension"""
        if abs(representation_dimension) == 3:
            return self.CF
        elif representation_dimension == 8:
            return self.CA


# The usual convention for QCD
SU3 = SU_N_group(nc=3,tr=1./2.)
