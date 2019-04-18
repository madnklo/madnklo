# =============================================
#               MPL CLASS
# =============================================
import ctypes
import logging
import math
from madgraph import MG5DIR
from madgraph.various import misc


logger=logging.getLogger("SubtractionCurrents")

import madgraph.various.math_tools.mpl as MPL

# =============================================
#          Harcoded Expression Class (HE)
# =============================================

class HE(object):
    """
    This is a class whose static methods are the hardcoded expressions for the integrated currents
    """

    ############################################################################################################
    # Beam-soft integrated counterterms
    # These are (Eikonal - its collinear limits) integrated over the unresolved phase space of Vittorio's
    # soft-recoiling-against-initial-state mapping. This is defined in sec 5.1.2 of his notes
    # Each of these should multiply a color-correlated amplitude and when summed over colors reconstruct the S+CS counterterms
    # Detailed information in Nicolas Deutschmann's handwritten notes from 03.10.2018
    # calculation in ndeutsch/workdir/research/Active/Subtractions/ReverseUnitarityCT/\
    # Vittorio_Soft/ExpandedSbars
    ############################################################################################################
    @staticmethod
    def integrated_bs_endpoint_pole(dipole_invariant):
        """Pole of the endpoint contribution of the collinear-subtracted eikonal integrated over the SoftVsInitial mapping unresolved PS (BS)"""
        return 2.*math.log(dipole_invariant)

    @staticmethod
    def integrated_bs_endpoint_finite(dipole_invariant):
        """Finite part of the endpoint contribution of the collinear-subtracted eikonal integrated over the SoftVsInitial mapping unresolved PS (BS)"""
        return  - math.log(16.)*math.log(dipole_invariant) - math.log(dipole_invariant)**2 + 2.*(MPL.G([1,0],dipole_invariant)-math.pi**2/6.)

    @staticmethod
    def integrated_bs_bulk_finite(dipole_invariant,xi):
        """Finite part of the bulk contribution of the collinear-subtracted eikonal integrated over the SoftVsInitial mapping unresolved PS (BS)"""
        return (8. * xi * math.log(dipole_invariant)) / ((-1. + xi) * (1. + xi))

    @staticmethod
    def integrated_bs_counterterm_finite(dipole_invariant,xi):
        """Finite part of the counterterm contribution of the collinear-subtracted eikonal integrated over the SoftVsInitial mapping unresolved PS (BS)"""
        return (4. * math.log(dipole_invariant) ) / (-1. + xi)

