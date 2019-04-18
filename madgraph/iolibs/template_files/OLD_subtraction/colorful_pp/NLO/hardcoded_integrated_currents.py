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
    @staticmethod
    def CqqFF_Finite_Gabor_EEJJ(a0,y12in):
        '''Final-final collinear integrated counterterm for qq for m=2'''
        if abs(y12in-1.) < 1e-6:
            y12 = 1-1e-6
        else:
            y12 = y12in
        #misc.sprint("In CqqFF")
        #misc.sprint("y12 = " + str(y12))
        return (-10/9 - (8*a0)/(3*(a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)) +
   (4*a0*y12)/((a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)) - 
   (4*a0*y12**2)/(3*(a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)) + (2*MPL.G([0], a0))/3 + (2*MPL.G([0], y12))/3 + 
   (2*y12*MPL.G([y12/(a0*(-2 + y12))], 1))/(3*(-2 + y12)**2) - (16*a0*MPL.G([y12/(a0*(-1 + y12))], 1))/
    (3*(a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)) - (8*y12*MPL.G([y12/(a0*(-1 + y12))], 1))/
    (3*(a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)) + (8*a0*y12*MPL.G([y12/(a0*(-1 + y12))], 1))/
    ((a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)) + (8*y12**2*MPL.G([y12/(a0*(-1 + y12))], 1))/
    (3*(a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)) - (4*a0*y12**2*MPL.G([y12/(a0*(-1 + y12))], 1))/
    ((a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)) - (2*y12**3*MPL.G([y12/(a0*(-1 + y12))], 1))/
    (3*(a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)) + (2*a0*y12**3*MPL.G([y12/(a0*(-1 + y12))], 1))/
    (3*(a0*(-2 + y12) - y12)*(-2 + y12)**2*(-1 + y12)))

    @staticmethod
    def CggFF_Finite_Gabor_EEJJ(a0,y12in):
        '''Final-final collinear integrated counterterm for gg for m=2'''
        if abs(y12in-1.) < 1e-6:
            y12 = 1-1e-6
        else:
            y12 = y12in
        #misc.sprint("In CggFF")
        #misc.sprint("y12 = " + str(y12))
        return (100/9 - (2*math.pi**2)/3 + (40 - 42*y12 + 13*y12**2)/(3*(-2 + y12)**2*(-1 + y12)) -
   (y12*(4 - 84*a0 - 46*y12 + 126*a0*y12 + 42*y12**2 - 64*a0*y12**2 - 11*y12**3 + 11*a0*y12**3))/
    (3*(-2 + y12)**2*(-1 + y12)*(-2*a0 - y12 + a0*y12)) - (11*MPL.G([0], a0))/3 - (11*MPL.G([0], y12))/3 + 
   2*MPL.G([0], a0)*MPL.G([0], y12) + 2*MPL.G([0], a0)*MPL.G([-a0], y12) - 2*MPL.G([0], y12)*MPL.G([-y12], a0) - 
   (2*y12*(2 + y12)*MPL.G([y12/(-2 + y12)], a0))/(3*a0*(-2 + y12)**3) + 
   (2*y12*(2 - a0*(-2 + y12) + y12)*MPL.G([y12/(-2 + y12)], a0))/(3*a0*(-2 + y12)**3) - 
   (11*y12*MPL.G([y12/(-1 + y12)], a0))/(3*a0*(-1 + y12)**2) + (11*(a0 + y12 - a0*y12)*MPL.G([y12/(-1 + y12)], a0))/
    (3*a0*(-1 + y12)**2) - (4*MPL.G([0], y12)*MPL.G([y12/(-1 + y12)], a0))/(-1 + y12) + 2*MPL.G([0, 0], a0) + 
   2*MPL.G([0, 0], y12) - 4*MPL.G([0, -y12], a0) - (4*MPL.G([0, y12/(-1 + y12)], a0))/(-1 + y12) + 
   (4*y12*MPL.G([0, y12/(-1 + y12)], a0))/(-1 + y12) + 2*MPL.G([-a0, 0], y12) - 2*MPL.G([-y12, 0], a0) + 
   (4*MPL.G([y12/(-1 + y12), 0], a0))/(-1 + y12) - (4*MPL.G([y12/(-1 + y12), y12/(-1 + y12)], a0))/(-1 + y12))

    @staticmethod
    def CqgFF_Finite_Gabor_EEJJ(a0,y12in):
        '''Final-final collinear integrated counterterm for qg for m=2'''
        if abs(y12in-1.) < 1e-6:
            y12 = 1-1e-6
        else:
            y12 = y12in
        #misc.sprint("In CqgFF")
        #misc.sprint("y12 = " + str(y12))
        #misc.sprint("a0 = " + str(a0))
        return (5 - math.pi**2/3 + 3/(2*(-1 + y12)) - (3*y12)/(2*(-1 + y12)) - (3*MPL.G([0], a0))/2 - (3*MPL.G([0], y12))/2 +
   MPL.G([0], a0)*MPL.G([0], y12) + MPL.G([0], a0)*MPL.G([-a0], y12) - MPL.G([0], y12)*MPL.G([-y12], a0) + 
   (3*(a0 + y12 - a0*y12)*MPL.G([y12/(-1 + y12)], a0))/(2*a0*(-1 + y12)**2) - 
   ((3*y12 - 4*a0*MPL.G([0], y12) + 4*a0*y12*MPL.G([0], y12))*MPL.G([y12/(-1 + y12)], a0))/(2*a0*(-1 + y12)**2) + 
   MPL.G([0, 0], a0) + MPL.G([0, 0], y12) - 2*MPL.G([0, -y12], a0) - (2*MPL.G([0, y12/(-1 + y12)], a0))/(-1 + y12) + 
   (2*y12*MPL.G([0, y12/(-1 + y12)], a0))/(-1 + y12) + MPL.G([-a0, 0], y12) - MPL.G([-y12, 0], a0) + 
   (2*MPL.G([y12/(-1 + y12), 0], a0))/(-1 + y12) - (2*MPL.G([y12/(-1 + y12), y12/(-1 + y12)], a0))/(-1 + y12))


    @staticmethod
    def SoftFF_Finite_Gabor_EEJJ(y0,Yin):
        '''Final-final soft+soft-colinear integrated counterterm for qg for m=2'''
        if abs(Yin-1.) < 1e-6:
            Y = 1-1e-6
        else:
            Y = Yin
        return (-math.pi**2/6 + 2*(y0 + MPL.G([0], Y)*(y0 - MPL.G([0], y0))) - MPL.G([0, 0], Y) + MPL.G([1, 0], Y))


    @staticmethod
    def CqqFF_Finite_Gabor_DIVJAC_NOD0(a0,y12in):
        '''Final-final collinear integrated counterterm for qq canonically normalized for the rescaling mapping'''
        if abs(y12in-1.) < 1e-6:
            y12 = 1-1e-6
        else:
            y12 = y12in                
        return ((10*(-2 + y12)*y12 - 2*a0*(20 - 17*y12 + 5*y12**2))/(9*(a0*(-2 + y12) - y12)*(-2 + y12)) + 
   (4*(-1 + y12)*MPL.G([0], 2))/(3*(-2 + y12)**2) + (4*(3 - 3*y12 + y12**2)*MPL.G([0], a0))/(3*(-2 + y12)**2) - 
   (4*(-1 + y12)*MPL.G([0], y12))/(3*(-2 + y12)**2) + (2*MPL.G([a0/(-1 + a0)], y12))/3 + 
   (4*(-1 + y12)*MPL.G([(2*a0)/(-1 + a0)], y12))/(3*(-2 + y12)**2))

    @staticmethod
    def CggFF_Finite_Gabor_DIVJAC_NOD0(a0,y12in):
        '''Final-final collinear integrated counterterm for gg canonically normalized for the rescaling mapping'''
        if abs(y12in-1.) < 1e-6:
            y12 = 1-1e-6
        else:
            y12 = y12in        
        return (((-67 + 9*math.pi**2)*(-2 + y12)*y12 + a0*(268 - 9*math.pi**2*(-2 + y12)**2 - 262*y12 + 67*y12**2))/
    (9*(a0*(-2 + y12) - y12)*(-2 + y12)) - (4*(-1 + y12)*MPL.G([0], 2))/(3*(-2 + y12)**2) + 6*MPL.G([0], y12)**2 + 
   MPL.G([0], a0)*((-2*(42 - 42*y12 + 11*y12**2))/(3*(-2 + y12)**2) + (8*(-2 + y12)*MPL.G([1], y12))/y12) + 
   MPL.G([0], y12)*((4*(-1 + y12))/(3*(-2 + y12)**2) - 2*MPL.G([-a0], y12) - 4*MPL.G([a0/(-1 + a0)], y12)) - 
   (11*MPL.G([a0/(-1 + a0)], y12))/3 - (4*(-1 + y12)*MPL.G([(2*a0)/(-1 + a0)], y12))/(3*(-2 + y12)**2) - 
   4*MPL.G([0, 0], y12) + 2*MPL.G([0, -a0], y12) + 4*MPL.G([0, a0/(-1 + a0)], y12) + (-8 + 16/y12)*MPL.G([1, 0], y12) + 
   (8*(-2 + y12)*MPL.G([1, a0/(-1 + a0)], y12))/y12 + 2*MPL.G([-a0, 0], y12) + 4*MPL.G([a0/(-1 + a0), 0], y12) - 
   4*MPL.G([a0/(-1 + a0), a0/(-1 + a0)], y12))

    @staticmethod
    def CqgFF_Finite_Gabor_DIVJAC_NOD0(a0,y12in):
        '''Final-final collinear integrated counterterm for qg canonically normalized for the rescaling mapping'''
        if abs(y12in-1.) < 1e-6:
            y12 = 1-1e-6
        else:
            y12 = y12in        
        return ((7 - math.pi**2)/2 + 3*MPL.G([0], y12)**2 + MPL.G([0], a0)*(-3 + (4 - 8/y12)*MPL.G([1], y12)) + 
   MPL.G([0], y12)*(-MPL.G([-a0], y12) - 2*MPL.G([a0/(-1 + a0)], y12)) - (3*MPL.G([a0/(-1 + a0)], y12))/2 - 
   2*MPL.G([0, 0], y12) + MPL.G([0, -a0], y12) + 2*MPL.G([0, a0/(-1 + a0)], y12) + (-4 + 8/y12)*MPL.G([1, 0], y12) + 
   (4 - 8/y12)*MPL.G([1, a0/(-1 + a0)], y12) + MPL.G([-a0, 0], y12) + 2*MPL.G([a0/(-1 + a0), 0], y12) - 
   2*MPL.G([a0/(-1 + a0), a0/(-1 + a0)], y12))

    @staticmethod
    def SoftFF_Finite_Gabor_DIVJAC_NOD0(y0,Yin):
        '''Final-final soft+soft-colinear integrated counterterm canonically normalized for the soft rescaling mapping'''
        if abs(Yin-1.) < 1e-6:
            Y = 1-1e-6
        else:
            Y = Yin        
        return (-math.pi**2/6 + 2*(y0 + MPL.G([0], Y)*(y0 - MPL.G([0], y0))) - MPL.G([0, 0], Y) + MPL.G([1, 0], Y))

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

