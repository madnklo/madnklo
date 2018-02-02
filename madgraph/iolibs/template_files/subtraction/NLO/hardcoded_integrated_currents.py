# =============================================
#               MPL CLASS
# =============================================
import ctypes
import logging
import math

logger=logging.getLogger("SubtractionCurrents")

class MPL(object):
    """ MPL is a class to numerically evaluate multiple polylogarithms
Usage: MPL.G([a_1,...,a_n],x)
Output: float
The output is equal to G(a_1,...,a_n,x) evaluated withour linked Ginac code """

    # establish the interface with ginacg
    _ginacG = ctypes.CDLL("ginacg.so")
    _ginacG.GinacG.argtypes = (ctypes.POINTER(ctypes.c_float), ctypes.c_float, ctypes.c_int)
    _ginacG.GinacG.restype = ctypes.c_float
    text= "Ginac MPL interface established"
    logger.info(text)

    @classmethod
    def G(cls,l, x, *args, **opts):
        w = len(l)
        array_type = ctypes.c_float * w
        return cls._ginacG.GinacG(array_type(*l), ctypes.c_float(x), ctypes.c_int(w))


# =============================================
#          Harcoded Expression Class (HE)
# =============================================

class HE(object):
    """
    This is a class whose static methods are the hardcoded expressions for the integrated currents
    """
    @staticmethod
    def CqqFF_Finite_Gabor_EEJJ(a0,y12):
        '''Final-final collinear integrated counterterm for qq for m=2'''
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
    def CggFF_Finite_Gabor_EEJJ(a0,y12):
        '''Final-final collinear integrated counterterm for gg for m=2'''
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
    def CqgFF_Finite_Gabor_EEJJ(a0,y12):
        '''Final-final collinear integrated counterterm for qg for m=2'''
        return (5 - math.pi**2/3 + 3/(2*(-1 + y12)) - (3*y12)/(2*(-1 + y12)) - (3*MPL.G([0], a0))/2 - (3*MPL.G([0], y12))/2 + 
   MPL.G([0], a0)*MPL.G([0], y12) + MPL.G([0], a0)*MPL.G([-a0], y12) - MPL.G([0], y12)*MPL.G([-y12], a0) + 
   (3*(a0 + y12 - a0*y12)*MPL.G([y12/(-1 + y12)], a0))/(2*a0*(-1 + y12)**2) - 
   ((3*y12 - 4*a0*MPL.G([0], y12) + 4*a0*y12*MPL.G([0], y12))*MPL.G([y12/(-1 + y12)], a0))/(2*a0*(-1 + y12)**2) + 
   MPL.G([0, 0], a0) + MPL.G([0, 0], y12) - 2*MPL.G([0, -y12], a0) - (2*MPL.G([0, y12/(-1 + y12)], a0))/(-1 + y12) + 
   (2*y12*MPL.G([0, y12/(-1 + y12)], a0))/(-1 + y12) + MPL.G([-a0, 0], y12) - MPL.G([-y12, 0], a0) + 
   (2*MPL.G([y12/(-1 + y12), 0], a0))/(-1 + y12) - (2*MPL.G([y12/(-1 + y12), y12/(-1 + y12)], a0))/(-1 + y12))

    @staticmethod
    def SoftFF_Finite_Gabor_EEJJ(y0,Y):
        '''Final-final soft+soft-colinear integrated counterterm for qg for m=2'''
        return (-math.pi**2/6 + 2*(y0 + MPL.G([0], Y)*(y0 - MPL.G([0], y0))) - MPL.G([0, 0], Y) + MPL.G([1, 0], Y))