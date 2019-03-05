from madgraph import MG5DIR
import ctypes
import logging
import math

logger=logging.getLogger(__name__)

class GinacMPL(object):
    """ MPL is a class to numerically evaluate multiple polylogarithms
Usage: MPL.G([a_1,...,a_n],x)
Output: float
The output is equal to G(a_1,...,a_n,x) evaluated withour linked Ginac code """

    # establish the interface with ginacg
    # try:
    #     os.environ["LD_LIBRARY_PATH"]=MG5DIR + "/HEPTools/lib:"+os.environ["LD_LIBRARY_PATH"]
    # except KeyError:
    #     os.environ["LD_LIBRARY_PATH"] = MG5DIR + "/HEPTools/lib"
    # print os.environ["LD_LIBRARY_PATH"]

    _interface_established=False
    _interface_so_path = MG5DIR + "/HEPTools/lib/ginac_mg5_interface.so"
    @classmethod
    def initialize(cls):
        if cls._interface_established:
            text = "Ginac MPL interface already established"
            logger.info(text)
            return

        cls._ginacG = ctypes.CDLL(cls._interface_so_path)
        cls._ginacG.GinacG.argtypes = (ctypes.POINTER(ctypes.c_float), ctypes.c_float, ctypes.c_int)
        cls._ginacG.GinacG.restype = ctypes.c_float
        text= "Ginac MPL interface established"
        logger.info(text)
        cls._interface_established = True

    @classmethod
    def G(cls,l, x, *args, **opts):
        if not cls._interface_established:
            cls.initialize()
        w = len(l)
        array_type = ctypes.c_float * w
        return cls._ginacG.GinacG(array_type(*l), ctypes.c_float(x), ctypes.c_int(w))

    @staticmethod
    def G_conditions(l,x):
        return True
