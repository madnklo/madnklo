################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################

from __future__ import division
import os
import logging
import math
import shutil
import numpy as np

try:
    import madgraph
except ImportError:
    MADEVENT = True
    import internal.banner as bannermod
    import internal.misc as misc
    import internal.files as files
    import internal.cluster as cluster
    import internal.lhe_parser as lhe_parser
    import internal.misc_integrator as misc_integrator
    import internal.histograms as histograms
else:
    MADEVENT= False
    import madgraph.various.misc as misc
    import madgraph.iolibs.files as files
    import madgraph.various.cluster as cluster
    import madgraph.various.lhe_parser as lhe_parser
    import madgraph.integrator.misc_integrator as misc_integrator
    import madgraph.various.histograms as histograms

logger = logging.getLogger('madgraph.integrator')
pjoin = os.path.join

class VirtualObservable(object):
    """A mother base class that specified mandatory feature that any observable should implement."""

    def __init__(self, name = 'default'):
        self.name = name
        pass    

    def __call__(self, wgt, *args, **opts):
        """ Integrand function call, with list of continuous and discrete input values for all dimensions."""
        assert(isinstance(wgt, float))
        
        # A unique float must be returned, but additional return values can 
        # be specified in the body of the function or via options.
        return True

class ObservableList(list):
    """A mother base class that specified mandatory feature that any observable list should implement."""

    def __init__(self, *args, **opts):
        super(ObservableList, self).__init__(*args, **opts)

    def apply_observables(self, wgt, *args, **opts):
        """ Apply observables."""
        for obs in self:
            obs(wgt,*args, **opts)

    def append(self, arg, **opts):
        """ Type-checking. """
        assert(isinstance(arg, VirtualObservable))
        super(ObservableList, self).append(arg, **opts)


##########################################################
#                     HwU filling tools
##########################################################

class ObservableFunctions(object):
    """A compendium of basic observables implemented as static methods
    Make your class inherit from this to have access to the library
    """
    @staticmethod
    def inclusive_xsec(*args,**kwargs):
        """Total cross section. No need to bin data by kinematics"""
        return 1

    @staticmethod
    def scalar_pt_sum(data_for_observables,*args,**kwargs):
        """Sum of the transverse momentum of all particles in the final state """
        PS_point = data_for_observables['PS_point'].to_list()
        flavors = data_for_observables['flavors']
        pt = 0
        for p in PS_point[len(flavors[0]):]:
            pt+=p.pt()
        return pt



class HwUObservable(VirtualObservable,ObservableFunctions):
    """Class that creates and fills in a HwU histogram for an observable given as a function"""

    def __init__(self, name='default', observable_function = None, range = [0,2], nbins = 1):
        super(HwUObservable, self).__init__(name)
        if observable_function is None:
            self.observable_function = self.inclusive_xsec
        else:
            self.observable_function = observable_function
        self.range = [float(range[0]),float(range[1])]
        assert self.range[0] < self.range[1]
        self.nbins = int(nbins)
        self.bin_size = (self.range[1]-self.range[0])/float(nbins)
        assert self.bin_size > 0

        self.create_HwU()

    def create_HwU(self):
        bins = histograms.BinList(bin_range=self.range+[self.bin_size])
        self.HwU = histograms.HwU(title=self.name,bins=bins)

    def __call__(self, wgt, *args, **kwargs):
        assert super(HwUObservable, self).__call__(wgt, *args, **kwargs)
        observable = self.observable_function(*args, **kwargs)
        self.HwU.addEvent(observable, wgt)




class HwUObservableList(ObservableList):
    def __init__(self,*args,**opts):
        super(HwUObservableList, self).__init__(*args,**opts)
        assert all([isinstance(x,HwUObservable) for x in self])
        self.create_HwUList()

    def create_HwUList(self):
        HwUs=[]
        for obs in self:
            HwUs.append(obs.HwU)
        self.HwUList = histograms.HwUList(HwUs)

    def append(self, arg, **opts):
        """ Type-checking. """
        assert(isinstance(arg, HwUObservable))
        super(ObservableList, self).append(arg, **opts)

    def output(self,path,*args,**opts):
        self.HwUList.output(path,*args,**opts)
