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
    import internal.integrands as integrands
else:
    MADEVENT= False
    import madgraph.various.misc as misc
    import madgraph.iolibs.files as files
    import madgraph.various.cluster as cluster
    import madgraph.various.lhe_parser as lhe_parser
    import madgraph.integrator.integrands as integrands


logger = logging.getLogger('madgraph.integrator')
pjoin = os.path.join

class VirtualObservable(object):
    """A mother base class that specified mandatory feature that any observable should implement."""

    def __init__(self):
        self.name = 'default'
        pass    

    def __call__(self, continuous_inputs, discrete_inputs, wgt, **opts):
        """ Integrand function call, with list of continuous and discrete input values for all dimensions."""
        integrands.VirtualIntegrand.check_input_types(continuous_inputs, discrete_inputs)
        assert(isinstance(wgt, float))
        
        # A unique float must be returned, but additional return values can 
        # be specified in the body of the function or via options.
        return True

class ObservableList(list):
    """A mother base class that specified mandatory feature that any observable list should implement."""

    def __init__(self, *args, **opts):
        super(list, self).__init__(self, *args, **opts)

    def apply_observables(self, continuous_inputs, discrete_inputs, wgt, **opts):
        """ Apply observables."""
        integrands.VirtualIntegrand.check_input_types(continuous_inputs, discrete_inputs)
        for obs in self:
            obs(continuous_inputs, discrete_inputs, wgt)

    def append(self, arg, **opts):
        """ Type-checking. """
        assert(isinstance(arg, VirtualObservable))
        super(list, self).append(self, arg, **opts)
    
    
    