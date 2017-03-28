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
import random

try:
    import madgraph
except ImportError:
    MADEVENT = True
    import internal.banner as bannermod
    import internal.misc as misc
    import internal.files as files
    import internal.cluster as cluster
    import internal.lhe_parser as lhe_parser
    import internal.functions as functions
    import internal.observables as observables
else:
    MADEVENT= False
    import madgraph.various.misc as misc
    import madgraph.iolibs.files as files
    import madgraph.various.cluster as cluster
    import madgraph.various.lhe_parser as lhe_parser
    import madgraph.integrator.functions as functions
    import madgraph.integrator.observables as observables

logger = logging.getLogger('madgraph.integrator')
pjoin = os.path.join

class VirtualIntegrand(object):
    """A mother base class that specifies the feature that any integrand should implement."""
    
    @classmethod
    def check_input_types(continuous_inputs, discrete_inputs):
        """ Wrapper to easily check input types."""
        assert(isinstance(discrete_inputs, np.ndarray))
        assert(isinstance(i, np.int64) for i in discrete_inputs)
        assert(isinstance(continous_inputs, np.ndarray))
        assert(isinstance(f, numpy.float64) for f in continuous_inputs)        
        return True

    def __init__(self, dimensions=DimensionList()):
        self.dimensions                 = dimensions
        self.apply_observables          = True
        self.observable_list            = VirtualObservableList()
        self.function_list              = VirtualFunctionList()
        pass
    
    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Integrand function call, with list of continuous and discrete input values for all dimensions."""
        assert(check_input_types(continuous_inputs, discrete_inputs))
        assert(len(discrete_inputs)==
               len([_ for d in self.dimensions.get_discrete_dimensions() if not d.folded]))
        assert(len(continuous_inputs)==
               len([_ for d in self.dimensions.get_continuous_dimensions() if not d.folded]))
        
        # A unique float must be returned
        wgt = 0.0
        
        # This is the hearth of where the complication for building the integrand lie.
        # This virtual integrator implements the simplest possible case-use
        wgt = sum(f(continuous_inputs, discrete_inputs) for f in self.function_list)
        
        if self.apply_observables:
            self.observable_list.apply_observables(continuous_inputs, discrete_inputs, wgt)

        return wgt

class Dimension(object):
    """ A dimension object specifying a specific integration dimension."""
    
    def __init__(self, name, folded=False):
        self.name   = name
        self.folded = folded
    
    def length(self):
        raise NotImplemented
    
    def random_sample(self):        
        raise NotImplemented

class DiscreteDimension(Dimension):
    """ A dimension object specifying a specific discrete integration dimension."""
    
    def __init__(self, name, values, **opts):
        try:
            self.normalized = opts.pop('normalized')
        except:
            self.normalized = False
        super(Dimension, self).__init__(name, **opts)
        assert(isinstance(values, list))
        self.values = values
    
    def length(self):
        if normalized:
            return 1.0/float(len(values))
        else:
            return 1.0
    
    def random_sample(self):
        return np.int64(random.choice(values))
        
class ContinuousDimension(Dimension):
    """ A dimension object specifying a specific discrete integration dimension."""
    
    def __init__(self, name, lower_bound=0.0, upper_bound=1.0):
        super(Dimension, self).__init__(name, **opts)
        assert(upper_bound>lower_bound)
        self.lower_bound  = lower_bound
        self.upper_bound  = upper_bound 

    def length(self):
        return (self.upper_bound-self.lower_bound)

    def random_sample(self):
        return np.float64(lower_bound+random.random()*(upper_bound-lower_bound))

class DimensionList(list):
    """A DimensionList."""

    def __init__(self, *args, **opts):
        super(list, self).__init__(self, *args, **opts)

    def volume(self):
        """ Returns the volue of the complete list of dimensions."""
        vol = 1.0
        for d in self:
            vol *= d.length()
        return vol
    
    def append(self, arg, **opts):
        """ Type-checking. """
        assert(isinstance(arg, Dimension))
        super(list, self).append(self, arg, **opts)
        
    def get_discrete_dimensions(self):
        """ Access all discrete dimensions. """
        return DimensionList(d for d in self if isinstance(d, DiscreteDimension))
    
    def get_continuous_dimensions(self):
        """ Access all discrete dimensions. """
        return DimensionList(d for d in self if isinstance(d, ContinuousDimension))
    
    def random_sample(self):
        return np.array(d.random_sample() for d in self)