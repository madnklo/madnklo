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
    import internal.functions as functions

else:
    MADEVENT= False
    import madgraph.various.misc as misc
    import madgraph.iolibs.files as files
    import madgraph.various.cluster as cluster
    import madgraph.various.lhe_parser as lhe_parser
    import madgraph.integrator.integrands as integrands
    import madgraph.integrator.functions as functions

logger = logging.getLogger('madgraph.integrator')
pjoin = os.path.join

class VirtualIntegrator(object):
    """A mother base class that specifies the feature that any integrator should implement."""
    
    def __init__(self, integrand, **opts):
        self.integrand = integrand
        
    def integrate(self):
        """ Return the final integral and error estimators. """
        raise NotImplementedError
    
class SimpleMonteCarloIntegrator(VirtualIntegrator):
    """ Implements the simplest Monte-Carlo integrator. """
    
    def __init__(self, integrand, **opts):
        """ Initialize the simplest MC integrator."""
        try:
            self.accuracy_target = opts.pop('accuracy_target')
        except KeyError:
            self.accuracy_target = 0.01
        try:
            self.n_iterations = opts.pop('n_iterations')
        except KeyError:
            self.n_iterations = None
        try:
            self.n_points_per_iterations = opts.pop('n_points_per_iterations')
        except KeyError:
            self.n_points_per_iterations = 100
        super(VirtualIntegratorm,self).__init__(integrand, **opts)
        
    def integrate(self):
        """ Return the final integral and error estimates."""
        
        iteration_number   = 0
        integral_estimate  = 0.0
        error_estimate     = 0.0
           
        while (self.n_iterations is None or iteration_number <= self.n_iterations ) and \
              (self.accuracy_target is None or self.error_estimate < self.accuracy_target):   
            n_points = 0
            
            while n_points <= self.n_points_per_iterations:
                n_points += 1
                discrete_dimensions = self.integrand.dimensions.get_discrete_dimensions().randome_sample()
                continuous_dimensions = self.integrand.dimensions.get_continuous_dimensions().randome_sample()
                new_wgt = self.integrand(discrete_dimensions,continuous_dimensions)
                self.integral_estimate += new_wgt
                self.error_estimate += new_wgt**2
            
            integral_estimate /= n_points
            error_estimate = (math.sqrt(self.error_estimate) / n_points) 
        
        phase_space_volume = self.integrand.dimensions.volume()
        return phase_space_volume*integral_estimate, phase_space_volume*error_estimate
            
                
if __name__ == "__main__":

    # Example usage of the new integrator framework
    
    # First build the integrand
    my_integrand = integrands.VirtualIntegrand()
    # Define its dimensions
    my_integrand.dimensions = integrands.DimensionList([
            integrands.ContinuousDimension('x'),
            integrands.ContinuousDimension('y' ,lower_bound=0.0, upper_bound=3.0)
             ])
                    
    # Define its constituting functions
    my_integrand.function_list = functions.FunctionList([
                functions.FunctionFromPythonExpression('math.sin(x/y)', dimensions=my_integrand.dimensions),
                functions.FunctionFromPythonExpression('math.tan(x*y)', dimensions=my_integrand.dimensions)
            ])
    
    # Then the integrator
    my_integrator = SimpleMonteCarloIntegrator(my_integrand)
    
    # Finally integrate
    print my_integrator.integrate()