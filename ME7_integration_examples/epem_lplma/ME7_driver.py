#! /usr/bin/env python

################################################################################
#
# This is an example of how one can use the new MadEvent7 integration framework
# to integrate the process e+ e- > l+ l- a. Everything is tailored-hardcoded
# for that process for now.
#
################################################################################

import sys
import os
import logging

root_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(root_path, os.path.pardir, os.path.pardir))

try:
    import madgraph
except ImportError:
    MADEVENT = True
    import internal.misc as misc
    import internal.integrators as integrators
    import internal.pyCubaIntegrator as pyCubaIntegrator
    import internal.integrands as integrands
    import internal.functions as functions
    import internal.phase_space_generators as phase_space_generators

else:
    MADEVENT= False
    import madgraph.various.misc as misc
    import madgraph.integrator.integrators as integrators
    import madgraph.integrator.pyCubaIntegrator as pyCubaIntegrator
    import madgraph.integrator.integrands as integrands
    import madgraph.integrator.functions as functions
    import madgraph.integrator.phase_space_generators as phase_space_generators

logger = logging.getLogger(' MadEvent7 @ e+ e- > l+ l- a ')
logger.setLevel(logging.INFO)
logging.basicConfig()

pjoin = os.path.join

class MadEvent7_driver(object):

    def __init__(self, E_cm = 1000.0, integrator_class=integrators.SimpleMonteCarloIntegrator, integrator_options={}):

        self.E_cm = E_cm
        self.ps_generator = phase_space_generators.FlatInvertiblePhasespace([0.]*2, [0.]*3)
        self.build_integrand()
        self.set_integrator(integrator_class, integrator_options)

    def set_integrator(self, integrator_class=integrators.SimpleMonteCarloIntegrator, integrator_options={}):
        """ Assign a particular integrator."""
        self.integrator = integrator_class(self.integrands, **integrator_options)

    def build_integrand(self):
        """ Build the integrand for the process e+ e- > l+ l- a """

        my_integrand = integrands.VirtualIntegrand()

        # Define its dimensions
        my_integrand.set_dimensions( integrands.DimensionList([
            integrands.ContinuousDimension('x',lower_bound=2.0, upper_bound=5.0),
            integrands.ContinuousDimension('y' ,lower_bound=1.0, upper_bound=3.0)
             ]) )
                    
        # Define its constituting functions
        my_integrand.function_list = functions.FunctionList([
                functions.FunctionFromPythonExpression('math.sin(x/y)', dimensions=my_integrand.get_dimensions()),
                functions.FunctionFromPythonExpression('math.cos(x*y)', dimensions=my_integrand.get_dimensions())
            ])

        self.integrands = [my_integrand, my_integrand]

    def run(self):
        logger.info("Starting integration of 'e+ e- > l+ l- a' with %s"%self.integrator.get_name())

        xsec, error = self.integrator.integrate()
        logger.info("Inclusive cross-section: %f +/- %f"%(xsec, error))

#####################
# Steering MadEvent #
#####################

# Try with a naive MonteCarloFirst
verbosity  = 0
parameters = {
    'E_cm'               : 1000.0,
    'integrator_class'   : integrators.SimpleMonteCarloIntegrator,
    'integrator_options' : {  'n_iterations'            : 10,
                              'n_points_per_iterations' : 1000,
                              'accuracy_target'         : None,
                              'verbosity'               : verbosity }
    }

# Initialize a MadEvent driver
ME7 = MadEvent7_driver( **parameters )

# Now launch the integration
ME7.run()

# Change to Vegas integrator for example 
ME7.set_integrator(pyCubaIntegrator.pyCubaIntegrator, {'algorithm':'Vegas', 'verbosity': verbosity})
# And now run again.
ME7.run()

# Finally change to Cuhre for good measure 
ME7.set_integrator(pyCubaIntegrator.pyCubaIntegrator, {'algorithm':'Cuhre', 'verbosity': verbosity})
# And now run again.
ME7.run()
