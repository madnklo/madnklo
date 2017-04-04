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
import numpy as np

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

class ME_function(functions.VirtualFunction):

    def __init__(self, process, *args, **opts):
        assert(process in ['e+ e- > e+ e- a', 'e+ e- > mu+ mu- a'])

        ## Setup access to f2py'ed MEs

        super(ME_function, self).__init__(*args, **opts)

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Integrand function call, with list of continuous and discrete input values for all dimensions."""
        self.check_input_types(continuous_inputs, discrete_inputs)

        ## Call the f2py'ed MEs        

        # Return a data structure with 'weight' in it.
        return {'weight': 0.0}

class epem_lplma_integrand(integrands.VirtualIntegrand):

    def __init__(self, process, phase_space_generator, *args, **opts):
        assert(process in ['e+ e- > e+ e- a', 'e+ e- > mu+ mu- a'])        
        self.process = process
        self.phase_space_generator = phase_space_generator

        default_opts = {'E_cm': 1000.0}
        # Set instance attributes options
        for opt in default_opts:
            if opt in opts:
                setattr(self,opt,opts.pop(opt))
            else:
                setattr(self,opt,default_opts[opt])

        super(epem_lplma_integrand, self).__init__(*args, **opts)

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Integrand function call, with list of continuous and discrete input values for all dimensions."""
        
        assert(self.check_input_types(continuous_inputs, discrete_inputs))
        assert(len(discrete_inputs)==
               len([1 for d in self.discrete_dimensions if not d.folded]))
        assert(len(continuous_inputs)==
               len([1 for d in self.continuous_dimensions if not d.folded]))
        
        # A unique float must be returned
        wgt = 0.0
        
        # Call the f2py'ed MEs
        data = [f(continuous_inputs, discrete_inputs) for f in self.function_list]
        wgt = sum(d['weight'] for d in data)
       
        # Post-process call with flux factors, PS weight , etc..

        if self.apply_observables:
            self.observable_list.apply_observables(continuous_inputs, discrete_inputs, data)

        return wgt

class MadEvent7_driver(object):

    def __init__(self, E_cm = 1000.0, integrator_class=integrators.SimpleMonteCarloIntegrator, integrator_options={}):

        self.E_cm = E_cm
        self.phase_space_generator = phase_space_generators.FlatInvertiblePhasespace([0.]*2, [0.]*3)
        self.build_integrand()
        self.set_integrator(integrator_class, integrator_options)

    def set_integrator(self, integrator_class=integrators.SimpleMonteCarloIntegrator, integrator_options={}):
        """ Assign a particular integrator."""
        self.integrator = integrator_class(self.integrands, **integrator_options)

    def build_integrand(self):
        """ Build the integrand for the process e+ e- > l+ l- a """
        
        PS_dimensions = integrands.DimensionList(
            [ integrands.ContinuousDimension('x_%d'%i,lower_bound=0.0, upper_bound=1.0) 
                            for i in range(1,self.phase_space_generator.nDimPhaseSpace()+1) ] )

        epem_mupmuma_integrand = epem_lplma_integrand('e+ e- > mu+ mu- a', self.phase_space_generator, PS_dimensions, E_cm = self.E_cm)
        epem_epema_integrand = epem_lplma_integrand('e+ e- > e+ e- a', self.phase_space_generator, PS_dimensions, E_cm = self.E_cm)
        
        self.integrands = [epem_mupmuma_integrand, epem_epema_integrand]

    def build_dummy_integrand(self):
        """ Build the dummy integrand for testing purposes. """

        dummy_integrand = integrands.VirtualIntegrand()

        # Define its dimensions
        dummy_integrand.set_dimensions( integrands.DimensionList([
            integrands.ContinuousDimension('x',lower_bound=2.0, upper_bound=5.0),
            integrands.ContinuousDimension('y' ,lower_bound=1.0, upper_bound=3.0)
             ]) )
                    
        # Define its constituting functions
        dummy_integrand.function_list = functions.FunctionList([
                functions.FunctionFromPythonExpression('math.sin(x/y)', dimensions=my_integrand.get_dimensions()),
                functions.FunctionFromPythonExpression('math.cos(x*y)', dimensions=my_integrand.get_dimensions())
            ])

        self.integrands = [dummy_integrand, dummy_integrand]

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
