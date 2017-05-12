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
import math
import logging
import numpy as np
import subprocess

root_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(root_path, os.path.pardir, os.path.pardir))

logger = logging.getLogger(' MadEvent7 ')
logger.setLevel(logging.INFO)
logging.basicConfig()

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

pjoin = os.path.join

# Import the matrix elements
sys.path.append(pjoin(root_path,'MEs_for_epem_lplma','SubProcesses'))
try:
    import P1_epem_mupmuma.matrix2py as epem_mupmuma
    import P2_epem_epema.matrix2py as epem_epema
    import P3_epem_mupmum_no_z.matrix2py as epem_mupmum_no_z
except ImportError:
    logger.info("Attempting to recompile the matrix element python modules with f2py...")
    devnull = open(os.devnull,'w')
    subprocess.call(['./compile_MEs.sh'],shell=True,cwd=root_path, stdout=devnull)
    try:
        import P1_epem_mupmuma.matrix2py as epem_mupmuma
        import P2_epem_epema.matrix2py as epem_epema
    except:
        msg = ["Could not import the following modules:"]
        msg.append("  -> %s"%os.path.join(root_path,'MEs_for_epem_lplma','SubProcesses','P1_epem_mupmuma','matrix2py.so'))
        msg.append("  -> %s"%os.path.join(root_path,'MEs_for_epem_lplma','SubProcesses','P2_epem_epema','matrix2py.so'))
        msg.append("Try again after having compiled these libraries with f2py manually.")
        logger.critical() 
        sys.exit(0)

class ME_function(functions.VirtualFunction):

    def __init__(self, process, *args, **opts):
        assert(process in ['e+ e- > e+ e- a', 'e+ e- > mu+ mu- a', 'e+ e- > mu+ mu- / z'])

        ## Setup access to f2py'ed MEs
        if process == 'e+ e- > e+ e- a':
            self.ME_accessor = epem_epema
        elif process == 'e+ e- > mu+ mu- a':
            self.ME_accessor = epem_mupmuma
        elif process == 'e+ e- > mu+ mu- / z':
            self.ME_accessor = epem_mupmum_no_z
        else:
            logger.critical('Should not happend.')
            sys.exit(0)

        self.ME_accessor.initialise('./Cards/param_card.dat')        

        super(ME_function, self).__init__(*args, **opts)

    @staticmethod
    def format_momenta_for_f2py(p):
        """ fortran/C-python do not order table in the same order.
        Also, remove the mass component of the momenta. """
        new_p = []
        for i in range(4):
            new_p.append([0]*len(p))
        for i, onep in enumerate(p):
            for j, x in enumerate(onep):
                if j==4: continue
                new_p[j][i] = x
        return new_p

    def __call__(self, PS_point, **opts):
        """ ME_function call , for now with only the PS point specified """

        out_data = {'weight': 0.0}

        # Alpha_s is irrelevant here
        alpha_s = 0.118
        # We sum explicitely over helicities for now
        nhel = 0
         
        # Call the f2py'ed MEs
        out_data['weight'] = self.ME_accessor.get_me(
            self.format_momenta_for_f2py(PS_point), alpha_s, nhel)

        # Return a data structure with 'weight' in it.
        return out_data

class epem_lplma_integrand(integrands.VirtualIntegrand):

    def __init__(self, process, phase_space_generator, *args, **opts):
        """ Instantiate the integrand for a specific process."""

        assert(process in ['e+ e- > e+ e- a', 'e+ e- > mu+ mu- a', 'e+ e- > mu+ mu- / z'])
        self.process = process

        default_opts = {'E_cm': 1000.0}
        # Set instance attributes options
        for opt in default_opts:
            if opt in opts:
                setattr(self,opt,opts.pop(opt))
            else:
                setattr(self,opt,default_opts[opt])

        super(epem_lplma_integrand, self).__init__(*args, **opts)
        
        self.function_list.append(ME_function(process))
        
        self.phase_space_generator = phase_space_generator

    @staticmethod
    def Lambda(s,sqrMA,sqrMB):
        """ Kahlen function."""

        return s**2 + sqrMA**2 + sqrMB**2 - \
             2.*s*sqrMA - 2.*sqrMB*sqrMA - 2.*s*sqrMB

    def pass_cuts(self, PS_point):
        """ Implementation of the isolation cuts. """
       
#       Uncomment below to debug cuts
##        misc.sprint( "\nProcessings cuts for process '%s' and PS point:\n\n%s"%( self.process,
##                phase_space_generators.FlatInvertiblePhasespace.nice_momenta_string(
##                            PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial) ))

        pt_cut = 100.0
        dr_cut = 0.4
        
        # Apply a uniform P_t cut first
        for p in PS_point[2:]:
            # misc.sprint(p.pt())            
            if p.pt() < pt_cut:
                return False

        for i, p1 in enumerate(PS_point[2:]):
            for j, p2 in enumerate(PS_point[2+i+1:]):
                # misc.sprint(2+i,2+i+1+j,p1.deltaR(p2))
                if p1.deltaR(p2) < dr_cut:
                    return False

        return True
    
    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Integrand function call, with list of continuous and discrete input values for all dimensions."""
        
        assert(self.check_input_types(continuous_inputs, discrete_inputs))
        assert(len(discrete_inputs)==
               len([1 for d in self.discrete_dimensions if not d.folded]))
        assert(len(continuous_inputs)==
               len([1 for d in self.continuous_dimensions if not d.folded]))
         
        # A unique float must be returned
        wgt = 1.0
         
        # Random variables sent
        random_variables = list(continuous_inputs)
        # Avoid extremum since the phase-space generation algorithm doesn't like it
        epsilon = 1e-10
        random_variables = [min(max(rv,epsilon),1.-epsilon) for rv in random_variables]

        # Generate a PS point
        PS_point, PS_weight = self.phase_space_generator.generateKinematics(self.E_cm, random_variables)
       
        # Apply cuts
        if not self.pass_cuts(PS_point):
            return 0.0

        # Account for PS weight
        wgt *= PS_weight

        # Call the f2py'ed MEs
        ME_outputs = [ME_function(PS_point) for ME_function in self.function_list]
        wgt *= sum(o['weight'] for o in ME_outputs)
       
#       Uncomment below to debug MEs
##        misc.sprint(random_variables)
##        misc.sprint( "\nEvaluation for process '%s' and following PS point:\n\n%s"%( self.process,
##                phase_space_generators.FlatInvertiblePhasespace.nice_momenta_string(
##                            PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial) ))
##        misc.sprint("\n ME weights = [ %s ]\n"%(', '.join('%.16e'%o['weight'] for o in ME_outputs)))
##        sys.exit(0)

        ##################################
        # Post-processing of the weight
        ##################################

        # Include flux factor
        flux = 1.
        if self.phase_space_generator.n_initial == 2:
            flux = 1. / (2.*math.sqrt(self.Lambda(self.E_cm**2, 
                            self.phase_space_generator.initial_masses[0],
                            self.phase_space_generator.initial_masses[1])))
        elif self.phase_space_generator.n_initial == 1:
            flux = 1. / (2.*self.E_cm)
        flux /= math.pow(2.*math.pi, 3*self.phase_space_generator.n_final - 4)
        wgt *= flux

        # Convert GeV^-2 to picobarns
        wgt *= 0.389379304e9

        if self.apply_observables:
            data_for_observables = {'PS_point': PS_point, 'ME_outputs' : ME_outputs}
            self.observable_list.apply_observables(continuous_inputs, discrete_inputs, data_for_observables)
        
        return wgt

class MadEvent7_driver(object):

    def __init__(self, E_cm = 1000.0,
                 process = 'e+ e- > l+ l- a',
                 integrator_class=integrators.SimpleMonteCarloIntegrator, 
                 integrator_options={}):
        
        assert(process in ['e+ e- > l+ l- a','e+ e- > mu+ mu- / z'])
        self.process = process
        self.n_final = None
        if self.process == 'e+ e- > l+ l- a':
            self.n_final = 3
        elif self.process == 'e+ e- > mu+ mu- / z':
            self.n_final = 2
        self.E_cm = E_cm
        self.phase_space_generator = phase_space_generators.FlatInvertiblePhasespace([0.]*2, [0.]*self.n_final)
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

        if self.process == 'e+ e- > l+ l- a':
            epem_mupmuma_integrand = epem_lplma_integrand('e+ e- > mu+ mu- a', self.phase_space_generator, PS_dimensions, E_cm = self.E_cm)
            epem_epema_integrand = epem_lplma_integrand('e+ e- > e+ e- a', self.phase_space_generator, PS_dimensions, E_cm = self.E_cm)
            self.integrands = [epem_mupmuma_integrand, epem_epema_integrand]
        
        else:
            integrand_for_process = epem_lplma_integrand(self.process, self.phase_space_generator, PS_dimensions, E_cm = self.E_cm)
            self.integrands = [integrand_for_process, ]

    def build_dummy_integrand(self):
        """ Build the dummy integrand for testing purposes. """

        dummy_integrand = integrands.VirtualIntegrand()

        # Define its dimensions
        dummy_integrand.set_dimensions( integrands.DimensionList([
            integrands.ContinuousDimension('x' ,lower_bound=2.0, upper_bound=5.0),
            integrands.ContinuousDimension('y' ,lower_bound=1.0, upper_bound=3.0)
             ]) )
                    
        # Define its constituting functions
        dummy_integrand.function_list = functions.FunctionList([
                functions.FunctionFromPythonExpression('math.sin(x/y)', dimensions=my_integrand.get_dimensions()),
                functions.FunctionFromPythonExpression('math.cos(x*y)', dimensions=my_integrand.get_dimensions())
            ])

        self.integrands = [dummy_integrand, dummy_integrand]

    def run(self):
        logger.info("Now Integrating '%s' with %s"%(self.process, self.integrator.get_name()))

        xsec, error = self.integrator.integrate()
        print ''
        logger.info("="*100)
        logger.info('{:^100}'.format("\033[92mCross-section for process '%s' and integrator '%s':\033[0m"
                                                             %(self.process, self.integrator.get_name())))
        logger.info('{:^100}'.format("\033[94m%.5e +/- %.2e [pb]\033[0m"%(xsec, error)))
        logger.info("="*100+"\n")


# ----------------------------------------------------------------------------------------------
#                                  #########################
#                                  # Steering of MadEvent7 #
#                                  #########################
# ----------------------------------------------------------------------------------------------


# Use 0 for minimal verbosity and 2 for talkative integrators
verbosity  = 2

if verbosity == 0:
    logger.info("Set verbosity to a value larger than 0 to watch more closely the integration progress.")

# Try with a naive MonteCarloFirst
parameters = {
    'process'            : 'e+ e- > l+ l- a',
#    'process'            : 'e+ e- > mu+ mu- / z',
    'E_cm'               : 1000.0
}

ME_references_values = { 'e+ e- > l+ l- a' : (0.07329, 9.583e-05),
                         'e+ e- > mu+ mu- / z': (0.09005, 4.84e-05)}

if parameters['process'] in ME_references_values:
    xsec, err = ME_references_values[parameters['process']]
    logger.info("MadEvent6 reference result for process '%s' = \033[92m%.5e +/- %.2e [pb]\033[0m "%(parameters['process'], xsec, err) )

# Initialize a MadEvent driver
ME7 = MadEvent7_driver( **parameters )

# Now launch the integration (or skip it by commenting below if not interested in the stupid integrator.)
## ME7.run()

# Try several integrators
integrators = {

   'Naive' : (integrators.SimpleMonteCarloIntegrator, {  'n_iterations'            : 10,
                                                         'n_points_per_iterations' : 1000,
                                                         'accuracy_target'         : None,
                                                         'verbosity'               : verbosity } ),

   'VEGAS' : (pyCubaIntegrator.pyCubaIntegrator, { 'algorithm' : 'Vegas', 
                                                   'verbosity' : verbosity,
                                                   'seed'      : 3,
                                                   'target_accuracy' : 1.0e-3,
                                                   'n_start'   : 10000,
                                                   'n_increase': 5000,
                                                   'n_batch'   : 1000,
                                                   'max_eval'  : 100000,
                                                   'min_eval'  : 0}),
   
   'SUAVE'   : (pyCubaIntegrator.pyCubaIntegrator, { 'algorithm' :'Suave', 
                                                     'verbosity' : verbosity,
                                                     'target_accuracy' : 1.0e-3,
                                                     'max_eval'  : 100000,
                                                     'min_eval'  : 0 } ),
  
   'DIVONNE' : (pyCubaIntegrator.pyCubaIntegrator, { 'algorithm' : 'Divonne', 
                                                     'verbosity': verbosity,
                                                     'target_accuracy' : 1.0e-5,
                                                     'max_eval'  : 100000000,
                                                     'min_eval'  : 0 } ),

   'CUHRE'   : (pyCubaIntegrator.pyCubaIntegrator, { 'algorithm' : 'Cuhre',
                                                     'verbosity' : verbosity,
                                                     'target_accuracy' : 1.0e-3,
                                                     'max_eval'  : 100000,
                                                     'min_eval'  : 0 } ),
}

# Now choose which integrator(s) to try. 
#for name in ['Naive','VEGAS','SUAVE','DIVONNE','CUHRE']:
for name in ['DIVONNE']:
    # Set the integrator
    ME7.set_integrator(*integrators[name])
    # Run one last time
    ME7.run()

