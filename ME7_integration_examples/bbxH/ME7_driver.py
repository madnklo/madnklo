#! /usr/bin/env python

################################################################################
#
# This is an example of how one can use the new MadEvent7 integration framework
# to integrate the process b b~ > h (possibly @NLO). Everything is 
# tailored-hardcoded for that process for now.
#
################################################################################

import sys
import os
import time
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

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

############################################################################################
# This whole part about ME generation and loading is temporary and will be made more 
# generic in due time... 
############################################################################################

# Import the matrix elements
sys.path.append(pjoin(root_path,'bbxH_tree_MEs','SubProcesses'))
sys.path.append(pjoin(root_path,'bbxH_loop_ME','SubProcesses'))

devnull = open(os.devnull,'w')
if not os.path.isdir(pjoin(root_path,'bbxH_tree_MEs')):
    logger.info(bcolors.OKBLUE + "Generation of the b b~ > h +X tree and loop matrix elements with MG5_aMC..." + bcolors.ENDC)
    subprocess.call(['./generate_MEs.sh'], shell=True, cwd=root_path, stdout=devnull)

try:
    import P1_bbx_h.matrix2py as bbx_h 
    import P2_bbx_hg.matrix2py as bbx_hg
    import P3_gbx_hbx.matrix2py as gbx_hbx
    import P4_bg_hb.matrix2py as bg_hb
    import P5_bbx_h.matrix2py as bbx_h_loop
except ImportError:    
    logger.info(bcolors.OKBLUE + "Compilation of the matrix element python modules with f2py..." + bcolors.ENDC)
    subprocess.call(['./compile_MEs.sh'], shell=True, cwd=root_path, stdout=devnull)
    try:
        import P1_bbx_h.matrix2py as bbx_h 
        import P2_bbx_hg.matrix2py as bbx_hg
        import P3_gbx_hbx.matrix2py as gbx_hbx
        import P4_bg_hb.matrix2py as bg_hb
        import P5_bbx_h.matrix2py as bbx_h_loop
    except ZeroDivisionError:
        msg = ["Could not import some of the following modules:"]
        msg.append("  -> %s"%pjoin(root_path,'bbxH_loop_ME','SubProcesses','P5_bbx_h','matrix2py.so'))
        msg.append("  -> %s"%pjoin(root_path,'bbxH_tree_MEs','SubProcesses','P1_bbx_h','matrix2py.so'))
        msg.append("  -> %s"%pjoin(root_path,'bbxH_tree_MEs','SubProcesses','P2_bbx_hg','matrix2py.so'))
        msg.append("  -> %s"%pjoin(root_path,'bbxH_tree_MEs','SubProcesses','P3_gbx_hbx','matrix2py.so'))
        msg.append("  -> %s"%pjoin(root_path,'bbxH_tree_MEs','SubProcesses','P4_bg_hb','matrix2py.so')) 
        msg.append("Try again after having compiled these libraries with f2py manually.")
        logger.critical('\n'.join(msg)) 
        sys.exit(0)

############################################################################################

class ME_function(functions.VirtualFunction):

    def __init__(self, process, *args, **opts):

        processes = {'b b~ > h'         : bbx_h,
                     'b b~ > h @loop'   : bbx_h_loop,
                     'b b~ > h g'       : bbx_hg,
                     'g b~ > h b~'      : gbx_hbx,
                     'b g > h b'        : bg_hb}

        assert(process in processes.keys())
        
        self.ME_accessor = processes[process] 

        self.ME_accessor.initialise(pjoin(root_path,'Cards','param_card.dat'))      

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

    def __call__(self, PS_point, alpha_s, nhel, **opts):
        """ ME_function call , for now with only the PS point specified.
        nhel corresponds to the chosen helicity, and 0 means summing over them."""

        out_data = {'weight': 0.0}
 
        # Call the f2py'ed MEs
        out_data['weight'] = self.ME_accessor.get_me(
            self.format_momenta_for_f2py(PS_point), alpha_s, nhel)

        # Return a data structure with 'weight' in it.
        return out_data

class bbxH_integrand_LO(integrands.VirtualIntegrand):

    def __init__(self, phase_space_generator, options, *args, **opts):
        """ Instantiate the integrand for a specific process."""

        self.options               = options
        self.phase_space_generator = phase_space_generator

        self.n_final    = 1
        self.n_initial  = 2
        # Load the PDFs
        lhapdf_config = misc.which('lhapdf-config')
        lhapdf        = misc.import_python_lhapdf(lhapdf_config)
        if not lhapdf:
            raise Exception('The python lhapdf API could not be loaded.')
        lhapdf.setVerbosity(self.options['verbosity'])
        self.pdfsets = lhapdf.getPDFSet(self.options['PDF'])
        # Pick the central PDF for now
        self.pdf = self.pdfsets.mkPDF(0)
        # Setup the PDF cache (warning, implement a max size to this buffer if memory becomes an issue)
        self.pdfQ2 = {}

        super(bbxH_integrand_LO, self).__init__(*args, **opts)
        
        self.function_list.append(ME_function('b b~ > h'))
 
    @staticmethod
    def Lambda(s,sqrMA,sqrMB):
        """ Kahlen function."""

        return s**2 + sqrMA**2 + sqrMB**2 - \
             2.*s*sqrMA - 2.*sqrMB*sqrMA - 2.*s*sqrMB

    def pass_cuts(self, PS_point):
        """ Implementation of the isolation cuts. """
     
        # Not cuts for b b~ > h @ LO
        return True

        logger.debug( "Processings cuts for process '%s' and PS point:\n%s"%( self.options['process'],
                phase_space_generators.FlatInvertiblePhasespace.nice_momenta_string(
                            PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial) ))

        pt_cut = -100.0
        dr_cut = -0.4
        
        # Apply a uniform P_t cut first
        for i, p in enumerate(PS_point[2:]):
            logger.debug('p_%i.pt()=%.5e'%((self.n_initial+i),p.pt()))
            if p.pt() < pt_cut:
                return False

        for i, p1 in enumerate(PS_point[2:]):
            for j, p2 in enumerate(PS_point[2+i+1:]):
                logger.debug('deltaR(p_%i,p_%i)=%.5e'%(
                     self.n_initial+i, self.n_final+j, p1.deltaR(p2)))
                if p1.deltaR(p2) < dr_cut:
                    return False

        return True
    
    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Integrand function call, with list of continuous and discrete input values for all dimensions."""
       
        logger.debug("="*80)       
        logger.debug('Starting a new evaluation of the integrand %s'%self.__class__.__name__)

        # Just sanity checks, can be disregarded for now
        assert(self.check_input_types(continuous_inputs, discrete_inputs))
        assert(len(discrete_inputs)==
               len([1 for d in self.discrete_dimensions if not d.folded]))
        assert(len(continuous_inputs)==
               len([1 for d in self.continuous_dimensions if not d.folded]))
         
        # A unique float must be returned
        wgt = 1.0
         
        # Random variables sent
        random_variables    = list(continuous_inputs)
        logger.debug('Random variables received: %s'%str(random_variables))

        # Avoid extrema since the phase-space generation algorithm doesn't like it
        epsilon = 1e-10
        random_variables = [min(max(rv,epsilon),1.-epsilon) for rv in random_variables]

        PDF_random_variables = random_variables[:1]
        PS_random_variables  = random_variables[1:]
        
        # We generate the PDF from two variables \tau = x1*x2 and ycm = 1/2 * log(x1/x2)
        # So that:
        #  x_1 = sqrt(tau) * exp(ycm)
        #  x_2 = sqrt(tau) * exp(-ycm)
        # The jacobian of this transformation is 1.
        # Here tau is fixed by the \delta(xb_1*xb_2*s - m_h**2) which sets tau to 
        E_cm = self.options['MH'] 
        tau = (E_cm**2/self.options['collider_energy']**2)
        tau_min = tau
        # And we can now rescale ycm appropriately
        ycm = (-0.5 + PDF_random_variables[0]) * math.log(1./tau_min)
        # We need correct the phase-space volume appropriately
        wgt *= math.log(1./tau_min)

        xb_1 = math.sqrt(tau) * math.exp(ycm)
        xb_2 = math.sqrt(tau) * math.exp(-ycm)

        # Account for the \delta(xb_1*xb_2*s - m_h**2) and corresponding y_cm matching to unit volume
        wgt *= (1./self.options['collider_energy']**2)

        # Make sure the Bjorken x's are physical:
        if xb_1>1. or xb_2>1.:
            logger.debug('Unphysical configuration: x1, x2 = %.5e, %.5e'%(xb_1, xb_2))
            logger.debug(bcolors.OKGREEN + 'Returning a weight of 0. for this integrand evaluation.' + bcolors.ENDC)
            return 0.0

        # Generate a PS point
        PS_point, PS_weight = self.phase_space_generator.generateKinematics(E_cm, PS_random_variables)

        # We must now boost the PS point in the lab frame
        ref_lab = PS_point[0]*xb_1 + PS_point[2]*xb_2
        ref_lab.setMass(ref_lab.calculateMass())
        for p in PS_point:
            p.boost(ref_lab.boostVector())

        logger.debug("Evaluation for process '%s' and following PS point:\n%s"%( self.options['process'],
                phase_space_generators.FlatInvertiblePhasespace.nice_momenta_string(
                            PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial) ))

        # Apply cuts
        if not self.pass_cuts(PS_point):
            logger.debug('Event failing the generation-level cuts.')
            logger.debug(bcolors.OKGREEN + 'Returning a weight of 0. for this integrand evaluation.' + bcolors.ENDC)            
            return 0.0

        # Account for PS weight
        wgt *= PS_weight
        logger.debug("PS_weight: %.5e"%PS_weight)

        # Call the f2py'ed MEs
        helicity = 0 # Means that we sum over all helicity configuration
        alpha_s  = self.pdf.alphasQ(self.options['mu_r'])
        logger.debug("mu_r, mu_f, alpha_s = %.5e, %.5e, %.5e"%(self.options['mu_r'], self.options['mu_f'], alpha_s))

        ME_outputs = [ME_function(PS_point, alpha_s, helicity) for ME_function in self.function_list]
        wgt *= sum(o['weight'] for o in ME_outputs)

        logger.debug("ME weights = [ %s ]"%(', '.join('%.16e'%o['weight'] for o in ME_outputs)))

        # Apply PDF weight
        logger.debug("Bjorken x's x1, x2, sqrt(x1*x2*s): %.5e, %.5e. %.5e"%(
            xb_1, xb_2, math.sqrt(self.options['collider_energy']**2*xb_1*xb_2)))
        PDF1 = self.get_pdfQ2(self.pdf, 5, xb_1, self.options['mu_f']**2)
        PDF2 = self.get_pdfQ2(self.pdf, 5, xb_2, self.options['mu_f']**2)
        wgt *= PDF1*PDF2

        logger.debug("PDF(x1), PDF(x2) = %.5e, %.5e"%(PDF1, PDF2))

        # Include flux factor
        flux = 1.
        if self.phase_space_generator.n_initial == 2:
            flux = 1. / (2.*math.sqrt(self.Lambda(E_cm**2, 
                            self.phase_space_generator.initial_masses[0],
                            self.phase_space_generator.initial_masses[1])))
        elif self.phase_space_generator.n_initial == 1:
            flux = 1. / (2.*E_cm)
        flux /= math.pow(2.*math.pi, 3*self.phase_space_generator.n_final - 4)

        wgt *= flux
        logger.debug("Flux factor: %.5e"%flux)

        # Convert GeV^-2 to picobarns
        wgt *= 0.389379304e9

        if self.apply_observables:
            data_for_observables = {'PS_point': PS_point, 'ME_outputs' : ME_outputs}
            self.observable_list.apply_observables(continuous_inputs, discrete_inputs, data_for_observables)
        
        logger.debug(bcolors.OKGREEN + "Final weight returned: %.5e"%wgt + bcolors.ENDC)
        logger.debug("="*80)

        return wgt

    def get_pdfQ2(self, pdf, pdg, x, scale2):
       
        if pdg not in [21,22] and abs(pdg) not in range(1,7):
            return 1.
                
        if (pdf, pdg, x, scale2) in self.pdfQ2:
            return self.pdfQ2[(pdf, pdg, x, scale2)]
        
        # Call to lhapdf API
        f = pdf.xfxQ2(pdg, x, scale2)/x

        self.pdfQ2[(pdf, pdg,x,scale2)] = f
        return f 

class MadEvent7_driver(object):

    def __init__(self,
                 options = { 'collider_energy'         : 14000.0,
                             'process'                 : 'b b~ > h @LO',
                             'mu_r'                    : 91.188,
                             'mu_f'                    : 91.188,                            
                             'PDF'                     : 'NNPDF23_nlo_as_0119',
                             'MH'                      : 125.0,
                             'verbosity'               : 0
                           },
                 integrator_class       = integrators.SimpleMonteCarloIntegrator, 
                 integrator_options     = {}):
                
        self.options = options
        self.process = self.options['process']
        
        # List of allowed processe for now.
        assert(self.process in ['b b~ > h @LO','b b~ > h @NLO'])
        
        self.n_final = None
        if self.process in ['b b~ > h @LO','b b~ > h @NLO']:
            self.n_final = 1
        self.collider_energy = self.options['collider_energy']

        if self.process in ['b b~ > h @LO','b b~ > h @NLO']:                      
            self.born_phase_space_generator = phase_space_generators.FlatInvertiblePhasespace(
                                                      [0.]*2, [self.options['MH']]+[0.]*(self.n_final-1))

        self.build_integrand()

        self.set_integrator(integrator_class, integrator_options)

    def set_integrator(self, integrator_class=integrators.SimpleMonteCarloIntegrator, integrator_options={}):
        """ Assign a particular integrator."""
        self.integrator = integrator_class(self.integrands, **integrator_options)

    def build_integrand(self):
        """ Build the integrand for the process b b~ > h (j) @(N)LO """
        
        if self.process in ['b b~ > h @LO','b b~ > h @NLO']:
            bjorken_dimensions = [ integrands.ContinuousDimension('ycms',lower_bound=0.0, upper_bound=1.0) ]

        PS_dimensions = integrands.DimensionList(
            bjorken_dimensions+
            [ integrands.ContinuousDimension('x_%d'%i,lower_bound=0.0, upper_bound=1.0) 
                            for i in range(1,self.born_phase_space_generator.nDimPhaseSpace()+1) ] )

        if self.process in ['b b~ > h @LO'] :
            bbxH_LO_integrand = bbxH_integrand_LO(self.born_phase_space_generator, self.options, PS_dimensions)
            self.integrands = [bbxH_LO_integrand]

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


# Available verbosity levels:
#   0  -> Nothing
#   1  -> Periodic reports from the integrator
#   10 -> Debug printouts from the integration framework (mainly from __call__ of the integrands)
verbosity  = 1

if verbosity == 0:
    logger.info("Set verbosity to a value larger than 0 to watch more closely the integration progress.")

if verbosity >= 10:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

# Try with a naive MonteCarloFirst
parameters = {
    'process'               : 'b b~ > h @LO',
    'collider_energy'       : 13000.0,
    'mu_r'                  : 91.188,
    'mu_f'                  : 91.188,
    # Corresponds to lhaid = 230000
    'PDF'                   : 'NNPDF23_nlo_as_0119',
    'MH'                    : 125.0,
    'verbosity'             : verbosity,
}

# This is for a Higgs mass of 125.0 and a bottom yukawa of 4.3 GeV
ME_references_values = { 'b b~ > h @LO'  : (8.294e-01, 6.1e-04),
                         'b b~ > h @NLO' : (7.308e-01, 1.3e-03)  }

# Initialize a MadEvent driver
ME7 = MadEvent7_driver( parameters )

# Now launch the integration (or skip it by commenting below if not interested in the stupid integrator.)
## ME7.run()

# Try several integrators
integrators = {

   'Naive' : (integrators.SimpleMonteCarloIntegrator, {  'n_iterations'            : 10,
                                                         'n_points_per_iterations' : 100,
                                                         'accuracy_target'         : None,
                                                         'verbosity'               : verbosity } ),

   'VEGAS' : (pyCubaIntegrator.pyCubaIntegrator, { 'algorithm' : 'Vegas', 
                                                   'verbosity' : verbosity,
                                                   'seed'      : 3,
                                                   'target_accuracy' : 1.0e-3,
                                                   'n_start'   : 1000,
                                                   'n_increase': 500,
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
for name in ['VEGAS']:
    # Set the integrator
    ME7.set_integrator(*integrators[name])
    # Run one last time
    ME7.run()

# Print out the reference value if specified
if parameters['process'] in ME_references_values:
    xsec, err = ME_references_values[parameters['process']]
    logger.info("MadEvent6 reference result for process '%s' = \033[92m%.5e +/- %.2e [pb]\033[0m "%(parameters['process'], xsec, err)) 
