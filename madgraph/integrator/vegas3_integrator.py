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

import sys
import os
import logging

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir))

try:
    import vegas 
except ImportError:
    print "WARNING :: Vegas 3 could not be found, install it if needed."

try:
    import numpy as np
except ImportError:
    print "ERROR :: Vegas 3 requires nunmpy, install it first."
    sys.exit(1);

try:
    import gvar as gv
except ImportError:
    print "WARNING :: Vegas 3 requires gvar (gaussian variable handling package), install it if needed."


try:
    import madgraph
except ImportError:
    MADEVENT = True
    import internal.banner as bannermod
    import internal.misc as misc
    import internal.files as files
    import internal.integrator.integrators as integrators
    from internal.integrator.integrators import IntegratorError
    import internal.integrator.integrands as integrands
    import internal.integrator.functions as functions
    from internal import InvalidCmd, MadGraph5Error

else:
    MADEVENT= False
    import madgraph.various.misc as misc
    import madgraph.iolibs.files as files
    import madgraph.integrator.integrators as integrators
    from madgraph.integrator.integrators import IntegratorError
    import madgraph.integrator.integrands as integrands
    import madgraph.integrator.functions as functions
    from madgraph import InvalidCmd, MadGraph5Error, MG5DIR

logger = logging.getLogger('madgraph.vegas3')

class Vegas3Integrator(integrators.VirtualIntegrator):

    def __init__(self, integrands, **opts):
        """ Initialize the Vegas3 integrator."""
       
        # General steering parameters
        default_opts= {
           'verbosity':0,
           'target_accuracy':1.0e-3}

        # Parameter for controlling VEGAS3
        # ------------------------------------------------
        
        # the number of iterations for the training session. 
        default_opts['survey_n_iterations'] = 10
        # the number of points per iteration in the training session. 
        default_opts['survey_n_points'] = 2000

        # the number of iterations for the production session. 
        default_opts['refine_n_iterations'] = 10
        # the number of points per iteration in the production session. 
        default_opts['refine_n_points'] = 10000

        # Set instance attributes options
        for opt in default_opts:
            if opt in opts:
                setattr(self,opt,opts.pop(opt))
            else:
                setattr(self,opt,default_opts[opt])

        super(Vegas3Integrator,self).__init__(integrands, **opts)

        # Save here the main instance for vegas3
        self.vegas3_integrator = None 
        
        if self.verbosity > 0:
            logger.level = logging.DEBUG
        else:
            logger.level = logging.INFO

    def wrapped_integrand(self, x_inputs):
        """ Function to wrap the integrand to the Vegas3 standards."""

        dims = self.integrands[0].continuous_dimensions
        fct_inputs = np.array([dims[i].lower_bound + x*(dims[i].upper_bound - dims[i].lower_bound) 
                                                                         for i, x in enumerate(x_inputs)])
        all_results = [] 
        for integrand in self.integrands:
            res = integrand(fct_inputs, np.array([], dtype=int))
            res *= integrand.get_dimensions().volume()
            all_results.append(res)
            
        self.n_function_evals +=1
        if self.n_function_evals%(self.curr_n_evals_per_iterations/4)==0 and \
                                                      self.vegas3_integrator.mpi_rank == 0:
            logger.debug('Evaluation #%d / %d*%d  (%.3g%%)'%(self.n_function_evals, self.curr_n_iterations, self.curr_n_evals_per_iterations, 
                                            (100.0*self.n_function_evals / (self.curr_n_iterations*self.curr_n_evals_per_iterations))))
        
        return dict( ('I%d'%i, res) for i, res in enumerate(all_results) )

    def get_name(self):
        """ Returns the integrator name."""
        return 'Vegas3' 

    def integrate(self):
        """ Return the final integral and error estimators. """

        for integrand in self.integrands:
            if len(integrand.discrete_dimensions)>0:
                raise IntegratorError("Vegas3 integrator does not support discrete dimensions for now.")

        if len(set([len(integrand.continuous_dimensions) for integrand in self.integrands]))!=1:
            raise IntegratorError("Vegas3 only supports multiple integrands with all the same number of dimensions.")
        
        n_dimensions = len(self.integrands[0].continuous_dimensions)
        self.vegas3_integrator = vegas.Integrator(n_dimensions * [[0., 1.]])

        self.tot_func_evals = 0
        # Train grid
        if self.vegas3_integrator.mpi_rank == 0:
            logger.debug("=================================")
            logger.debug("Vegas3 starting the survey stage.")
            logger.debug("=================================")
        self.n_function_evals = 0
        self.curr_n_iterations = self.survey_n_iterations
        self.curr_n_evals_per_iterations = self.survey_n_points
        training = self.vegas3_integrator(self.wrapped_integrand, 
                                          nitn=self.survey_n_iterations, neval=self.survey_n_points)

        if self.vegas3_integrator.mpi_rank == 0:
            logger.debug('\n'+training.summary())

        self.tot_func_evals += self.n_function_evals         
        # Final integration
        if self.vegas3_integrator.mpi_rank == 0:
            logger.debug("=============================================")
            logger.debug("Vegas3 starting the refined integration stage")
            logger.debug("=============================================")

        self.n_function_evals = 0
        self.curr_n_iterations = self.refine_n_iterations
        self.curr_n_evals_per_iterations = self.refine_n_points
        result = self.vegas3_integrator(self.wrapped_integrand, 
                                        nitn=self.refine_n_iterations, neval=self.refine_n_points) 
        if self.vegas3_integrator.mpi_rank == 0:
            logger.debug('\n'+result.summary())
        
        self.tot_func_evals += self.n_function_evals
        # Result is a list of instances of vegas.RAvg variablesr, with the following attributes: 
        # RAvg.mean : mean 
        # RAvg.sdev : standard dev.
        # RAvg.chi2 : Chi^2 of the weighted average
        # RAvg.dof  : number of degreeas of freedom
        # RAvg.Q    : p-value of the weighted average 
        # RAvg.itn_results : list of the integral estimates for each iteration
        # RAvg.summary() : summary of the integration    
        summed_result = sum(result.values())
        
        if self.vegas3_integrator.mpi_rank == 0:
            logger.debug("===============================================================")
            logger.debug('Vegas3 used a total of %d function evaluations.'%self.tot_func_evals)
            logger.debug('Vegas3 returned final results : %s'%summed_result)
            logger.debug("===============================================================")

        return summed_result.mean, summed_result.sdev

    def show_grid(self, n_grid=40, shrink=False, axes=None):
        """ Shows the integration grid of Vegas3 in a neat matplotlib window.
        n_grid is the number of grid nodes in each direction to plot, 
        shrink is whether to show the entire range or just the n_grid nodes and
        axes can be a list of tuples specifying the list of directions to show."""

        if not self.vegas3_integrator:
            raise InvalidCmd("Vegas3 can only show grids after one integration has been performed.") 
        if not axes and not shrink:
            self.vegas3_integrator.map.show_grid(n_grid)
        else:
            raise IntegratorError("Options for show_grid not correctly handled yet.")
            self.vegas3_integrator.map.show_grid(n_grid, shrink, axes=axes)

if __name__ == '__main__':

    # Example usage of the new integrator framework
    
    # First build the integrand
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

    # Then the integrators
    # For a typical verbosity, chose 2
    verbosity = 1
    my_integrators = {}

    my_integrators['Naive MC'] = integrators.SimpleMonteCarloIntegrator(my_integrand, 
        n_iterations=10, n_points_per_iterations=1000, accuracy_target=None, verbosity = verbosity)

    my_integrators['VEGAS3'] = Vegas3Integrator(my_integrand, verbosity = verbosity)

    chosen_integrators = ['Naive MC', 'VEGAS3']

    for name in chosen_integrators:
        res = my_integrators[name].integrate()
        # Finally integrate
        print '\n'+'='*70
        print '%-30s = %.4e +/- %.2e'%('%s :: Final result'%name, res[0], res[1])
        print '='*70

