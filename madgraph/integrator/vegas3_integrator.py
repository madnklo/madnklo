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
import time
import os
import logging
import random
from multiprocessing import Process

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir))

try:
    import vegas 
except ImportError:
    print "WARNING :: Vegas 3 could not be found, install it if needed."

try:
    import numpy as np
except ImportError:
    print "ERROR :: Vegas 3 requires numpy, install it first."
    sys.exit(1)

try:
    import gvar as gv
except ImportError:
    print "WARNING :: Vegas 3 requires gvar (gaussian variable handling package), install it if needed."

import madgraph.various.misc as misc
import madgraph.iolibs.files as files
import madgraph.integrator.integrators as integrators
from madgraph.integrator.integrators import IntegratorError
import madgraph.integrator.integrands as integrands
import madgraph.integrator.functions as functions
import madgraph.various.cluster as cluster
import madgraph.iolibs.save_load_object as save_load_object

from madgraph import InvalidCmd, MadGraph5Error, MG5DIR, MPI_ACTIVE, MPI_SIZE, MPI_RANK
from multiprocessing import Value, Array

if MPI_ACTIVE:
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD

numpy = np

logger = logging.getLogger('madgraph.vegas3')

class ParallelWrappedIntegrand(vegas.BatchIntegrand):
    """ Class to wrap the integrand to the Vegas3 standards with parallelization."""
    
    def __init__(self, integrator, integrands, cluster_instance, start_time=None):
        " Save integrand; registers a cluster instance."
        self.integrator = integrator
        self.integrands = integrands
        self.cluster = cluster_instance
        self.start_time = start_time
        if not isinstance(cluster_instance,(cluster.MultiCore,cluster.MPICluster)):
            raise IntegratorError("Vegas3 integrator does not have cluster support "+
                                                                     "yet, only multicore")
        
    def __del__(self):
        " Standard cleanup. "
        pass

    @staticmethod
    def batch_integrand_process(integrands, inputs_list, result):
        """ Function to call the wrapped integrand for a batch of inputs"""
        
        for input_position, inputs in enumerate(inputs_list):
            sum = 0.
            dims = integrands[0].continuous_dimensions
            fct_inputs = np.array(
                [dims[i].lower_bound + x*(dims[i].upper_bound - dims[i].lower_bound) 
                                                          for i, x in enumerate(inputs)])
            for integrand in integrands:
                res = integrand(fct_inputs, np.array([], dtype=int))
                res *= integrand.get_dimensions().volume()
                sum += res

            result[input_position] = sum

        return 0

    def batch_integrand(self, inputs_list, result):
        """ Integrand wrapper spawning an independent process for the computation of
        the batch of points listed in inputs_list."""
        p = Process(target=ParallelWrappedIntegrand.batch_integrand_process, args=(
                                                    self.integrands, inputs_list, result))
        p.start()
        p.join()
        return 0
        
    def wait_monitoring(self, Idle, Running, Done):
        if Idle+Running+Done == 0:
            return
        logger.debug('MadEvent7 integration: %d Idle, %d Running, %d Done'%(Idle, Running, Done))

    def __call__(self, x):
        " Divide x into self.nproc chunks, feeding one to each process. "
        batch_size = x.shape[0]
        nx = batch_size // self.cluster.nb_core + 1
        inputs  = [x[i*nx : (i+1)*nx] for i in range(self.cluster.nb_core)]
        # The results enttried will be accessed from independent processes and must
        # therefore be defined as shared memory
        results = [Array('d', [0.]*len(input)) for input in inputs]
        # launch evaluation of self.integrands for each chunk, in parallel
        if not self.start_time is None:
            logger.debug('Dispatching batch of %d points on %d cores: [%s]'%(
            batch_size,self.cluster.nb_core, misc.format_time(time.time()-self.start_time)))
        else:
            logger.debug('Dispatching batch of %d points on %d cores:'%(
                                                          batch_size,self.cluster.nb_core))            
        for position, input in enumerate(inputs):
            self.cluster.submit(self.batch_integrand, [input, results[position]])            

        # Wait for all jobs to finish.
        self.cluster.wait('', self.wait_monitoring,update_first=self.wait_monitoring)
        
        self.integrator.n_function_evals += batch_size
        return np.concatenate([ np.array(result[:]) for result in results ])

class Vegas3Integrator(integrators.VirtualIntegrator):

    def __init__(self, integrands, **opts):
        """ Initialize the Vegas3 integrator."""
       
        # General steering parameters
        default_opts= {
           'verbosity'          :   0,
           'target_accuracy'    :   1.0e-3,
           'seed'               :   None,
           'save_grids'         :   None,
           'load_grids'         :   None,
        }

        # Parameter for controlling VEGAS3
        # ------------------------------------------------
        
        # the number of iterations for the training session. 
        default_opts['n_iterations_survey'] = 10
        # the number of points per iteration in the training session. 
        default_opts['n_points_survey'] = 2000

        # the number of iterations for the production session. 
        default_opts['n_iterations_refine'] = 10
        # the number of points per iteration in the production session. 
        default_opts['n_points_refine'] = 10000

        # Steering parallelization
        if not MPI_ACTIVE:
            default_opts['cluster'] = cluster.onecore
        else:
            default_opts['cluster'] = cluster.MPICluster(MPI_RANK,MPI_SIZE)

        # Default batch size per node
        default_opts['batch_size'] = 1000

        # Set instance attributes options
        for opt in default_opts:
            if opt in opts:
                setattr(self,opt,opts.pop(opt))
            else:
                setattr(self,opt,default_opts[opt])
                
        

        super(Vegas3Integrator,self).__init__(integrands, **opts)

        if self.seed is not None:
            random.seed(self.seed)
            np.random.seed(self.seed)

        # Adjust batch size proportionally to the number of nodes
        self.batch_size *= self.cluster.nb_core

        # Save here the main instance for vegas3
        self.vegas3_integrator = None 
        
        if self.verbosity > 0:
            logger.level = logging.DEBUG
        else:
            logger.level = logging.INFO

    def observables_normalization(self, n_integrand_calls):
        """ Given the number of integrand calls, return the appropriate overall normalization
        to apply to the observables."""
        return 1.0/float(self.n_iterations_refine)

    def wrapped_integrand(self, x_inputs):
        """ Function to wrap the integrand to the Vegas3 standards."""

        dims = self.integrands[0].continuous_dimensions
        if len(list(x_inputs))>len(dims):
            jacobian = x_inputs[-1]*dims.volume()
            x_inputs = list(x_inputs)[:-1]
        else:
            x_inputs = list(x_inputs)
            jacobian = 1.0*dims.volume()

        fct_inputs = np.array([dims[i].lower_bound + x*(dims[i].upper_bound - dims[i].lower_bound) 
                                                                         for i, x in enumerate(x_inputs)])
        
        all_results = []
        for integrand in self.integrands:
            res = integrand(fct_inputs, np.array([], dtype=int), integrator_jacobian=jacobian)
            res *= integrand.get_dimensions().volume()
            all_results.append(res)
            
        self.n_function_evals +=1
        if self.n_function_evals%(max(self.curr_n_evals_per_iterations/2,1))==0 and MPI_RANK == 0:
            logger.debug('Evaluation #%d / %d*%d  (%.3g%%) [%s]'%(
                self.n_function_evals*MPI_SIZE, self.curr_n_iterations, self.curr_n_evals_per_iterations, 
                (100.0*self.n_function_evals*MPI_SIZE / (self.curr_n_iterations*self.curr_n_evals_per_iterations)),
                                            misc.format_time(time.time()-self.start_time)))
        
        return dict( ('I%d'%i, res) for i, res in enumerate(all_results) )

    def get_name(self):
        """ Returns the integrator name."""
        return 'Vegas3' 

    def integrate(self):
        """ Return the final integral and error estimators. """

        start_time = time.time()

        for integrand in self.integrands:
            if len(integrand.discrete_dimensions)>0:
                raise IntegratorError("Vegas3 integrator does not support discrete dimensions for now.")

        if len(set([len(integrand.continuous_dimensions) for integrand in self.integrands]))!=1:
            raise IntegratorError("Vegas3 only supports multiple integrands with all the same number of dimensions.")

        # During the survey, it is sub-optimal to fill in histograms, and we therefore disable here
        # until the refine stage is reached.
        # We must of course keep a back-up copy of the 'apply_observable' attribute of the
        # integrand instances so as to be able to reinstate it before the refine.
        if MPI_RANK == 0:
            apply_observables_for_integrands_back_up = []
            for i_integrand, integrand in enumerate(self.integrands):
                apply_observables_for_integrands_back_up.append(integrand.apply_observables)
                integrand.apply_observables = False

        n_dimensions = len(self.integrands[0].continuous_dimensions)
        # sync_ran is to decide if VEGAS3 random number generator should produce the same
        # number on different processors.
        if any(apply_observables_for_integrands_back_up) and self.cluster.nb_core==1:
#            if self.n_iterations_refine > 1:
#                logger.warning("Vegas3 can only run a single refine iteration when a fixed-order analysis is active.\n"+
#                               "The parameter 'n_iterations_refine' will consequently be forced to take the value 1.")
#                self.n_iterations_refine = 1
            self.vegas3_integrator = VegasWithJacobianInFunctionInput(n_dimensions * [[0., 1.]],
                analyzer        = vegas.reporter() if self.verbosity>1 else None, 
                nhcube_batch    = self.batch_size,
                # We must disable the hypercube optimisation for now when plotting observables
                max_nhcube      = 1,
                beta            = 0.,
                alpha           = 0.5,
                sync_ran        = True)
        else:
            self.vegas3_integrator = vegas.Integrator(n_dimensions * [[0., 1.]],
                analyzer        = vegas.reporter() if self.verbosity>1 else None, 
                nhcube_batch    = self.batch_size,
                # We must disable the hypercube optimisation for now when plotting observables
                max_nhcube      = 1e9,
                beta            = 0.75,
                alpha           = 0.5,
                sync_ran        = True)
        
        self.start_time = time.time()
        # Access/generate the wrapped integrand
        if isinstance(self.cluster,cluster.MPICluster) or self.cluster.nb_core==1:
            wrapped_integrand = self.wrapped_integrand
        else:
            wrapped_integrand = ParallelWrappedIntegrand(
                                      self, self.integrands, self.cluster, self.start_time)          

        self.tot_func_evals = 0
        result = None
        
        # Load pre-exisitng grids if specified. Note that no information on Vegas' hypercubes is supplied
        # here, only the standard integration grids along each dimension.
        if self.load_grids is not None:
            try:
                vegas_map = save_load_object.load_from_file(self.load_grids)
            except Exception as e:
                raise IntegratorError("Vegas3 could not load its integration grid from '%s'. Error: %s"%(
                                                                    self.load_grids,str(e)))
            logger.info('Vegas3 reading integration grids from file: %s'%self.load_grids)
            self.vegas3_integrator.set(map=vegas_map)
        
        if self.n_iterations_survey > 0 and self.n_points_survey > 0:
            # Train grid
            if MPI_RANK == 0:
                logger.debug("=================================")
                logger.debug("Vegas3 starting the survey stage.")
                logger.debug("=================================")
            self.n_function_evals = 0
            self.curr_n_iterations = self.n_iterations_survey
            self.curr_n_evals_per_iterations = self.n_points_survey
            result = self.vegas3_integrator(wrapped_integrand, 
                        nitn=self.n_iterations_survey, neval=self.n_points_survey, adapt=True)
    
            if MPI_RANK == 0:
                logger.debug('\n'+result.summary())
                
            self.tot_func_evals += self.n_function_evals
        else:
            logger.info("The survey step of Vegas3 integration is being skipped as per user's request.")
           
        # Save generated grids if asked for. Notice that no information on Vegas' hypercubes is saved here,
        # but just the standard integration grids.
        if self.save_grids is not None:
            try:
                vegas_map = save_load_object.save_to_file(self.save_grids, self.vegas3_integrator.map) 
            except Exception as e:
                raise IntegratorError("Vegas3 could not save its integration grid to '%s'. Error: %s"%(
                                                                    self.load_grids,str(e)))
            logger.info('Vegas3 saved its integration grids to file: %s'%self.save_grids)
           
        if self.n_iterations_refine > 0 and self.n_points_refine > 0:
            # Final integration
            if MPI_RANK == 0:
                logger.debug("=============================================")
                logger.debug("Vegas3 starting the refined integration stage")
                logger.debug("=============================================")
    
                # Restore the filling of the histograms for the refine, only if not running in parallel
                if self.cluster.nb_core==1:
                    for i_integrand, integrand in enumerate(self.integrands):
                        integrand.apply_observables=apply_observables_for_integrands_back_up[i_integrand]
                else:
                    if any(apply_observables_for_integrands_back_up):
                        logger.warning('Filling of the histograms of the fixed-order analysis will'+
                            ' now be disabled as this functionality is not yet available for parallel integration.')
    
            self.n_function_evals = 0
            self.curr_n_iterations = self.n_iterations_refine
            self.curr_n_evals_per_iterations = self.n_points_refine
            result = self.vegas3_integrator(wrapped_integrand, 
                        nitn=self.n_iterations_refine, neval=self.n_points_refine, adapt=False) 
            if MPI_RANK == 0:
                logger.debug('\n'+result.summary())
        else:
            logger.info("The refine step of Vegas3 integration is being skipped as per user's request.")

        if MPI_RANK == 0:
            # Make sure to restore original settings for the 'apply_observable' attribute
            for i_integrand, integrand in enumerate(self.integrands):
                integrand.apply_observables=apply_observables_for_integrands_back_up[i_integrand]
       
        self.tot_func_evals += self.n_function_evals
        if MPI_ACTIVE:
            # Aggregate the number of function evaluation from each MPI rank
            if MPI_RANK == 0:
                for rank in range(1,MPI_SIZE):
                    n_points_in_slave = mpi_comm.recv(source=rank,tag=11)
                    self.tot_func_evals += n_points_in_slave
            else:
                mpi_comm.send(self.n_function_evals, dest=0, tag=11)

        # Result is a list of instances of vegas.RAvg variablesr, with the following attributes: 
        # RAvg.mean : mean 
        # RAvg.sdev : standard dev.
        # RAvg.chi2 : Chi^2 of the weighted average
        # RAvg.dof  : number of degreeas of freedom
        # RAvg.Q    : p-value of the weighted average 
        # RAvg.itn_results : list of the integral estimates for each iteration
        # RAvg.summary() : summary of the integration
        if isinstance(self.cluster,cluster.MPICluster) or self.cluster.nb_core==1:
            summed_result = sum(result.values())
        else:
            summed_result = result
           
        integration_time = time.time() - start_time

        if MPI_RANK == 0:
            logger.debug("===============================================================")
            logger.debug('Total integration time: %s'%misc.format_time(integration_time))
            logger.debug('Vegas3 used a total of %d function evaluations, on %d core(s).'%
                                                  (self.tot_func_evals,self.cluster.nb_core))
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


# Wrapper around the original Vegas.Integrator which adds the jacobian of the point as
# the last input of the 1-dimensional array passed to the integrand function in the 
# __call__ function of the integrator.
# WARNING: does not work with vegas.Integrator.beta > 0 because Cython somehow makes the
# float attribute 'sum_sigf' of the integrand read-only, making the line:
#     self.sum_sigf = sum_sigf
# crash.
class VegasWithJacobianInFunctionInput(vegas.Integrator):

    def __init__(self, *args, **opts):
        super(VegasWithJacobianInFunctionInput, self).__init__(*args, **opts)

    def __call__(self, fcn, **kargs):
        """ Integrate integrand ``fcn``.

        A typical integrand has the form, for example::

            def f(x):
                return x[0] ** 2 + x[1] ** 4

        The argument ``x[d]`` is an integration point, where
        index ``d=0...`` represents direction within the
        integration volume.

        Integrands can be array-valued, representing multiple
        integrands: e.g., ::

            def f(x):
                return [x[0] ** 2, x[0] / x[1]]

        The return arrays can have any shape. Dictionary-valued
        integrands are also supported: e.g., ::

            def f(x):
                return {'a':x[0] ** 2, 'b':[x[0] / x[1], x[1] / x[0]]}


        Integrand functions that return arrays or dictionaries
        are useful for multiple integrands that are closely related,
        and can lead to substantial reductions in the errors for
        ratios or differences of the results.

        It is usually much faster to use |vegas| in batch
        mode, where integration points are presented to the
        integrand in batches. A simple batch integrand might
        be, for example::

            @vegas.batchintegrand
            def f(x):
                return x[:, 0] ** 2 + x[:, 1] ** 4

        where decorator ``@vegas.batchintegrand`` tells
        |vegas| that the integrand processes integration
        points in batches. The array ``x[i, d]``
        represents a collection of different integration
        points labeled by ``i=0...``. (The number of points is controlled
        |Integrator| parameter ``nhcube_batch``.) The batch index
        is always first.

        Batch integrands can also be constructed from classes
        derived from :class:`vegas.BatchIntegrand`.

        Batch mode is particularly useful (and fast) when the class
        derived from :class:`vegas.BatchIntegrand` is coded
        in Cython. Then loops over the integration points
        can be coded explicitly, avoiding the need to use
        :mod:`numpy`'s whole-array operators if they are not
        well suited to the integrand.

        Any |vegas| parameter can also be reset: e.g.,
        ``self(fcn, nitn=20, neval=1e6)``.
        """

        wgt      = numpy.empty(1, numpy.float_)
        mean     = numpy.empty(1, numpy.float_)
        var      = numpy.empty((1, 1), numpy.float_)

        neval_batch = self.nhcube_batch * self.min_neval_hcube
        self.fdv2 = numpy.empty(neval_batch, numpy.float_)

        firsteval = True

        self.neval_hcube = (
            numpy.zeros(self.nhcube_batch, numpy.intp) + self.min_neval_hcube
            )

        if kargs:
            self.set(kargs)

        # synchronize random numbers across all processes (mpi)
        if self.sync_ran:
            self.synchronize_random()

        # Put integrand into standard form
        fcn = vegas._vegas.VegasIntegrand(fcn)

        sigf = self.sigf
        for itn in range(self.nitn):
            # if self.minimize_mem:
            #     self.set()
            if self.analyzer is not None:
                self.analyzer.begin(itn, self)

            # initalize arrays that accumulate results for a single iteration
            mean[:] = 0.0
            var[:, :] = 0.0
            sum_sigf = 0.0
                        
            # iterate batch-slices of integration points
            for x, y, wgt, hcube in self.random_batch(
                yield_hcube=True, yield_y=True, fcn=fcn
                ):
                fdv2 = self.fdv2        # must be inside loop
               
                ##############################################################################################
                # START OF MODIFICATION FROM THE ORIGINAL IMPLEMENTATION OF VEGAS3
                ##############################################################################################
                # Append the jacobian of the integrator as the last entry of the inputs passed to the function
                x = numpy.asarray( [ numpy.append(one_x,[wgt[i]]) for i, one_x in enumerate(x)], numpy.float_)
                ##############################################################################################
                # EMD OF MODIFICATION FROM THE ORIGINAL IMPLEMENTATION OF VEGAS3
                ##############################################################################################

                # evaluate integrand at all points in x                
                fx = fcn.eval(x)

                if firsteval:
                    # allocate work arrays on first pass through;
                    # (needed a sample fcn evaluation in order to do this)
                    firsteval = False
                    wf = numpy.empty(fcn.size, numpy.float_)
                    sum_wf = numpy.empty(fcn.size, numpy.float_)
                    sum_wf2 = numpy.empty((fcn.size, fcn.size), numpy.float_)
                    mean = numpy.empty(fcn.size, numpy.float_)
                    var = numpy.empty((fcn.size, fcn.size), numpy.float_)
                    mean[:] = 0.0
                    var[:, :] = 0.0
                    result = vegas._vegas.VegasResult(fcn, weighted=self.adapt)

                # compute integral and variance for each h-cube
                # j is index of hcube within batch, i is absolute index
                j = 0
                for i in range(hcube[0], hcube[-1] + 1):
                    # iterate over h-cubes
                    sum_wf[:] = 0.0
                    sum_wf2[:, :] = 0.0
                    neval = 0
                    while j < len(hcube) and hcube[j] == i:
                        for s in range(fcn.size):
                            wf[s] = wgt[j] * fx[j, s]
                            sum_wf[s] += wf[s]
                            for t in range(s + 1):
                                sum_wf2[s, t] += wf[s] * wf[t]
                        fdv2[j] = (wf[0] * self.neval_hcube[i - hcube[0]]) ** 2
                        j += 1
                        neval += 1
                    for s in range(fcn.size):
                        mean[s] += sum_wf[s]
                        for t in range(s + 1):
                            var[s, t] += (
                                sum_wf2[s, t] * neval - sum_wf[s] * sum_wf[t]
                                ) / (neval - 1.)
                        if var[s, s] <= 0:
                            var[s, s] = mean[s] ** 2 * 1e-15 + TINY
                    sigf2 = abs(sum_wf2[0, 0] * neval - sum_wf[0] * sum_wf[0])
                    if self.beta > 0 and self.adapt:
                        if not self.minimize_mem:
                            sigf[i] = sigf2 ** (self.beta / 2.)
                            sum_sigf += sigf[i]
                        else:
                            sum_sigf += sigf2 ** (self.beta / 2.)
                    if self.adapt_to_errors and self.adapt:
                        # replace fdv2 with variance
                        fdv2[j - 1] = sigf2
                        self.map.add_training_data(
                            y[j - 1:, :], fdv2[j - 1:], 1
                            )
                if (not self.adapt_to_errors) and self.adapt and self.alpha > 0:
                    self.map.add_training_data(y, fdv2, y.shape[0])

            for s in range(var.shape[0]):
                for t in range(s):
                    var[t, s] = var[s, t]

            # accumulate result from this iteration
            result.update(mean, var)

            if self.beta > 0 and self.adapt:
                self.sum_sigf = sum_sigf
            if self.alpha > 0 and self.adapt:
                self.map.adapt(alpha=self.alpha)
            if self.analyzer is not None:
                result.update_analyzer(self.analyzer)

            if result.converged(self.rtol, self.atol):
                break
        return result.result

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

