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

logger = logging.getLogger('pyCubaIntegrator')

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir))

try:
    import pycuba
except ImportError:
    print "WARNING :: The pyCubaIntegrator requires the pyCuba python module. Install it first if needed."

try:
    import numpy as np
except ImportError:
    print "ERROR :: The pyCubaIntegrator requires the numpy python module. Install it first."
    sys.exit(1);


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

else:
    MADEVENT= False
    import madgraph.various.misc as misc
    import madgraph.iolibs.files as files
    import madgraph.integrator.integrators as integrators
    from madgraph.integrator.integrators import IntegratorError
    import madgraph.integrator.integrands as integrands
    import madgraph.integrator.functions as functions

class tmpGlobal(object):
    """Designed so as to temporary add a global variable
    """

    def __init__(self, name, value):
        self.global_name  = name
        self.global_value = value

    def __enter__(self):
        globals()[self.global_name] = self.global_value

    def __exit__(self, ctype, value, traceback ):
        try:
            del globals()[self.global_name]
        except KeyError:
            pass

class pyCubaIntegrator(integrators.VirtualIntegrator):

    def __init__(self, integrands, **opts):
        """ Initialize the pyCuba integrator."""
       
        # General steering parameters
        default_opts= {
           'algorithm':'Vegas',
           'verbosity':2,
           'target_accuracy':1.0e-3}

        # Number of cores to be used in Cuba.
        # Recover nb_core in the option'cluster'
        if 'cluster' in opts and hasattr(opts['cluster'], 'nb_core'):
            default_opts['nb_core'] = opts['cluster'].nb_core
        else:
            # By default turn off parallelization as it is potentially problematic
            # (see warning message later).
            default_opts['nb_core'] = 0

        # Maximum number of evaluation before returning an answer
        default_opts['max_eval'] = 50000
        # Minimum number of evaluation before returning an answer
        default_opts['min_eval'] = 0

        # Parameter relevant for Suave integration method
        # ------------------------------------------------
        
        # the number of integrand evaluations per iteration to start with. 
        default_opts['n_start'] = 1000
        # the increase in the number of integrand evaluations per iteration. 
        default_opts['n_increase'] = 500
        # the batch size for sampling.
        # Vegas samples points not all at once, but in batches of size nbatch, to avoid exces-
        # sive memory consumption. 1000 is a reasonable value, though it should not affect
        # performance too much. 
        default_opts['n_batch'] = 1000
        # Random seed
        default_opts['seed'] = None

        # Parameter relevant for Suave integration method
        # ------------------------------------------------
        
        # the number of new integrand evaluations in each subdivision.
        default_opts['n_new'] = 1000
        # the minimum number of samples a former pass must contribute
        # to a subregion to be considered in that region's compound integral value. Increasing
        # nmin may reduce jumps in the chi^2 value.
        default_opts['n_min'] = 2

        # the parameter p in Eq. (1), i.e. the type of norm
        # used to compute the fluctuation of a sample. This determines how prominently 'out-liers'
        # . i.e. individual samples with a large fluctuation, figure in the total fluctuation,
        # which in turn determines how a region is split up. As suggested by its name, flatness
        # should be chosen large for 'flat' integrands and small for 'volatile' integrands with
        # high peaks. Note that since flatness appears in the exponent, one should not use
        # too large values (say, no more than a few hundred) lest terms be truncated internally
        # to prevent overflow.
        default_opts['flatness'] = 50.

        # Parameter relevant for Divonne integration method
        # ------------------------------------------------

        default_opts['key1'] = 47
        # *key1*: determines sampling in the partitioning phase:
        # key1 = 7, 9, 11, 13 selects the cubature rule of degree key1. Note that the degree-11
        # rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
        # For other values of key1, a quasi-random sample of n1 = |key1| points is used, where
        # the sign of key1 determines the type of sample,
        # - key1 > 0, use a Korobov quasi-random sample,
        # - key1 < 0, use a "standard" sample (a Sobol quasi-random sample if seed = 0,
        # otherwise a pseudo-random sample).

        default_opts['key2'] = 1
        # *key2*: determines sampling in the final integration phase:
        # key2 = 7, 9, 11, 13 selects the cubature rule of degree key2. Note that the degree-11
        # rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
        # For other values of key2, a quasi-random sample is used, where the sign of key2
        # determines the type of sample,
        # - key2 > 0, use a Korobov quasi-random sample,
        # - key2 < 0, use a 'standard' sample (see description of key1 above),
        # and n2 = |key2| determines the number of points,
        # - n2 >= 40, sample n2 points,
        # - n2 < 40, sample n2 need points, where nneed is the number of points needed to
        # reach the prescribed accuracy, as estimated by Divonne from the results of the
        # partitioning phase.

        default_opts['key3'] = 1
        # *key3* sets the strategy for the refinement phase:
        # key3 = 0, do not treat the subregion any further.
        # key3 = 1, split the subregion up once more.
        # Otherwise, the subregion is sampled a third time with key3 specifying the sampling
        # parameters exactly as key2 above.

        default_opts['max_pass'] = 5
        # *maxpass*: controls the thoroughness of the partitioning phase: The
        # partitioning phase terminates when the estimated total number of integrand evalu-
        # ations (partitioning plus final integration) does not decrease for maxpass successive
        # iterations.
        # A decrease in points generally indicates that Divonne discovered new structures of
        # the integrand and was able to find a more effective partitioning. maxpass can be
        # understood as the number of 'safety' iterations that are performed before the par-
        # tition is accepted as final and counting consequently restarts at zero whenever new
        # structures are found.

        default_opts['border'] = 0.
        # *border*: the width of the border of the integration region.
        # Points falling into this border region will not be sampled directly, but will be extrap-
        # olated from two samples from the interior. Use a non-zero border if the integrand
        # subroutine cannot produce values directly on the integration boundary.
        
        default_opts['maxchisq'] = 10.
        # *maxchisq*: the maximum chisquare value a single subregion is al-
        # lowed to have in the final integration phase. Regions which fail this chisquare test and whose
        # sample averages differ by more than mindeviation move on to the refinement phase.
 
        default_opts['min_deviation'] = .25 
        # *mindeviation*: a bound, given as the fraction of the re-
        #  quested error of the entire integral, which determines whether it is worthwhile fur-
        # ther examining a region that failed the chisquare test. Only if the two sampling averages
        # obtained for the region differ by more than this bound is the region further treated.

    
        # *ldxgiven*: the leading dimension of xgiven, i.e. the offset between one
        # point and the next in memory.
        default_opts['ldxgiven'] = None
        default_opts['x_given'] = None
        # *xgiven(ldxgiven,n_given)*: a list of points where the inte-
        # grand might have peaks. Divonne will consider these points when partitioning the
        # integration region. The idea here is to help the integrator find the extrema of the in-
        # tegrand in the presence of very narrow peaks. Even if only the approximate location
        # of such peaks is known, this can considerably speed up convergence.

        default_opts['n_extra'] = 0
        # *nextra*: the maximum number of extra points the peak-finder subrou-
        # tine will return. If nextra is zero, peakfinder is not called and an arbitrary object
        # may be passed in its place, e.g. just 0.

        default_opts['peak_finder'] = None
        # *peakfinder*: the peak-finder subroutine. This subroutine is called
        # whenever a region is up for subdivision and is supposed to point out possible peaks
        # lying in the region, thus acting as the dynamic counterpart of the static list of points
        # supplied in xgiven. It is expected to be declared as
        # def peakfinder(ndim, b, n, x):
        # The bounds of the subregion are passed in the array b, where b(1,d ) is the lower and
        # b(2,d ) the upper bound in dimension d . On entry, n specifies the maximum number
        # of points that may be written to x. On exit, n must contain the actual number of
        # points in x.


        # Parameter relevant for Cuhre integration method
        # ------------------------------------------------

        default_opts['key'] = 0
        # *key* chooses the basic integration rule:
        # key = 7, 9, 11, 13 selects the cubature rule of degree key. Note that the degree-11
        # rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
        # For other values, the default rule is taken, which is the degree-13 rule in 2 dimensions,
        # the degree-11 rule in 3 dimensions, and the degree-9 rule otherwise.

        # Number of min and max evaluations for the integration with Cuhre
        default_opts['min_eval'] = 0
        default_opts['max_eval'] = 50000

        # Set instance attributes options
        for opt in default_opts:
            if opt in opts:
                setattr(self,opt,opts.pop(opt))
            else:
                setattr(self,opt,default_opts[opt])

        super(pyCubaIntegrator,self).__init__(integrands, **opts)

        supported_algorithms = ['Vegas', 'Suave', 'Divonne', 'Cuhre']
        if self.algorithm not in supported_algorithms:
            raise IntegratorError("The pyCubaIntegrator only supports the following algorithms: [ %s ]"%\
                                                                      (', '.join(supported_algorithms)))

        # After each integrations, results will be available from here
        self.last_integration_results = None
    
    @staticmethod
    def wrapped_integrand(ndim, xx, ncomp, ff, userdata):
        """ Function to wrap the integrand to the pyCuba standards.
        IntegratorInstance is taken from a global which is set only during the
        time of the integration with the low-level language integrator implementation."""
        try:
            n_integrands_to_consider = ncomp.contents.value
        
            for integrand_index, integrand in enumerate(IntegratorInstance.integrands[:n_integrands_to_consider]):
                dimensions = integrand.continuous_dimensions
                input_point = np.array([dimensions[i].lower_bound + xx[i]*(dimensions[i].upper_bound - dimensions[i].lower_bound) 
                                                                                         for i in range(ndim.contents.value)])
                ff[integrand_index] = IntegratorInstance.integrands[integrand_index](input_point, np.array([], dtype=int))
                # All integration methods here assume unit phase-space volume, so we need to rescale our evaluations here.
                ff[integrand_index] *= IntegratorInstance.integrands[integrand_index].get_dimensions().volume()

            return 0
        except KeyboardInterrupt:
            logger.warning('Integration with pyCuba aborted by user.')
            sys.exit(1)

    @staticmethod
    def wrapped_integrand_divonne(ndim, xx, ncomp, ff, phase):
        """ Wrapper for divonne which can specify the argument "phase".
            0, sampling of the points in xgiven,
            1, partitioning phase,
            2, final integration phase,
            3, refinement phase.
        This information might be useful if the integrand takes long to compute and a sufficiently
        accurate approximation of the integrand is available. The actual value of the integral is only
        of minor importance in the partitioning phase, which is instead much more dependent on
        the peak structure of the integrand to find an appropriate tessellation. An approximation
        which reproduces the peak structure while leaving out the fine details might hence be a
        perfectly viable and much faster substitute when phase < 2."""
        try:        
            return pyCubaIntegrator.wrapped_integrand(ndim, xx, ncomp, ff, None)
        except KeyboardInterrupt:
            print 'Integration with pyCuba aborted by user.'
            sys.exit(1)

    def get_name(self):
        """ Returns the integrator name."""
        return self.__class__.__name__+'@%s'%self.algorithm

    def aggregate_results(self):
        """ Aggregate results from a MC run into a global summed contribution and error."""
        
        integral = 0.0
        error    = 0.0
        for res in self.last_integration_results['results']:
            integral += res['integral']
            error    += res['error']
        return (integral, error)

    def integrate(self):
        """ Return the final integral and error estimators. """

        for integrand in self.integrands:
            if len(integrand.discrete_dimensions)>0:
                raise IntegratorError("pyCuba integrator does not support discrete dimensions for now.")

        if len(set([len(integrand.continuous_dimensions) for integrand in self.integrands]))!=1:
            raise IntegratorError("pyCuba only supports multiple integrands with all the same number of dimensions.")

        if self.nb_core > 1:
            logger.warning(
"""\n\n-------
Parallelization with pyCuba is handled directly by the Cuba C++ library.
This is potentially dangerous because python does not know that the wrapped integrands are called
in parallel and it may overwrite memory. In practice, this seems however not to be a problem for a wide
range of applications. To be safe however, consider making sure that your result match what is obtained
with less statistics using one core only.
-------\n""")
            os.environ['CUBACORES'] = '%d'%self.nb_core
        else:
            os.environ['CUBACORES'] = '0'

        n_dimensions = len(self.integrands[0].continuous_dimensions)

        if self.algorithm == 'Vegas':
            with tmpGlobal('IntegratorInstance',self):
                self.last_integration_results = pycuba.Vegas(
                    pyCubaIntegrator.wrapped_integrand, 
                    n_dimensions,
                    ncomp       = len(self.integrands),
                    verbose     = self.verbosity,
                    userdata    = 0,
                    epsrel      = self.target_accuracy,                    
                    nstart      = self.n_start,
                    nincrease   = self.n_increase,
                    nbatch      = self.n_batch,
                    seed        = self.seed,
                    mineval     = self.min_eval,
                    maxeval     = self.max_eval
                )

            return self.aggregate_results()

        elif self.algorithm == 'Suave':
            with tmpGlobal('IntegratorInstance',self):
                self.last_integration_results = pycuba.Suave(
                    pyCubaIntegrator.wrapped_integrand, 
                    n_dimensions,
                    ncomp       = len(self.integrands),
                    verbose     = self.verbosity,
                    userdata    = 0,
                    epsrel      = self.target_accuracy,                    
                    nnew        = self.n_new,
                    nmin        = self.n_min,
                    flatness    = self.flatness,
                    mineval     = self.min_eval,
                    maxeval     = self.max_eval
                )
            return self.aggregate_results()

        elif self.algorithm == 'Divonne':
            with tmpGlobal('IntegratorInstance',self):
                self.last_integration_results = pycuba.Divonne(
                    pyCubaIntegrator.wrapped_integrand_divonne, 
                    n_dimensions,
                    self.key1,
                    self.key2,
                    self.key3,
                    self.max_pass,
                    self.border,
                    self.maxchisq,
                    self.min_deviation,
                    epsrel          = self.target_accuracy,
                    ldxgiven        = self.ldxgiven,
                    xgiven          = self.x_given,
                    nextra          = self.n_extra,
                    peakfinder      = self.peak_finder,
                    ncomp           = len(self.integrands),
                    verbose         = self.verbosity,
                    userdata        = 0,
                    mineval     = self.min_eval,
                    maxeval     = self.max_eval)
            return self.aggregate_results()

        elif self.algorithm == 'Cuhre':
            with tmpGlobal('IntegratorInstance',self):
                self.last_integration_results = pycuba.Cuhre(
                    pyCubaIntegrator.wrapped_integrand, 
                    n_dimensions,
                    key         = self.key, 
                    mineval     = self.min_eval,
                    maxeval     = self.max_eval,
                    epsrel      = self.target_accuracy,                
                    ncomp       = len(self.integrands),
                    verbose     = self.verbosity,
                    userdata    = 0)

            return self.aggregate_results()

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
    verbosity = 0
    my_integrators = {}

    my_integrators['Naive MC'] = integrators.SimpleMonteCarloIntegrator(my_integrand, 
        n_iterations=10, n_points_per_iterations=1000, accuracy_target=None, verbosity = verbosity)

    my_integrators['VEGAS'] = pyCubaIntegrator(my_integrand, algorithm='Vegas', verbosity = verbosity)

    my_integrators['SUAVE'] = pyCubaIntegrator(my_integrand, algorithm='Suave', verbosity = verbosity)

    my_integrators['DIVONNE'] = pyCubaIntegrator(my_integrand, algorithm='Divonne', verbosity = verbosity)

    my_integrators['CUHRE'] = pyCubaIntegrator(my_integrand, algorithm='Cuhre', verbosity = verbosity)

    chosen_integrators = ['Naive MC', 'VEGAS', 'SUAVE', 'DIVONNE', 'CUHRE']

    for name in chosen_integrators:
        res = my_integrators[name].integrate()
        # Finally integrate
        print '\n'+'='*70
        print '%-30s = %.4e +/- %.2e'%('%s :: Final result'%name, res[0], res[1])
        print '='*70

