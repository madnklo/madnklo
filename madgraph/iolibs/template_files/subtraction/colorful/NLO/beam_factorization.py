##########################################################################################
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
##########################################################################################
"""Implementation of NLO beam_factorisation currents. These are the PDF counterterms as well
as the integrated initial state collinear counterterms."""

import os
import math
from madgraph.core.base_objects import EpsilonExpansion
import madgraph.various.misc as misc

try:
    # First try to import this in the context of the exported currents
    import SubtractionCurrents.subtraction_current_implementations_utils as utils
    import SubtractionCurrents.QCD_local_currents as currents
except ImportError:
    # If not working, then it must be within MG5_aMC context:
    import madgraph.iolibs.template_files.\
                   subtraction.subtraction_current_implementations_utils as utils
    import madgraph.iolibs.template_files.\
                   subtraction.QCD_local_currents as currents

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError

log = math.log
pi = math.pi

# All counterterms here adopt a xi-dependent distribution of the following form:
#
#        Counterterm(xi) = F_+(xi) + [F] \delta(xi-1)
#    (which can also be explicitely written)
#        Counterterm(xi) = F(xi) + {F(xi)} \delta(xi-1) + [F] \delta(xi-1)
#
#  where 'F' can either be a PDF counterterm or an interated collinear ISR counterterm.
#  Then each piece of the distribution is assigned a different value for its attribute
#  'distribution_type' as follows:
#
#     F(xi)   --> distribution_type = 'bulk'
#     {F(xi)} --> distribution_type = 'counterterm'
#     [F(xi)] --> distribution_type = 'endpoint'

#=========================================================================================
# PDF Counterterm
#=========================================================================================
class QCD_beam_factorization_F0(currents.QCDBeamFactorizationCurrent):
    """Implements the NLO QCD PDF counterterm of type F(xi)"""

    distribution_types_implemented_in_this_class = ['bulk','counterterm','endpoint']
    
    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None

        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that it involves exactly one F structure with one leg. 
        if len(ss.substructures)==0:
            factorization_structure = ss
        elif len(ss.substructures)==1 and len(ss.substructures[0].substructures)==0:
            factorization_structure = ss.substructures[0]
        else:      
            return None
        
        if factorization_structure.name() != 'F':
            return None
        if len(factorization_structure.legs) != 1:
            return None
        
        # Make sure the one leg of the F structure is initial-state
        if not cls.is_initial(factorization_structure.legs[0]): 
            return None

        # The current is valid (remember that this implements the PDF counterterm of 
        # all possible incoming flavors.
        return init_vars

    def evaluate_kernel(self, PS_point, process, xi, mu_r, mu_f, normalization, 
                                                  allowed_backward_evolved_flavors='ALL'):
        """ Return an instance of BeamFactorizationCurrentEvaluation, whose 'values' entry
        are dictionaries specifying the counterterm in flavor space, for the value of xi 
        specified in argument."""

        if allowed_backward_evolved_flavors != 'ALL':
            raise CurrentImplementationError('The current %s must always be called with'%self.__class__.__name__+
                "allowed_backward_evolved_flavors='ALL', not %s"%str(allowed_backward_evolved_flavors))

        # Only the order epsilon of the scales pre-factor matters here.
        prefactor = EpsilonExpansion({
                0 : 1.,
                1 : log(mu_r**2 / mu_f**2)
        })
        prefactor *= EpsilonExpansion({-1:1.})*normalization

        # Assign a fake xi for now if the distribution type is 'endpoint'
        # TODO: this is not optimal, eventually we should put each of these three pieces in
        # separate currents
        if self.distribution_type == 'endpoint':
            xi = 0.5

        # Define the NLO QCD PDF counterterms kernels
        kernel_gg = { 
            'bulk'          :     prefactor*( 
                2.*self.CA*( 1./ (1.-xi) + (1.-xi)/xi -1. + xi*(1-xi) ) 
            ),
            'counterterm'   :     prefactor*( 2.*self.CA / (1.-xi) ),
            'endpoint'      :     prefactor*( 11./6.*self.CA - 2./3.*self.NF*self.TR)
        }
                 
        kernel_gq = { 
            'bulk'          :     prefactor*( self.CF*(1.+(1.-xi)**2)/xi ),
            'counterterm'   :     None,
            'endpoint'      :     None
        }
        
        kernel_qg = { 
            'bulk'          :     prefactor*( self.TR*(xi**2 + (1.-xi)**2) ),
            'counterterm'   :     None,
            'endpoint'      :     None
        }
        
        kernel_qq = { 
            'bulk'          :     prefactor*( self.CF*((1.+xi**2)/(1.-xi)) ),
            'counterterm'   :     prefactor*( self.CF*((1.+xi**2)/(1.-xi)) ),
            'endpoint'      :     None
        }

        active_quark_PDGs = tuple([pdg for pdg in range(1,7)+range(-1,-7,-1) 
                                                                 if pdg in self.beam_PDGs])
        
        # Build the NLO flavor matrix
        flavor_matrix = {}
        for reduced_flavor in self.beam_PDGs:
            # Gluon backward evolution
            if reduced_flavor==21:
                gluon_dict = {}
                if kernel_gg[self.distribution_type] is not None:
                    gluon_dict[(21,)] = kernel_gg[self.distribution_type]
                if active_quark_PDGs and kernel_gq[self.distribution_type] is not None:
                    gluon_dict[active_quark_PDGs] = kernel_gq[self.distribution_type]
                if gluon_dict:
                    flavor_matrix[21] = gluon_dict
            
            # Quark backward evolution            
            if reduced_flavor in active_quark_PDGs:
                quark_dict = {}
                if kernel_qg[self.distribution_type] is not None:
                    quark_dict[(21,)] = kernel_qg[self.distribution_type]
                if kernel_qq[self.distribution_type] is not None:
                    quark_dict[(reduced_flavor,)] = kernel_qq[self.distribution_type]
                if quark_dict:                    
                    flavor_matrix[reduced_flavor] = quark_dict

        # Truncate all entries of the flavor matrix so as to remove irrelevant O(\eps) terms
        for flav_in, flav_outs in flavor_matrix.items():
            for flav_out, eps_expansion in flav_outs.items():
                eps_expansion.truncate(max_power=0)

        # Now assign the flavor matrix in the BeamFactorizationCurrentEvaluation instance
        evaluation = utils.BeamFactorizationCurrentEvaluation({
            'spin_correlations'         : [None,],
            'color_correlations'        : [None,],
            'values'                    : { (0,0) : flavor_matrix }
        })

        return evaluation

#=========================================================================================
# PDF integrated initial-state single collinear counterterm
#=========================================================================================
class QCD_beam_factorization_single_collinear(currents.QCDBeamFactorizationCurrent):
    """Implements the NLO QCD initial-state single collinear integgratated counterterm of type F(xi)"""

    distribution_types_implemented_in_this_class = ['bulk','counterterm','endpoint']
    
    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None

        # Retrieve singular structure
        ss = current.get('singular_structure')

        # Check that it involves exactly one collinear structure with two legs. 
        if len(ss.substructures)!=1:
            return None
        
        collinear_structure = ss.substructures[0]
        if collinear_structure.name() != 'C':
            return None
        if len(collinear_structure.legs) != 2:
            return None
        
        # Make sure that one of the two legs of the C structure is initial-state
        if not any(cls.is_initial(leg) for leg in collinear_structure.legs):
            return None

        # The current is valid (remember that this implements the integrated
        # initial state collinear counterterm of all possible incoming flavors.
        return init_vars

    def evaluate_kernel(self, PS_point, process, xi, mu_r, mu_f, normalization,
                                                   allowed_backward_evolved_flavors='ALL'):
        """ Return an instance of BeamFactorizationCurrentEvaluation, whose 'values' entry
        are dictionaries specifying the counterterm in flavor space, for the value of xi 
        specified in argument."""

        # Obtain Q_square
        Q        = sum([PS_point[l.get('number')] for l in process.get_initial_legs()])
        Q_square = Q.square()

        # Only up to the order epsilon^2 of the scales prefactor matters here.
        logMuQ = log(mu_r**2/Q_square)
        prefactor = EpsilonExpansion({ 0 : 1., 1 : logMuQ, 2 : 0.5*logMuQ**2 })
        prefactor *= self.SEpsilon*normalization

        # The additional 1/x part of the prefactor is included later during the PDF
        # convolution of the event (using its 'Bjorken rescaling' attribute) because
        # we must make sure that the plus distribution hits on it.
        # Also, the same 1/x appears in the PDF counterterms as a result of the change
        # of variable necessary to bring them in the form where the plus distribution 
        # only acts on the PDF. So it makes sense to keep it completely factorised.

        # Input variables
        y_0 = currents.SomogyiChoices.y_0
        logy0 = log(y_0)
        # Assign a fake x for now if the distribution type is 'endpoint'
        # TODO: this is not optimal, eventually we should put each of these three pieces in
        # separate currents
        if self.distribution_type == 'endpoint':
            x = 0.5
        else:
            x  = xi

        log1mx = log(1.-x)
        
        # Heaviside
        theta_x_1my0 = 1. if (x-(1-y_0)) >= 0. else 0.
        theta_1my0_x = 1. if ((1-y_0)-x) >= 0. else 0.

        # Define the NLO QCD integrate initial-state single collinear counterterms kernels
        color_factor = self.CA
        kernel_gg = { 
            'bulk'          :     prefactor*color_factor*(EpsilonExpansion({ 
                -1 : -2.*( 1./(1.-x) + (1.-x)/x - 1 + x*(1-x) ),
                 0 : (2.*log1mx / (1.-x))*(1.+theta_x_1my0) + (2.*logy0/(1.-x))*theta_1my0_x
                     + 2.*( ((1.-x)/x) -1. + x*(1.-x) )*( log1mx*(1.+theta_x_1my0) + logy0*theta_1my0_x )
            })),
            'counterterm'   :     prefactor*color_factor*(EpsilonExpansion({ 
                -1 : -2.* ( 1./(1.-x) )  ,
                 0 : (2.*log1mx / (1.-x))*(1.+theta_x_1my0)  ,
            })),
            'endpoint'      :     prefactor*color_factor*(EpsilonExpansion({
                -2 : 1.  ,
                -1 : 0.  ,
                 0 : -(math.pi**2/6.) + logy0**2
            }))
        }
        
        color_factor = self.CA        
        kernel_gq = { 
            'bulk'          :     prefactor*color_factor*(EpsilonExpansion({ 
                -1 : -(self.CF/self.CA)*(1.+(1.-x)**2) / x  ,
                 0 : (self.CF/self.CA)*( ((1.+(1.-x)**2)/x)*( log1mx*(1.+theta_x_1my0) + logy0*theta_1my0_x ) + x )
            })),
            'counterterm'   :     None,
            'endpoint'      :     None
        }
        
        color_factor = self.CF
        kernel_qg = { 
            'bulk'          :     prefactor*color_factor*(EpsilonExpansion({ 
                -1 : -(self.TR/self.CF)*(x**2+(1.-x)**2)  ,
                 0 : (self.TR/self.CF)*( (x**2 + (1.-x)**2)*( log1mx*(1.+theta_x_1my0) + logy0*theta_1my0_x ) + 2.*x*(1.-x) )
            })),
            'counterterm'   :     None,
            'endpoint'      :     None
        }
        
        color_factor = self.CF
        kernel_qq = { 
            'bulk'          :     prefactor*color_factor*(EpsilonExpansion({ 
                -1 : -((1.+x**2)/(1.-x))  ,
                 0 : (2.*log1mx / (1.-x))*(1.+theta_x_1my0) + (2.*logy0/(1.-x))*theta_1my0_x
                     - ( (1.+x)*( log1mx*(1.+theta_x_1my0)+logy0*theta_1my0_x ) -1.+x )
            })),
            'counterterm'   :     prefactor*color_factor*(EpsilonExpansion({ 
                -1 : -((1.+x**2)/(1.-x))  ,
                 0 : (2.*log1mx / (1.-x))*(1.+theta_x_1my0)  ,
            })),
            'endpoint'      :     prefactor*color_factor*(EpsilonExpansion({ 
                -2 : 1.     ,
                -1 : 3./2.  ,
                 0 : -(math.pi**2/6.) + logy0**2
            }))
        }

        active_quark_PDGs = tuple([pdg for pdg in range(1,7)+range(-1,-7,-1) 
                                                                 if pdg in self.beam_PDGs])

        # Build the NLO flavor matrix
        flavor_matrix = {}
        for reduced_flavor in self.beam_PDGs:
            # Gluon backward evolution
            if reduced_flavor==21:
                gluon_dict = {}
                if kernel_gg[self.distribution_type] is not None:
                    gluon_dict[(21,)] = kernel_gg[self.distribution_type]
                if active_quark_PDGs and kernel_gq[self.distribution_type] is not None:
                    gluon_dict[active_quark_PDGs] = kernel_gq[self.distribution_type]
                if gluon_dict:
                    flavor_matrix[21] = gluon_dict
            
            # Quark backward evolution            
            if reduced_flavor in active_quark_PDGs:
                quark_dict = {}
                if kernel_qg[self.distribution_type] is not None:
                    quark_dict[(21,)] = kernel_qg[self.distribution_type]
                if kernel_qq[self.distribution_type] is not None:
                    quark_dict[(reduced_flavor,)] = kernel_qq[self.distribution_type]
                if quark_dict:                    
                    flavor_matrix[reduced_flavor] = quark_dict

        # Truncate all entries of the flavor matrix so as to remove irrelevant O(\eps) terms
        for flav_in, flav_outs in flavor_matrix.items():
            for flav_out, eps_expansion in flav_outs.items():
                eps_expansion.truncate(max_power=0)

        # Now apply the mask 'allowed_backward_evolved_flavors' if not set to 'ALL'
        if allowed_backward_evolved_flavors != 'ALL':
            filtered_flavor_matrix = { reduced_flavor: 
                { end_flavor : wgt for end_flavor, wgt in  flavor_matrix[reduced_flavor].items() if 
                  end_flavor in allowed_backward_evolved_flavors
                }  for reduced_flavor in flavor_matrix }

        # Now assign the flavor matrix in the BeamFactorizationCurrentEvaluation instance
        evaluation = utils.BeamFactorizationCurrentEvaluation({
            'spin_correlations'         : [None,],
            'color_correlations'        : [None,],
            'values'                    : { (0,0) : filtered_flavor_matrix }
        })

        return evaluation

#=========================================================================================
# Integrated single soft counterterm when recoiling equally against the two incoming beams
# This implementation is only used for the colorful currents scheme when used together with
# mapping like ppToOneWalker where the soft momentum recoils democratically against the two
# incoming beams
#=========================================================================================
class QCD_beam_factorization_single_soft(currents.QCDBeamFactorizationCurrent):
    """Implements the NLO QCD initial-state single soft integgratated counterterm."""

    distribution_types_implemented_in_this_class = ['bulk','counterterm','endpoint']

    # These integrated contributions are not really directly related to the physical
    # properties of beam factorization (for instance they don't act on the flavor space) and
    # therefore apply independely of it.
    beam_types_implemented_in_this_class = 'ALL'
    beam_PDGs_implemented_in_this_class = 'ALL'

    def __init__(self, *args, **opts):
        super(QCD_beam_factorization_single_soft, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None:
            return None

        # Check the structure is an integrated simple soft
        if len(current.get('singular_structure').substructures)==0:
            return None
        singular_structure = current.get('singular_structure').substructures[0]
        if singular_structure.name() != 'S': 
            return None
        if len(singular_structure.substructures)>0:
            return None
        
        # All checks passed
        return init_vars

    def evaluate_kernel(self, PS_point, process, xi, mu_r, mu_f, normalization,
                                                    allowed_backward_evolved_flavors='ALL'):
        """ Return an instance of SubtractionCurrentEvaluation, whose 'values' entry
        are simple EpsilonExpansions since soft-integrated counterterms convoluted in a 
        correlated fashion with the initial state beams *cannot* act in flavor space."""

        if allowed_backward_evolved_flavors != 'ALL':
            raise CurrentImplementationError('The current %s must always be called with'%self.__class__.__name__+
                "allowed_backward_evolved_flavors='ALL', not %s"%str(allowed_backward_evolved_flavors))

        if process is None:
            raise CurrentImplementationError(self.name() + " requires a reduced_process.")

        # Now find all colored leg numbers in the reduced process
        all_colored_parton_numbers = []
        for leg in process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color')==1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [],
            'values'              : {}
        })

        # Obtain Q_square
        Q        = sum([PS_point[l.get('number')] for l in process.get_initial_legs()])
        Q_square = Q.square()

        # Only up to the order epsilon^2 of the scales prefactor matters here.
        logMuQ = log(mu_r**2/Q_square)
        prefactor = EpsilonExpansion({ 0 : 1., 1 : logMuQ, 2 : 0.5*logMuQ**2 })
        prefactor *= self.SEpsilon*normalization

        #TODO Look at how this works beyond Initial-initial
        # For II: Soft +CS(beam 1) + BS(beam 2) = -1/2 Soft
        prefactor *= -1./2.


        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        for i, a in enumerate(all_colored_parton_numbers):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in all_colored_parton_numbers[i:]:
                # Write the integrated eikonal for that pair
                if a!=b:
                    mult_factor = 2.
                else:
                    mult_factor = 1.
                evaluation['color_correlations'].append( ((a, b), ) )
                
                pa = PS_point[a]
                pb = PS_point[b]
                if self.distribution_type == 'bulk':
                    # TODO Place-holder example
                    kernel = EpsilonExpansion({
                                0 : - 16.*xi * log(1.-xi**2) /(1.-xi**2),
                                -1 : 8. * xi / (1.-xi**2),
                                -2 : 0.
                    })
                elif self.distribution_type == 'counterterm':
                    # TODO Place-holder example, but it indeed regulates the above example
                    kernel = EpsilonExpansion({
                                0 : -8.*log(2.*(1.-xi))/(1.-xi),
                                -1 : 4./(1.-xi),
                                -2 : 0.
                    })
                elif self.distribution_type == 'endpoint':
                    # TODO Place-holder example
                    kernel = EpsilonExpansion({
                                0 : pi**2./3.-4.*log(2.)**2,
                                -1 : 4.*log(2.),
                                -2 : -2.
                    })
                else:
                    raise CurrentImplementationError("Distribution type '%s' not supported."
                                                                    %self.distribution_type)

                evaluation['values'][(0, color_correlation_index)] = kernel*normalization
                color_correlation_index += 1

        return evaluation

#=========================================================================================
# PDF integrated initial-state single soft-collinear counterterm
#=========================================================================================
class QCD_beam_factorization_single_softcollinear(currents.QCDBeamFactorizationCurrent):
    """Implements the NLO QCD initial-state single soft-collinear integgratated counterterm
    of type F(xi). These are zero here since they have already been accounted for
    in the soft counterterms."""

    distribution_types_implemented_in_this_class = ['bulk','counterterm','endpoint']

    # These integrated contributions are not really directly related to the physical
    # properties of beam factorization (for instance they don't act on the flavor space) and
    # therefore apply independely of it.
    beam_types_implemented_in_this_class = 'ALL'
    beam_PDGs_implemented_in_this_class = 'ALL'

    # The soft-collinear integrated counterterm has been accounted for completely in the 
    # soft integrated counterterm
    is_zero = True

    def __init__(self, *args, **opts):
        
        # Make sure it is initialized with the proper set of options and remove them
        # before calling the mother constructor
        if 'color_charge' not in opts:
            raise CurrentImplementationError(
                "The current '%s' must be instantiated with "%self.__class__.__name__+
                                                " a 'color_charge' option specified.")
        color_charge = opts.pop('color_charge')
        
        super(QCD_beam_factorization_single_softcollinear, self).__init__(*args, **opts)
        self.supports_helicity_assignment = False
        # At this state color_charge is the string of the argument to retrieve ('CA' or 'CF')
        # And now that the mother constructor is called, the group factors have been initialized
        # and we can retrieve them.
        self.color_charge = getattr(self, color_charge)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None:
            return None

        # Retrieve singular structure
        ss = current.get('singular_structure')

        # Check that it involves exactly one collinear structure with two legs. 
        if len(ss.substructures)!=1:
            return None
        
        # Finally check that the singular structure and PDG matches
        singular_structure = current.get('singular_structure').substructures[0]
        
        # It main structure should be of collinear type
        if singular_structure.name()!='C':
            return None        
        
        # It should have only one leg left, the other one being in the nested soft structure
        # It must be an initial-state leg.
        if len(singular_structure.legs)!=1:
            return None

        # The leg not soft must be quark or a gluon      
        if not abs(singular_structure.legs[0].pdg) in [21,]+range(1,7):
            return None

        # It should have exactly one nested structures
        if len(singular_structure.substructures)!=1:
            return None
        
        sub_singular_structure = singular_structure.substructures[0]
        
        # Make sure this substructure is soft
        if sub_singular_structure.name()!='S':
            return None
        
        # Make sure it contains a single soft leg
        if len(sub_singular_structure.legs)!=1:
            return None
        
        soft_leg = sub_singular_structure.legs[0]
        
        # Make sure the soft leg is massless final and a gluon
        if model.get_particle(soft_leg.pdg).get('mass').upper()!='ZERO':
            return None
        if soft_leg.pdg != 21:
            return None
        
        # We now know that this current is implemented here. We return
        # the specific color charge to instantiate this kernel with,
        # in the form of a the name of the group factor to retrieve upon
        # initialization.
        if singular_structure.legs[0].pdg == 21:
            # This is a 'g > g g' soft-collinear splitting
            init_vars['color_charge'] = 'CA'
        else:
            # This is a 'q > g g' soft-collinear splitting
            init_vars['color_charge'] = 'CA'

        return init_vars
