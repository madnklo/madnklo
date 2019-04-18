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
"""Implementation of NLO beam_factorization currents. These are the PDF counterterms as well
as the integrated initial state collinear counterterms."""

import os
import math
from madgraph.core.base_objects import EpsilonExpansion
import madgraph.various.misc as misc

import commons.utils as utils
import commons.QCD_local_currents as currents

from commons.integrated_current_expressions import HE

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

        # If this is a BF current it will not have substructures
        singular_structure = current.get('singular_structure')
        if len(singular_structure.substructures)==0:
            return None
        singular_structure = singular_structure.substructures[0]

        # Check the structure is an integrated simple soft
        if len(singular_structure.substructures)==0:
            return None
        singular_structure = singular_structure.substructures[0]

        if singular_structure.name() != 'S': 
            return None
        if len(singular_structure.substructures)>0:
            return None
        
        # All checks passed
        return init_vars

    def evaluate_kernel(self, PS_point, process, xi, mu_r, mu_f, Q, normalization,
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
        # also find the initial state colored partons to tag initial/final eikonals
        all_colored_parton_numbers = []
        colored_initial_parton_numbers = []
        all_initial_numbers = [l.get('number') for l in process.get_initial_legs()]
        for leg in process.get('legs'):
            leg_number =leg.get('number')
            if self.model.get_particle(leg.get('id')).get('color')==1:
                continue
            all_colored_parton_numbers.append(leg_number)
            if leg_number in all_initial_numbers:
                colored_initial_parton_numbers.append(leg.get('number'))

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [],
            'values'              : {}
        })

        # Obtain Q_square
        #Q        = sum([PS_point[l.get('number')] for l in process.get_initial_legs()])
        Q_square = Q.square()

        # Only up to the order epsilon^2 of the scales prefactor matters here.
        logMuQ = log(mu_r**2/Q_square)
        # Correction for the counterterm: in BS (bulk+counterterm), the variable Q_square corresponds to that
        # of the real event. However the counterterm corresponds to the residue of the bulk at xi=1.
        # This is effectively obtained by multiplying by xi: Q_residue = Q_real * xi.
        # Note for future dumb-me: log(mu_r**2/(Q_square*xi**2)) = logMuQ - log(xi**2)
        if self.distribution_type == 'counterterm':
            logMuQ-=log(xi**2)        
        prefactor = EpsilonExpansion({ 0 : 1., 1 : logMuQ, 2 : 0.5*logMuQ**2 })
        prefactor *= normalization

        color_correlation_index = 0
        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        for i, a in enumerate(all_colored_parton_numbers):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in all_colored_parton_numbers[i+1:]:
                # Write the integrated eikonal for that pair

                evaluation['color_correlations'].append( ((a, b), ) )

                pa = PS_point[a]
                pb = PS_point[b]

                # We can only handle massless particles
                try:
                    assert pa.square()/Q_square < 1.e-09
                except AssertionError:
                    misc.sprint("No massive particles in soft currents for now")
                    raise

                # Assign the type of dipole
                # dipole_type = [bool_a, bool_b] tags the dipole tag, with True indicating initial state and False final
                dipole_type = [a in colored_initial_parton_numbers, b in colored_initial_parton_numbers]

                if all(dipole_type): # Initial-initial
                    # Initial-initial: S+CS = 0
                    if self.distribution_type == 'bulk':
                        kernel = EpsilonExpansion({0:0})
                    elif self.distribution_type == 'counterterm':
                        kernel = EpsilonExpansion({0:0})
                    elif self.distribution_type == 'endpoint':
                        kernel = EpsilonExpansion({0:0})
                    else:
                        raise CurrentImplementationError("Distribution type '%s' not supported."
                                                                        %self.distribution_type)
                else: # At least one leg final
                    # The integrated counterterms are evaluated in terms of
                    # dipole_invariant = 1-cos(angle between the dipole momenta)
                    dipole_invariant = 0.5*pa.dot(pb)*Q.square()/(pa.dot(Q)*pb.dot(Q))
                    if self.distribution_type == 'bulk':
                        #The factor xi^2 below corrects the flux factor used in the bulk BS which has a 1/xi^2 too many
                        #A more permanent change is warranted after testing.
                        #See github issue #9 for reference 
                        kernel = EpsilonExpansion({0:xi**2*HE.integrated_bs_bulk_finite(dipole_invariant,xi)})
                    elif self.distribution_type == 'counterterm':
                        kernel = EpsilonExpansion({0:HE.integrated_bs_counterterm_finite(dipole_invariant,xi)})
                    elif self.distribution_type == 'endpoint':
                        kernel = EpsilonExpansion({-1:HE.integrated_bs_endpoint_pole(dipole_invariant),
                                                    0:HE.integrated_bs_endpoint_finite(dipole_invariant)})
                    else:
                        raise CurrentImplementationError("Distribution type '%s' not supported."
                                                                        %self.distribution_type)

                    # Former implementation of the II soft+SC. Commented by Nicolas
                    # While no longer useful, this is kept for now to remember how non-zero integrated soft shoud be implemented
                    # if self.distribution_type == 'bulk':
                    #     kernel = EpsilonExpansion({
                    #                 0 : - 16.*xi * log(1.-xi**2) /(1.-xi**2),
                    #                 -1 : 8. * xi / (1.-xi**2),
                    #                 -2 : 0.
                    #     })
                    # elif self.distribution_type == 'counterterm':
                    #     kernel = EpsilonExpansion({
                    #                 0 : -8.*log(2.*(1.-xi))/(1.-xi),
                    #                 -1 : 4./(1.-xi),
                    #                 -2 : 0.
                    #     })
                    # elif self.distribution_type == 'endpoint':
                    #     kernel = EpsilonExpansion({
                    #                 0 : pi**2./3.-4.*log(2.)**2,
                    #                 -1 : 4.*log(2.),
                    #                 -2 : -2.
                    #     })
                    # else:
                    #     raise CurrentImplementationError("Distribution type '%s' not supported."
                    #                                                     %self.distribution_type)

                evaluation['values'][(0, color_correlation_index)] = kernel*prefactor
                color_correlation_index += 1

        return evaluation

