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

# All counterterms here adopt a chsi-dependent distribution of the following form:
#
#        Counterterm(\chsi) = F_+(\chsi) + [F] \delta(\chsi-1)
#    (which can also be explicitely written)
#        Counterterm(\chsi) = F(\chsi) + {F(\chsi)} \delta(\chsi-1) + [F] \delta(\chsi-1)
#
#  where 'F' can either be a PDF counterterm or an interated collinear ISR counterterm.
#  Then each piece of the distribution is assigned a different value for its attribute
#  'distribution_type' as follows:
#
#     F(\chsi)   --> distribution_type = 'bulk'
#     {F(\chsi)} --> distribution_type = 'counterterm'
#     [F(\chsi)] --> distribution_type = 'endpoint'

#=========================================================================================
# PDF Counterterm
#=========================================================================================
class QCD_beam_factorization_F0(currents.QCDBeamFactorizationCurrent):
    """Implements the NLO QCD PDF counterterm of type F(\chsi)"""

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

    def evaluate_kernel(self, chsi, normalization):
        """ Return an instance of BeamFactorizationCurrentEvaluation, whose 'values' entry
        are dictionaries specifying the counterterm in flavor space, for the value of chsi 
        specified in argument."""

        # Define the NLO QCD PDF counterterms kernels
        kernel_gg = { 
            'bulk' :  EpsilonExpansion({
                0 : 1.
            })*normalization,
            'counterterm' :  EpsilonExpansion({
                0 : 1.
            })*normalization,
            'endpoint' :  EpsilonExpansion({
                0 : 1.
            })*normalization
        }
                 
        kernel_gq = { 
            'bulk' :   EpsilonExpansion({
                0 : 1.
            })*normalization,
            'counterterm' :  EpsilonExpansion({
                0 : 1.
            })*normalization,
            'endpoint' :  EpsilonExpansion({
                0 : 1.
            })*normalization
        }
        
        kernel_qg = { 
            'bulk' :   EpsilonExpansion({
                0 : 1.
            })*normalization,
            'counterterm' :  EpsilonExpansion({
                0 : 1.
            })*normalization,
            'endpoint' :  EpsilonExpansion({
                0 : 1.
            })*normalization
        }
        
        kernel_qq = { 
            'bulk' :   EpsilonExpansion({
                0 : 1.
            })*normalization,
            'counterterm' :  EpsilonExpansion({
                0 : 1.
            })*normalization,
            'endpoint' :  EpsilonExpansion({
                0 : 1.
            })*normalization
        }

        active_quark_PDGs = tuple([pdg for pdg in range(1,7)+range(-1,-7,-1) 
                                                                 if pdg in self.beam_PDGs])
        
        # Build the NLO flavor matrix
        flavor_matrix = {}
        for reduced_flavor in self.beam_PDGs:
            # Gluon backward evolution
            if reduced_flavor==21:
                flavor_matrix[21] = { (21,) : kernel_gg[self.distribution_type] }
                if active_quark_PDGs:
                    flavor_matrix[21][active_quark_PDGs] = kernel_gq[self.distribution_type]
            
            # Quark backward evolution            
            if reduced_flavor in active_quark_PDGs:
                flavor_matrix[reduced_flavor] = { 
                    (21,) : kernel_qg[self.distribution_type],
                    (reduced_flavor,) : kernel_qq[self.distribution_type]
                }

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
    """Implements the NLO QCD initial-state single collinear integgratated counterterm of type F(\chsi)"""

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

    def evaluate_kernel(self, chsi, normalization):
        """ Return an instance of BeamFactorizationCurrentEvaluation, whose 'values' entry
        are dictionaries specifying the counterterm in flavor space, for the value of chsi 
        specified in argument."""

        # Define the NLO QCD integrate initial-state single collinear counterterms kernels
        kernel_gg = { 
            'bulk' :   EpsilonExpansion({
                0 : 1.
            })*normalization,
            'counterterm' :  EpsilonExpansion({
                0 : 1.
            })*normalization,
            'endpoint' :  EpsilonExpansion({
                0 : 1.
            })*normalization
        }
                 
        kernel_gq = { 
            'bulk' :   EpsilonExpansion({
                0 : 1.
            })*normalization,
            'counterterm' :  EpsilonExpansion({
                0 : 1.
            })*normalization,
            'endpoint' :  EpsilonExpansion({
                0 : 1.
            })*normalization
        }
        
        kernel_qg = { 
            'bulk' :   EpsilonExpansion({
                0 : 1.
            })*normalization,
            'counterterm' :  EpsilonExpansion({
                0 : 1.
            })*normalization,
            'endpoint' :  EpsilonExpansion({
                0 : 1.
            })*normalization
        }
        
        kernel_qq = { 
            'bulk' :   EpsilonExpansion({
                0 : 1.
            })*normalization,
            'counterterm' :  EpsilonExpansion({
                0 : 1.
            })*normalization,
            'endpoint' :  EpsilonExpansion({
                0 : 1.
            })*normalization
        }

        active_quark_PDGs = tuple([pdg for pdg in range(1,7)+range(-1,-7,-1) 
                                                                 if pdg in self.beam_PDGs])
        
        # Build the NLO flavor matrix
        flavor_matrix = {}
        for reduced_flavor in self.beam_PDGs:
            # Gluon backward evolution
            if reduced_flavor==21:
                flavor_matrix[21] = { (21,) : kernel_gg[self.distribution_type] }
                if active_quark_PDGs:
                    flavor_matrix[21][active_quark_PDGs] = kernel_gq[self.distribution_type]
                                
            # Quark backward evolution            
            if reduced_flavor in active_quark_PDGs:
                flavor_matrix[reduced_flavor] = { 
                    (21,) : kernel_qg[self.distribution_type],
                    (reduced_flavor,) : kernel_qq[self.distribution_type]
                }

        # Now assign the flavor matrix in the BeamFactorizationCurrentEvaluation instance
        evaluation = utils.BeamFactorizationCurrentEvaluation({
            'spin_correlations'         : [None,],
            'color_correlations'        : [None,],
            'values'                    : { (0,0) : flavor_matrix }
        })

        return evaluation
