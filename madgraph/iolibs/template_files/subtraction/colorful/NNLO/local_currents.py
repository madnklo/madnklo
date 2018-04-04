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
"""Implementation of NNLO type of currents."""

import os
import math

import madgraph.core.subtraction as subtraction
import madgraph.integrator.mappings as mappings
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

#=========================================================================================
# Auxiliary functions
#=========================================================================================

def kTdiff(i, j, zs, kTs):
    """Compute the difference vector (kTi/zi - kTj/zj)."""

    return kTs[i-1]/zs[i-1] - kTs[j-1]/zs[j-1]

def sij(i, j, zs, kTs, kTdiffij=None):
    """Compute the invariant mass sij from zs and kTs."""

    if kTdiffij: kTdiffij2 = kTdiffij.square()
    else: kTdiffij2 = kTdiff(i, j, zs, kTs).square()
    return -zs[i-1]*zs[j-1]*kTdiffij2

def tijk(i, j, k, zs, kTs, s_ij=None, s_ik=None, s_jk=None):
    """Compute the invariant t_{ij,k} from zs and kTs."""

    if not s_ij: s_ij = sij(i, j, zs, kTs)
    if not s_ik: s_ik = sij(i, k, zs, kTs)
    if not s_jk: s_jk = sij(j, k, zs, kTs)
    return (2*(zs[i-1]*s_jk-zs[j-1]*s_ik)+(zs[i-1]-zs[j-1])*s_ij)/(zs[i-1]+zs[j-1])

#=========================================================================================
# NNLO final-collinear currents
#=========================================================================================

class QCD_final_collinear_0_QQxq(currents.QCDLocalCollinearCurrent):
    """Q Q~ q collinear tree-level current."""

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NNLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 4, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that there are 3 massless final state quarks or antiquarks
        if len(ss.legs) != 3: return None
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
            if not cls.is_quark(leg, model): return None
        # Look for a quark/antiquark pair
        pair = None
        for i in range(len(ss.legs)):
            for j in range(i+1, len(ss.legs)):
                if cls.are_antiparticles(ss.legs[i], ss.legs[j]):
                    pair = (ss.legs[i], ss.legs[j])
                    continue
            if pair is not None: continue
        if pair is None: return None
        # Identify the remaining quark
        other_quarks = [leg for leg in ss.legs if leg not in pair]
        # Since leg numbers have been discarded, equal legs will not appear here
        # Thus if the quark species were the same, other_quarks = []
        if len(other_quarks) != 1: return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        if cls.are_antiparticles(legs[0], legs[1]):
            return (legs[0].n, legs[1].n, legs[2].n)
        elif cls.are_antiparticles(legs[0], legs[2]):
            return (legs[0].n, legs[2].n, legs[1].n)
        else:
            return (legs[1].n, legs[2].n, legs[0].n)

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables and compute basic quantities
        z1, z2, z3 = zs
        s12 = sij(1, 2, zs, kTs)
        s13 = sij(1, 3, zs, kTs)
        s23 = sij(2, 3, zs, kTs)
        s123 = s12 + s13 + s23
        t123 = tijk(1, 2, 3, zs, kTs, s_ij=s12, s_ik=s13, s_jk=s23)
        # Assemble kernel
        sqrbrk  = -(t123 ** 2)/(s12*s123)
        sqrbrk += (4*z3 + (z1-z2)**2) / (z1+z2)
        sqrbrk += z1 + z2 - s12/s123
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None}}
        })
        evaluation['values'][(0, 0)]['finite'] = self.CF*self.TR * s123 / (2*s12) * sqrbrk
        return evaluation
    
class QCD_final_soft_0_gg(currents.QCDLocalSoftCurrent):
    """Double soft tree-level current."""
    
    is_cut = staticmethod(currents.SomogyiChoices.cut_soft)
    factor = staticmethod(currents.SomogyiChoices.factor_soft)
    
    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NNLO QCD soft tree-level current
        init_vars = cls.common_does_implement_this_current(current, QCD_squared_order=4, n_loops=0)
        if init_vars is None: return None
        
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that there are 2 massless gluons
        if len(ss.legs) != 2: return None
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
            if not cls.is_gluon(leg, model): return None

        # The current is valid
        return init_vars
   
    @staticmethod
    def eikonal(PS_point, i, j, r):
        """Eikonal factor for soft particle with number 'r'
        emitted from 'i' and reconnecting to 'j'.
        """

        pipj = PS_point[i].dot(PS_point[j])
        pipr = PS_point[i].dot(PS_point[r])
        pjpr = PS_point[j].dot(PS_point[r])
        return pipj/(pipr*pjpr)

    @staticmethod
    def strongly_ordered_double_soft_kernel(PS_point, i, j, r1, r2):
        """Non abelian piece of the double soft current with first soft particle with 
        number 'r1' and second with number 'r2' with emitting index 'i'  and reconnecting 
        index 'j'. See Eq.(109-111) of hep-ph/9908523v1 (Catani-Grazzini)
        """

        pipj = PS_point[i].dot(PS_point[j])
        pr1pr2 = PS_point[r1].dot(PS_point[r2])
        pipr1 = PS_point[i].dot(PS_point[r1])
        pipr2 = PS_point[i].dot(PS_point[r2])
        pjpr1 = PS_point[j].dot(PS_point[r1])
        pjpr2 = PS_point[j].dot(PS_point[r2])
        
        return ( 
            (pipj / pr1pr2) * ( 1. / (pipr1*pjpr2) + 1. / (pjpr1*pipr2) ) - 
                                                    (pipj)**2 / (pipr1*pjpr2*pipr2*pjpr1) 
        )

    @staticmethod
    def non_abelian_eikonal(PS_point, i, j, r1, r2):
        """Non abelian piece of the double soft current with first soft particle with 
        number 'r1' and second with number 'r2' with emitting index 'i'  and reconnecting 
        index 'j'. See Eq.(109-111) of hep-ph/9908523v1 (Catani-Grazzini)
        """

        pipj = PS_point[i].dot(PS_point[j])
        pr1pr2 = PS_point[r1].dot(PS_point[r2])
        pipr1 = PS_point[i].dot(PS_point[r1])
        pipr2 = PS_point[i].dot(PS_point[r2])
        pjpr1 = PS_point[j].dot(PS_point[r1])
        pjpr2 = PS_point[j].dot(PS_point[r2])
        
        fully_expanded_expression = False 
        if not fully_expanded_expression:
            so_double_soft = QCD_final_soft_0_gg.strongly_ordered_double_soft_kernel(
                                                                  PS_point, i, j, r1, r2 )
            return (
                so_double_soft + ((pipr1*pjpr2+pipr2*pjpr1)/((pipr1+pipr2)*(pjpr1+pjpr2)))*(
                    (1./(pr1pr2**2)) - 0.5*so_double_soft
                ) - ((2.*pipj)/(pr1pr2*(pipr1+pipr2)*(pjpr1+pjpr2)))
            )

        else:
            return (
                 (1.0/(pr1pr2**2)) * ( (pipr1*pjpr2 + pipr2*pjpr1) / ( (pipr1+pipr2)*(pjpr1+pjpr2) ) )
               - ( ((pipj)**2) / (2.*pipr1*pjpr2*pipr2*pjpr1) )*( 2. - (pipr1*pjpr2+pipr2*pjpr1)/( (pipr1+pipr2)*(pjpr1+pjpr2) ) )
               + ( ( pipj / (2.*pr1pr2) ) * (
                   2./(pipr1*pjpr2) + 2./(pjpr1*pipr2) - (1./((pipr1+pipr2)*(pjpr1+pjpr2)))*(4.+ ((pipr1*pjpr2+pipr2*pjpr1)**2)/(pipr1*pjpr2*pipr2*pjpr1))
                   )
                 )
            )        
   
    def create_CataniGrazzinni_correlator(self, (i,j),(k,l)):
        """ Returns the correlator of Catani-Grazzini (Eq.113 of hep-ph/9908523v1)
                <M| ( T^-1_i \dot T^-1_j ) * ( T^-1_k \dot T^-1_l ) | M > 
                
            converted into MadGraph's conventions:
            
              ( (a,-1,a),(b,-2,b) ) , ( (c,-1,c),(d,-2,d) ) --> T^-1_a T^-2_b T^-1_c T^-2_d
            """

        # It is important to never commute two color operators acting on the same index, so we must chose carefully which
        # index to pick to carry the gluon index '-2' of the first connection. This can be either 'k' or 'l'.
        if j!=k and j!=l:
            # If all indices are different, we can pick either k or l, it is irrelevant
            index1, index2, index3, index4 = i, k, j, l
        elif j==k and j!=l:
            # If j is equal to k, we must pick l
            index1, index2, index3, index4 = i, l, j, k
        elif j==l and j!=k:
            # If j is equal to l, we must pick k
            index1, index2, index3, index4 = i, k, j, l
        elif j==l and j==k:
            # If j is equal to both l and k, then agin it doesn't matter and we can pick k
            index1, index2, index3, index4 = i, k, j, l

        # The sorting according to the first index of each tuple of each of the two convention is to match
        # Madgraph's convention for sorting color connection in the color correlators definition
        return (
            tuple(sorted([ (index1,-1,index1), (index2,-2,index2) ], key = lambda el: el[0], reverse=True)),
            tuple(sorted([ (index3,-1,index3), (index4,-2,index4) ], key = lambda el: el[0], reverse=True))
        )

    def evaluate_subtraction_current(
        self, current, PS_point,
        reduced_process=None, hel_config=None,
        mapping_variables=None, leg_numbers_map=None ):

        if not hel_config is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s" % self.__class__.__name__ +
                " does not support helicity assignment.")
        if reduced_process is None:
            raise CurrentImplementationError(
                "Subtraction current implementation %s" % self.__class__.__name__ +
                " requires a reduced_process.")
        
        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r    = model_param_dict['MU_R']

        # Now find all colored leg numbers in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color')==1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        # Identify the soft leg numbers
        soft_leg_number_A = current.get('singular_structure').legs[0].n
        soft_leg_number_B = current.get('singular_structure').legs[1].n

        # Include the counterterm only in a part of the phase space
        if self.is_cut(mapping_variables, None):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Now instantiate what the result will be
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'   : [ None ],
            'color_correlations'  : [],
            'values'              : {}
        })
        
        # Normalization factors
        couplings_factors = 4. * math.pi * alpha_s
        norm = (couplings_factors)**2*self.factor(mapping_variables, None)

        # Keep track of the color correlators added
        color_correlators_added = {}
        color_correlation_max_index = 0
        # Now loop over the colored parton numbers to implement the squared double-soft
        # current flow. Note that significant improvement can be obtained by taking advantage
        # of the symmetries, as well as writing fully expanded hard-coded expressions with
        # all dot-products cached.
        # For efficiency, we choose here to define coefficients of each single product of
        # correlators. So, typically in the abelian piece, the color correlator reads:
        #       ( '{}' denotes an anti-commutator below and '.' denotes a dot-product )
        #    { T_i . T_j, T_k . T_l }
        # And we will register *two* entries:
        #    evaluation['color_correlations'].append(
        #      ( 
        #         ( ( (i,-1,i), (k,-2,k) ), ( (j,-1,j), (l,-2,l) ) ),
        #      )
        #    )
        # And:
        #    evaluation['color_correlations'].append(
        #      (     
        #         ( ( (i,-2,i), (k,-1,k) ), ( (j,-2,j), (l,-1,l) ) ),
        #      )
        #    )
        # As opposed to directly defining their sum:
        #    evaluation['color_correlations'].append(
        #            ( 
        #                    ( ( (i,-1,i), (k,-2,k) ), ( (j,-1,j), (l,-2,l) ) ),
        #                    ( ( (i,-2,i), (k,-1,k) ), ( (j,-2,j), (l,-1,l) ) )
        #            )
        #    )
        for i in all_colored_parton_numbers:
            for j in all_colored_parton_numbers:
                non_abelian_kernel = - self.CA * norm * self.non_abelian_eikonal(
                                      PS_point, i, j, soft_leg_number_A, soft_leg_number_B)
                
                # Implement the non-abelian piece
                non_abelian_correlator = ( 
                    ( ( (i,-1,i), ), ( (j,-1,j), ) ), 
                )
                if non_abelian_correlator in color_correlators_added:
                    color_correlation_index = color_correlators_added[non_abelian_correlator]
                    evaluation['values'][(0, color_correlation_index)]['finite'] += non_abelian_kernel
                else:
                    evaluation['color_correlations'].append(non_abelian_correlator)
                    color_correlation_index = color_correlation_max_index
                    color_correlators_added[non_abelian_correlator] = color_correlation_max_index
                    color_correlation_max_index += 1
                    evaluation['values'][(0, color_correlation_index)] = {
                                                             'finite': non_abelian_kernel }
    
                for k in all_colored_parton_numbers:
                    for l in all_colored_parton_numbers:
                    
                        # Implement the abelian piece
                        abelian_kernel = 0.5 * norm * \
                                        self.eikonal(PS_point, i, j, soft_leg_number_A) * \
                                        self.eikonal(PS_point, k, l, soft_leg_number_B)
                        
                        abelian_correlator_A = ( self.create_CataniGrazzinni_correlator((i,j),(k,l)), )
                        abelian_correlator_B = ( self.create_CataniGrazzinni_correlator((k,l),(i,j)), )

                        for correlator in [abelian_correlator_A, abelian_correlator_B]:
                            if correlator in color_correlators_added:
                                color_correlation_index = color_correlators_added[correlator]
                                #misc.sprint('Adding %f ((%d,%d,%d)->%f, (%d,%d,%d)->%f) to CC: %d, %s'%\
                                #            (abelian_kernel,
                                #             i,j,soft_leg_number_A,self.eikonal(PS_point, i, j, soft_leg_number_A),
                                #             k,l,soft_leg_number_B,self.eikonal(PS_point, k, l, soft_leg_number_B),
                                #             color_correlation_index, str(correlator)
                                #            ))
                                evaluation['values'][(0, color_correlation_index)]['finite'] += abelian_kernel
                            else:
                                evaluation['color_correlations'].append(correlator)
                                color_correlation_index = color_correlation_max_index
                                color_correlators_added[correlator] = color_correlation_max_index
                                color_correlation_max_index += 1
                                evaluation['values'][(0, color_correlation_index)] = {
                                                                 'finite': abelian_kernel }
                                #misc.sprint('Adding %f ((%d,%d,%d)->%f, (%d,%d,%d)->%f) to CC: %d, %s'%\
                                #            (abelian_kernel,
                                #             i,j,soft_leg_number_A,self.eikonal(PS_point, i, j, soft_leg_number_A),
                                #             k,l,soft_leg_number_B,self.eikonal(PS_point, k, l, soft_leg_number_B),
                                #             color_correlation_index, str(correlator)
                                #            ))
        
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())) )
        
        #misc.sprint('==BELOW LIST CC==')
        #for i, cc in enumerate(evaluation['color_correlations']):
        #    misc.sprint('Color correlator: %d : %s = %f'%(i, cc, evaluation['values'][(0,i)]['finite']))
        #misc.sprint('==ABOVE LIST CC==')

        return result



    
