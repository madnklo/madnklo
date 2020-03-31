from madgraph.core.base_objects import EpsilonExpansion
import madgraph.various.misc as misc

#=========================================================================================
# Class listing collinear Altarelli-Parisi Kernels for distributed softs
#=========================================================================================

class AltarelliParisiKernels_soft:
    """ Implementation of AltarelliParisiKernels. Notice that the the first argument color_factors is always
    any object with the float attributes CF, CA, etc... to access color factors."""

    # None  : -g^{mu,nu}
    # kT    : k_T^mu k_T^nu

    @staticmethod
    def P_qg_averaged(color_factors,z):
        return EpsilonExpansion({
            0: color_factors.CF*(1-z), # ((1.+z**2)/(1.-z))
            1: color_factors.CF*(-(1-z)) # same
        })

    @staticmethod
    def P_qg(color_factors, z, kT):
        return [
            ( None, AltarelliParisiKernels_soft.P_qg_averaged(color_factors, z) ), #same
        ]

    @staticmethod
    def P_qq(color_factors, z, kT):
        return [
            ( None, EpsilonExpansion({
                    0: color_factors.TR, #same
            })),
            ( (kT,) , EpsilonExpansion({
                0: color_factors.TR * (4.*z*(1.-z)*(1./kT.square())), #(4.*z*(1.-z)*(1./kT.square())) same
            })),
        ]

    @staticmethod
    def P_gg(color_factors, z, kT):

        return [
            # The line below implements the -g_{\mu\nu} part of the splitting kernel.
            # Notice that the extra longitudinal terms included in the spin-correlation 'None'
            # from the relation:
            #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
            #    = g^{\mu\nu} + longitudinal terms
            # are irrelevant because Ward identities evaluate them to zero anyway.
            ( None, EpsilonExpansion({
                0: 0, #2 * color_factors.CA * (z/(1-z) + (1-z)/z)
            })),
            ( (kT,) , EpsilonExpansion({
                0: -4 * color_factors.CA * z * (1-z) / kT.square(), # same
                1: 4 * color_factors.CA * z * (1-z) / kT.square(), # no term
            })),
        ]

    @staticmethod
    def P_gq(color_factors, z, kT):

        return [
            ( None, AltarelliParisiKernels_soft.P_qg_averaged(color_factors, 1.-z) ), #same
        ]



#=========================================================================================
# Class listing soft kernels
#=========================================================================================

class SoftKernels_soft:
    """ Implementation of the universal soft kernels for distributed softs."""

    @staticmethod
    def eikonal_dipole_soft(pi, pj, ps): # in Lionetti's thesis (pi = p_j, pj = p_k, ps = p_i)
        """Eikonal factor for soft particle with momentum ps
        emitted from the dipole with momenta pi and pj.
        Modified for the subtraction scheme distributed softs.
        """
        pipj = pi.dot(pj)
        pips = pi.dot(ps)
        p_ij = pi + pj
        return pipj / pips * (1.0/ps.dot(p_ij))













