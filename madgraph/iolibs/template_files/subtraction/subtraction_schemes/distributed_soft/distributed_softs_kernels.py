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

    # NNLO splitting functions

    @staticmethod
    def P_qqx_sub(color_factors, z, kT):
        # splinting function for the colinear quark antiquark counterterm with the soft colinear counterterm subtracted

        return [
            ( None, EpsilonExpansion({
                    0: color_factors.TR,
            })),
            ( (kT,) , EpsilonExpansion({
                0: color_factors.TR * (4.*z*(1.-z)*(1./kT.square())),
            })),
        ]

    @staticmethod
    def P_qxqq (color_factors, z_1, z_2, z_3,s_12, s_13, s_23, kT_1, kT_2, kT_3):
        # 1,2 same quark flavour :: 2 quark, 1 antiquark
        # 3 different quark flavour
        # subtract soft colinear counterterm
        # for pure colinear counterterm
        s_123 = s_12 + s_13 + s_23

        t_123 = 2.0 * (z_1 * s_23 - z_2 * s_13)/(z_1 + z_2) + (z_1 - z_2)/(z_1 +z_2) * s_12
        return [
            (None, EpsilonExpansion({
                0 : 0.5 * color_factors.CF * color_factors.TR * (s_123/s_12) *
                    (-t_123**2/(s_12 * s_123) + (4.0 * z_3 + (z_1 - z_2)**2)/(z_1 + z_2)+(z_1 + z_2 - s_12/s_123)),
                1 : 0.5 * color_factors.CF * color_factors.TR * (-2.0 * (z_1 + z_2 - s_12/s_123)),
            })),
        ]

    @staticmethod
    def P_qqxq_sub (color_factors, z_1, z_2, z_3,s_12, s_13, s_23, kT_1, kT_2, kT_3):
        # 1,2 same quark flavour :: 1 quark, 2 antiquark
        # 3 different quark flavour
        # subtract soft colinear counterterm
        # subtract colinear colinear counterterm Todo ?
        # subtract colinear soft colinear counterterm Todo?

        s_123 = s_12 + s_13 + s_23
        z_12 = z_1 + z_2
        # print(z_12, 1. - z_3)
        print (s_123/s_12, s_13/s_12, s_23/s_12)
        print ((z_1**2 + z_2**2) / z_12, z_1, z_2, z_3, z_12)
        print ((s_123/s_12) * ((z_1**2 + z_2**2) / z_12) - 1.0)

        return [
            (None, EpsilonExpansion({
                0 : color_factors.CF * color_factors.TR * ((s_123/s_12) * ((z_1**2 + z_2**2) / z_12) - 1.0),
                1 : color_factors.CF * color_factors.TR * (1.0 + s_123/s_12*z_12),
            })),
        ]

        # return [
        #     (None, EpsilonExpansion({
        #         0 : color_factors.CF * color_factors.TR * 0.5 * (s_123/s_12 * (1.0 + (z_1**2 + z_2**2) / z_12) - 1.0),
        #         1 : color_factors.CF * color_factors.TR * (1.0 - s_123/s_12*z_12),
        #     })),
        # ]

    @staticmethod
    def P_qpqp_q (color_factors, z_1, z_2, z_hat_12, z_hat_3, kT_1, kT_2, kT_hat_12, kT_hat_3):
         # for C12 C123 counterterm

        return [(
            None, EpsilonExpansion({
                0 : color_factors.CF * color_factors.TR *
                    ((2.0 * z_hat_3 / z_hat_12) + z_hat_12
                    - 2.0 * z_1 * z_2 * (z_hat_12 + z_hat_3 / z_hat_12 * (2.0 * kT_1.dot(kT_hat_3))**2/(kT_1.square() * kT_hat_3.square()))),
                1 : - color_factors.CF * color_factors.TR * z_hat_12,
            })
        )]

    @staticmethod
    def P_qpqp_q_sub(color_factors, z_1, z_2, z_hat_12, z_hat_3, kT_1, kT_2, kT_hat_12, kT_hat_3):
        # for C12 C123 counterterm - S12 C12 C123 countertem

        return [(
            None, EpsilonExpansion({
                0 : color_factors.CF * color_factors.TR * (z_hat_12 - 2.0 * z_1 * z_2 * z_hat_12),
                1 : - color_factors.CF * color_factors.TR * z_hat_12,
            })
        )]

    @staticmethod
    def P_S_qpqp_q (color_factors, z_1, z_2, z_hat_12, z_hat_3, kT_1, kT_2, kT_hat_12, kT_hat_3):
         # for S 12 C12 C123 counterterm

        return [(
            None, EpsilonExpansion({
                0 : color_factors.CF * color_factors.TR *
                    (2.0 * z_hat_3 / z_hat_12) * (
                            1.0 - z_1 * z_2 * (2.0 * kT_1.dot(kT_hat_3))**2/(kT_1.square() * kT_hat_3.square())),
            })
        )]



#=========================================================================================
# Class listing soft kernels
#=========================================================================================

class SoftKernels_soft:
    """ Implementation of the universal soft kernels for distributed softs."""

    @staticmethod
    def eikonal_dipole_soft(pi, pj, ps): # in Lionetti's thesis (pi = p_j, pj = p_k, ps = p_i)
        """Eikonal factor for soft particle with momentum ps
        emitted from the dipole with momenta pi and pj.
        Modified for the subtraction scheme distributed soft.
        """
        pipj = pi.dot(pj)
        pips = pi.dot(ps)
        p_ij = pi + pj
        return (pipj / pips) * (1.0/ps.dot(p_ij))

    @staticmethod
    def eikonal_dipole_soft_mod(pi, pj, ps, s_tilde): # in Lionetti's thesis (pi = p_j, pj = p_k, ps = p_i)
        """Eikonal factor for soft particle with momentum ps
        emitted from the dipole with momenta pi and pj.
        Modified for the subtraction scheme distributed soft.
        With term (s_tilde/s_ijk) for easier integration
        """
        pipj = pi.dot(pj)
        s_ijk = 2.0 * pj.dot(pi + ps)
        pips = pi.dot(ps)
        return (pipj / pips) * (1.0/(ps.dot(pi)+ps.dot(pj) * (s_tilde/s_ijk))) * (s_tilde/s_ijk)

    @staticmethod
    def eikonal_qqx_soft_colinear(color_factors, z_1, z_2, s_12, kT_1, kT_2, p_hat_12, p_hat_i, p_hat_j):
        """distributed soft eikonal factor"""

        # p_1, p_2 : q, qx :: symmetric in exchange of p_1 and p_2
        # p_hat_i : other colinear momentum mapped for colinear q qx
        # p_hat_j is other momentum mapped

        s_hat_ij = 2.0 * p_hat_i.dot(p_hat_j)
        s_hat_12i = 2.0 * p_hat_12.dot(p_hat_i)
        s_hat_12j = 2.0 * p_hat_12.dot(p_hat_j)

        # print (s_hat_12i)
        # print (s_hat_12j)
        # print (s_hat_12j/s_hat_12i)
        # print (kT_1.dot(p_hat_i))
        # print (kT_1.dot(p_hat_j))

        # print (1.0 - 0.5 * s_hat_12j/s_hat_12i * (kT_1.dot(p_hat_i))/(kT_1.dot(p_hat_j)) - 0.5 * s_hat_12i/s_hat_12j * (kT_1.dot(p_hat_j))/(kT_1.dot(p_hat_i)))

        return ((color_factors.TR / s_12) * (1.0/(s_hat_12i + s_hat_12j)) *
                (
                    - s_hat_ij/(s_hat_12i) + 2.0 * z_1 * z_2 * (((2.0 * kT_1.dot(p_hat_i)) * (2.0 * kT_1.dot(p_hat_j)))/(s_hat_12i * kT_1.square())) *
                    (1.0 - 0.5 * s_hat_12j/s_hat_12i * (kT_1.dot(p_hat_i))/(kT_1.dot(p_hat_j)) - 0.5 * s_hat_12i/s_hat_12j * (kT_1.dot(p_hat_j))/(kT_1.dot(p_hat_i)))
                ))


        # return ((color_factors.TR / s_12) * # not distributed among the colinears
        #         (
        #             - s_hat_ij/(s_hat_12i * s_hat_12j) + 2.0 * z_1 * z_2 * (((2.0 * kT_1.dot(p_hat_i)) * (2.0 * kT_1.dot(p_hat_j)))/(s_hat_12i * s_hat_12j * kT_1.square())) *
        #             (1.0 - 0.5 * s_hat_12j/s_hat_12i * (kT_1.dot(p_hat_i))/(kT_1.dot(p_hat_j)) - 0.5 * s_hat_12i/s_hat_12j * (kT_1.dot(p_hat_j))/(kT_1.dot(p_hat_i)))
        #         ))

    @staticmethod
    def eikonal_qqx(color_factors,p_1, p_2, p_i, p_j):
        """distributed soft eikonal factor"""

        # p_1, p_2 : q, qx :: symmetric in exchange of p_1 and p_2
        # p_i : other colinear momentum
        # p_j is other momentum

        s_12 = 2.0 * p_1.dot(p_2)
        s_1i = 2.0 * p_1.dot(p_i)
        s_1j = 2.0 * p_1.dot(p_j)
        s_2i = 2.0 * p_2.dot(p_i)
        s_2j = 2.0 * p_2.dot(p_j)
        s_ij = 2.0 * p_i.dot(p_j)
        p_12 = p_1 + p_2
        s_12i = 2.0 * p_12.dot(p_i)
        s_12j = 2.0 * p_12.dot(p_j)

        return ((color_factors.TR / s_12**2) * (s_12j/(s_12i + s_12j)) *
                ((s_1i * s_2j + s_1j * s_2i - s_12 * s_ij)/(s_12i * s_12j) - s_1i * s_2i / s_12i**2 - s_1j * s_2j / s_12j**2))















