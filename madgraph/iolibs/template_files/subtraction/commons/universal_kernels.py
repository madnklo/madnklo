from madgraph.core.base_objects import EpsilonExpansion
import madgraph.various.misc as misc

#=========================================================================================
# Class listing collinear Altarelli-Parisi Kernels
#=========================================================================================

class AltarelliParisiKernels:
    """ Implementation of AltarelliParisiKernels. Notice that the the first argument color_factors is always
    any object with the float attributes CF, CA, etc... to access color factors."""

    @staticmethod
    def P_qg_averaged(color_factors,z):
        return EpsilonExpansion({
            0: color_factors.CF*((1.+z**2)/(1.-z)),
            1: color_factors.CF*(-(1-z))
        })

    @staticmethod
    def P_qq(color_factors, z, kT):
        return [
            ( None, EpsilonExpansion({
                    0: color_factors.TR,
            })),
            ( kT, EpsilonExpansion({
                0: color_factors.TR * (4.*z*(1.-z)*(1./kT.square())),
            })),
        ]

    @staticmethod
    def P_qqpqp(color_factors, z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s):
        """ Kernel for the q -> q qp' qp' splitting. The return value is not a float but a list of tuples:
                ( spin_correlation_vectors_with_parent, weight )
            where spin_correlation_vector_with_parent can be None if None is required.
        """

        # Compute handy term and variables
        s_irs = s_ir+s_is+s_rs
        t_rs_i = 2.*(z_r*s_is-z_s*s_ir)/(z_r+z_s) + ((z_r-z_s)/(z_r+z_s))*s_rs
        dimensional_term = z_r + z_s - s_rs/s_irs

        # Overall prefactor
        prefactor = (1./2.)*color_factors.CF*color_factors.TR*(s_irs/s_rs)

        return [
            ( None,
                EpsilonExpansion({
                    0 : -(t_rs_i**2/(s_rs*s_irs)) + (4.*z_i + (z_r - z_s)**2)/(z_r+z_s) + dimensional_term,
                    1 : -2.*dimensional_term
                })*prefactor
            )
        ]

#TZaddition
    @staticmethod
    def P_qqq(color_factors, z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s):
        """ Kernel for the q -> q q q splitting. The return value is not a float but a list of tuples:
                ( spin_correlation_vectors_with_parent, weight )
            where spin_correlation_vector_with_parent can be None if None is required.
        """

        # Compute handy term and variables
        s_irs = s_ir+s_is+s_rs

        fac1 = s_irs/s_rs
        fac2 = s_irs/s_ir
        fac3 = s_irs**2/s_rs/s_ir

        P1id = EpsilonExpansion({
                    0 : 2.*s_is/s_rs + fac1*( (1.+z_r**2)/(1.-z_s) - 2.*z_s/(1.-z_i) ) - fac3*z_r/2.*( (1.+z_r**2)/(1.-z_s)/(1.-z_i) ),
                    1 : -1. - 2.*s_is/s_rs - fac1*( (1.-z_i)**2/(1.-z_s) + 1. + z_r - 2.*z_s/(1.-z_i) ) + fac3*z_r/2.*(1.+2.*(1.-z_s)/(1.-z_i)),
                    2 : 1. - fac1*(1.-z_i) + fac3*z_r/2.
                })

        P2id = EpsilonExpansion({
                    0 : 2.*s_is/s_ir + fac2*( (1.+z_r**2)/(1.-z_i) - 2.*z_i/(1.-z_s) ) - fac3*z_r/2.*( (1.+z_r**2)/(1.-z_i)/(1.-z_s) ),
                    1 : -1. - 2.*s_is/s_ir - fac2*( (1.-z_s)**2/(1.-z_i) + 1. + z_r - 2.*z_i/(1.-z_s) ) + fac3*z_r/2.*(1.+2.*(1.-z_i)/(1.-z_s)),
                    2 : 1. - fac2*(1.-z_s) + fac3*z_r/2.
                })


        # Prefactor for the id term
        prefactor_id = color_factors.CF*(color_factors.CF - (1./2.)*color_factors.CA)

        pqqpqp1=AltarelliParisiKernels.P_qqpqp(color_factors, z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s)[0][1]
        pqqpqp2=AltarelliParisiKernels.P_qqpqp(color_factors, z_s, z_r, z_i, s_rs, s_is, s_ir, kT_s, kT_r, kT_i)[0][1]

        return [
            ( None,
                pqqpqp1 + pqqpqp2 +(P1id + P2id)*prefactor_id
            )
        ]

    #TZaddition
    @staticmethod
    def P_qgg_ab(color_factors, z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s):
        """ Abelian part of the kernel for the q -> q g g splitting. The return value is not a float but a list of tuples:
                ( spin_correlation_vectors_with_parent, weight )
            where spin_correlation_vector_with_parent can be None if None is required.
        """

        # Compute handy term and variables
        s_irs = s_ir+s_is+s_rs
        fac_1 = s_irs**2/2./s_ir/s_is*z_i
        fac_2 = s_irs/s_ir
        fac_3 = s_is/s_ir

        return [
            ( None,
                EpsilonExpansion({
                    0 : fac_1*(1.+z_i**2)/z_r/z_s
                        + fac_2*(z_i*(1.-z_r)+(1.-z_s)**3)/z_r/z_s
                        - fac_3,
                    1 : - fac_1*( (z_r**2+z_s**2)/z_r/z_s + 1. )
                        - fac_2*(z_r**2+z_r*z_s+z_s**2)*(1.-z_s)
                                /z_r/z_s
                        + 1. + 2.*fac_3,
                    2 : - fac_1 + fac_2*(1.+z_i) - 1. - fac_3
                })
            )
        ]

    #TZaddition
    @staticmethod
    def P_qgg_nab(color_factors, z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s):
        """ Non-abelian part of the kernel for the q -> q g g splitting. The return value is not a float but a list of tuples:
                ( spin_correlation_vectors_with_parent, weight )
            where spin_correlation_vector_with_parent can be None if None is required.
        """

        # Compute handy term and variables
        s_irs = s_ir+s_is+s_rs
        t_rs_i = 2.*(z_r*s_is-z_s*s_ir)/(z_r+z_s) + ((z_r-z_s)/(z_r+z_s))*s_rs
        fac_1 = s_irs**2/2./s_rs/s_ir
        fac_2 = s_irs**2/4./s_ir/s_is*z_i
        fac_3 = s_irs/2./s_rs
        fac_4 = s_irs/2./s_ir

        long_term_1 = ( z_r*(2.-2.*z_r+z_r**2) - z_s*(6.-6.*z_s+z_s**2) )/z_s/(1.-z_i)
        long_term_2 = ((1.-z_s)**3+z_i**2-z_s)/z_s/(1.-z_i)


        return [
            ( None,
                EpsilonExpansion({
                    0 : t_rs_i**2/4./s_rs**2 + 1./4.
                        + fac_1*( ((1.-z_i)**2 + 2.*z_i)/z_s 
                                + (z_s**2 + 2.*(1.-z_s))/(1.-z_i) )
                        - fac_2*( ((1.-z_i)**2 + 2.*z_i)/z_r/z_s  )
                        + fac_3*long_term_1
                        + fac_4*( long_term_2 - (z_i*(1.-z_r)+(1.-z_s)**3)/z_r/z_s ),
                    1 : - t_rs_i**2/4./s_rs**2 - 1./4. - 1./2.
                        + fac_1*( -(1.-z_i)**2/z_s - z_s**2/(1.-z_i) )
                        - fac_2*( -(1.-z_i)**2/z_r/z_s + 1. )
                        + fac_3*( -long_term_1 + 2.*(z_i*(z_r-2.*z_s)-z_s)/z_s/(1.-z_i) )
                        + fac_4*( -long_term_2 - 2.*(1.-z_s)*(z_s-z_i)/z_s/(1.-z_i) + z_r - z_s + (1.-z_s)*(z_r**2+z_s**2)/z_r/z_s ),
                    2 : 1./2. + fac_2 - fac_4*(1.-z_s)
                })
            )
        ]

    #TZaddition
    @staticmethod
    def P_qgg(color_factors, z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s):
        """ Kernel for the q -> q g g splitting. The return value is not a float but a list of tuples:
                ( spin_correlation_vectors_with_parent, weight )
            where spin_correlation_vector_with_parent can be None if None is required.
        """

        # Compute handy term and variables
        prefactor_ab  = color_factors.CF**2
        prefactor_nab = color_factors.CF * color_factors.CA

        pqgg_ab_1 = AltarelliParisiKernels.P_qgg_ab(color_factors, z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s)[0][1]
        pqgg_ab_2 = AltarelliParisiKernels.P_qgg_ab(color_factors, z_i, z_s, z_r, s_is, s_ir, s_rs, kT_i, kT_s, kT_r)[0][1]

        pqgg_nab_1 = AltarelliParisiKernels.P_qgg_nab(color_factors, z_i, z_r, z_s, s_ir, s_is, s_rs, kT_i, kT_r, kT_s)[0][1]
        pqgg_nab_2 = AltarelliParisiKernels.P_qgg_nab(color_factors, z_i, z_s, z_r, s_is, s_ir, s_rs, kT_i, kT_s, kT_r)[0][1]

        return [
            ( None,
                ( pqgg_ab_1 + pqgg_ab_2 )*prefactor_ab
                + ( pqgg_nab_1 + pqgg_nab_2 )*prefactor_nab
            )
        ]


    @staticmethod
    def P_q_qpqp(color_factors, z_rs, z_i_rs, kT_rs, pi_hat, p_rs_hat):
        """ Kernel for the q -> q (qp' qp') strongly ordered splitting. The return value is not a float but a list of tuples:
                ( spin_correlation_vectors_with_parent, weight )
            where spin_correlation_vector_with_parent can be None if None is required.
        """

        # Overall prefactor
        result = ( AltarelliParisiKernels.P_qg_averaged(color_factors,z_i_rs) +
                   EpsilonExpansion({0:
                        -2.*color_factors.CF*z_rs*(1.-z_rs)*(1-z_i_rs-((2.*kT_rs.dot(pi_hat))**2)/(kT_rs.square()*(2.*pi_hat.dot(p_rs_hat))))
                    }))*color_factors.TR

        return [ ( None, result ) ]


#=========================================================================================
# Class listing soft kernels
#=========================================================================================

class SoftKernels:
    """ Implementation of the universal soft kernels."""

    @staticmethod
    def eikonal_g(color_factors, pi, pk, pr, spin_corr_vector=None):
        """ Gluon eikonal, with p_i^\mu p_i^\nu in the numerator, dotted into each other if spin_corr_vector=None
        otherwise dotted with the spin_corr_vector."""

        s_ir = pi.dot(pr)
        s_kr = pk.dot(pr)
        numerator = pi.dot(pk) if spin_corr_vector is None else pi.dot(spin_corr_vector)*pk.dot(spin_corr_vector)

        return EpsilonExpansion({'finite': numerator/(s_ir*s_kr)})

    @staticmethod
    def eikonal_qqx(color_factors, pi, pk, pr, ps):
        """ Taken from Gabor's notes on colorful ISR@NNLO."""

        s_ir = pi.dot(pr)
        s_ks = pk.dot(ps)
        s_is = pi.dot(ps)
        s_kr = pk.dot(pr)
        s_ik = pi.dot(pk)
        s_rs = pr.dot(ps)
        s_i_rs = pi.dot(pr+ps)
        s_k_rs = pk.dot(pr+ps)

        return EpsilonExpansion({'finite':
           color_factors.TR*(((s_ir*s_ks + s_is*s_kr - s_ik*s_rs)/(s_i_rs*s_k_rs)) - ((s_ir*s_is)/(s_i_rs**2)) - ((s_kr*s_ks)/(s_k_rs**2)))
        })

    @staticmethod
    def qqx(color_factors, colored_partons_momenta, soft_momenta, all_colored_parton_numbers):
        """ Computes the soft current for a quark-antiquark pair."""

        # A list of 2-tuple of the form (color_correlator, epsilon_expansion_weight)
        result = []

        # Now loop over the colored parton number pairs (a,b)
        # and add the corresponding contributions to this current
        for i, a in enumerate(all_colored_parton_numbers):
            # Use the symmetry of the color correlation and soft current (a,b) <-> (b,a)
            for b in all_colored_parton_numbers[i:]:
                # Write the eikonal for that pair
                if a != b:
                    mult_factor = 2.
                else:
                    mult_factor = 1.
                    # For massless partons we can simply skip this entry
                    continue
                pa, pb = colored_partons_momenta[a], colored_partons_momenta[b]
                eikonal = SoftKernels.eikonal_qqx(color_factors, pa, pb, soft_momenta[0], soft_momenta[1])
                result.append( (
                        ((a, b), ),
                        eikonal*mult_factor
                    )
                )

        return result

    #TZaddition
    @staticmethod
    def eikonal_2g(color_factors, pi, pk, pr, ps):
        """two-gluon eikonal"""

        strongly_ordered = SoftKernels.eikonal_g(color_factors, pi, pk, ps, spin_corr_vector=None) * (
                    SoftKernels.eikonal_g(color_factors, pi, ps, pr, spin_corr_vector=None) +
                    SoftKernels.eikonal_g(color_factors, pk, ps, pr, spin_corr_vector=None) -
                    SoftKernels.eikonal_g(color_factors, pi, pk, pr, spin_corr_vector=None)
                )

        S_ik_rs = SoftKernels.eikonal_g(color_factors, pi, pk, pr + ps, spin_corr_vector=None)

        s_ir = 2.*pi.dot(pr)
        s_is = 2.*pi.dot(ps)
        s_kr = 2.*pk.dot(pr)
        s_ks = 2.*pk.dot(ps)
        s_rs = 2.*pr.dot(ps)
        s_i_rs = s_ir + s_is
        s_k_rs = s_kr + s_ks

        psrsm2 = EpsilonExpansion({ 0 : 1./s_rs**2, 1 : -1./s_rs**2 })

        return (strongly_ordered + (psrsm2 - strongly_ordered/8.)*4.*(s_ir*s_ks+s_is*s_kr)/s_i_rs/s_k_rs - S_ik_rs*4./s_rs)#*color_factors.CA/4.
