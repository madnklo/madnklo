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

    @staticmethod
    def P_q_qpqp(color_factors, z_rs, z_i_rs, kT_rs, s_i_rs, p_i):
        """ Kernel for the q -> q (qp' qp') strongly ordered splitting. The return value is not a float but a list of tuples:
                ( spin_correlation_vectors_with_parent, weight )
            where spin_correlation_vector_with_parent can be None if None is required.
        """

        # Overall prefactor
        result = ( AltarelliParisiKernels.P_qg_averaged(color_factors,z_i_rs) +
                   EpsilonExpansion({0:
                        -2.*color_factors.CF*z_rs*(1.-z_rs)*(1-z_i_rs-((2.*kT_rs.dot(p_i))**2)/(kT_rs.square()*s_i_rs))
                    }))*color_factors.TR

        return [ ( None, result ) ]

#=========================================================================================
# Class listing soft kernels
#=========================================================================================

class SoftKernels:
    """ Implementation of the universal soft kernels."""

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
