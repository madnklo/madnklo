from madgraph.core.base_objects import EpsilonExpansion

#=========================================================================================
# Class listing Altarelli-Parisi Kernels
#=========================================================================================
class AltarelliParisiKernels:
    """ Implementation of AltarelliParisiKernels. Notice that the the first argument color_factors is always
    any object with the float attributes CF, CA, etc... to access color factors."""

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