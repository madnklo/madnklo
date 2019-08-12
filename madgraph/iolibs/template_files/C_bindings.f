C This file has been automatically generated for %(output_name)s
C It sets up all C bindings for interfacing to the fortran Matrix Element

C Additional interface subroutines to add:

C
C get_max_n_spin_corr_legs
C get_max_n_spin_corr_vectors
C get_n_color_correlators
C get_color_correlator_for_id
C get_color_connection_max_order
C get_n_color_connections_(nn...)lo
C index_to_color_connection_(nn..)lo
C get_nsqso_born
C get_split_order_names
C get_squared_orders_for_soindex
C

      SUBROUTINE %(binding_prefix)s%(proc_prefix)sME_ACCESSOR_HOOK(P,HEL,USER_ALPHAS,ANS) bind(c)
        USE iso_c_binding, ONLY: c_int, c_double
        INCLUDE 'nexternal.inc'
        INCLUDE 'nsqso_born.inc'
        INTEGER I
        real(c_double), intent(in)  :: P(0:3,NEXTERNAL)
        integer(c_int), intent(in)  :: HEL
        real(c_double), intent(in)  :: USER_ALPHAS
        real(c_double), intent(out) :: ANS(0:NSQSO_BORN)

        CALL %(proc_prefix)sME_ACCESSOR_HOOK(P,HEL,USER_ALPHAS,ANS)
      ENDSUBROUTINE %(binding_prefix)s%(proc_prefix)sME_ACCESSOR_HOOK

      SUBROUTINE %(binding_prefix)s%(proc_prefix)sINITIALISE(PATH) bind(c)
        USE iso_c_binding
        CHARACTER(c_char) :: PATH(512)
        CHARACTER(512) :: PATH_IN
        DO i=1,512
          PATH_IN(i:i) = PATH(i)
        ENDDO
        CALL %(proc_prefix)sINITIALISE(PATH_IN)
      ENDSUBROUTINE %(binding_prefix)s%(proc_prefix)sINITIALISE

## if (spin_correlation) {
      SUBROUTINE %(binding_prefix)s%(proc_prefix)sSET_SPIN_CORRELATION_VECTORS(LEG_INDEX, N_VECTORS, VECTORS) bind(c)
        USE iso_c_binding, ONLY: c_int, c_double

        INCLUDE 'spin_correlations.inc'

        integer(c_int), intent(in)    :: LEG_INDEX
        integer(c_int), intent(in)    :: N_VECTORS
        integer(c_double), intent(in) :: VECTORS(MAX_N_SPIN_CORR_VECTORS,4)

        CALL %(proc_prefix)sSET_SPIN_CORRELATION_VECTORS(LEG_INDEX, N_VECTORS, VECTORS)

      END SUBROUTINE %(binding_prefix)s%(proc_prefix)sSET_SPIN_CORRELATION_VECTORS

      SUBROUTINE %(binding_prefix)s%(proc_prefix)sRESET_SPIN_CORRELATION_VECTORS() bind(c)

        CALL %(proc_prefix)sRESET_SPIN_CORRELATION_VECTORS()

      END SUBROUTINE %(binding_prefix)s%(proc_prefix)sSET_SPIN_CORRELATION_VECTORS
## }

## if (color_correlation) {
      SUBROUTINE %(binding_prefix)s%(proc_prefix)sSET_COLOR_CORRELATORS_TO_CONSIDER(FIRST_CONNECTION, SECOND_CONNECTION) bind(c)
        USE iso_c_binding, ONLY: c_int
        integer(c_int), intent(in)  :: FIRST_CONNECTION
        integer(c_int), intent(in)  :: SECOND_CONNECTION

        CALL %(proc_prefix)sSET_COLOR_CORRELATORS_TO_CONSIDER(FIRST_CONNECTION, SECOND_CONNECTION)

      END SUBROUTINE %(binding_prefix)s%(proc_prefix)sSET_COLOR_CORRELATORS_TO_CONSIDER

      SUBROUTINE %(binding_prefix)s%(proc_prefix)sADD_COLOR_CORRELATORS_TO_CONSIDER(FIRST_CONNECTION, SECOND_CONNECTION) bind(c)
        USE iso_c_binding, ONLY: c_int
        integer(c_int), intent(in)  :: FIRST_CONNECTION
        integer(c_int), intent(in)  :: SECOND_CONNECTION

        CALL %(proc_prefix)sADD_COLOR_CORRELATORS_TO_CONSIDER(FIRST_CONNECTION, SECOND_CONNECTION)

      END SUBROUTINE %(binding_prefix)s%(proc_prefix)sADD_COLOR_CORRELATORS_TO_CONSIDER
## }