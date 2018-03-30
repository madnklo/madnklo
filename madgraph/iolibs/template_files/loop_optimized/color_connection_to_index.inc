C    =========================================================
C    Mapping of %(order)s color connections to a unique index
C    =========================================================

      SUBROUTINE %(proc_prefix)sGET_N_COLOR_CONNECTIONS_%(order)s(N_CC)
        IMPLICIT NONE
        INTEGER N_CC_FOR_THAT_ORDER
        PARAMETER (N_CC_FOR_THAT_ORDER=%(n_cc_for_that_order)d)
        INTEGER N_CC
CF2PY INTENT(OUT) :: N_CC
        N_CC = N_CC_FOR_THAT_ORDER
      END SUBROUTINE %(proc_prefix)sGET_N_COLOR_CONNECTIONS_%(order)s
         
      SUBROUTINE %(proc_prefix)sAPPLY_COLOR_CONNECTION_DEFINITIONS_%(order)s()
C     This routine initialise the mapping of the color connections for order %(order)s
      
      INTEGER    PERT_ORDER
      PARAMETER (PERT_ORDER=%(pert_order)d)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    N_CC_DEFINITION_ENTRIES
      PARAMETER (N_CC_DEFINITION_ENTRIES=%(n_cc_definition_entries)s)
      
      %(color_connection_definition)s
      LOGICAL COLOR_CONNECTION_INITIALIZED
C     If the color connection doesn't exist, the index returned will be zero
      DATA COLOR_CONNECTIONS/N_CC_DEFINITION_ENTRIES*0/
      DATA COLOR_CONNECTION_INITIALIZED/.False./
      COMMON/%(proc_prefix)sCOLOR_CONNECTION_DEFINITIONS_%(order)s/COLOR_CONNECTIONS, COLOR_CONNECTION_INITIALIZED

%(color_connection_entries)s
      
      END SUBROUTINE %(proc_prefix)sAPPLY_COLOR_CONNECTION_DEFINITIONS_%(order)s
      
      SUBROUTINE %(proc_prefix)sCOLOR_CONNECTION_TO_INDEX_%(order)s(%(color_connection_arguments)s, INDEX_OUT)
	  IMPLICIT NONE
C
C CONSTANTS
C
      INTEGER    PERT_ORDER
      PARAMETER (PERT_ORDER=%(pert_order)d)
      INTEGER    CC_SPECIFIER_DIMENSION
      PARAMETER (CC_SPECIFIER_DIMENSION=2*PERT_ORDER)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)      
C  
C ARGUMENTS 
C  
	  %(color_connection_arguments_definition)s
%(color_connection_arguments_definition_f2py)s
	  
	  INTEGER INDEX_OUT
CF2PY INTENT(OUT) :: INDEX_OUT
	  
C
C LOCAL VARIABLES
C
      INTEGER I, J
      LOGICAL MATCHED
      INTEGER CC_SPECIFIER(CC_SPECIFIER_DIMENSION)
C
C GLOBAL VARIABLES
C
      %(color_connection_definition)s
      LOGICAL COLOR_CONNECTION_INITIALIZED
      COMMON/%(proc_prefix)sCOLOR_CONNECTION_DEFINITIONS_%(order)s/COLOR_CONNECTIONS, COLOR_CONNECTION_INITIALIZED

C ----------
C BEGIN CODE
C ----------

      IF (.NOT.COLOR_CONNECTION_INITIALIZED) THEN
        COLOR_CONNECTION_INITIALIZED = .True.
	    CALL %(proc_prefix)sAPPLY_COLOR_CONNECTION_DEFINITIONS_%(order)s()
      ENDIF
     
C     Set the cc_specifier
      %(set_cc_specifier)s
      
C     And return it now
      %(return_connection_index)s

	  END SUBROUTINE %(proc_prefix)sCOLOR_CONNECTION_TO_INDEX_%(order)s
	  
      SUBROUTINE %(proc_prefix)sINDEX_TO_COLOR_CONNECTION_%(order)s(INDEX_IN,%(color_connection_arguments)s)
	  IMPLICIT NONE
C
C CONSTANTS
C
      INTEGER    PERT_ORDER
      PARAMETER (PERT_ORDER=%(pert_order)d)
      INTEGER NCOLORCONNECTIONS_%(order)s
      PARAMETER (NCOLORCONNECTIONS_%(order)s=%(n_color_connections_for_that_order)d)
      INTEGER INDEX_%(order)s_OFFSET
      PARAMETER (INDEX_%(order)s_OFFSET=%(index_offset_for_that_order)d)      
      
C  
C ARGUMENTS 
C  
	  %(color_connection_arguments_definition)s
%(color_connection_arguments_definition_f2py_out)s
	  
	  INTEGER INDEX_IN
CF2PY INTENT(IN) :: INDEX_IN
	  
C
C LOCAL VARIABLES
C
      INTEGER I, J, K
C
C GLOBAL VARIABLES
C
      INTEGER CC_INDEX_TO_DEFINITION_%(order)s(NCOLORCONNECTIONS_%(order)s,PERT_ORDER,3)
%(cc_index_to_definition_data)s
      COMMON/%(proc_prefix)sCC_INDEX_TO_DEFINITION_%(order)s/CC_INDEX_TO_DEFINITION_%(order)s

C ----------
C BEGIN CODE
C ----------

C     Return the color connection for this index
      %(set_returned_cc_for_index)s

	  END SUBROUTINE %(proc_prefix)sINDEX_TO_COLOR_CONNECTION_%(order)s

      FUNCTION %(proc_prefix)sGET_COLOR_CONNECTION_INDEX_%(order)s(%(color_connection_arguments)s)
	  IMPLICIT NONE     
C  
C ARGUMENTS 
C  
	  %(color_connection_arguments_definition)s
      
      INTEGER %(proc_prefix)sGET_COLOR_CONNECTION_INDEX_%(order)s
C
C LOCAL VARIABLES
C
      INTEGER CC_INDEX
 
C ----------
C BEGIN CODE
C ----------

      CALL %(proc_prefix)sCOLOR_CONNECTION_TO_INDEX_%(order)s(%(color_connection_arguments)s, CC_INDEX)
      %(proc_prefix)sGET_COLOR_CONNECTION_INDEX_%(order)s = CC_INDEX
      
      END FUNCTION %(proc_prefix)sGET_COLOR_CONNECTION_INDEX_%(order)s