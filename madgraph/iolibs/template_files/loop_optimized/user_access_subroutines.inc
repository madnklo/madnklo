C --=========================================--
C   Access subroutines that can be used by the
C   user to steer the computation.
C --=========================================--

C --=========================================--
C   Functions for dealing with the ordering
C   and indexing of split order contributions
C --=========================================--

      SUBROUTINE %(proc_prefix)sGET_NSQSO_LOOP(NSQSO)
C
C     Simple subroutine returning the number of squared split order
C     contributions returned in ANS when calling sloopmatrix 
C
      INTEGER    NSQUAREDSO
	  PARAMETER (NSQUAREDSO=%(nSquaredSO)d)

	  INTEGER NSQSO

	  NSQSO=NSQUAREDSO

      END

      SUBROUTINE %(proc_prefix)sGET_ANSWER_DIMENSION(ANS_DIM)
C
C     MadLoop subroutines return an array of dimension ANS(0:3,0:ANS_DIM)
C     In order for the user program to be able to correctly declare this
c     array when calling MadLoop, this subroutine returns its dimension
C
      INTEGER    NSQUAREDSO
	  PARAMETER (NSQUAREDSO=%(nSquaredSO)d)
	  INTEGER ANS_DIM 

	  %(get_nsqso_born)s 

	  ANS_DIM=MAX(NSQSO_BORN,NSQUAREDSO)

      END

      INTEGER FUNCTION %(proc_prefix)sML5SOINDEX_FOR_SQUARED_ORDERS(ORDERS)
C
C This functions returns the integer index identifying the split orders list passed in argument which correspond to the values of the following list of couplings (and in this order):
C %(split_order_str_list)s
C
C CONSTANTS
C
      INTEGER    NSO, NSQSO
	  PARAMETER (NSO=%(nSO)d, NSQSO=%(nSquaredSO)d)
C
C ARGUMENTS
C
      INTEGER ORDERS(NSO)
C
C LOCAL VARIABLES
C
      INTEGER I,J
	  INTEGER SQPLITORDERS(NSQSO,NSO)
	  %(SquaredSO)s
	  COMMON/%(proc_prefix)sML5SQPLITORDERS/SQPLITORDERS
C
C BEGIN CODE
C
      DO I=1,NSQSO
	    DO J=1,NSO
		  IF (ORDERS(J).NE.SQPLITORDERS(I,J)) GOTO 1009
		ENDDO
		%(proc_prefix)sML5SOINDEX_FOR_SQUARED_ORDERS = I
		RETURN
1009    CONTINUE
	  ENDDO

	  WRITE(*,*) 'ERROR:: Stopping function %(proc_prefix)sML5SOINDEX_FOR_SQUARED_ORDERS'
	  WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
	  STOP

      END

      INTEGER FUNCTION %(proc_prefix)sML5SOINDEX_FOR_BORN_AMP(AMPID)
C
C For a given born amplitude number, it returns the ID of the split orders it has
C
C CONSTANTS
C
	  INTEGER NBORNAMPS
	  PARAMETER (NBORNAMPS=%(nBornAmps)d)
C
C ARGUMENTS
C
      INTEGER AMPID
C
C LOCAL VARIABLES
C
	  INTEGER BORNAMPORDERS(NBORNAMPS)
	  %(BornAmpSO)s
C -----------
C BEGIN CODE
C -----------
	  IF (AMPID.gt.NBORNAMPS) THEN
	    WRITE(*,*) 'ERROR:: Born amplitude ID ',AMPID,' above the maximum ',NBORNAMPS	
      ENDIF
	  %(proc_prefix)sML5SOINDEX_FOR_BORN_AMP = BORNAMPORDERS(AMPID)

      END

      INTEGER FUNCTION %(proc_prefix)sML5SOINDEX_FOR_LOOP_AMP(AMPID)
C
C For a given loop amplitude number, it returns the ID of the split orders it has
C
C CONSTANTS
C
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=%(nloopamps)d)
C
C ARGUMENTS
C
      INTEGER AMPID
C
C LOCAL VARIABLES
C
      INTEGER LOOPAMPORDERS(NLOOPAMPS)
	  %(loopAmpSO)s
C -----------
C BEGIN CODE
C -----------
	  IF (AMPID.gt.NLOOPAMPS) THEN
	    WRITE(*,*) 'ERROR:: Loop amplitude ID ',AMPID,' above the maximum ',NLOOPAMPS	
	  ENDIF
	  %(proc_prefix)sML5SOINDEX_FOR_LOOP_AMP = LOOPAMPORDERS(AMPID)

      END


      INTEGER FUNCTION %(proc_prefix)sML5SQSOINDEX(ORDERINDEXA, ORDERINDEXB)
C
C This functions plays the role of the interference matrix. It can be hardcoded or 
C made more elegant using hashtables if its execution speed ever becomes a relevant
C factor. From two split order indices, it return the corresponding index in the squared 
c order canonical ordering.
C
C CONSTANTS
C      
      INTEGER    NSO, NSQUAREDSO, NAMPSO
	  PARAMETER (NSO=%(nSO)d, NSQUAREDSO=%(nSquaredSO)d, NAMPSO=%(nAmpSO)d)
C
C ARGUMENTS
C
	  INTEGER ORDERINDEXA, ORDERINDEXB
C
C LOCAL VARIABLES
C
      INTEGER I, SQORDERS(NSO)
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
	  %(ampsplitorders)s
	  COMMON/%(proc_prefix)sML5AMPSPLITORDERS/AMPSPLITORDERS
C
C FUNCTION
C
      INTEGER %(proc_prefix)sML5SOINDEX_FOR_SQUARED_ORDERS
C
C BEGIN CODE
C
      DO I=1,NSO
	    SQORDERS(I)=AMPSPLITORDERS(ORDERINDEXA,I)+AMPSPLITORDERS(ORDERINDEXB,I)
	  ENDDO
	  %(proc_prefix)sML5SQSOINDEX=%(proc_prefix)sML5SOINDEX_FOR_SQUARED_ORDERS(SQORDERS)
	  END

C This is the inverse subroutine of ML5SOINDEX_FOR_SQUARED_ORDERS. Not directly useful, but provided nonetheless.
      SUBROUTINE %(proc_prefix)sML5GET_SQUARED_ORDERS_FOR_SOINDEX(SOINDEX,ORDERS)
C
C This functions returns the orders identified by the squared split order index in argument. Order values correspond to following list of couplings (and in this order):
C %(split_order_str_list)s
C
C CONSTANTS
C
      INTEGER    NSO, NSQSO
	  PARAMETER (NSO=%(nSO)d, NSQSO=%(nSquaredSO)d)
C
C ARGUMENTS
C
      INTEGER SOINDEX, ORDERS(NSO)
C
C LOCAL VARIABLES
C
      INTEGER I
	  INTEGER SQPLITORDERS(NSQSO,NSO)
	  COMMON/%(proc_prefix)sML5SQPLITORDERS/SQPLITORDERS	  
C
C BEGIN CODE
C
      IF (SOINDEX.gt.0.and.SOINDEX.le.NSQSO) THEN
        DO I=1,NSO
          ORDERS(I) =  SQPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

	  WRITE(*,*) 'ERROR:: Stopping function %(proc_prefix)sML5GET_SQUARED_ORDERS_FOR_SOINDEX'
	  WRITE(*,*) 'Could not find squared orders index ',SOINDEX
	  STOP

      END SUBROUTINE

C This is the inverse subroutine of getting amplitude SO orders. Not directly useful, but provided nonetheless.
      SUBROUTINE %(proc_prefix)sML5GET_ORDERS_FOR_AMPSOINDEX(SOINDEX,ORDERS)
C
C This functions returns the orders identified by the split order index in argument. Order values correspond to following list of couplings (and in this order):
C %(split_order_str_list)s
C
C CONSTANTS
C
      INTEGER    NSO, NAMPSO
	  PARAMETER (NSO=%(nSO)d, NAMPSO=%(nAmpSO)d)
C
C ARGUMENTS
C
      INTEGER SOINDEX, ORDERS(NSO)
C
C LOCAL VARIABLES
C
      INTEGER I
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
	  COMMON/%(proc_prefix)sML5AMPSPLITORDERS/AMPSPLITORDERS
C
C BEGIN CODE
C
      IF (SOINDEX.gt.0.and.SOINDEX.le.NAMPSO) THEN
        DO I=1,NSO
          ORDERS(I) =  AMPSPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

	  WRITE(*,*) 'ERROR:: Stopping function %(proc_prefix)sML5GET_ORDERS_FOR_AMPSOINDEX'
	  WRITE(*,*) 'Could not find amplitude split orders index ',SOINDEX
	  STOP

      END SUBROUTINE


C This function is not directly useful, but included for completeness
      INTEGER FUNCTION %(proc_prefix)sML5SOINDEX_FOR_AMPORDERS(ORDERS)
C
C This functions returns the integer index identifying the amplitude split orders passed in argument which correspond to the values of the following list of couplings (and in this order):
C %(split_order_str_list)s
C
C CONSTANTS
C
      INTEGER    NSO, NAMPSO
	  PARAMETER (NSO=%(nSO)d, NAMPSO=%(nAmpSO)d)
C
C ARGUMENTS
C
      INTEGER ORDERS(NSO)
C
C LOCAL VARIABLES
C
      INTEGER I,J
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
	  COMMON/%(proc_prefix)sML5AMPSPLITORDERS/AMPSPLITORDERS
C
C BEGIN CODE
C
      DO I=1,NAMPSO
	    DO J=1,NSO
		  IF (ORDERS(J).NE.AMPSPLITORDERS(I,J)) GOTO 1009
		ENDDO
		%(proc_prefix)sML5SOINDEX_FOR_AMPORDERS = I
		RETURN
1009    CONTINUE
	  ENDDO

	  WRITE(*,*) 'ERROR:: Stopping function %(proc_prefix)sML5SOINDEX_FOR_AMPORDERS'
	  WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
	  STOP

      END

C --=========================================--
C   Definition of additional access routines
C --=========================================--

	  SUBROUTINE %(proc_prefix)sCOLLIER_COMPUTE_UV_POLES(ONOFF)
C
C This function can be called by the MadLoop user so as to chose to have COLLIER
C compute the UV pole or not (it costs more time).
C
      LOGICAL ONOFF

	  include 'MadLoopParams.inc'

      LOGICAL FORCED_CHOICE_OF_COLLIER_UV_POLE_COMPUTATION, FORCED_CHOICE_OF_COLLIER_IR_POLE_COMPUTATION
      LOGICAL COLLIER_UV_POLE_COMPUTATION_CHOICE, COLLIER_IR_POLE_COMPUTATION_CHOICE
	  COMMON/%(proc_prefix)sCOLLIERPOLESFORCEDCHOICE/FORCED_CHOICE_OF_COLLIER_UV_POLE_COMPUTATION, FORCED_CHOICE_OF_COLLIER_IR_POLE_COMPUTATION,COLLIER_UV_POLE_COMPUTATION_CHOICE,COLLIER_IR_POLE_COMPUTATION_CHOICE

	  COLLIERComputeUVpoles                        = ONOFF
C     This is just so that if we read the param again, we don't overwrite the choice made here
	  FORCED_CHOICE_OF_COLLIER_UV_POLE_COMPUTATION = .TRUE.
      COLLIER_UV_POLE_COMPUTATION_CHOICE           = ONOFF

	  END SUBROUTINE

	  SUBROUTINE %(proc_prefix)sCOLLIER_COMPUTE_IR_POLES(ONOFF)
C
C This function can be called by the MadLoop user so as to chose to have COLLIER
C compute the IR pole or not (it costs more time).
C
      LOGICAL ONOFF

	  include 'MadLoopParams.inc'

      LOGICAL FORCED_CHOICE_OF_COLLIER_UV_POLE_COMPUTATION, FORCED_CHOICE_OF_COLLIER_IR_POLE_COMPUTATION
      LOGICAL COLLIER_UV_POLE_COMPUTATION_CHOICE, COLLIER_IR_POLE_COMPUTATION_CHOICE
	  COMMON/%(proc_prefix)sCOLLIERPOLESFORCEDCHOICE/FORCED_CHOICE_OF_COLLIER_UV_POLE_COMPUTATION, FORCED_CHOICE_OF_COLLIER_IR_POLE_COMPUTATION,COLLIER_UV_POLE_COMPUTATION_CHOICE,COLLIER_IR_POLE_COMPUTATION_CHOICE

	  COLLIERComputeIRpoles         = ONOFF
C     This is just so that if we read the param again, we don't overwrite the choice made here
	  FORCED_CHOICE_OF_COLLIER_IR_POLE_COMPUTATION = .TRUE.
      COLLIER_IR_POLE_COMPUTATION_CHOICE           = ONOFF

	  END SUBROUTINE

	  SUBROUTINE %(proc_prefix)sFORCE_STABILITY_CHECK(ONOFF)
C
C This function can be called by the MadLoop user so as to always have stability
C checked, even during initialisation, when calling the *_thres routines.
C
      LOGICAL ONOFF

	  LOGICAL BYPASS_CHECK, ALWAYS_TEST_STABILITY
	  DATA BYPASS_CHECK, ALWAYS_TEST_STABILITY /.FALSE.,.FALSE./
	  COMMON/%(proc_prefix)sBYPASS_CHECK/BYPASS_CHECK, ALWAYS_TEST_STABILITY

	  ALWAYS_TEST_STABILITY = ONOFF

	  END SUBROUTINE

	  SUBROUTINE %(proc_prefix)sSET_AUTOMATIC_CACHE_CLEARING(ONOFF)
C
C This function can be called by the MadLoop user so as to manually chose when
C to reset the TIR cache.
C
	  	IMPLICIT NONE

	    include 'MadLoopParams.inc'

		LOGICAL ONOFF

      	LOGICAL AUTOMATIC_CACHE_CLEARING
      	DATA AUTOMATIC_CACHE_CLEARING/.TRUE./
      	COMMON/%(proc_prefix)sRUNTIME_OPTIONS/AUTOMATIC_CACHE_CLEARING
      	
	    INTEGER N_DP_EVAL, N_QP_EVAL
	    COMMON/%(proc_prefix)sN_EVALS/N_DP_EVAL,N_QP_EVAL

## if(not TIRCaching){
       WRITE(*,*) 'Warning: No TIR caching implemented. Call to SET_AUTOMATIC_CACHE_CLEARING did nothing.'
## } else {

		AUTOMATIC_CACHE_CLEARING = ONOFF
		
		IF (NRotations_DP.ne.0.or.NRotations_QP.ne.0) THEN
		  WRITE(*,*) 'Warning: One cannot remove the TIR cache automatic clearing while at the same time keeping Lorentz rotations for stability tests.'
		  WRITE(*,*) 'MadLoop will therefore automatically set NRotations_DP and NRotations_QP to 0.'
		  NRotations_DP = 0
		  NRotations_QP = 0
  		  CALL %(proc_prefix)sSET_N_EVALS(N_DP_EVAL,N_QP_EVAL)
		ENDIF
## }
	  END SUBROUTINE
	  
	  SUBROUTINE %(proc_prefix)sSET_COUPLINGORDERS_TARGET(SOTARGET)
	  IMPLICIT NONE
C
C     This routine can be accessed by an external user to set the squared split order target.
C     If set to a value different than -1, the code will try to avoid computing anything which
C     does not contribute to contributions of squared split orders SQSO_TARGET and below.
C     This can considerably speed up the code. However, keep in mind that any contribution of
C     'squared order index' larger than SQSO_TARGET cannot be trust.
C
C     ARGUMENTS
C
      INTEGER SOTARGET
C
C     GLOBAL
C
	  INTEGER SQSO_TARGET
	  common/%(proc_prefix)sSOCHOICE/SQSO_TARGET
C ----------
C BEGIN CODE
C ----------
      SQSO_TARGET = SOTARGET
	  END

      SUBROUTINE %(proc_prefix)sSET_LEG_POLARIZATION(LEG_ID, LEG_POLARIZATION)
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER LEG_ID
      INTEGER LEG_POLARIZATION
C
C     LOCALS
C
      INTEGER I
      INTEGER LEG_POLARIZATIONS(0:5)
C     ----------
C     BEGIN CODE
C     ----------

      IF (LEG_POLARIZATION.eq.-10000) THEN
        LEG_POLARIZATIONS(0)=-1
        DO I=1,5
          LEG_POLARIZATIONS(I)=-10000
        ENDDO      
      ELSE
        LEG_POLARIZATIONS(0)=1
        LEG_POLARIZATIONS(1)=LEG_POLARIZATION
        DO I=2,5
          LEG_POLARIZATIONS(I)=-10000
        ENDDO
      ENDIF
      CALL %(proc_prefix)sSET_LEG_POLARIZATIONS(LEG_ID,LEG_POLARIZATIONS)

      END

      SUBROUTINE %(proc_prefix)sSET_LEG_POLARIZATIONS(LEG_ID, LEG_POLARIZATIONS)
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER NPOLENTRIES
      PARAMETER (NPOLENTRIES=(NEXTERNAL+1)*6)
      INTEGER    NCOMB
      PARAMETER (NCOMB=%(ncomb)d)
C
C     ARGUMENTS
C
      INTEGER LEG_ID
      INTEGER LEG_POLARIZATIONS(0:5)
C
C     LOCALS
C
      INTEGER I,J
      LOGICAL ALL_SUMMED_OVER
C
C     GLOBALS
C
C     Entry 0 of the first dimension is all -1 if there is no polarization requirement.
C     Then for each leg with ID legID, it is either summed over if
C     POLARIZATIONS(legID,0) is -1, or the list of helicity considered for that
C     leg is POLARIZATIONS(legID,1: POLARIZATIONS(legID,0)   ).
      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      DATA ((POLARIZATIONS(I,J),I=0,NEXTERNAL),J=0,5)/NPOLENTRIES*-1/
      COMMON/%(proc_prefix)sBEAM_POL/POLARIZATIONS

      INTEGER BORN_POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/%(proc_prefix)sBORN_BEAM_POL/BORN_POLARIZATIONS

C     ----------
C     BEGIN CODE
C     ----------

      IF (LEG_POLARIZATIONS(0).eq.-1) THEN
        DO I=0,5
          POLARIZATIONS(LEG_ID,I)=-1
        ENDDO
      ELSE
        DO I=0,LEG_POLARIZATIONS(0)
          POLARIZATIONS(LEG_ID,I)=LEG_POLARIZATIONS(I)
        ENDDO
        DO I=LEG_POLARIZATIONS(0)+1,5
          POLARIZATIONS(LEG_ID,I)=-10000        
        ENDDO
      ENDIF

      ALL_SUMMED_OVER = .True.
      DO I=1,NEXTERNAL
        IF (POLARIZATIONS(I,0).NE.-1) THEN
          ALL_SUMMED_OVER = .False.
          EXIT
        ENDIF
      ENDDO
      IF (ALL_SUMMED_OVER) THEN
        DO I=0,5
          POLARIZATIONS(0,I)=-1
        ENDDO
      ELSE
        DO I=0,5
          POLARIZATIONS(0,I)=0
        ENDDO
      ENDIF

      DO I=0,NEXTERNAL
        DO J=0,5
          BORN_POLARIZATIONS(I,J) = POLARIZATIONS(I,J)
        ENDDO
      ENDDO


      RETURN

      END

      SUBROUTINE %(proc_prefix)sSLOOPMATRIXHEL(P,HEL,ANS)
      IMPLICIT NONE
C  
C CONSTANTS
C
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NSQUAREDSO
      PARAMETER (NSQUAREDSO=%(nSquaredSO)d)
C  
C ARGUMENTS 
C  
      %(real_dp_format)s P(0:3,NEXTERNAL)
      %(real_dp_format)s ANS(0:3,0:NSQUAREDSO)
	  INTEGER HEL, USERHEL
	  common/%(proc_prefix)sUSERCHOICE/USERHEL
C ----------
C BEGIN CODE
C ----------
	  USERHEL=HEL
      CALL %(proc_prefix)sSLOOPMATRIX(P,ANS)
	  END

      SUBROUTINE %(proc_prefix)sSLOOPMATRIXHEL_THRES(P,HEL,ANS,PREC_ASKED,PREC_FOUND,RET_CODE)
	  IMPLICIT NONE
C  
C CONSTANTS
C
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NSQUAREDSO
      PARAMETER (NSQUAREDSO=%(nSquaredSO)d)
C  
C ARGUMENTS 
C  
      %(real_dp_format)s P(0:3,NEXTERNAL)
      %(real_dp_format)s ANS(0:3,0:NSQUAREDSO)
	  INTEGER HEL, RET_CODE
	  %(real_dp_format)s PREC_ASKED,PREC_FOUND(0:NSQUAREDSO)
C
C LOCAL VARIABLES
C
      INTEGER I
C
C GLOBAL VARIABLES
C
	  %(real_dp_format)s USER_STAB_PREC
	  COMMON/%(proc_prefix)sUSER_STAB_PREC/USER_STAB_PREC

	  INTEGER H,T,U
	  %(real_dp_format)s ACCURACY(0:NSQUAREDSO)
	  common/%(proc_prefix)sACC/ACCURACY,H,T,U

	  LOGICAL BYPASS_CHECK, ALWAYS_TEST_STABILITY
	  COMMON/%(proc_prefix)sBYPASS_CHECK/BYPASS_CHECK, ALWAYS_TEST_STABILITY

C ----------
C BEGIN CODE
C ----------
      USER_STAB_PREC = PREC_ASKED

      CALL %(proc_prefix)sSLOOPMATRIXHEL(P,HEL,ANS)
	  IF(ALWAYS_TEST_STABILITY.AND.(H.eq.1.OR.ACCURACY(0).lt.0.0d0)) THEN
	    BYPASS_CHECK = .TRUE.
        CALL %(proc_prefix)sSLOOPMATRIXHEL(P,HEL,ANS) 
	    BYPASS_CHECK = .FALSE.
C       Make sure we correctly return an initialization-type T code
		IF (T.eq.2) T=4
		IF (T.eq.1) T=3
	  ENDIF

C Reset it to default value not to affect next runs
	  USER_STAB_PREC = -1.0d0
	  
      DO I=0,NSQUAREDSO
	    PREC_FOUND(I)=ACCURACY(I)
	  ENDDO
      RET_CODE=100*H+10*T+U

	  END

      SUBROUTINE %(proc_prefix)sSLOOPMATRIX_THRES(P,ANS,PREC_ASKED,PREC_FOUND,RET_CODE)
C
C     Inputs are:
C     P(0:3, Nexternal)  double  :: Kinematic configuration (E,px,py,pz)
C     PEC_ASKED          double  :: Target relative accuracy, -1 for default
C
C     Outputs are:
C     ANS(3)             double  :: Result (finite, single pole, double pole) 
C     PREC_FOUND         double  :: Relative accuracy estimated for the result
C                                   Returns -1 if no stab test could be performed.
C	  RET_CODE			 integer :: Return code. See below for details
C
C     Return code conventions: RET_CODE = H*100 + T*10 + U
C
C     H == 1
C         Stability unknown.
C     H == 2
C         Stable PS (SPS) point.
C         No stability rescue was necessary.
C     H == 3
C         Unstable PS (UPS) point.
C         Stability rescue necessary, and successful.
C     H == 4
C         Exceptional PS (EPS) point.
C         Stability rescue attempted, but unsuccessful.
C
C     T == 1
C         Default computation (double prec.) was performed.
C     T == 2
C         Quadruple precision was used for this PS point.
C     T == 3
C         MadLoop in initialization phase. Only double precision used.
C     T == 4
C         MadLoop in initialization phase. Quadruple precision used.
C
C     U == 0
C         Not stable.
C     U == 1
C         Stable with CutTools in double precision.
C     U == 2
C         Stable with PJFry++.
C     U == 3
C         Stable with IREGI.
C     U == 4
C         Stable with Golem95
C     U == 5
C         Stable with Samurai
C     U == 6
C         Stable with Ninja in double precision
C     U == 8
C         Stable with Ninja in quadruple precision
C     U == 9
C         Stable with CutTools in quadruple precision.         
C
      IMPLICIT NONE
C  
C CONSTANTS
C
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NSQUAREDSO
      PARAMETER (NSQUAREDSO=%(nSquaredSO)d)
C  
C ARGUMENTS 
C  
      %(real_dp_format)s P(0:3,NEXTERNAL)
      %(real_dp_format)s ANS(0:3,0:NSQUAREDSO)
	  %(real_dp_format)s PREC_ASKED,PREC_FOUND(0:NSQUAREDSO)
	  INTEGER RET_CODE	  
C
C LOCAL VARIABLES
C
      INTEGER I
C
C GLOBAL VARIABLES
C
	  %(real_dp_format)s USER_STAB_PREC
	  COMMON/%(proc_prefix)sUSER_STAB_PREC/USER_STAB_PREC

	  INTEGER H,T,U
	  %(real_dp_format)s ACCURACY(0:NSQUAREDSO)
	  common/%(proc_prefix)sACC/ACCURACY,H,T,U

	  LOGICAL BYPASS_CHECK, ALWAYS_TEST_STABILITY
	  COMMON/%(proc_prefix)sBYPASS_CHECK/BYPASS_CHECK, ALWAYS_TEST_STABILITY

C ----------
C BEGIN CODE
C ----------
      USER_STAB_PREC = PREC_ASKED
      CALL %(proc_prefix)sSLOOPMATRIX(P,ANS)
	  IF(ALWAYS_TEST_STABILITY.AND.(H.eq.1.OR.ACCURACY(0).lt.0.0d0)) THEN
	    BYPASS_CHECK = .TRUE.	  
        CALL %(proc_prefix)sSLOOPMATRIX(P,ANS)
		BYPASS_CHECK = .FALSE.
C     Make sure we correctly return an initialization-type T code
		IF (T.eq.2) T=4
		IF (T.eq.1) T=3
	  ENDIF

C Reset it to default value not to affect next runs
	  USER_STAB_PREC = -1.0d0
      DO I=0,NSQUAREDSO	  
	    PREC_FOUND(I)=ACCURACY(I)
	  ENDDO
	  RET_CODE=100*H+10*T+U

	  END
	  
C The subroutine below perform clean-up duties for MadLoop like de-allocating
c arrays
	  SUBROUTINE %(proc_prefix)sEXIT_MADLOOP()
## if(ComputeColorFlows) {
	    CALL %(proc_prefix)sDEALLOCATE_COLOR_FLOWS()
## }
        CONTINUE
	  END

## if (spin_correlation) {

C --------------------------------------------------------------
C Beginning of a series of subroutine related to the evaluation
C and steering of the computation of spin-correlated MEs.
C --------------------------------------------------------------
      
      SUBROUTINE %(proc_prefix)sSET_SPIN_CORRELATION_VECTORS(LEG_INDEX, N_VECTORS, VECTORS)
        IMPLICIT NONE
C       Sets a list of up to MAX_N_VECTORS to use instead of the polarization vectors for leg number LEG_INDEX
C       Notice that of course only fermions and vectors can be assigned correlations here.

C       
C       CONSTANTS AND GLOBALS
C       

        INTEGER    NEXTERNAL
        PARAMETER (NEXTERNAL=%(nexternal)d)
        include 'spin_correlations.inc'     

C       
C       ARGUMENTS
C       
		INTEGER LEG_INDEX
		INTEGER N_VECTORS
		REAL*8 VECTORS(MAX_N_SPIN_CORR_VECTORS,4)
CF2PY INTENT(IN) :: LEG_INDEX
CF2PY INTENT(IN) :: N_VECTORS
CF2PY INTENT(IN) :: VECTORS

C       
C       LOCAL VARIABLES
C  
        
        INTEGER I,J,N
        INTEGER OLD_N_SPIN_CORR_RUNS
        INTEGER N_LEGS_WITH_SPIN_CORRELATORS

C
C       GLOBAL ARRAYS INITIALIZATION
C

        INTEGER SPIN_CORR_VECTORS_N_ENTRIES
        PARAMETER (SPIN_CORR_VECTORS_N_ENTRIES=(4*NEXTERNAL*MAX_N_SPIN_CORR_VECTORS))
        DATA SPIN_CORR_VECTORS/SPIN_CORR_VECTORS_N_ENTRIES*0.0d0/
        DATA SYSTEM_SPIN_CORR_VECTORS/SPIN_CORR_VECTORS_N_ENTRIES*0.0d0/
        DATA N_SPIN_CORR_VECTORS/NEXTERNAL*0/
        INTEGER SPIN_CORR_RUNS_N_ENTRIES
        PARAMETER (SPIN_CORR_RUNS_N_ENTRIES=(MAX_SPIN_CORR_RUNS+1)*NEXTERNAL)  
C       The integer -1000 will indicate to simply keep the original helicity polarization vector assigned to that leg      
        DATA SPIN_CORR_RUNS/SPIN_CORR_RUNS_N_ENTRIES*-1000/
        DATA N_SPIN_CORR_RUNS/0/

		INTEGER BYPASS_CHECK_VALUE_BACKUP
		COMMON/%(proc_prefix)sSPIN_CORR_PARAMS_BACKUP/BYPASS_CHECK_VALUE_BACKUP
	    LOGICAL BYPASS_CHECK, ALWAYS_TEST_STABILITY
	    COMMON/%(proc_prefix)sBYPASS_CHECK/BYPASS_CHECK, ALWAYS_TEST_STABILITY
C       ----------
C       BEGIN CODE
C       ----------

C       Turn off the helicity filter and loop checks as one should not use it in conjunction with spin-correlations.
C       Back-up its value however so as to be able to recover it when disabling again spin-correlations
	    BYPASS_CHECK_VALUE_BACKUP = BYPASS_CHECK
	    BYPASS_CHECK = .True.

C       Make sure the user does not try to set the spin-correlation vectors of a leg again
        IF (N_SPIN_CORR_VECTORS(LEG_INDEX).ne.0) THEN
          WRITE(*,*) 'Cannot define for a second time the spin-correlation vectors of a lef. Call the subroutine %(proc_prefix)sRESET_SPIN_CORRELATION_VECTORS() first.'
          STOP
        ENDIF

C       Make sure it does not exceed the maximum number of legs that can get spin correlators assigned.
		N_LEGS_WITH_SPIN_CORRELATORS = 0
		DO I=1,NEXTERNAL
		  IF (N_SPIN_CORR_VECTORS(I).gt.0) THEN
		    N_LEGS_WITH_SPIN_CORRELATORS = N_LEGS_WITH_SPIN_CORRELATORS + 1
		  ENDIF
		ENDDO
		IF (N_LEGS_WITH_SPIN_CORRELATORS.ge.MAX_LEGS_WITH_SPIN_CORR) THEN
		  WRITE(*,*) 'In the current output, only a maximum of ',MAX_LEGS_WITH_SPIN_CORR,' can be assigned spin-correlations.'
		  WRITE(*,*) 'Change the parameter MAX_LEGS_WITH_SPIN_CORR in spin_correlations.inc or re-generate the process with a higher power of N for the option --spin_correlators=N^iLO'
		  STOP
		ENDIF
		
C       First add the vectors to the list of spin-correlation vectors defined
        N_SPIN_CORR_VECTORS(LEG_INDEX) = N_VECTORS
        DO I=1,N_VECTORS
          DO J=1,4
C           Leave the flexibility of having these vectors complex, although we take them real here.
            SPIN_CORR_VECTORS(LEG_INDEX,I,J) = DCMPLX(VECTORS(I,J),0.0d0)
          ENDDO
        ENDDO
		
C       Now assign the additional runs to be performed within the loop over helicity configuration so as to capture these spin-correlations.
        OLD_N_SPIN_CORR_RUNS = N_SPIN_CORR_RUNS
        DO I=1,N_VECTORS
C         Add a run for each of the existing ones. This means that if user specifies vector VEC1A and VEC1B for leg 1, and VEC2A and VEC2B for leg 2, we will include all 4 combinations. (VEC1A, VEC2A), (VEC1A, VEC2B), etc...
C         The MIN(1,OLD_N_SPIN_CORR_RUNS) is just here because the first time we add spin correlation we must start from the blank entries in SPIN_CORR_RUNS(0,*)
          DO J=MIN(1,OLD_N_SPIN_CORR_RUNS),OLD_N_SPIN_CORR_RUNS
            N_SPIN_CORR_RUNS = N_SPIN_CORR_RUNS + 1
            DO N=1,NEXTERNAL
              IF (N.eq.LEG_INDEX) THEN
C               If the leg is the one selected, enforce the vectors specified by the user
                SPIN_CORR_RUNS(N_SPIN_CORR_RUNS,N) = I 
              ELSE
C               Otherwise, just copy the ones of the previous runs.
                SPIN_CORR_RUNS(N_SPIN_CORR_RUNS,N) = SPIN_CORR_RUNS(J,N)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        
C       That's it, the corresponding additional runs for the new spin correlations asked for have been added to the list SPIN_CORR_RUNS

      END SUBROUTINE %(proc_prefix)sSET_SPIN_CORRELATION_VECTORS
      
      SUBROUTINE %(proc_prefix)sRESET_SPIN_CORRELATION_VECTORS()
        IMPLICIT NONE
C       Resets the spin correlation vectors definition

C       
C       CONSTANTS AND GLOBALS
C       

        INTEGER    NEXTERNAL
        PARAMETER (NEXTERNAL=%(nexternal)d)
        include 'spin_correlations.inc'     

C       
C       LOCAL VARIABLES
C  
        
        INTEGER I,J,N

C       
C       GLOBAL VARIABLES
C
		INTEGER BYPASS_CHECK_VALUE_BACKUP
		COMMON/%(proc_prefix)sSPIN_CORR_PARAMS_BACKUP/BYPASS_CHECK_VALUE_BACKUP
	    LOGICAL BYPASS_CHECK, ALWAYS_TEST_STABILITY
	    COMMON/%(proc_prefix)sBYPASS_CHECK/BYPASS_CHECK, ALWAYS_TEST_STABILITY

C       ----------
C       BEGIN CODE
C       ----------
		
C       Restore the old value for bypass_check
		BYPASS_CHECK = BYPASS_CHECK_VALUE_BACKUP
		
		DO N=1,NEXTERNAL
		  N_SPIN_CORR_VECTORS(N) = 0
		  DO I=1,MAX_N_SPIN_CORR_VECTORS
		    DO J=1,4
		      SPIN_CORR_VECTORS(N,I,J) = 0.0d0
		    ENDDO
		  ENDDO
		ENDDO
		
		N_SPIN_CORR_RUNS = 0
		DO I=1,MAX_SPIN_CORR_RUNS
		  DO J=1,NEXTERNAL
		    SPIN_CORR_RUNS(I,J) = -1000
		  ENDDO
		ENDDO
		
		END SUBROUTINE %(proc_prefix)sRESET_SPIN_CORRELATION_VECTORS

      SUBROUTINE %(proc_prefix)sGET_MAX_N_SPIN_CORR_VECTORS(MAX_N_SC_VECTORS)
        IMPLICIT NONE
        INTEGER    NEXTERNAL
        PARAMETER (NEXTERNAL=%(nexternal)d)
        include 'spin_correlations.inc'
        INTEGER MAX_N_SC_VECTORS
CF2PY INTENT(OUT) :: MAX_N_SC_VECTORS
		MAX_N_SC_VECTORS = MAX_N_SPIN_CORR_VECTORS
      END SUBROUTINE %(proc_prefix)sGET_MAX_N_SPIN_CORR_VECTORS
      
      SUBROUTINE %(proc_prefix)sGET_MAX_N_SPIN_CORR_LEGS(MAX_N_SC_LEGS)
        IMPLICIT NONE
        INTEGER    NEXTERNAL
        PARAMETER (NEXTERNAL=%(nexternal)d)
        include 'spin_correlations.inc'
        INTEGER MAX_N_SC_LEGS
CF2PY INTENT(OUT) :: MAX_N_SC_LEGS
		MAX_N_SC_LEGS = MAX_LEGS_WITH_SPIN_CORR
      END SUBROUTINE %(proc_prefix)sGET_MAX_N_SPIN_CORR_LEGS
## }
