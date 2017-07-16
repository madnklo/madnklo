      SUBROUTINE %(proc_prefix)sSMATRIX(P,ANS_SUMMED)
C
C     Simple routine wrapper to provide the same interface for
C     backward compatibility for usage without split orders.
C
C  
C CONSTANTS
C 
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER NSQAMPSO
      PARAMETER (NSQAMPSO=%(nSqAmpSplitOrders)d)
C  
C ARGUMENTS 
C 
      REAL*8 P(0:3,NEXTERNAL), ANS_SUMMED
C
C VARIABLES
C
      INTEGER I
	  REAL*8 ANS(0:NSQAMPSO)
C
C BEGIN CODE
C
	  CALL %(proc_prefix)sSMATRIX_SPLITORDERS(P,ANS)
	  ANS_SUMMED=ANS(0)

	  END

      SUBROUTINE %(proc_prefix)sSMATRIXHEL(P,HEL,ANS)
      IMPLICIT NONE
C
C CONSTANT
C
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=%(ncomb)d)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),ANS
	  INTEGER HEL
C
C GLOBAL VARIABLES
C
      INTEGER USERHEL
	  COMMON/%(proc_prefix)sHELUSERCHOICE/USERHEL
C ----------
C BEGIN CODE
C ----------
      USERHEL=HEL
	  CALL %(proc_prefix)sSMATRIX(P,ANS)
	  USERHEL=-1

	  END

      SUBROUTINE %(proc_prefix)sSMATRIX_SPLITORDERS(P,ANS)
C  
%(info_lines)s
C 
C MadGraph StandAlone Version
C 
C Returns amplitude squared summed/avg over colors
c and helicities
c for the point in phase space P(0:3,NEXTERNAL)
C  
%(process_lines)s
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
	  INTEGER    NINITIAL 
      PARAMETER (NINITIAL=%(nincoming)d)
      INTEGER NPOLENTRIES
      PARAMETER (NPOLENTRIES=(NEXTERNAL+1)*6)
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=%(ncomb)d)
      INTEGER NSQAMPSO
      PARAMETER (NSQAMPSO=%(nSqAmpSplitOrders)d)
	  INTEGER HELAVGFACTOR
	  PARAMETER (HELAVGFACTOR=%(hel_avg_factor)d)
	  LOGICAL CHOSEN_SO_CONFIGS(NSQAMPSO)
	  DATA CHOSEN_SO_CONFIGS/%(chosen_so_configs)s/
	  COMMON/%(proc_prefix)sCHOSEN_BORN_SQSO/CHOSEN_SO_CONFIGS
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),ANS(0:NSQAMPSO)
C  
C LOCAL VARIABLES 
C
      LOGICAL DO_INCLUDE_HEL_CONTRIB
      INTEGER NTRY
      REAL*8 T(NSQAMPSO), BUFF
      INTEGER IHEL,IDEN, I, J
C     For a 1>N process, them BEAMTWO_HELAVGFACTOR would be set to 1.
      INTEGER BEAMS_HELAVGFACTOR(2)
	  DATA (BEAMS_HELAVGFACTOR(I),I=1,2)/%(beamone_helavgfactor)d,%(beamtwo_helavgfactor)d/
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./
      
      INTEGER HELICITIES(NEXTERNAL)
      
%(den_factor_line)s
C
C GLOBAL VARIABLES
C
      INTEGER NHEL(NEXTERNAL,NCOMB)
%(helicity_lines)s
      COMMON/%(proc_prefix)sBORN_HEL_CONFIGS/NHEL
      
      INTEGER USERHEL
	  DATA USERHEL/-1/	  	  
	  COMMON/%(proc_prefix)sHELUSERCHOICE/USERHEL
      
	  INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      DATA ((POLARIZATIONS(I,J),I=0,NEXTERNAL),J=0,5)/NPOLENTRIES*-1/
      COMMON/%(proc_prefix)sBORN_BEAM_POL/POLARIZATIONS
      
## if (color_correlation) {
C
C COLOR CORRELATORS DATA
C
C     CCI will be the general color correlated index in loops
C     and CCIA is the index of the loop over COLOR_CORRELATORS_TO_CONSIDER
      INTEGER CCI, CCIA
      
      INTEGER NCOLORCORRELATORS
      PARAMETER (NCOLORCORRELATORS=%(n_color_correlators)d)
      INTEGER NCC_TIMES_NSQAMPSO
      PARAMETER (NCC_TIMES_NSQAMPSO=NCOLORCORRELATORS*(NSQAMPSO+1))
      
C     Index 0 is the number of correlators to consider and the next indices are which one to consider
	  INTEGER COLOR_CORRELATORS_TO_CONSIDER(0:NCOLORCORRELATORS)   
      DATA COLOR_CORRELATORS_TO_CONSIDER(0)/NCOLORCORRELATORS/
%(color_correlators_to_consider_initialization)s

      REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS,0:NSQAMPSO)
      DATA ((COLOR_CORRELATED_EVALS(I,J), I=1,NCOLORCORRELATORS),J=0,NSQAMPSO) / NCC_TIMES_NSQAMPSO*0.0d0 /       
      COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER, COLOR_CORRELATED_EVALS
## }

## if (spin_correlation) {
      include 'spin_correlations.inc'
      INTEGER SPINCORR_HEL_INDICES(NEXTERNAL)
      LOGICAL SKIP_THIS_HELICITY
      INTEGER ISC
      COMPLEX*16 WEXT(20,NEXTERNAL)
## }
 
C
C     FUNCTIONS
C
      LOGICAL %(proc_prefix)sIS_BORN_HEL_SELECTED

C ----------
C BEGIN CODE
C ----------

## if (color_correlation) {
C
C Reinitialize color correlated values
C
     DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
       CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
       DO J=1,NSQAMPSO
         COLOR_CORRELATED_EVALS(CCI,J) = 0.0d0
       ENDDO
     ENDDO
## }

      NTRY=NTRY+1
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
	  DO I=1,NSQAMPSO
      	ANS(I) = 0D0	  
	  ENDDO
C     When spin-2 particles are involved, the Helicity filtering is dangerous for the 2->1 topology.
C     This is because depending on the MC setup the initial PS points have back-to-back initial states
C     for which some of the spin-2 helicity configurations are zero. But they are no longer zero
C     if the point is boosted on the z-axis. Remember that HELAS helicity amplitudes are no longer
C     lorentz invariant with expternal spin-2 particles (only the helicity sum is).
C     For this reason, we simply remove the filterin when there is only three external particles.
      IF (NEXTERNAL.LE.3) THEN
        DO IHEL=1,NCOMB
		  GOODHEL(IHEL)=.TRUE.
		ENDDO
	  ENDIF
	  
## if (spin_correlation) {
C     If spin correlations are active, then remove the helicity filtering.
      IF (N_SPIN_CORR_RUNS.GT.0) THEN
        DO IHEL=1,NCOMB
		  GOODHEL(IHEL)=.TRUE.
		ENDDO   
      ENDIF
## }

## if (spin_correlation) {
C     Save helicity indices of the spin-correlated legs that we will use to insert the extra spin. corr. runs.
      DO I=1,NEXTERNAL
        SPINCORR_HEL_INDICES(I) = -1
      ENDDO

## }
      DO IHEL=1,NCOMB
	     IF (USERHEL.EQ.-1.OR.USERHEL.EQ.IHEL) THEN
           IF (GOODHEL(IHEL) .OR. NTRY .LT. 2 .OR.USERHEL.NE.-1) THEN
              IF(NTRY.GE.2.AND.POLARIZATIONS(0,0).ne.-1.and.(.not.%(proc_prefix)sIS_BORN_HEL_SELECTED(IHEL))) THEN
                CYCLE
              ENDIF
## if (spin_correlation) {

              IF (N_SPIN_CORR_RUNS.GT.0) THEN
C               Now we know that we are in an helicity configuration that was selected by the user
C               Record the helicity index of legs that are spin-correlated
                DO I=1,NEXTERNAL
                  IF (N_SPIN_CORR_VECTORS(I).gt.0.and.SPINCORR_HEL_INDICES(I).eq.-1) THEN
                    SPINCORR_HEL_INDICES(I) = NHEL(I,IHEL)
                  ENDIF
                ENDDO

C               Now make sure to skip all helicity configurations with helicity indices for the spin-correlated legs that
C               does not match the first one encountered
                SKIP_THIS_HELICITY = .FALSE.
                DO I=1,NEXTERNAL
                  IF (SPINCORR_HEL_INDICES(I).ne.-1.and.SPINCORR_HEL_INDICES(I).ne.NHEL(I,IHEL)) THEN
                    SKIP_THIS_HELICITY = .TRUE.
                    EXIT
                  ENDIF
                ENDDO
                IF (SKIP_THIS_HELICITY) THEN
                  CYCLE
                ENDIF
              ENDIF
              
C             Now insert the extra runs for the spin-correlations
C             The MIN statement allows to make sure that we proceed with a normal run if not spin-correlations are active
              DO ISC=MIN(N_SPIN_CORR_RUNS,1),N_SPIN_CORR_RUNS              
## }	     
              DO_INCLUDE_HEL_CONTRIB = (POLARIZATIONS(0,0).eq.-1.or.%(proc_prefix)sIS_BORN_HEL_SELECTED(IHEL))
              DO I=1,NEXTERNAL
## if (not spin_correlation) {
                HELICITIES(I) = NHEL(1,IHEL)
## } else {
                IF (SPIN_CORR_RUNS(ISC, I).eq.-1000) THEN
C                  No spin correlation, leave it as it is
                   HELICITIES(I) = NHEL(1,IHEL)
                   DO J=1,20
                     WEXT(J,I) = DCMPLX(0.0d0,0.0d0)
                   ENDDO
                ELSE
C                  Spin correlation, enforce the polarization vector
                   HELICITIES(I) = -1000
                   DO J=1,4
                     WEXT(J,I) = SPIN_CORR_VECTORS(I,SPIN_CORR_RUNS(ISC, I),J)
                   ENDDO
                   DO J=5,20
                     WEXT(J,I) = DCMPLX(0.0d0,0.0d0)
                   ENDDO
                ENDIF
## }
              ENDDO
## if (not color_correlation) {
## if (not spin_correlation) {
              CALL %(proc_prefix)sMATRIX(P ,HELICITIES(1),JC(1), T)
## } else {
              CALL %(proc_prefix)sMATRIX(P ,HELICITIES(1),JC(1), T, WEXT(1,1))
## }
## } else {
## if (not spin_correlation) {
              CALL %(proc_prefix)sMATRIX(P ,HELICITIES(1),JC(1), T, DO_INCLUDE_HEL_CONTRIB)
## } else {
              CALL %(proc_prefix)sMATRIX(P ,HELICITIES(1),JC(1), T, DO_INCLUDE_HEL_CONTRIB, WEXT(1,1))
## }
## }
	 	      BUFF=0D0
			  DO I=1,NSQAMPSO
                IF(DO_INCLUDE_HEL_CONTRIB) THEN
                  ANS(I)=ANS(I)+T(I)
				ENDIF
			    BUFF=BUFF+T(I)
		      ENDDO
              IF (BUFF .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
                GOODHEL(IHEL)=.TRUE.
              ENDIF
## if (spin_correlation) {
C             End of the do-loop for the insertion of the extra spin-correlation runs
              ENDDO
##}           
            ENDIF
		  ENDIF
       ENDDO
	  ANS(0)=0.0d0
## if (color_correlation) {
     DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
       CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
       COLOR_CORRELATED_EVALS(CCI,0) = 0.0d0
     ENDDO
## }

	  DO I=1,NSQAMPSO
        ANS(I)=ANS(I)/DBLE(IDEN)
## if (color_correlation) {
        DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
           CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
           COLOR_CORRELATED_EVALS(CCI,I) = COLOR_CORRELATED_EVALS(CCI,I)/DBLE(IDEN)
        ENDDO
## }
		IF (CHOSEN_SO_CONFIGS(I)) THEN
		  ANS(0)=ANS(0)+ANS(I)
## if (color_correlation) {
          DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
             CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
             COLOR_CORRELATED_EVALS(CCI,0) = COLOR_CORRELATED_EVALS(CCI,0) + COLOR_CORRELATED_EVALS(CCI,I)
          ENDDO
## }
		ENDIF
	  ENDDO
	  IF(USERHEL.NE.-1) THEN
		DO I=0,NSQAMPSO
          ANS(I)=ANS(I)*HELAVGFACTOR
## if (color_correlation) {
          DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
             CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
             COLOR_CORRELATED_EVALS(CCI,I) = COLOR_CORRELATED_EVALS(CCI,I)*HELAVGFACTOR
          ENDDO
## }
		ENDDO
	  ELSE
		DO J=1,NINITIAL
		  IF (POLARIZATIONS(J,0).ne.-1) THEN
		    DO I=0,NSQAMPSO
              ANS(I)=ANS(I)*BEAMS_HELAVGFACTOR(J)
              ANS(I)=ANS(I)/POLARIZATIONS(J,0)
## if (color_correlation) {
              DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
                CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
                COLOR_CORRELATED_EVALS(CCI,I) = COLOR_CORRELATED_EVALS(CCI,I)*BEAMS_HELAVGFACTOR(J)
                COLOR_CORRELATED_EVALS(CCI,I) = COLOR_CORRELATED_EVALS(CCI,I)/POLARIZATIONS(J,0)  
              ENDDO
## }
		    ENDDO
		  ENDIF
		ENDDO
	  ENDIF
      END

      SUBROUTINE %(proc_prefix)sSMATRIXHEL_SPLITORDERS(P,HEL,ANS)
      IMPLICIT NONE
C
C CONSTANT
C
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=%(ncomb)d)
      INTEGER NSQAMPSO
      PARAMETER (NSQAMPSO=%(nSqAmpSplitOrders)d)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),ANS(0:NSQAMPSO)
	  INTEGER HEL
C
C GLOBAL VARIABLES
C
      INTEGER USERHEL
	  COMMON/%(proc_prefix)sHELUSERCHOICE/USERHEL
C ----------
C BEGIN CODE
C ----------
      USERHEL=HEL
	  CALL %(proc_prefix)sSMATRIX_SPLITORDERS(P,ANS)
	  USERHEL=-1

	  END
 
## if (not color_correlation) {
## if (not spin_correlation) {
      SUBROUTINE %(proc_prefix)sMATRIX(P,NHEL,IC,RES)
## } else {
      SUBROUTINE %(proc_prefix)sMATRIX(P,NHEL,IC,RES,WEXT)
## }
## } else {
## if (not spin_correlation) {
      SUBROUTINE %(proc_prefix)sMATRIX(P,NHEL,IC,RES, DO_INCLUDE_HEL_CONTRIB)
## } else {
      SUBROUTINE %(proc_prefix)sMATRIX(P,NHEL,IC,RES, DO_INCLUDE_HEL_CONTRIB,WEXT)
## }
## }

C  
%(info_lines)s
C
C Returns amplitude squared summed/avg over colors
c for the point with external lines W(0:6,NEXTERNAL)
C  
%(process_lines)s
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=%(ngraphs)d) 
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=%(nwavefuncs)d, NCOLOR=%(ncolor)d) 
      INTEGER NAMPSO, NSQAMPSO
      PARAMETER (NAMPSO=%(nAmpSplitOrders)d, NSQAMPSO=%(nSqAmpSplitOrders)d)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
	  REAL*8 RES(NSQAMPSO)
## if (color_correlation) {
      LOGICAL DO_INCLUDE_HEL_CONTRIB
## }
## if (spin_correlation) {
      COMPLEX*16 WEXT(20,NEXTERNAL)
## }
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,M,N, SQSOIND
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS)
	  COMPLEX*16 JAMP(NCOLOR,NAMPSO)
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0d0, 0d0), (1d0, 0d0)/
C
C FUNCTION
C
      INTEGER %(proc_prefix)sSQSOINDEX
C  
C GLOBAL VARIABLES
C  
      include 'coupl.inc'
## if (spin_correlation) {
      include 'spin_correlations.inc'
## }
C  
C COLOR DATA
C  
%(color_data_lines)s
## if (color_correlation) {

C
C COLOR CORRELATORS DATA
C
C     CCI will be the general color correlated index in loops
C     and CCIA is the index of the loop over COLOR_CORRELATORS_TO_CONSIDER
INTEGER CCI, CCIA
 
INTEGER NCOLORCORRELATORS
PARAMETER (NCOLORCORRELATORS=%(n_color_correlators)d)

COMPLEX*16 CCZTEMP(NCOLORCORRELATORS)

C Index 0 is the number of correlators to consider and the next indices are which one to consider
INTEGER COLOR_CORRELATORS_TO_CONSIDER(0:NCOLORCORRELATORS)
REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS, 0:NSQAMPSO)
COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER, COLOR_CORRELATED_EVALS

INTEGER CC_DENOM(NCOLORCORRELATORS,NCOLOR), CC_CF(NCOLORCORRELATORS,NCOLOR,NCOLOR)
%(color_correlators_data_lines)s

## }
C ----------
C BEGIN CODE
C ----------

## if (spin_correlation) {
IF (N_SPIN_CORR_RUNS.gt.0) THEN
C First copy the external wavefunction with those provided
  DO I=1,NEXTERNAL
    DO J=1,20
     W(J,I) = WEXT(J,I)
    ENDDO
  ENDDO
ENDIF
## }

%(helas_calls)s
%(jamp_lines)s

      RES = 0.D0 
	  DO M = 1, NAMPSO
        DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
## if (color_correlation) {
     DO CCI=1,NCOLORCORRELATORS
       CCZTEMP(CCI) = (0.D0,0.D0)
     ENDDO
## }
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J,M)
## if (color_correlation) {
     IF (DO_INCLUDE_HEL_CONTRIB) THEN
       DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
         CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
         CCZTEMP(CCI) = CCZTEMP(CCI) + CC_CF(CCI,J,I)*JAMP(J,M)
       ENDDO 
     ENDIF
## }
          ENDDO
		  DO N = 1, NAMPSO
		    SQSOIND = %(proc_prefix)sSQSOINDEX(M,N)
            RES(SQSOIND) = RES(SQSOIND) + ZTEMP*DCONJG(JAMP(I,N))/DENOM(I)
## if (color_correlation) {
     IF (DO_INCLUDE_HEL_CONTRIB) THEN
       DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
         CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
         COLOR_CORRELATED_EVALS(CCI,SQSOIND) = COLOR_CORRELATED_EVALS(CCI,SQSOIND) + CCZTEMP(CCI)*DCONJG(JAMP(I,N))/CC_DENOM(CCI,I)
       ENDDO
     ENDIF 
## }
		  ENDDO
        ENDDO
	  ENDDO

      END

      SUBROUTINE %(proc_prefix)sGET_ME(P, ALPHAS, NHEL ,ANS)
      IMPLICIT NONE   
C
C CONSTANT
C
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),ANS
      INTEGER NHEL
      DOUBLE PRECISION ALPHAS 
      real*8 pi     
CF2PY INTENT(OUT) :: ANS  
CF2PY INTENT(IN) :: NHEL   
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL) 
CF2PY INTENT(IN) :: ALPHAS
C     ROUTINE FOR F2PY to read the benchmark point.    
C     the include file with the values of the parameters and masses 
      include "coupl.inc"
      
      pi = 3.141592653589793d0
      G = 2* DSQRT(ALPHAS*pi)
      call UPDATE_AS_PARAM()
      if (NHEL.ne.0) then
        CALL %(proc_prefix)sSMATRIXHEL(P, NHEL, ANS)
      else
        CALL %(proc_prefix)sSMATRIX(P, ANS)
      endif
      return 
      end

      SUBROUTINE %(proc_prefix)sINITIALISE(PATH)
C     ROUTINE FOR F2PY to read the benchmark point.    
      IMPLICIT NONE   
      CHARACTER*512 PATH
CF2PY INTENT(IN) :: PATH
C     USE SETPARA2 and not SETPARA so that ident_card.dat can be automatically found
      call setpara2(PATH)  !first call to setup the paramaters    
      return 
      end

      LOGICAL FUNCTION %(proc_prefix)sIS_BORN_HEL_SELECTED(HELID)
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NCOMB
      PARAMETER (NCOMB=%(ncomb)d)
C
C     ARGUMENTS
C
      INTEGER HELID
C
C     LOCALS
C
      INTEGER I,J
      LOGICAL FOUNDIT
C
C     GLOBALS
C
      INTEGER HELC(NEXTERNAL,NCOMB)
      COMMON/%(proc_prefix)sBORN_HEL_CONFIGS/HELC

      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/%(proc_prefix)sBORN_BEAM_POL/POLARIZATIONS
C     ----------
C     BEGIN CODE
C     ----------
      
      %(proc_prefix)sIS_BORN_HEL_SELECTED = .True.
      if (POLARIZATIONS(0,0).eq.-1) THEN
        RETURN
      ENDIF

      DO I=1,NEXTERNAL
        IF (POLARIZATIONS(I,0).eq.-1) THEN
          CYCLE
        ENDIF
        FOUNDIT = .FALSE.
        DO J=1,POLARIZATIONS(I,0)
          IF (HELC(I,HELID).eq.POLARIZATIONS(I,J)) THEN
            FOUNDIT = .True.
            EXIT
          ENDIF
        ENDDO
        IF(.not.FOUNDIT) THEN
          %(proc_prefix)sIS_BORN_HEL_SELECTED = .False.
          RETURN
        ENDIF
      ENDDO

      RETURN
      END
      
## if (color_correlation) {

C --------------------------------------------------------------
C Beginning of a series of subroutine related to the evaluation
C and steering of the computation of color-correlated MEs.
C --------------------------------------------------------------

      SUBROUTINE %(proc_prefix)sGET_COLOR_CORRELATED_ME(CORRELATED_MES)
        IMPLICIT NONE
C       Provide access to the color correlated evaluations
		INTEGER CCI, I, J
        INTEGER    NEXTERNAL
        PARAMETER (NEXTERNAL=%(nexternal)d)
        INTEGER NAMPSO, NSQAMPSO
        PARAMETER (NAMPSO=%(nAmpSplitOrders)d, NSQAMPSO=%(nSqAmpSplitOrders)d)
        INTEGER NCOLORCORRELATORS
        PARAMETER (NCOLORCORRELATORS=%(n_color_correlators)d)

		REAL*8 CORRELATED_MES(NCOLORCORRELATORS, 0:NSQAMPSO)
CF2PY INTENT(OUT) :: CORRELATED_MES
  
C Index 0 is the number of correlators to consider and the next indices are which one to consider
        INTEGER COLOR_CORRELATORS_TO_CONSIDER(0:NCOLORCORRELATORS)        
        REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS, 0:NSQAMPSO)                
        COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER, COLOR_CORRELATED_EVALS
        
        DO CCI=1,NCOLORCORRELATORS
          DO J=0,NSQAMPSO
            CORRELATED_MES(CCI, J) = COLOR_CORRELATED_EVALS(CCI, J)
          ENDDO
        ENDDO

      END
      
      SUBROUTINE %(proc_prefix)sGET_COLOR_CORRELATOR_ID(CORRELATED_LEG_A, CORRELATED_LEG_B, ID)
        IMPLICIT NONE
        
        INTEGER    NEXTERNAL
        PARAMETER (NEXTERNAL=%(nexternal)d)
        INTEGER NCOLORCORRELATORS
        PARAMETER (NCOLORCORRELATORS=%(n_color_correlators)d)
       
        INTEGER CORRELATED_LEG_A, CORRELATED_LEG_B, ID, IDX
CF2PY INTENT(IN) :: CORRELATED_LEG_A
CF2PY INTENT(IN) :: CORRELATED_LEG_B
CF2PY INTENT(OUT) :: ID

        INTEGER COLOR_CORRELATOR_TO_INDEX(NEXTERNAL,NEXTERNAL)
        INTEGER INDEX_TO_COLOR_CORRELATOR(NCOLORCORRELATORS, 2)
        COMMON/%(proc_prefix)sCOLOR_CORRELATION_MAPS/COLOR_CORRELATOR_TO_INDEX, INDEX_TO_COLOR_CORRELATOR

		IDX = COLOR_CORRELATOR_TO_INDEX(CORRELATED_LEG_A, CORRELATED_LEG_B)
		
		IF (IDX.LE.0) THEN
		  ID = -1
		ELSE
		  ID = IDX
		ENDIF

      END
      
      SUBROUTINE %(proc_prefix)sGET_COLOR_CORRELATOR_FOR_ID(IDX, CORRELATED_LEG_A, CORRELATED_LEG_B)
        IMPLICIT NONE
        
        INTEGER    NEXTERNAL
        PARAMETER (NEXTERNAL=%(nexternal)d)
        INTEGER NCOLORCORRELATORS
        PARAMETER (NCOLORCORRELATORS=%(n_color_correlators)d)

        INTEGER CORRELATED_LEG_A, CORRELATED_LEG_B, IDX
CF2PY INTENT(IN) :: IDX
CF2PY INTENT(OUT) :: CORRELATED_LEG_A
CF2PY INTENT(OUT) :: CORRELATED_LEG_B

        INTEGER COLOR_CORRELATOR_TO_INDEX(NEXTERNAL,NEXTERNAL)
        INTEGER INDEX_TO_COLOR_CORRELATOR(NCOLORCORRELATORS, 2)
        COMMON/%(proc_prefix)sCOLOR_CORRELATION_MAPS/COLOR_CORRELATOR_TO_INDEX, INDEX_TO_COLOR_CORRELATOR

		CORRELATED_LEG_A = INDEX_TO_COLOR_CORRELATOR(IDX, 1)
		CORRELATED_LEG_B = INDEX_TO_COLOR_CORRELATOR(IDX, 2)

      END
      
      SUBROUTINE %(proc_prefix)sGET_N_COLOR_CORRELATORS(N_COLOR_CORRELATORS)
C       Provide access to the number of color correlators

		INTEGER N_COLOR_CORRELATORS
CF2PY INTENT(OUT) :: N_COLOR_CORRELATORS  
		
        INTEGER NCOLORCORRELATORS
        PARAMETER (NCOLORCORRELATORS=%(n_color_correlators)d)
        
        N_COLOR_CORRELATORS = NCOLORCORRELATORS

      END
      
      SUBROUTINE %(proc_prefix)sSET_COLOR_CORRELATORS_TO_CONSIDER(FIRST_LEG_INDEX, SECOND_LEG_INDEX)
        IMPLICIT NONE
C       Sets which color correlators to compute. If both leg_index are specified, this will select only one.
C       For each leg index set to -1, we will span the entire set of correlators for that leg index.
		
		INTEGER FIRST_LEG_INDEX, SECOND_LEG_INDEX
CF2PY INTENT(IN) :: FIRST_LEG_INDEX
CF2PY INTENT(IN) :: SECOND_LEG_INDEX

        INTEGER    NEXTERNAL
        PARAMETER (NEXTERNAL=%(nexternal)d)
        INTEGER NAMPSO, NSQAMPSO
        PARAMETER (NAMPSO=%(nAmpSplitOrders)d, NSQAMPSO=%(nSqAmpSplitOrders)d)
        INTEGER NCOLORCORRELATORS
        PARAMETER (NCOLORCORRELATORS=%(n_color_correlators)d)
        
		INTEGER I, J, K, NCC_TO_CONSIDER, CC_TO_ADD
        LOGICAL FOUNDIT

        INTEGER COLOR_CORRELATOR_TO_INDEX(NEXTERNAL,NEXTERNAL)
%(color_correlator_to_index_data)s
        
        INTEGER INDEX_TO_COLOR_CORRELATOR(NCOLORCORRELATORS, 2)
%(index_to_color_correlator_data)s
        
        COMMON/%(proc_prefix)sCOLOR_CORRELATION_MAPS/COLOR_CORRELATOR_TO_INDEX, INDEX_TO_COLOR_CORRELATOR
        
C       Index 0 is the number of correlators to consider and the next indices are which one to consider
		INTEGER COLOR_CORRELATORS_TO_CONSIDER(0:NCOLORCORRELATORS)   
        REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS, 0:NSQAMPSO)                
        COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER, COLOR_CORRELATED_EVALS
        
        NCC_TO_CONSIDER = 0
        
C       First reset the array
	 	COLOR_CORRELATORS_TO_CONSIDER(0) = 0
        DO I=1,NCOLORCORRELATORS
          COLOR_CORRELATORS_TO_CONSIDER(I) = -1
        ENDDO
        
C       Index 0 on the first leg simply indicate the desire to reset the selection
        IF (FIRST_LEG_INDEX.eq.0) THEN
          RETURN
        ENDIF
        
        DO I=1, NEXTERNAL
          IF (FIRST_LEG_INDEX.ne.-1.and.FIRST_LEG_INDEX.ne.I) THEN
            CYCLE
          ENDIF
          DO J=1, NEXTERNAL
            IF (SECOND_LEG_INDEX.ne.-1.and.SECOND_LEG_INDEX.ne.J) THEN
              CYCLE
            ENDIF
C           If one of those legs is colorless and does not have a correlator, skip it.
            CC_TO_ADD = COLOR_CORRELATOR_TO_INDEX(I,J)
            IF (CC_TO_ADD.le.0) THEN
              CYCLE
            ENDIF
C           Now make sure we have not considered it yet, since (i, j) maps to (j, i) in these correlators
            FOUNDIT = .FALSE.
            DO K=1,COLOR_CORRELATORS_TO_CONSIDER(0)
              IF (COLOR_CORRELATORS_TO_CONSIDER(K).EQ.CC_TO_ADD) THEN
                FOUNDIT = .TRUE.
                EXIT
              ENDIF
            ENDDO
            IF (FOUNDIT) THEN
              CYCLE
            ENDIF
            NCC_TO_CONSIDER = NCC_TO_CONSIDER+1
            COLOR_CORRELATORS_TO_CONSIDER(0) = NCC_TO_CONSIDER
            COLOR_CORRELATORS_TO_CONSIDER(NCC_TO_CONSIDER) = CC_TO_ADD
          ENDDO
        ENDDO

      END
## }

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
        DATA N_SPIN_CORR_VECTORS/NEXTERNAL*0/
        INTEGER SPIN_CORR_RUNS_N_ENTRIES
        PARAMETER (SPIN_CORR_RUNS_N_ENTRIES=(MAX_SPIN_CORR_RUNS+1)*NEXTERNAL)  
C       The integer -1000 will indicate to simply keep the original helicity polarization vector assigned to that leg      
        DATA SPIN_CORR_RUNS/SPIN_CORR_RUNS_N_ENTRIES*-1000/
        DATA N_SPIN_CORR_RUNS/0/

C       ----------
C       BEGIN CODE
C       ----------

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
            SPIN_CORR_VECTORS(LEG_INDEX,I,J) = DBLE(VECTORS(I,J))
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

C       ----------
C       BEGIN CODE
C       ----------

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
CF2PY INTENT(IN) :: MAX_N_SC_VECTORS
		MAX_N_SC_VECTORS = MAX_N_SPIN_CORR_VECTORS
      END SUBROUTINE %(proc_prefix)sGET_MAX_N_SPIN_CORR_VECTORS
## }