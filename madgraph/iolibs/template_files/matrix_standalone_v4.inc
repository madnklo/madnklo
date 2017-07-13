      SUBROUTINE %(proc_prefix)sSMATRIXHEL(P,HEL,ANS)
      IMPLICIT NONE
C
C CONSTANT
C
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=%(ncomb)d)
CF2PY INTENT(OUT) :: ANS                                                                                                                                                                                     
CF2PY INTENT(IN) :: HEL                                                                                                                                                                                    
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL)  

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

      SUBROUTINE %(proc_prefix)sSMATRIX(P,ANS)
C  
%(info_lines)s
C 
C MadGraph5_aMC@NLO StandAlone Version
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
	  INTEGER HELAVGFACTOR
	  PARAMETER (HELAVGFACTOR=%(hel_avg_factor)d)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),ANS
CF2PY INTENT(OUT) :: ANS
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      LOGICAL DO_INCLUDE_HEL_CONTRIB
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 %(proc_prefix)sMATRIX
      INTEGER IHEL,IDEN, I, J
C     For a 1>N process, them BEAMTWO_HELAVGFACTOR would be set to 1.
      INTEGER BEAMS_HELAVGFACTOR(2)
	  DATA (BEAMS_HELAVGFACTOR(I),I=1,2)/%(beamone_helavgfactor)d,%(beamtwo_helavgfactor)d/
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./

C
C GLOBAL VARIABLES
C
      INTEGER USERHEL
	  COMMON/%(proc_prefix)sHELUSERCHOICE/USERHEL
	  DATA USERHEL/-1/	  

%(helicity_lines)s
%(den_factor_line)s

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
      
C     Index 0 is the number of correlators to consider and the next indices are which one to consider
	  INTEGER COLOR_CORRELATORS_TO_CONSIDER(0:NCOLORCORRELATORS)   
      DATA COLOR_CORRELATORS_TO_CONSIDER(0)/NCOLORCORRELATORS/
%(color_correlators_to_consider_initialization)s

      REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS)
      DATA (COLOR_CORRELATED_EVALS(I), I=1,NCOLORCORRELATORS) / NCOLORCORRELATORS*0.0d0 /       
      COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER, COLOR_CORRELATED_EVALS
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
       COLOR_CORRELATED_EVALS(CCI) = 0.0d0
     ENDDO
## }

      IF(USERHEL.EQ.-1) NTRY=NTRY+1
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
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
      ANS = 0D0
          DO IHEL=1,NCOMB
	         IF (USERHEL.EQ.-1.OR.USERHEL.EQ.IHEL) THEN		  
                IF (GOODHEL(IHEL) .OR. NTRY .LT. 20.OR.USERHEL.NE.-1) THEN
				    IF(NTRY.GE.2.AND.POLARIZATIONS(0,0).ne.-1.and.(.not.%(proc_prefix)sIS_BORN_HEL_SELECTED(IHEL))) THEN
                        CYCLE
                    ENDIF
                    DO_INCLUDE_HEL_CONTRIB = (POLARIZATIONS(0,0).eq.-1.or.%(proc_prefix)sIS_BORN_HEL_SELECTED(IHEL))
## if (not color_correlation) {             
                    T=%(proc_prefix)sMATRIX(P ,NHEL(1,IHEL),JC(1))
## } else {
                    T=%(proc_prefix)sMATRIX(P ,NHEL(1,IHEL),JC(1), DO_INCLUDE_HEL_CONTRIB)
## }
                    IF(DO_INCLUDE_HEL_CONTRIB) THEN
                        ANS=ANS+T
                    ENDIF
				    IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
                        GOODHEL(IHEL)=.TRUE.
                    ENDIF
                ENDIF
		     ENDIF
          ENDDO
      ANS=ANS/DBLE(IDEN)
## if (color_correlation) {
      DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
         CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
         COLOR_CORRELATED_EVALS(CCI) = COLOR_CORRELATED_EVALS(CCI)/DBLE(IDEN)
      ENDDO
## }
	  IF(USERHEL.NE.-1) THEN
          ANS=ANS*HELAVGFACTOR
## if (color_correlation) {
          DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
            CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
            COLOR_CORRELATED_EVALS(CCI) = COLOR_CORRELATED_EVALS(CCI)*HELAVGFACTOR
          ENDDO
## }
	  ELSE
		DO J=1,NINITIAL
		  IF (POLARIZATIONS(J,0).ne.-1) THEN
            ANS=ANS*BEAMS_HELAVGFACTOR(J)
            ANS=ANS/POLARIZATIONS(J,0)
## if (color_correlation) {
          DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
            CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
            COLOR_CORRELATED_EVALS(CCI) = COLOR_CORRELATED_EVALS(CCI)*BEAMS_HELAVGFACTOR(J)
            COLOR_CORRELATED_EVALS(CCI) = COLOR_CORRELATED_EVALS(CCI)/POLARIZATIONS(J,0)
          ENDDO
## }
		  ENDIF
		ENDDO
	  ENDIF
	  END
       
## if (not color_correlation) {                 
      REAL*8 FUNCTION %(proc_prefix)sMATRIX(P,NHEL,IC)
## } else {
      REAL*8 FUNCTION %(proc_prefix)sMATRIX(P,NHEL,IC, DO_INCLUDE_HEL_CONTRIB)
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
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
## if (color_correlation) {
      LOGICAL DO_INCLUDE_HEL_CONTRIB
## }
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0d0, 0d0), (1d0, 0d0)/
C  
C GLOBAL VARIABLES
C  
      include 'coupl.inc'
%(global_variable)s
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
REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS)
COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER, COLOR_CORRELATED_EVALS

INTEGER CC_DENOM(NCOLORCORRELATORS,NCOLOR), CC_CF(NCOLORCORRELATORS,NCOLOR,NCOLOR)
%(color_correlators_data_lines)s

## }
C ----------
C BEGIN CODE
C ----------
%(helas_calls)s
%(jamp_lines)s

      %(proc_prefix)sMATRIX = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
## if (color_correlation) {
     DO CCI=1,NCOLORCORRELATORS
       CCZTEMP(CCI) = (0.D0,0.D0)
     ENDDO
## }
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
## if (color_correlation) {
     IF (DO_INCLUDE_HEL_CONTRIB) THEN
       DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
         CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
         CCZTEMP(CCI) = CCZTEMP(CCI) + CC_CF(CCI,J,I)*JAMP(J)
       ENDDO
     ENDIF 
## }
          ENDDO
          %(proc_prefix)sMATRIX = %(proc_prefix)sMATRIX+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
## if (color_correlation) {
     IF (DO_INCLUDE_HEL_CONTRIB) THEN
       DO CCIA=1,COLOR_CORRELATORS_TO_CONSIDER(0)
         CCI = COLOR_CORRELATORS_TO_CONSIDER(CCIA)
         COLOR_CORRELATED_EVALS(CCI) = COLOR_CORRELATED_EVALS(CCI) + CCZTEMP(CCI)*DCONJG(JAMP(I))/CC_DENOM(CCI,I)
       ENDDO
     ENDIF 
## }      
      ENDDO
%(amp2_lines)s
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
C     USE SETPARA2 and not SETPARA so that ident_card.dat can be automatically found.
      CALL SETPARA2(path)
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
      SUBROUTINE %(proc_prefix)sGET_COLOR_CORRELATED_ME(CORRELATED_MES)
        IMPLICIT NONE
C       Provide access to the color correlated evaluations
		INTEGER CCI, I, J
        INTEGER    NEXTERNAL
        PARAMETER (NEXTERNAL=%(nexternal)d)
        INTEGER NCOLORCORRELATORS
        PARAMETER (NCOLORCORRELATORS=%(n_color_correlators)d)

		REAL*8 CORRELATED_MES(NCOLORCORRELATORS, 0:1)
CF2PY INTENT(OUT) :: CORRELATED_MES
  
C Index 0 is the number of correlators to consider and the next indices are which one to consider
        INTEGER COLOR_CORRELATORS_TO_CONSIDER(0:NCOLORCORRELATORS)        
        REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS)                
        COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER, COLOR_CORRELATED_EVALS
        
        DO CCI=1,NCOLORCORRELATORS
            CORRELATED_MES(CCI, 0) = COLOR_CORRELATED_EVALS(CCI)
            CORRELATED_MES(CCI, 1) = COLOR_CORRELATED_EVALS(CCI)            
        ENDDO

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
        REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS)                
        COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER, COLOR_CORRELATED_EVALS
        
        NCC_TO_CONSIDER = 0
        
C       First reset the array
	 	COLOR_CORRELATORS_TO_CONSIDER(0) = 0
        DO I=1,NCOLORCORRELATORS
          COLOR_CORRELATORS_TO_CONSIDER(I) = -1
        ENDDO
        
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