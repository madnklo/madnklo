      Subroutine Initialise(path)

      CHARACTER(512) path
CF2PY intent(in)::path
 
C     INCLUDE FILES
C     
C     the include file with the values of the parameters and masses   
C        
      INCLUDE 'coupl.inc'
C     USE SETPARA2 and not SETPARA so that ident_card.dat can be automatically found.
      CALL SETPARA2(path)
      return 
      end

      SUBROUTINE INITIALISE_MADLOOP_PATH(PATH)

      CHARACTER(512) PATH
CF2PY intent(in)::path
      CALL SETMADLOOPPATH(PATH)
      RETURN
      END

      SUBROUTINE GET_ME(P, ALPHAS, SCALE2, NHEL , ANS,RETURNCODE)
      IMPLICIT NONE
C     
C     CONSTANTS  
C     
      REAL*8 ZERO
      PARAMETER (ZERO=0D0)

C     integer nexternal C number particles (incoming+outgoing) in the
C      me 
      include 'nexternal.inc'

C      CHARACTER(512) MADLOOPRESOURCEPATH
C     
C     INCLUDE FILES
C     
C     the include file with the values of the parameters and masses   
C        
      INCLUDE 'coupl.inc'
C     particle masses
      REAL*8 PMASS(NEXTERNAL)
C     integer    n_max_cg
      INCLUDE 'ngraphs.inc'
      INCLUDE 'nsquaredSO.inc'

C     LOCAL
C     
      INTEGER I
C     four momenta. Energy is the zeroth component.
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER MATELEM_ARRAY_DIM
      REAL*8 , ALLOCATABLE :: MATELEM(:,:)
      INTEGER RETURNCODE
      INTEGER NSQUAREDSO_LOOP
      REAL*8 , ALLOCATABLE :: PREC_FOUND(:)

      DOUBLE PRECISION ANS
      INTEGER NHEL
      DOUBLE PRECISION ALPHAS, SCALE2
CF2PY INTENT(OUT) :: ANS
CF2PY INTENT(OUT) :: RETURNCODE
CF2PY INTENT(IN) :: NHEL   
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL) 
CF2PY INTENT(IN) :: ALPHAS
CF2PY INTENT(IN) :: SCALE2


C     
C     GLOBAL VARIABLES
C     
C     This is from ML code for the list of split orders selected by
C     the process definition
C
      LOGICAL CHOSEN_LOOP_SO_CONFIGS(NSQUAREDSO)
      COMMON/%(proc_prefix)sCHOSEN_LOOP_SQSO/CHOSEN_LOOP_SO_CONFIGS
C     
C     BEGIN CODE
C     
      CALL %(proc_prefix)sFORCE_STABILITY_CHECK(.True.)
      CALL %(proc_prefix)sGET_ANSWER_DIMENSION(MATELEM_ARRAY_DIM)
      ALLOCATE(MATELEM(0:3,0:MATELEM_ARRAY_DIM))
      CALL %(proc_prefix)sGET_NSQSO_LOOP(NSQUAREDSO_LOOP)
      ALLOCATE(PREC_FOUND(0:NSQUAREDSO_LOOP))
      INCLUDE 'pmass.inc'

C       Update the couplings with the new ALPHAS
        CALL UPDATE_AS_PARAM2(SCALE2, ALPHAS)

C       
C       Now we can call the matrix element
C       
        if (NHEL.eq.0) then
          CALL %(proc_prefix)sSLOOPMATRIX_THRES(P,MATELEM,-1.0D0, PREC_FOUND, RETURNCODE)
       else
       CALL %(proc_prefix)sSLOOPMATRIXHEL_THRES(P,NHEL, MATELEM,-1.0D0, PREC_FOUND, RETURNCODE)
       endif
       
c      loop induce -> only finite part 
       ANS =  MATELEM(1,0)

      END







