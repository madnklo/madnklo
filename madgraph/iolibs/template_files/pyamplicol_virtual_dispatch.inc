      SUBROUTINE %(proc_prefix)sSLOOPMATRIX(P,ANS)
      IMPLICIT NONE
      INTEGER NEXTERNAL,NSQUAREDSO,USERHEL,RET_CODE,H,T,U
      PARAMETER(NEXTERNAL=%(pyamp_n)d,NSQUAREDSO=%(nSquaredSO)d)
      REAL*8 P(0:3,NEXTERNAL),ANS(0:3,0:NSQUAREDSO)
      REAL*8 USER_STAB_PREC,ACCURACY(0:NSQUAREDSO)
      DATA USERHEL/-1/
      DATA USER_STAB_PREC/-1.0D0/
      COMMON/%(proc_prefix)sUSERCHOICE/USERHEL
      COMMON/%(proc_prefix)sUSER_STAB_PREC/USER_STAB_PREC
      COMMON/%(proc_prefix)sACC/ACCURACY,H,T,U
      CALL %(proc_prefix)sPYAMP_VIRTUAL(P,USERHEL,ANS,USER_STAB_PREC,ACCURACY,RET_CODE)
      H=RET_CODE/100
      T=MOD(RET_CODE,100)/10
      U=MOD(RET_CODE,10)
      END

      SUBROUTINE %(proc_prefix)sGET_PYAMP_PROVIDER(ENABLED)
      INTEGER ENABLED
CF2PY INTENT(OUT) :: ENABLED
      ENABLED=1
      END

      SUBROUTINE %(proc_prefix)sPYAMP_VIRTUAL(P,HEL,ANS,PREC_ASKED,PREC_FOUND,RET_CODE)
      IMPLICIT NONE
      INTEGER HEL,RET_CODE
      REAL*8 P(0:3,%(pyamp_n)d),ANS(0:3,0:%(nSquaredSO)d)
      REAL*8 PREC_ASKED,PREC_FOUND(0:%(nSquaredSO)d)
      REAL*8 CC(MAX(1,%(pyamp_nc)d),0:3,0:%(nSquaredSO)d)
      CALL %(proc_prefix)sPYAMP_VIRTUAL_SELECTED(P,HEL,ANS,PREC_ASKED,PREC_FOUND,RET_CODE,CC)
      END

## if (color_correlation) {
      SUBROUTINE %(proc_prefix)sGET_VIRTUAL_WITH_COLOR(P,HEL,ALPHAS,MUR,PREC,ANS,FOUND,STATUS,CC)
      IMPLICIT NONE
      INTEGER HEL,STATUS
      REAL*8 P(0:3,%(pyamp_n)d),ALPHAS,MUR,PREC
      REAL*8 ANS(0:3,0:%(nSquaredSO)d),FOUND(0:%(nSquaredSO)d)
      REAL*8 CC(%(pyamp_nc)d,0:3,0:%(nSquaredSO)d)
CF2PY INTENT(IN) :: P, HEL, ALPHAS, MUR, PREC
CF2PY INTENT(OUT) :: ANS, FOUND, STATUS, CC
      CALL %(proc_prefix)sPYAMP_UPDATE_VIRTUAL(ALPHAS,MUR)
      CALL %(proc_prefix)sPYAMP_VIRTUAL_SELECTED(P,HEL,ANS,PREC,FOUND,STATUS,CC)
      END
## }

      SUBROUTINE %(proc_prefix)sPYAMP_VIRTUAL_SELECTED(P,HEL,ANS,PREC_ASKED,PREC_FOUND,RET_CODE,CC)
      IMPLICIT NONE
      INTEGER NEXTERNAL,HEL,RET_CODE
      PARAMETER(NEXTERNAL=%(pyamp_n)d)
      REAL*8 P(0:3,NEXTERNAL),ANS(0:3,0:%(nSquaredSO)d)
      REAL*8 PREC_ASKED,PREC_FOUND(0:%(nSquaredSO)d)
      REAL*8 CC(MAX(1,%(pyamp_nc)d),0:3,0:%(nSquaredSO)d)
      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/%(proc_prefix)sBEAM_POL/POLARIZATIONS
## if (spin_correlation) {
      INCLUDE '%(proc_prefix)sspin_correlations.inc'
## }
## if (not spin_correlation) {
      INTEGER MAX_N_SPIN_CORR_VECTORS,N_SPIN_CORR_VECTORS(NEXTERNAL)
      PARAMETER(MAX_N_SPIN_CORR_VECTORS=1)
      COMPLEX*16 SPIN_CORR_VECTORS(NEXTERNAL,1,4)
      DATA N_SPIN_CORR_VECTORS/NEXTERNAL*0/
## }
      IF(POLARIZATIONS(0,0).NE.-1) THEN
         WRITE(*,*) 'pyAmpliCol: selected beam polarizations are unsupported'
         STOP 1
      ENDIF
      CALL %(proc_prefix)sPYAMP_VIRTUAL_EVAL(P,HEL,ANS,PREC_ASKED,PREC_FOUND,RET_CODE,MAX_N_SPIN_CORR_VECTORS,N_SPIN_CORR_VECTORS,SPIN_CORR_VECTORS,CC)
      END

      SUBROUTINE %(proc_prefix)sPYAMP_UPDATE_VIRTUAL(ALPHAS,MUR)
      IMPLICIT NONE
      REAL*8 ALPHAS,MUR
      IF (.NOT.(ALPHAS.GT.0D0.AND.ALPHAS.LT.HUGE(1D0).AND.MUR.GT.0D0.AND.MUR.LT.HUGE(1D0))) THEN
         WRITE(*,*) 'pyAmpliCol: alpha_s and mu_r must be positive and finite'
         STOP 1
      ENDIF
      CALL UPDATE_AS_PARAM2(MUR,ALPHAS)
      END

      SUBROUTINE %(proc_prefix)sPYAMP_VIRTUAL_EVAL(P,HEL,ANS,PREC_ASKED,PREC_FOUND,RET_CODE,MAX_VECTORS,N_VECTORS,VECTORS,CC)
      IMPLICIT NONE
      INTEGER NEXTERNAL,NSQUAREDSO,MNK_NP,MNK_NI
      PARAMETER(NEXTERNAL=%(pyamp_n)d,NSQUAREDSO=%(nSquaredSO)d)
      PARAMETER(MNK_NP=%(pyamp_np)d,MNK_NI=%(pyamp_ni)d)
      INTEGER MNK_PDGS(NEXTERNAL),HEL,RET_CODE
      DATA MNK_PDGS/%(pyamp_pdgs)s/
      REAL*8 P(0:3,NEXTERNAL),ANS(0:3,0:NSQUAREDSO)
      REAL*8 PREC_ASKED,PREC_FOUND(0:NSQUAREDSO),MNK_PARAMS(MNK_NP)
      INTEGER NC
      PARAMETER(NC=%(pyamp_nc)d)
      INTEGER MAX_VECTORS,N_VECTORS(NEXTERNAL)
      COMPLEX*16 VECTORS(NEXTERNAL,MAX_VECTORS,4)
      REAL*8 CC(MAX(1,NC),0:3,0:NSQUAREDSO)
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
%(pyamp_parameters)s
      CALL MNK_PYAMP_VIRTUAL_CORRELATED_EVAL(NEXTERNAL,MNK_NI,MNK_PDGS,P,HEL,NSQUAREDSO,MU_R,PREC_ASKED,MNK_NP,MNK_PARAMS,MAX_VECTORS,N_VECTORS,VECTORS,NC,ANS,CC,PREC_FOUND,RET_CODE)
      END

      SUBROUTINE %(proc_prefix)sPYAMP_VIRTUAL_POINT(P,ALPHAS,MUR,LEG,KT,ANS,CC,STATUS)
      IMPLICIT NONE
      INTEGER N,LEG,STATUS
      PARAMETER(N=%(pyamp_n)d)
      INTEGER COUNTS(N)
      COMPLEX*16 VECTORS(N,1,4)
      REAL*8 P(0:3,N),ALPHAS,MUR,KT(0:3),ANS(0:3,0:%(nSquaredSO)d)
      REAL*8 CC(MAX(1,%(pyamp_nc)d),0:3,0:%(nSquaredSO)d)
      REAL*8 FOUND(0:%(nSquaredSO)d)
      CALL %(proc_prefix)sPYAMP_UPDATE_VIRTUAL(ALPHAS,MUR)
      COUNTS=0
      VECTORS=0D0
C A zero leg requests a spin sum, regardless of the legacy selectors.
      IF(LEG.NE.0) THEN
         IF(LEG.LT.1.OR.LEG.GT.N) THEN
            WRITE(*,*) 'pyAmpliCol: invalid virtual spin leg'
            STOP 1
         ENDIF
         COUNTS(LEG)=1
         VECTORS(LEG,1,:)=DCMPLX(KT,0D0)
      ENDIF
      CALL %(proc_prefix)sPYAMP_VIRTUAL_EVAL(P,-1,ANS,-1D0,FOUND,STATUS,1,COUNTS,VECTORS,CC)
      END

## if (spin_correlation) {
      SUBROUTINE %(proc_prefix)sGET_KKVNLO(P,ALPHAS,MUR,LEG,KT,RESULT)
      IMPLICIT NONE
      INTEGER LEG,STATUS
      REAL*8 P(0:3,%(pyamp_n)d),ALPHAS,MUR,KT(0:3),RESULT(-2:0)
      REAL*8 ANS(0:3,0:%(nSquaredSO)d),CC(MAX(1,%(pyamp_nc)d),0:3,0:%(nSquaredSO)d)
CF2PY INTENT(IN) :: P, ALPHAS, MUR, LEG, KT
CF2PY INTENT(OUT) :: RESULT
      IF(LEG.EQ.0) THEN
         WRITE(*,*) 'pyAmpliCol: invalid virtual spin leg'
         STOP 1
      ENDIF
      CALL %(proc_prefix)sPYAMP_VIRTUAL_POINT(P,ALPHAS,MUR,LEG,KT,ANS,CC,STATUS)
      CALL %(proc_prefix)sPYAMP_CHECK_VIRTUAL(STATUS)
      RESULT=ANS(3:1:-1,0)
      END
## }

## if (color_correlation) {
      SUBROUTINE %(proc_prefix)sGET_CCVNLO(P,ALPHAS,MUR,A,B,RESULT)
      IMPLICIT NONE
      INTEGER A,B,I,STATUS
      INTEGER DIPOLES(2,%(pyamp_nc)d)
      DATA DIPOLES/%(color_dipoles)s/
      REAL*8 P(0:3,%(pyamp_n)d),ALPHAS,MUR,RESULT(-2:0),KT(0:3)
      REAL*8 ANS(0:3,0:%(nSquaredSO)d),CC(%(pyamp_nc)d,0:3,0:%(nSquaredSO)d)
CF2PY INTENT(IN) :: P, ALPHAS, MUR, A, B
CF2PY INTENT(OUT) :: RESULT
      KT=0D0
      CALL %(proc_prefix)sPYAMP_VIRTUAL_POINT(P,ALPHAS,MUR,0,KT,ANS,CC,STATUS)
      CALL %(proc_prefix)sPYAMP_CHECK_VIRTUAL(STATUS)
      DO I=1,%(pyamp_nc)d
         IF(A.NE.DIPOLES(1,I).OR.B.NE.DIPOLES(2,I)) CYCLE
         RESULT=CC(I,3:1:-1,1)
         RETURN
      ENDDO
      WRITE(*,*) 'pyAmpliCol: unknown virtual colour dipole',A,B
      STOP 1
      END

C Return the ordinary and all dipole coefficients from one spin-summed point.
      SUBROUTINE %(proc_prefix)sGET_VIRTUAL_DIPOLES(P,ALPHAS,MUR,V,VC,STATUS)
      IMPLICIT NONE
      INTEGER N,I,A,B,STATUS
      PARAMETER(N=%(pyamp_n)d)
      INTEGER DIPOLES(2,%(pyamp_nc)d)
      DATA DIPOLES/%(color_dipoles)s/
      REAL*8 P(0:3,N),ALPHAS,MUR,V(-2:0),VC(-2:0,N,N),KT(0:3)
      REAL*8 ANS(0:3,0:%(nSquaredSO)d),CC(%(pyamp_nc)d,0:3,0:%(nSquaredSO)d)
      KT=0D0
      CALL %(proc_prefix)sPYAMP_VIRTUAL_POINT(P,ALPHAS,MUR,0,KT,ANS,CC,STATUS)
      V=ANS(3:1:-1,0)
      VC=HUGE(1D0)
      DO I=1,%(pyamp_nc)d
         A=DIPOLES(1,I)
         B=DIPOLES(2,I)
         IF(MIN(A,B).LT.1.OR.MAX(A,B).GT.N) CYCLE
         VC(:,A,B)=CC(I,3:1:-1,1)
      ENDDO
      END
## }

      SUBROUTINE %(proc_prefix)sPYAMP_CHECK_VIRTUAL(STATUS)
      IMPLICIT NONE
      INTEGER STATUS
      IF(STATUS/100.LT.1.OR.STATUS/100.GT.3) THEN
         WRITE(*,*) 'pyAmpliCol: virtual evaluation failed',STATUS
         STOP 1
      ENDIF
      END
