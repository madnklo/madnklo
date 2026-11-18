C Native provider wrappers keep the public selected-state Fortran interface.
## if (spin_correlation) {
C Two rank-one spin contractions, with no saved Born or selection state.
      DOUBLE PRECISION FUNCTION %(proc_prefix)sGET_KKKKBLO(P,ALPHAS,LEG_A,KT_A,LEG_B,KT_B)
      IMPLICIT NONE
      INTEGER LEG_A,LEG_B,LEGS(2)
      REAL*8 P(0:3,%(nexternal)d),ALPHAS,KT_A(0:3),KT_B(0:3),PAIR(2)
      COMPLEX*16 VECTOR(0:3,2)
%(pyamp_declarations)s
      IF(MIN(LEG_A,LEG_B).LT.1.OR.MAX(LEG_A,LEG_B).GT.MNK_N.OR.LEG_A.EQ.LEG_B) THEN
         WRITE(*,*) 'pyAmpliCol: invalid double-spin legs'
         STOP 1
      ENDIF
      LEGS=(/LEG_A,LEG_B/)
      VECTOR(:,1)=DCMPLX(KT_A,0D0)
      VECTOR(:,2)=DCMPLX(KT_B,0D0)
      CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,P,ALPHAS,'born',4,2,LEGS,VECTOR,PAIR)
      CALL %(proc_prefix)sPYAMP_REAL(PAIR,0D0,%(proc_prefix)sGET_KKKKBLO)
      END
## }

## if (color_correlation) {
C Symmetric colour insertion: d_ABC on LEG_A, ordinary T on LEG_B.
C Both public getters take a complete point and do not read saved Born state.
      DOUBLE PRECISION FUNCTION %(proc_prefix)sGET_DCCBLO(P,ALPHAS,LEG_A,LEG_B)
      IMPLICIT NONE
      INTEGER LEG_A,LEG_B,LEGS(1)
      REAL*8 P(0:3,%(nexternal)d),ALPHAS,PAIR(2),SCALE
      COMPLEX*16 VECTOR(0:3,1)
      LEGS=0
      VECTOR=DCMPLX(0D0,0D0)
      CALL %(proc_prefix)sPYAMP_SYMMETRIC(P,ALPHAS,LEG_A,LEG_B,0,LEGS,VECTOR,PAIR,SCALE)
      CALL %(proc_prefix)sPYAMP_REAL(PAIR,SCALE,%(proc_prefix)sGET_DCCBLO)
      END

      SUBROUTINE %(proc_prefix)sPYAMP_SYMMETRIC(P,ALPHAS,LEG_A,LEG_B,NSPIN,LEGS,VECTOR,PAIR,SCALE)
      IMPLICIT NONE
      INTEGER LEG_A,LEG_B,NSPIN,LEGS(MAX(1,NSPIN))
      REAL*8 P(0:3,%(nexternal)d),ALPHAS,PAIR(2),SCALE,QUARK(2),ANTIQUARK(2)
      COMPLEX*16 VECTOR(0:3,MAX(1,NSPIN))
      CHARACTER*64 ID
%(pyamp_declarations)s
      IF(MIN(LEG_A,LEG_B).LT.1.OR.MAX(LEG_A,LEG_B).GT.MNK_N) THEN
         WRITE(*,*) 'pyAmpliCol: invalid symmetric colour legs'
         STOP 1
      ENDIF
      WRITE(ID,'(A,I0,A,I0)') 'Dq_',LEG_A,'_',LEG_B
      CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,P,ALPHAS,ID,LEN_TRIM(ID),NSPIN,LEGS,VECTOR,QUARK)
      WRITE(ID,'(A,I0,A,I0)') 'Dqb_',LEG_A,'_',LEG_B
      CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,P,ALPHAS,ID,LEN_TRIM(ID),NSPIN,LEGS,VECTOR,ANTIQUARK)
C S^dagger (T_q - T_qbar) S = TR*d, with TR=1/2. Keep complex
C components until the final difference; SCALE belongs to this point only.
      PAIR=2D0*(QUARK-ANTIQUARK)
      SCALE=2D0*(ABS(QUARK(1))+ABS(ANTIQUARK(1)))
      IF(ANY(PAIR.NE.PAIR).OR.ANY(ABS(PAIR).GT.HUGE(1D0)).OR.SCALE.GT.HUGE(1D0)) THEN
         WRITE(*,*) 'pyAmpliCol: nonfinite symmetric colour result'
         STOP 1
      ENDIF
      END
## }

## if (color_correlation and spin_correlation) {
      DOUBLE PRECISION FUNCTION %(proc_prefix)sGET_KKDCCBLO(P,ALPHAS,PARENT_LEG,KT,LEG_A,LEG_B)
      IMPLICIT NONE
      INTEGER PARENT_LEG,LEG_A,LEG_B,LEGS(1)
      REAL*8 P(0:3,%(nexternal)d),ALPHAS,KT(0:3),PAIR(2),SCALE
      COMPLEX*16 VECTOR(0:3,1)
      LEGS(1)=PARENT_LEG
      VECTOR(:,1)=DCMPLX(KT,0D0)
      CALL %(proc_prefix)sPYAMP_SYMMETRIC(P,ALPHAS,LEG_A,LEG_B,1,LEGS,VECTOR,PAIR,SCALE)
      CALL %(proc_prefix)sPYAMP_REAL(PAIR,SCALE,%(proc_prefix)sGET_KKDCCBLO)
      END

C Explicit joint spin-dipole contraction; no Born hook or setters required.
      DOUBLE PRECISION FUNCTION %(proc_prefix)sGET_KKCCBLO(P,ALPHAS,PARENT_LEG,KT,LEG_A,LEG_B)
      IMPLICIT NONE
      INTEGER PARENT_LEG,LEG_A,LEG_B,LEGS(1)
      REAL*8 P(0:3,%(nexternal)d),ALPHAS,KT(0:3),PAIR(2)
      COMPLEX*16 VECTOR(0:3,1)
      CHARACTER*64 ID
%(pyamp_declarations)s
      IF(MIN(LEG_A,LEG_B).LT.1.OR.MAX(LEG_A,LEG_B).GT.MNK_N) THEN
         WRITE(*,*) 'pyAmpliCol: invalid joint spin-colour legs'
         STOP 1
      ENDIF
      LEGS(1)=PARENT_LEG
      VECTOR(:,1)=DCMPLX(KT,0D0)
      WRITE(ID,'(A,I0,A,I0)') 'T_',LEG_A,'_',LEG_B
      CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,P,ALPHAS,ID,LEN_TRIM(ID),1,LEGS,VECTOR,PAIR)
      CALL %(proc_prefix)sPYAMP_REAL(PAIR,0D0,%(proc_prefix)sGET_KKCCBLO)
      END
## }

## if (pyamp_quadruple_color) {
      DOUBLE PRECISION FUNCTION %(proc_prefix)sGET_QUADBLO(P,ALPHAS,A,B,C,D)
      IMPLICIT NONE
      INTEGER A,B,C,D,ORDERING,LEGS(1),IA(2),IB(2)
      REAL*8 P(0:3,%(nexternal)d),ALPHAS,PAIR(2),TOTAL(2)
      COMPLEX*16 VECTOR(0:3,1)
      CHARACTER*64 ID
%(pyamp_declarations)s
      IF(MIN(A,B,C,D).LT.1.OR.MAX(A,B,C,D).GT.MNK_N) THEN
         WRITE(*,*) 'pyAmpliCol: invalid quadruple colour legs'
         STOP 1
      ENDIF
C Sum both ordered dipole products before projecting onto the real part.
      IA=(/MIN(A,B),MIN(C,D)/)
      IB=(/MAX(A,B),MAX(C,D)/)
      LEGS=0
      VECTOR=0D0
      TOTAL=0D0
      DO ORDERING=1,2
         WRITE(ID,'(A,I0,A,I0,A,I0,A,I0)') 'mnk-extra/dipole-product/',IA(ORDERING),'_',IB(ORDERING),'/',IA(3-ORDERING),'_',IB(3-ORDERING)
         CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,P,ALPHAS,ID,LEN_TRIM(ID),0,LEGS,VECTOR,PAIR)
         TOTAL=TOTAL+PAIR
      ENDDO
      CALL %(proc_prefix)sPYAMP_REAL(TOTAL,0D0,%(proc_prefix)sGET_QUADBLO)
      END
## }

## if (color_correlation) {
      DOUBLE PRECISION FUNCTION %(proc_prefix)sGET_CCBLO(P,ALPHAS,A,B)
      IMPLICIT NONE
      INTEGER A,B,LEGS(1)
      REAL*8 P(0:3,%(nexternal)d),ALPHAS,PAIR(2)
      COMPLEX*16 VECTOR(0:3,1)
      CHARACTER*64 ID
%(pyamp_declarations)s
      LEGS=0
      VECTOR=0D0
      WRITE(ID,'(A,I0,A,I0)') 'T_',A,'_',B
      CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,P,ALPHAS,ID,LEN_TRIM(ID),0,LEGS,VECTOR,PAIR)
      CALL %(proc_prefix)sPYAMP_REAL(PAIR,0D0,%(proc_prefix)sGET_CCBLO)
      END
## }

## if (spin_correlation_general_resources) {
      DOUBLE PRECISION FUNCTION %(proc_prefix)sGET_KKBLO(P,ALPHAS,LEG,KT)
      IMPLICIT NONE
      INTEGER LEG,LEGS(1)
      REAL*8 P(0:3,%(nexternal)d),ALPHAS,KT(0:3),PAIR(2)
      COMPLEX*16 VECTOR(0:3,1)
%(pyamp_declarations)s
      LEGS(1)=LEG
      VECTOR(:,1)=DCMPLX(KT,0D0)
      CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,P,ALPHAS,'born',4,1,LEGS,VECTOR,PAIR)
      CALL %(proc_prefix)sPYAMP_REAL(PAIR,0D0,%(proc_prefix)sGET_KKBLO)
      END
## }

      SUBROUTINE %(proc_prefix)sPYAMP_EVALUATE(P,HEL,AS,ANS)
      IMPLICIT NONE
      INTEGER HEL
      REAL*8 P(0:3,%(nexternal)d),AS,ANS(0:1)
      REAL*8 CC(%(pyamp_ncorr)d,0:1)
      CALL %(proc_prefix)sPYAMP_DISPATCH(.FALSE.,P,HEL,AS,ANS,CC)
      END

## if (color_correlation) {
      SUBROUTINE %(proc_prefix)sGET_COLOR_CORRELATED_ME_AT_POINT(P,HEL,AS,CC)
      IMPLICIT NONE
      INTEGER HEL
      REAL*8 P(0:3,%(nexternal)d),AS,ANS(0:1)
      REAL*8 CC(%(pyamp_ncorr)d,0:1)
CF2PY INTENT(IN) :: P, HEL, AS
CF2PY INTENT(OUT) :: CC
      CALL %(proc_prefix)sPYAMP_DISPATCH(.TRUE.,P,HEL,AS,ANS,CC)
      END
## }

      SUBROUTINE %(proc_prefix)sPYAMP_DISPATCH(WANT_COLOR,P,HEL,AS,ANS,CC)
      IMPLICIT NONE
      INTEGER NEXTERNAL,NCORR
      PARAMETER(NEXTERNAL=%(nexternal)d,NCORR=%(pyamp_ncorr)d)
      INTEGER HEL,I,J,K,R,NRUN,COUNT,NSPIN
      INTEGER INDICES(NCORR),CA(NCORR),CB(NCORR)
      REAL*8 P(0:3,NEXTERNAL),AS,ANS(0:1),CC(NCORR,0:1)
      REAL*8 PAIR(2),BORN_VALUE,RAW(0:1),VALUE,OUTPUT(2,NCORR)
      LOGICAL WANT_COLOR
      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/%(proc_prefix)sBORN_BEAM_POL/POLARIZATIONS
      LOGICAL CHOSEN_SO_CONFIGS(1)
      COMMON/%(proc_prefix)sCHOSEN_BORN_SQSO/CHOSEN_SO_CONFIGS
%(pyamp_declarations)s
## if (spin_correlation) {
      INCLUDE '%(proc_prefix)sspin_correlations.inc'
      INTEGER SPIN_COUNTS(MAX_SPIN_CORR_RUNS)
      INTEGER SPIN_LEGS(NEXTERNAL,MAX_SPIN_CORR_RUNS)
      COMPLEX*16 SPIN_VECTORS(0:3,NEXTERNAL,MAX_SPIN_CORR_RUNS)
## } else {
      INTEGER SPIN_COUNTS(1),SPIN_LEGS(NEXTERNAL,1)
      COMPLEX*16 SPIN_VECTORS(0:3,NEXTERNAL,1)
## }
## if (color_correlation) {
      INTEGER SELECTED(0:NCORR)
      REAL*8 UNUSED_TABLE(NCORR,0:1)
      COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/SELECTED,UNUSED_TABLE
      INTEGER CONNECTION_MAP(%(pyamp_nconnections)d,%(pyamp_nconnections)d)
      INTEGER INDEX_MAP(NCORR,2)
      COMMON/%(proc_prefix)sCOLOR_CORRELATION_MAPS/CONNECTION_MAP,INDEX_MAP
## }
      IF(HEL.NE.-1.OR.POLARIZATIONS(0,0).NE.-1) THEN
         WRITE(*,*) 'pyAmpliCol: selected helicities/polarizations are unsupported'
         STOP 1
      ENDIF
      NRUN=1
      SPIN_COUNTS=0
      SPIN_LEGS=0
      SPIN_VECTORS=0D0
## if (spin_correlation) {
      IF(N_SPIN_CORR_RUNS.LT.0.OR.N_SPIN_CORR_RUNS.GT.MAX_SPIN_CORR_RUNS) THEN
         WRITE(*,*) 'pyAmpliCol: invalid spin run count'
         STOP 1
      ENDIF
      NRUN=MAX(1,N_SPIN_CORR_RUNS)
      DO R=1,N_SPIN_CORR_RUNS
         DO I=1,NEXTERNAL
            K=SPIN_CORR_RUNS(R,I)
            IF(K.EQ.-1000) CYCLE
            IF(K.LT.1.OR.K.GT.MAX_N_SPIN_CORR_VECTORS.OR.K.GT.N_SPIN_CORR_VECTORS(I)) THEN
               WRITE(*,*) 'pyAmpliCol: invalid spin run entry'
               STOP 1
            ENDIF
            J=SPIN_COUNTS(R)+1
            SPIN_COUNTS(R)=J
            SPIN_LEGS(J,R)=I
            SPIN_VECTORS(:,J,R)=SPIN_CORR_VECTORS(I,K,:)
         ENDDO
      ENDDO
## }

## if (color_correlation) {
      IF(WANT_COLOR) THEN
         COUNT=0
         IF(SELECTED(0).LT.0.OR.SELECTED(0).GT.NCORR) THEN
            WRITE(*,*) 'pyAmpliCol: invalid colour selection count'
            STOP 1
         ENDIF
         DO I=1,SELECTED(0)
            J=SELECTED(I)
            IF(J.LT.1.OR.J.GT.NCORR) THEN
               WRITE(*,*) 'pyAmpliCol: invalid colour selection ID'
               STOP 1
            ENDIF
            IF(ANY(INDICES(1:COUNT).EQ.J)) CYCLE
            COUNT=COUNT+1
            INDICES(COUNT)=J
            CA(COUNT)=INDEX_MAP(J,1)
            CB(COUNT)=INDEX_MAP(J,2)
         ENDDO
         CC=0D0
         DO R=1,NRUN
            IF(COUNT.EQ.0) EXIT
            NSPIN=SPIN_COUNTS(R)
            CALL MNK_PYAMP_LEGACY_MANY(MNK_N,MNK_NINITIAL,MNK_PDGS,P,AS,COUNT,CA,CB,NSPIN,SPIN_LEGS(:,R),SPIN_VECTORS(:,:,R),OUTPUT)
            DO I=1,COUNT
               CALL %(proc_prefix)sPYAMP_REAL(OUTPUT(:,I),0D0,VALUE)
               CC(INDICES(I),1)=CC(INDICES(I),1)+VALUE
            ENDDO
         ENDDO
         IF(CHOSEN_SO_CONFIGS(1)) CC(:,0)=CC(:,1)
         RETURN
      ENDIF
## }

      CALL MNK_PYAMP_EVAL(MNK_N,MNK_NINITIAL,MNK_PDGS,P,AS,HEL,1,RAW)
      BORN_VALUE=RAW(1)
      IF(NRUN.GT.1.OR.SPIN_COUNTS(1).GT.0) THEN
         BORN_VALUE=0D0
         DO R=1,NRUN
            NSPIN=SPIN_COUNTS(R)
            CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,P,AS,'born',4,NSPIN,SPIN_LEGS(:,R),SPIN_VECTORS(:,:,R),PAIR)
            CALL %(proc_prefix)sPYAMP_REAL(PAIR,ABS(RAW(1)),VALUE)
            BORN_VALUE=BORN_VALUE+VALUE
         ENDDO
      ENDIF
      ANS(1)=BORN_VALUE
      ANS(0)=0D0
      IF(CHOSEN_SO_CONFIGS(1)) ANS(0)=BORN_VALUE
      END

      SUBROUTINE %(proc_prefix)sPYAMP_REAL(PAIR,SCALE,VALUE)
      IMPLICIT NONE
      REAL*8 PAIR(2),SCALE,VALUE
      IF(ABS(PAIR(2)).GT.1D-10*MAX(ABS(PAIR(1)),SCALE,1D-100)) THEN
         WRITE(*,*) 'pyAmpliCol: complex result cannot enter a real legacy table; use the named API'
         STOP 1
      ENDIF
      VALUE=PAIR(1)
      END
