C Native provider wrappers keep the public selected-state Fortran interface.
## if (pyamp_quadruple_color) {
      DOUBLE PRECISION FUNCTION %(proc_prefix)sGET_QUADBLO(A,B,C,D)
      IMPLICIT NONE
      INTEGER A,B,C,D
      REAL*8 P(0:3,%(nexternal)d),AS,ANS(0:1),VALUE
      REAL*8 CC(%(pyamp_ncorr)d,0:1)
      CALL %(proc_prefix)sPYAMP_DISPATCH(4,P,-1,AS,ANS,CC,A,B,C,D,VALUE)
      %(proc_prefix)sGET_QUADBLO=VALUE
      END
## }

      SUBROUTINE %(proc_prefix)sPYAMP_EVALUATE(P,HEL,AS,ANS)
      IMPLICIT NONE
      INTEGER HEL
      REAL*8 P(0:3,%(nexternal)d),AS,ANS(0:1)
      REAL*8 CC(%(pyamp_ncorr)d,0:1),VALUE
      CALL %(proc_prefix)sPYAMP_DISPATCH(0,P,HEL,AS,ANS,CC,0,0,0,0,VALUE)
      END

      SUBROUTINE %(proc_prefix)sPYAMP_TABLE(CC)
      IMPLICIT NONE
      REAL*8 P(0:3,%(nexternal)d),AS,ANS(0:1),VALUE
      REAL*8 CC(%(pyamp_ncorr)d,0:1)
      CALL %(proc_prefix)sPYAMP_DISPATCH(1,P,-1,AS,ANS,CC,0,0,0,0,VALUE)
      END

      SUBROUTINE %(proc_prefix)sPYAMP_DIPOLE(A,B,VALUE)
      IMPLICIT NONE
      INTEGER A,B
      REAL*8 P(0:3,%(nexternal)d),AS,ANS(0:1),VALUE
      REAL*8 CC(%(pyamp_ncorr)d,0:1)
      CALL %(proc_prefix)sPYAMP_DISPATCH(2,P,-1,AS,ANS,CC,A,B,0,0,VALUE)
      END

      SUBROUTINE %(proc_prefix)sPYAMP_RESET()
      IMPLICIT NONE
      REAL*8 P(0:3,%(nexternal)d),AS,ANS(0:1),VALUE
      REAL*8 CC(%(pyamp_ncorr)d,0:1)
      CALL %(proc_prefix)sPYAMP_DISPATCH(-1,P,-1,AS,ANS,CC,0,0,0,0,VALUE)
      END

      SUBROUTINE %(proc_prefix)sPYAMP_SPIN(LEG,P,KT,VALUE)
      IMPLICIT NONE
      INTEGER LEG,LEGS(1)
      REAL*8 P(0:3,%(nexternal)d),KT(0:3),VALUE,AS,ANS(0:1)
      REAL*8 CC(%(pyamp_ncorr)d,0:1),PAIR(2)
      COMPLEX*16 VECTOR(0:3,1)
      CHARACTER*4 ID
%(pyamp_declarations)s
      CALL %(proc_prefix)sPYAMP_DISPATCH(3,P,-1,AS,ANS,CC,0,0,0,0,VALUE)
      LEGS(1)=LEG
      VECTOR(:,1)=DCMPLX(KT,0D0)
      ID='born'
      CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,P,AS,ID,4,1,LEGS,VECTOR,PAIR)
C Check reality against this contraction, not a previous Born magnitude.
      CALL %(proc_prefix)sPYAMP_REAL(PAIR,0D0,VALUE)
      END

      SUBROUTINE %(proc_prefix)sPYAMP_DISPATCH(MODE,P,HEL,AS,ANS,CC,A,B,C,D,VALUE)
      IMPLICIT NONE
      INTEGER NEXTERNAL,NCORR
      PARAMETER(NEXTERNAL=%(nexternal)d,NCORR=%(pyamp_ncorr)d)
      INTEGER MODE,HEL,A,B,C,D,I,J,K,R,NRUN,COUNT,NSPIN,ORDERING
      INTEGER FIRST_A,FIRST_B,SECOND_A,SECOND_B
      INTEGER LEGS(NEXTERNAL),INDICES(NCORR),CA(NCORR),CB(NCORR)
      REAL*8 P(0:3,NEXTERNAL),AS,ANS(0:1),CC(NCORR,0:1),VALUE
      REAL*8 LAST_P(0:3,NEXTERNAL),LAST_AS,BORN_SCALE
      REAL*8 TABLE(NCORR,0:1),PAIR(2),BORN_VALUE,RAW(0:1)
      REAL*8 OUTPUT(2,NCORR)
      COMPLEX*16 VECTORS(0:3,NEXTERNAL)
      CHARACTER*64 ID
      LOGICAL READY,TABLE_READY,CHOSEN
      DATA READY/.FALSE./,TABLE_READY/.FALSE./
      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/%(proc_prefix)sBORN_BEAM_POL/POLARIZATIONS
      LOGICAL CHOSEN_SO_CONFIGS(1)
      COMMON/%(proc_prefix)sCHOSEN_BORN_SQSO/CHOSEN_SO_CONFIGS
%(pyamp_declarations)s
## if (spin_correlation) {
      INCLUDE '%(proc_prefix)sspin_correlations.inc'
      INTEGER SAVED_NSPIN(MAX_SPIN_CORR_RUNS)
      INTEGER SAVED_LEGS(NEXTERNAL,MAX_SPIN_CORR_RUNS)
      COMPLEX*16 SAVED_VECTORS(0:3,NEXTERNAL,MAX_SPIN_CORR_RUNS)
## } else {
      INTEGER SAVED_NSPIN(1),SAVED_LEGS(NEXTERNAL,1)
      COMPLEX*16 SAVED_VECTORS(0:3,NEXTERNAL,1)
## }
## if (color_correlation) {
      INTEGER SELECTED(0:NCORR)
      REAL*8 CORRELATED_EVALS(NCORR,0:1)
      COMMON/%(proc_prefix)sCOLOR_CORRELATIONS/SELECTED,CORRELATED_EVALS
      INTEGER CONNECTION_MAP(%(pyamp_nconnections)d,%(pyamp_nconnections)d)
      INTEGER INDEX_MAP(NCORR,2)
      COMMON/%(proc_prefix)sCOLOR_CORRELATION_MAPS/CONNECTION_MAP,INDEX_MAP
## }
C Snapshot inputs so later selection/reset calls cannot change an evaluation.
      SAVE LAST_P,LAST_AS,BORN_SCALE,TABLE,READY,TABLE_READY,CHOSEN
      SAVE NRUN,COUNT,INDICES,CA,CB,SAVED_NSPIN,SAVED_LEGS,SAVED_VECTORS

      IF(MODE.EQ.-1) THEN
         READY=.FALSE.
         TABLE_READY=.FALSE.
## if (color_correlation) {
         CORRELATED_EVALS=0D0
## }
         RETURN
      ENDIF
      IF(MODE.EQ.0) THEN
         READY=.FALSE.
         TABLE_READY=.FALSE.
         IF(POLARIZATIONS(0,0).NE.-1) THEN
            WRITE(*,*) 'pyAmpliCol: selected helicities/polarizations are unsupported'
            STOP 1
         ENDIF
         CALL MNK_PYAMP_EVAL(MNK_N,MNK_NINITIAL,MNK_PDGS,P,AS,HEL,1,RAW)
         LAST_P=P
         LAST_AS=AS
         BORN_SCALE=ABS(RAW(1))
         CHOSEN=CHOSEN_SO_CONFIGS(1)
         NRUN=1
         SAVED_NSPIN(1)=0
## if (spin_correlation) {
         IF(N_SPIN_CORR_RUNS.LT.0.OR.N_SPIN_CORR_RUNS.GT.MAX_SPIN_CORR_RUNS) THEN
            WRITE(*,*) 'pyAmpliCol: invalid spin run count'
            STOP 1
         ENDIF
         NRUN=MAX(1,N_SPIN_CORR_RUNS)
         DO R=1,NRUN
            SAVED_NSPIN(R)=0
            IF(N_SPIN_CORR_RUNS.EQ.0) CYCLE
            DO I=1,NEXTERNAL
               K=SPIN_CORR_RUNS(R,I)
               IF(K.EQ.-1000) CYCLE
               IF(K.LT.1.OR.K.GT.MAX_N_SPIN_CORR_VECTORS.OR.K.GT.N_SPIN_CORR_VECTORS(I)) THEN
                  WRITE(*,*) 'pyAmpliCol: invalid spin run entry'
                  STOP 1
               ENDIF
               J=SAVED_NSPIN(R)+1
               SAVED_NSPIN(R)=J
               SAVED_LEGS(J,R)=I
               SAVED_VECTORS(:,J,R)=SPIN_CORR_VECTORS(I,K,:)
            ENDDO
         ENDDO
## }
         BORN_VALUE=RAW(1)
         IF(NRUN.GT.1.OR.SAVED_NSPIN(1).GT.0) THEN
            BORN_VALUE=0D0
            ID='born'
            DO R=1,NRUN
               NSPIN=SAVED_NSPIN(R)
               CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,LAST_P,LAST_AS,ID,LEN_TRIM(ID),NSPIN,SAVED_LEGS(:,R),SAVED_VECTORS(:,:,R),PAIR)
               CALL %(proc_prefix)sPYAMP_REAL(PAIR,BORN_SCALE,VALUE)
               BORN_VALUE=BORN_VALUE+VALUE
            ENDDO
         ENDIF
         ANS(1)=BORN_VALUE
         ANS(0)=0D0
         IF(CHOSEN) ANS(0)=BORN_VALUE
         COUNT=0
## if (color_correlation) {
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
         CORRELATED_EVALS=0D0
## }
         TABLE=0D0
         READY=.TRUE.
         RETURN
      ENDIF

      IF(.NOT.READY) THEN
         WRITE(*,*) 'pyAmpliCol: Born hook must precede correlated getter'
         STOP 1
      ENDIF
      IF(MODE.EQ.1) THEN
         IF(.NOT.TABLE_READY) THEN
            DO R=1,NRUN
               NSPIN=SAVED_NSPIN(R)
               IF(COUNT.EQ.0) EXIT
               CALL MNK_PYAMP_LEGACY_MANY(MNK_N,MNK_NINITIAL,MNK_PDGS,LAST_P,LAST_AS,COUNT,CA,CB,NSPIN,SAVED_LEGS(:,R),SAVED_VECTORS(:,:,R),OUTPUT)
               DO I=1,COUNT
                  CALL %(proc_prefix)sPYAMP_REAL(OUTPUT(:,I),BORN_SCALE,VALUE)
                  TABLE(INDICES(I),1)=TABLE(INDICES(I),1)+VALUE
               ENDDO
            ENDDO
            IF(CHOSEN) TABLE(:,0)=TABLE(:,1)
            TABLE_READY=.TRUE.
## if (color_correlation) {
            CORRELATED_EVALS=TABLE
## }
         ENDIF
         CC=TABLE
      ELSE IF(MODE.EQ.2) THEN
         WRITE(ID,'(A,I0,A,I0)') 'T_',A,'_',B
         BORN_VALUE=0D0
         DO R=1,NRUN
            NSPIN=SAVED_NSPIN(R)
            CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,LAST_P,LAST_AS,ID,LEN_TRIM(ID),NSPIN,SAVED_LEGS(:,R),SAVED_VECTORS(:,:,R),PAIR)
            CALL %(proc_prefix)sPYAMP_REAL(PAIR,BORN_SCALE,VALUE)
            BORN_VALUE=BORN_VALUE+VALUE
         ENDDO
         VALUE=BORN_VALUE
      ELSE IF(MODE.EQ.3) THEN
         AS=LAST_AS
## if (pyamp_quadruple_color) {
      ELSE IF(MODE.EQ.4) THEN
         IF(MIN(A,B,C,D).LT.1.OR.MAX(A,B,C,D).GT.NEXTERNAL) THEN
            WRITE(*,*) 'pyAmpliCol: invalid quadruple colour legs'
            STOP 1
         ENDIF
         BORN_VALUE=0D0
C B_abcd = <M|D_ab D_cd + D_cd D_ab|M>, with full intermediate colour.
         DO ORDERING=1,2
            IF(ORDERING.EQ.1) THEN
               FIRST_A=MIN(A,B)
               FIRST_B=MAX(A,B)
               SECOND_A=MIN(C,D)
               SECOND_B=MAX(C,D)
            ELSE
               FIRST_A=MIN(C,D)
               FIRST_B=MAX(C,D)
               SECOND_A=MIN(A,B)
               SECOND_B=MAX(A,B)
            ENDIF
            WRITE(ID,'(A,I0,A,I0,A,I0,A,I0)') 'mnk-extra/dipole-product/',FIRST_A,'_',FIRST_B,'/',SECOND_A,'_',SECOND_B
            DO R=1,NRUN
               NSPIN=SAVED_NSPIN(R)
               CALL MNK_PYAMP_CORRELATED_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,LAST_P,LAST_AS,ID,LEN_TRIM(ID),NSPIN,SAVED_LEGS(:,R),SAVED_VECTORS(:,:,R),PAIR)
C Keep complex parts until the Hermitian combination has been assembled.
               IF(ORDERING.EQ.1.AND.R.EQ.1) OUTPUT(:,1)=0D0
               OUTPUT(:,1)=OUTPUT(:,1)+PAIR
            ENDDO
         ENDDO
         CALL %(proc_prefix)sPYAMP_REAL(OUTPUT(:,1),BORN_SCALE,VALUE)
## }
      ELSE
         WRITE(*,*) 'pyAmpliCol: invalid dispatch mode'
         STOP 1
      ENDIF
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
