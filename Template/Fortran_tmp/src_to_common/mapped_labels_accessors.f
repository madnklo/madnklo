
      SUBROUTINE FILL_REAL_MAPPED_LABELS(P1,P2,LEG_PDGS,
     $  REAL_LEG_PDGS)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'mapped_labels_common.inc'
      INTEGER P1,P2,LEG_PDGS(NEXTERNAL),REAL_LEG_PDGS(NEXTERNAL-1)
      INTEGER OUT_DUMMY(NEXTERNAL)
      INTEGER K

C     if the combination P1-P2 is new, compute and save the result in
C     REAL_LABELS; if already done, do nothing (cache already filled)
      IF (.NOT. REAL_PAIR_DONE(P1,P2)) THEN
        CALL GET_MAPPED_LABELS(NEXTERNAL,P1,P2,LEG_PDGS
     $   ,REAL_LEG_PDGS,OUT_DUMMY)
        DO K=1,NEXTERNAL
          REAL_LABELS(P1,P2,K)=OUT_DUMMY(K)
        ENDDO
        REAL_PAIR_DONE(P1,P2)=.TRUE.
      ENDIF

      RETURN
      END

      
      SUBROUTINE FILL_BORN_MAPPED_LABELS(P1,P2,REAL_LEG_PDGS,
     $  BORN_LEG_PDGS)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'mapped_labels_common.inc'
      INTEGER P1,P2,REAL_LEG_PDGS(NEXTERNAL-1),BORN_LEG_PDGS(NEXTERNAL-2)
      INTEGER OUT_DUMMY(NEXTERNAL-1)
      INTEGER K

C     if the combination P1-P2 is new, compute and save the result in
C     BORN_LABELS; if already done, do nothing (cache already filled)
      IF (.NOT. BORN_PAIR_DONE(P1,P2)) THEN
        CALL GET_MAPPED_LABELS(NEXTERNAL-1,P1,P2,REAL_LEG_PDGS
     $   ,BORN_LEG_PDGS,OUT_DUMMY)
        DO K=1,NEXTERNAL-1
          BORN_LABELS(P1,P2,K)=OUT_DUMMY(K)
        ENDDO
        BORN_PAIR_DONE(P1,P2)=.TRUE.
      ENDIF

      RETURN
      END


c$$$      SUBROUTINE GET_BORN_MAPPED_LABELS(P1,P2,REAL_LEG_PDGS,
c$$$     $  BORN_LEG_PDGS,OUT)
c$$$      IMPLICIT NONE
c$$$      INCLUDE 'nexternal.inc'
c$$$      INCLUDE 'mapped_labels_common.inc'
c$$$      INTEGER P1,P2,REAL_LEG_PDGS(MAXLEG),BORN_LEG_PDGS(MAXLEG)
c$$$      INTEGER OUT(NEXTERNAL-1)
c$$$      INTEGER K
c$$$
c$$$C     if the combination P1-P2 is new, compute and save the result OUT in BORN_LABELS
c$$$C     if the combination has already been called, read the stored results from BORN_LABELS
c$$$      IF (.NOT. BORN_PAIR_DONE(P1,P2)) THEN
c$$$        CALL GET_MAPPED_LABELS(NEXTERNAL-1,P1,P2,REAL_LEG_PDGS
c$$$     $   ,BORN_LEG_PDGS,OUT)
c$$$        DO K=1,NEXTERNAL-1
c$$$          BORN_LABELS(P1,P2,K)=OUT(K)
c$$$        ENDDO
c$$$        BORN_PAIR_DONE(P1,P2)=.TRUE.
c$$$      ELSE
c$$$        DO K=1,NEXTERNAL-1
c$$$          OUT(K)=BORN_LABELS(P1,P2,K)
c$$$        ENDDO
c$$$      ENDIF
c$$$      
c$$$      RETURN
c$$$      END

    
c$$$      SUBROUTINE GET_REAL_MAPPED_LABELS(P1,P2,LEG_PDGS,
c$$$     $  REAL_LEG_PDGS,OUT)
c$$$      IMPLICIT NONE
c$$$      INCLUDE 'nexternal.inc'
c$$$      INCLUDE 'mapped_labels_common.inc'
c$$$      INTEGER P1,P2,LEG_PDGS(NEXTERNAL),REAL_LEG_PDGS(NEXTERNAL-1)
c$$$      INTEGER OUT(NEXTERNAL)
c$$$      INTEGER K
c$$$
c$$$C     if the combination ICONF/P1-P2 is new, compute and save
c$$$C     the result OUT in REAL_LABELS,
c$$$C     otherwise read the stored results from REAL_LABELS
c$$$      IF (.NOT. REAL_PAIR_DONE(P1,P2)) THEN
c$$$        CALL GET_MAPPED_LABELS(NEXTERNAL,P1,P2,LEG_PDGS
c$$$     $   ,REAL_LEG_PDGS,OUT)
c$$$        DO K=1,NEXTERNAL
c$$$          REAL_LABELS(P1,P2,K)=OUT(K)
c$$$        ENDDO
c$$$        REAL_PAIR_DONE(P1,P2)=.TRUE.
c$$$      ELSE
c$$$        DO K=1,NEXTERNAL
c$$$          OUT(K)=REAL_LABELS(P1,P2,K)
c$$$        ENDDO
c$$$      ENDIF
c$$$      
c$$$      RETURN
c$$$      END
