
      SUBROUTINE FILL_MAPPED_LABELS(P1,P2,REAL_LEG_PDGS,
     $  BORN_LEG_PDGS)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'rv_mapped_labels_common.inc'
      INTEGER P1,P2,REAL_LEG_PDGS(NEXTERNAL),BORN_LEG_PDGS(NEXTERNAL-1)
      INTEGER OUT_DUMMY(NEXTERNAL)
      INTEGER K

C     if the combination P1-P2 is new, compute and save the result in
C     MAPPED_LABELS; if already done, do nothing (cache already filled)
      IF (.NOT. MAPPED_PAIR_DONE(P1,P2)) THEN
        CALL GET_MAPPED_LABELS(NEXTERNAL,P1,P2,REAL_LEG_PDGS
     $   ,BORN_LEG_PDGS,OUT_DUMMY)
        DO K=1,NEXTERNAL
          MAPPED_LABELS(K,P1,P2)=OUT_DUMMY(K)
        ENDDO
        MAPPED_PAIR_DONE(P1,P2)=.TRUE.
      ENDIF

      RETURN
      END
