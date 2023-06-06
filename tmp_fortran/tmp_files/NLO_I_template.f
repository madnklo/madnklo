      subroutine int_counter_NLO(p,sLO,INLO,ierr)
c     MSbar finite part of the integrated counterterm
      implicit none
      INCLUDE 'nexternal.inc'
      INCLUDE 'math.inc'
      INCLUDE 'damping_factors.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'virtual_recoilers.inc'
      INCLUDE 'leg_PDGs_%(proc_prefix)s.inc'
      INCLUDE 'colored_partons.inc'
      integer i,j,r
      integer ierr
      double precision p(0:3,nexternal)
      double precision sLO(nexternal,nexternal)
      double precision INLO,pref
      double precision BLO,ccBLO
      double precision A20a,A21a,A20b,A20,A21
      DOUBLE PRECISION ALPHAS,ANS(0:NSQSO_BORN)
      DOUBLE PRECISION ALPHA_QCD
      INTEGER, PARAMETER :: HEL = - 1
      DOUBLE PRECISION  GET_CCBLO
      integer iref1(nexternal)
c
c     initialise
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
      pref=alphas/(2d0*pi)
      INLO=0d0
      iref1 = 0
      CCBLO = 0d0
      BLO = 0d0

      CALL ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BLO = ANS(0)
c      TODO: modify ierr      
c      if(ierr.eq.1)goto 999

c     TODO: add check 
      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo

c
c     Born contribution
      do i=1,nexternal
         if(leg_pdgs_%(proc_prefix)s(i).eq.21) then
            INLO=INLO+(CA/6d0+2*TR*Nf/3d0)*(log(sLO(i,iref1(i))/MU_R**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c           factor for conversion torino to ML, namely gamma[1-eps] to exp[- eps eulergamma]      
            INLO = INLO + pi**2/12d0 * CA
         elseif(leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6) then
            INLO=INLO+(CF/2d0)*(10d0-7d0*zeta2+log(sLO(i,iref1(i))/MU_R**2))
c           factor for conversion torino to ML, namely gamma[1-eps] to exp[- eps eulergamma]  
            INLO = INLO + pi**2/12d0 * CF
         endif
      enddo

c
c     Include damping factors
      A20a=A20(alpha)
      A21a=A21(alpha)
      A20b=A20(beta_FF)
      do i=1,nexternal
         if(leg_pdgs_%(proc_prefix)s(i).eq.21)INLO=INLO+CA*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_g-2d0*CA)*A20b
         if(leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6)INLO=INLO+CF*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_q-2d0*CF)*A20b
      enddo

      INLO=INLO*BLO
c
c     Colour-linked-Born contribution
      do i=1,nexternal
         do j=1,nexternal
            if(ISLOQCDPARTON(i) .and. ISLOQCDPARTON(j)) then
               CCBLO = GET_CCBLO(i,j)
            else
               cycle
            endif
            if(j.eq.i)cycle
            INLO=INLO+ccBLO*log(sLO(i,j)/MU_R**2)*(2d0-log(sLO(i,j)/MU_R**2)/2d0)
         enddo
      enddo
      INLO=INLO*pref
c
      if(abs(INLO).ge.huge(1d0).or.isnan(INLO))then
         write(77,*)'Exception caught in int_counter_NLO',INLO
         goto 999
      endif

c
      return
 999  ierr=1
      return
      end


      FUNCTION A10(W)
C                                                                                                                                                                                 
C     A10(w) = Psi0(w+1) + eulergamma                                                                                                                                             
      IMPLICIT NONE
      DOUBLE PRECISION A10,W
C                                                                                                                                                                                 
      IF(W.NE.0D0.AND.W.NE.1D0.AND.W.NE.2D0.AND.W.NE.3D0.AND.W.NE.4D0.AND.W.NE.5D0)THEN                                                                                           
        WRITE(*,*)'Value not coded in A10',W
        STOP
      ENDIF
C                                                                                                                                                                                 
      IF(W.EQ.0D0)A10=0D0
      IF(W.EQ.1D0)A10=1D0
      IF(W.EQ.2D0)A10=3D0/2D0
      IF(W.EQ.3D0)A10=11D0/6D0
      IF(W.EQ.4D0)A10=25D0/12D0
      IF(W.EQ.5D0)A10=137/60D0
C                                                                                                                                                                                 
      RETURN
      END


      FUNCTION A20(W)
C                                                                                                                                                                                 
C     A20(w) = Psi0(w+2) - 1 + eulergamma                                                                                                                                         
      IMPLICIT NONE
      DOUBLE PRECISION A20,W
C                                                                                                                                                                                 
      IF(W.NE.0D0.AND.W.NE.1D0.AND.W.NE.2D0.AND.W.NE.3D0.AND.W.NE.4D0.AND.W.NE.5D0)THEN                                                                                           
         WRITE(*,*)'Value not coded in A20',W
        STOP
      ENDIF
C                                                                                                                                                                                 
      IF(W.EQ.0D0)A20=0D0
      IF(W.EQ.1D0)A20=1D0/2D0
      IF(W.EQ.2D0)A20=5D0/6D0
      IF(W.EQ.3D0)A20=13D0/12D0
      IF(W.EQ.4D0)A20=77D0/60D0
      IF(W.EQ.5D0)A20=29D0/20D0
C                                                                                                                                                                                 
      RETURN
      END


      FUNCTION A21(W)
C                                                                                                                                                                                 
C     A21(w) = Psi1(w+2) + 1 - Zeta2                                                                                                                                              
      IMPLICIT NONE
      DOUBLE PRECISION A21,W
C                                                                                                                                                                                 
      IF(W.NE.0D0.AND.W.NE.1D0.AND.W.NE.2D0.AND.W.NE.3D0.AND.W.NE.4D0.AND.W.NE.5D0)THEN                                                                                           
        WRITE(*,*)'Value not coded in A21',W
        STOP
      ENDIF
C                                                                                                                                                                                 
      IF(W.EQ.0D0)A21=0D0
      IF(W.EQ.1D0)A21=-1D0/4D0
      IF(W.EQ.2D0)A21=-13D0/36D0
      IF(W.EQ.3D0)A21=-61D0/144D0
      IF(W.EQ.4D0)A21=-1669D0/3600D0
      IF(W.EQ.5D0)A21=-1769D0/3600D0
C                                                                                                                                                                                 
      RETURN
      END

