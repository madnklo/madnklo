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
      INCLUDE 'leg_PDGs_epem_ddx.inc'
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

      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo
      
c
c     Born contribution
      do i=1,nexternal
         if(leg_pdgs_epem_ddx(i).eq.21) then
            INLO=INLO+
     &           (CA/6d0+2*TR*Nf/3d0)*(log(sLO(i,iref1(i))/MU_R**2)-8d0/3d0)+
     &           CA*(6d0-7d0/2d0*zeta2)
c           factor for conversion torino to ML, namely gamma[1-eps] to exp[- eps eulergamma]      
            INLO = INLO + pi**2/12d0 * CA
         elseif(leg_pdgs_epem_ddx(i).ne.0 .and.
     &           abs(leg_pdgs_epem_ddx(i)).le.6) then
            INLO=INLO+
     &           (CF/2d0)*(10d0-7d0*zeta2+log(sLO(i,iref1(i))/MU_R**2))
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
         if(leg_pdgs_epem_ddx(i).eq.21)INLO=INLO+
     &   CA*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_g-2d0*CA)*A20b
         if(leg_pdgs_epem_ddx(i).ne.0 .and.
     &        abs(leg_pdgs_epem_ddx(i)).le.6)INLO=INLO+
     &   CF*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_q-2d0*CF)*A20b
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
            INLO=INLO+ccBLO*log(sLO(i,j)/MU_R**2)*
     &      (2d0-log(sLO(i,j)/MU_R**2)/2d0)
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


      function A10(w)
c 
c     A10(w) = Psi0(w+1) + eulergamma
      implicit none
      double precision A10,w
c
      if(w.ne.0d0.and.w.ne.1d0.and.
     &   w.ne.2d0.and.w.ne.3d0)then
         write(*,*)'Value not coded in A10',w
         stop
      endif
c
      if(w.eq.0d0)A10=0d0
      if(w.eq.1d0)A10=1d0
      if(w.eq.2d0)A10=3d0/2d0
      if(w.eq.3d0)A10=11d0/6d0
c 
      return
      end

      function A20(w)
c 
c     A20(w) = Psi0(w+2) - 1 + eulergamma
      implicit none
      double precision A20,w
c
      if(w.ne.0d0.and.w.ne.1d0.and.
     &   w.ne.2d0.and.w.ne.3d0)then
         write(*,*)'Value not coded in A20',w
         stop
      endif
c
      if(w.eq.0d0)A20=0d0
      if(w.eq.1d0)A20=1d0/2d0
      if(w.eq.2d0)A20=5d0/6d0
      if(w.eq.3d0)A20=13d0/12d0
c 
      return
      end


      function A21(w)
c 
c     A21(w) = Psi1(w+2) + 1 - Zeta2
      implicit none
      double precision A21,w
c
      if(w.ne.0d0.and.w.ne.1d0.and.
     &   w.ne.2d0.and.w.ne.3d0)then
         write(*,*)'Value not coded in A21',w
         stop
      endif
c
      if(w.eq.0d0)A21=0d0
      if(w.eq.1d0)A21=-1d0/4d0
      if(w.eq.2d0)A21=-13d0/36d0
      if(w.eq.3d0)A21=-61d0/144d0
c 
      return
      end
