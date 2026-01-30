      subroutine int_counter_I12_NNLO_%(isec)d_%(jsec)d(p,sNLO,sLO,I12NNLO,ierr)
c     MSbar integrated counterterm
c     FINITE_PART = I12NNLO(0)
c     SINGLE_POLE = I12NNLO(-1)
c     DOUBLE_POLE = I12NNLO(-2)
      use sectors2_module
      implicit none
      INCLUDE 'nexternal.inc'
      INCLUDE 'damping_factors.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'virtual_recoilers.inc'
      INCLUDE 'leg_PDGs_%(NLO_process)s.inc'
      INCLUDE 'colored_partons.inc'
      integer i,j,r
      integer l,m,n,q
      integer isec, jsec
      integer ierr
      double precision p(0:3,nexternal)
      double precision sNLO(nexternal,nexternal),sLO(nexternal-1,nexternal-1)
      double precision I12NNLO(-2:0),pref
      double precision BLO,ccBLO
      double precision A20a,A21a,A20b,A20,A21
      DOUBLE PRECISION ALPHAS,ANS(0:NSQSO_BORN)
      DOUBLE PRECISION ALPHA_QCD
      INTEGER, PARAMETER :: HEL = - 1
      DOUBLE PRECISION  GET_CCBLO
      DOUBLE PRECISION TRIBLO,QUADBLO
      integer iref1(nexternal)
      double precision vv,ypl,Q2,ddilog
      double precision pmass(nexternal)
      DOUBLE PRECISION SS,MK2,ML2
      DOUBLE PRECISION FF1,FF2,FF3
      PARAMETER(FF1=1D0,FF2=1D0,FF3=0D0)
      double precision res
      include 'pmass.inc'
c
c     initialise
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
      isec = %(isec)d
      jsec = %(jsec)d
      pref=alphas/(2d0*pi)
      I12NNLO = 0d0
      iref1 = 0
      CCBLO = 0d0
      BLO = 0d0
      res = 0d0
c
c     Simple Born
c      CALL ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BLO = ANS(0)
c
c     TODO: add check
      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo
c     Tripole Born
c      TRIBLO=GET_TRIBLO(L,M,Q)
c     Quadruple Born
c      QUADBLO=GET_QUADBLO(L,M,Q,N)
c
c     Soft contribution
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(NLO_process)s(i).eq.21) then
            I12NNLO(0) = I12NNLO(0) + (CA/6d0+2*TR*Nf/3d0)*(log(sNLO(i,iref1(i))/MU_R**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CA
            I12NNLO(-1) = I12NNLO(-1) + gamma_g
            I12NNLO(-2) = I12NNLO(-2) + CA
         elseif(leg_pdgs_%(NLO_process)s(i).ne.0 .and.abs(leg_pdgs_%(NLO_process)s(i)).le.6) then
            I12NNLO(0) = I12NNLO(0) + (CF/2d0)*(10d0-7d0*zeta2+log(sNLO(i,iref1(i))/MU_R**2))
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CF
            I12NNLO(-1) = I12NNLO(-1) + gamma_q
            I12NNLO(-2) = I12NNLO(-2) + CF
         endif
      enddo
c
c     Collinear contribution
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(NLO_process)s(i).eq.21) then
            I12NNLO(0) = I12NNLO(0) + (CA/6d0+2*TR*Nf/3d0)*(log(sNLO(i,iref1(i))/MU_R**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CA
            I12NNLO(-1) = I12NNLO(-1) + gamma_g
            I12NNLO(-2) = I12NNLO(-2) + CA
         elseif(leg_pdgs_%(NLO_process)s(i).ne.0 .and.abs(leg_pdgs_%(NLO_process)s(i)).le.6) then
            I12NNLO(0) = I12NNLO(0) + (CF/2d0)*(10d0-7d0*zeta2+log(sNLO(i,iref1(i))/MU_R**2))
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CF
            I12NNLO(-1) = I12NNLO(-1) + gamma_q
            I12NNLO(-2) = I12NNLO(-2) + CF
         endif
      enddo

c
c     Soft-Collinear  contribution
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(NLO_process)s(i).eq.21) then
            I12NNLO(0) = I12NNLO(0) + (CA/6d0+2*TR*Nf/3d0)*(log(sNLO(i,iref1(i))/MU_R**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CA
            I12NNLO(-1) = I12NNLO(-1) + gamma_g
            I12NNLO(-2) = I12NNLO(-2) + CA
         elseif(leg_pdgs_%(NLO_process)s(i).ne.0 .and.abs(leg_pdgs_%(NLO_process)s(i)).le.6) then
            I12NNLO(0) = I12NNLO(0) + (CF/2d0)*(10d0-7d0*zeta2+log(sNLO(i,iref1(i))/MU_R**2))
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CF
            I12NNLO(-1) = I12NNLO(-1) + gamma_q
            I12NNLO(-2) = I12NNLO(-2) + CF
         endif
      enddo

c
c     Include damping factors
      A20a=A20(alpha)
      A21a=A21(alpha)
      A20b=A20(beta_FF)
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(NLO_process)s(i).eq.21)I12NNLO(0) = I12NNLO(0) + CA*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_g-2d0*CA)*A20b
         if(leg_pdgs_%(NLO_process)s(i).ne.0.and.abs(leg_pdgs_%(NLO_process)s(i)).le.6)I12NNLO(0) = I12NNLO(0) + CF*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_q-2d0*CF)*A20b
      enddo
c
      I12NNLO=I12NNLO*BLO
c
c     Colour-linked-Born contribution
      do i=1,nexternal
         if(.not.ISNLOQCDPARTON(i))cycle
         do j=1,nexternal
            if(.not.ISNLOQCDPARTON(j))cycle
            if(j.eq.i)cycle
            CCBLO = GET_CCBLO(i,j)
            if(pmass(i).eq.0d0.and.pmass(j).eq.0d0)then
               I12NNLO(0) = I12NNLO(0) +ccBLO*log(sNLO(i,j)/MU_R**2)*(2d0-log(sNLO(i,j)/MU_R**2)/2d0)
               I12NNLO(-1) = I12NNLO(-1) + ccBLO*log(sNLO(i,j)/MU_R**2)
            elseif(pmass(i).eq.0d0.and.pmass(j).ne.0d0)then
               continue
            elseif(pmass(i).ne.0d0.and.pmass(j).eq.0d0)then
               continue
            elseif(pmass(i).ne.0d0.and.pmass(j).ne.0d0)then
               ML2=PMASS(I)**2
               MK2=PMASS(J)**2
               SS=SLO(I,J)
               VV=DSQRT(SS**2-4D0*ML2*MK2)/SS
               Q2=SS+ML2+MK2
               YPL=1D0+(DSQRT(ML2)-DSQRT(Q2))*2D0*DSQRT(ML2)/SS
C
c           this file was removed in the end - now there is I1 and I12; do we mean rvsub file here? (_y)
            call nnlo_irv_sub(ss,vv,mk2,ml2,mu_r,ccBLO,res)

            I12NNLO(0) = I12NNLO(0) + res

               I12NNLO(-1) = I12NNLO(-1) + CCBLO*(-1D0/2D0)*(2D0 - 1D0/VV*DLOG((1D0+VV)/(1D0-VV)) )
            endif
         enddo
      enddo
      I12NNLO = I12NNLO*pref
c
      if(abs(I12NNLO(0)).ge.huge(1d0).or.isnan(I12NNLO(0)))then
         write(77,*)'Exception caught in int_counter_I12_NNLO',I12NNLO(0)
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end




      FUNCTION A10(W)
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
