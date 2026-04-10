      subroutine int_counter_I1_NNLO_%(isec)d_%(jsec)d(p,sNLO,sLO,I1NNLO,ierr)
c     MSbar integrated counterterm
c     FINITE_PART = I1NNLO(0)
c     SINGLE_POLE = I1NNLO(-1)
c     DOUBLE_POLE = I1NNLO(-2)
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
      integer i,j,k,r,iref1(nexternal),ierr
      double precision p(0:3,nexternal)
      double precision sNLO(nexternal,nexternal),sLO(nexternal-1,nexternal-1)
      double precision I1NNLO(-2:0),I1NNLOS(-2:0),I1S(-2:0),pref
      double precision RNLO,CCRNLO
      double precision Lij,Lkr
      double precision A20a,A21a,A20b,A20,A21
      DOUBLE PRECISION ALPHAS,ANS(0:NSQSO_BORN)
      DOUBLE PRECISION ALPHA_QCD
      INTEGER, PARAMETER :: HEL = - 1
      double precision  %(NLO_proc_str)sGET_CCBLO
      double precision vv,ypl,Q2,ddilog
      DOUBLE PRECISION SS,MK2,ML2
      double precision res
      integer isec,jsec,ksec,lsec,iref0
      common/csecindices/isec,jsec,ksec,lsec,iref0
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
      isec = %(isec)d
      jsec = %(jsec)d
      pref=alphas/(2d0*pi)
      I1NNLO = 0d0
      I1NNLOS = 0d0
      iref1 = 0
      CCRNLO = 0d0
      RNLO = 0d0
      res = 0d0
c
c     recoilers
      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo
c
c     Hard-collinear contribution
      do k=1,nexternal
         if(pmass(k).ne.0d0)cycle
         if(iref1(k).eq.0)cycle
         if(iref1(k).eq.isec.or.iref1(k).eq.jsec.or.iref1(k).eq.k)cycle
         Lkr = log(sNLO(k,iref1(k))/MU_R**2)
         if(abs(leg_pdgs_%(NLO_process)s(k)).le.6) then
            I1NNLO(-2) = 0d0
            I1NNLO(-1) = I1NNLO(-1) + gamma_hc_q
            I1NNLO( 0) = I1NNLO( 0) - gamma_hc_q*Lkr + phi_hc_q
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I1NNLO( 0) = I1NNLO( 0) + pi**2/12d0 * CF
         elseif(leg_pdgs_%(NLO_process)s(k).eq.21) then
            I1NNLO(-2) = 0d0
            I1NNLO(-1) = I1NNLO(-1) + gamma_hc_g
            I1NNLO( 0) = I1NNLO( 0) - gamma_hc_g*Lkr + phi_hc_g
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I1NNLO( 0) = I1NNLO( 0) + pi**2/12d0 * CF
         endif
      enddo
c
c     TODO: check (gamma/gamma_hc)
c     Include damping factors
      A20a=A20(alpha)
      A21a=A21(alpha)
      A20b=A20(beta_FF)
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(NLO_process)s(i).eq.21)I1NNLO(0) = I1NNLO(0) + CA*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_g-2d0*CA)*A20b
         if(leg_pdgs_%(NLO_process)s(i).ne.0.and.abs(leg_pdgs_%(NLO_process)s(i)).le.6)I1NNLO(0) = I1NNLO(0) + CF*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_q-2d0*CF)*A20b
      enddo
c
c     real
      call %(NLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      RNLO = ANS(0)
      if(RNLO.lt.0d0.or.abs(RNLO).ge.huge(1d0).or.isnan(RNLO))then
         write(77,*) 'int_real: '
         write(77,*) 'Wrong RNLO', RNLO
         goto 999
      endif
      I1NNLO=I1NNLO*RNLO
c
c     Soft contribution
      do i=1,nexternal
         if(.not.ISNLOQCDPARTON(i))cycle
         do j=1,nexternal
            if(.not.ISNLOQCDPARTON(j))cycle
            if(j.eq.i)cycle
            if(pmass(i).eq.0d0.and.pmass(j).eq.0d0)then
               Lij = log(sNLO(i,j)/MU_R**2)
               I1S(-2) = 1d0
               I1S(-1) = 2d0 - Lij
               I1S( 0) = 6d0-7d0/2d0*zeta2 - 2d0*Lij + Lij**2d0/2d0
               call %(NLO_proc_str)sME_ACCESSOR_HOOK(p,hel,alphas,ans)
               CCRNLO = %(NLO_proc_str)sGET_CCBLO(I,J)
               I1S = -I1S*CCRNLO
               I1NNLOS = I1NNLOS + I1S
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
               call nnlo_rv_sub(ss,vv,mk2,ml2,mu_r,ccRNLO,res)
               I1NNLO(0) = res
               I1NNLO(-1) = I1NNLO(-1) + CCRNLO*(-1D0/2D0)*(2D0 - 1D0/VV*DLOG((1D0+VV)/(1D0-VV)) )
            endif
         enddo
      enddo
c
      I1NNLO=I1NNLO+I1NNLOS
      I1NNLO=pref*I1NNLO
c
      if(abs(I1NNLO(0)).ge.huge(1d0).or.isnan(I1NNLO(0)))then
         write(77,*)'Exception caught in int_counter_I1_NNLO',I1NNLO(0)
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
