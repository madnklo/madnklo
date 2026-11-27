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
      double precision Lij,Lkr,sigma_c,sigma_gamma
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
      ierr=0
c     renormalisation scale for the kinematics of this integrated counterterm
      call set_mur_from_momenta(p,nexternal,leg_pdgs_%(NLO_process)s)
      alphas=alphas_current
      isec = %(isec)d
      jsec = %(jsec)d
      pref=alphas/(2d0*pi)
      I1NNLO = 0d0
      I1NNLOS = 0d0
      iref1 = 0
      CCRNLO = 0d0
      RNLO = 0d0
      res = 0d0

c     TODO: only eps poles (indep. of kinematics)
      sigma_c = 2d0*cf+ca
      call epem_guux_me_accessor_hook(p,hel,alphas,ans)
      rnlo = ans(0)
      i1nnlo(-2) = i1nnlo(-2) + sigma_c * rnlo

      sigma_gamma = 2d0*gamma_q+gamma_g
      i1nnlo(-1) = i1nnlo(-1) + sigma_gamma*rnlo

      do i=1,nexternal
        if(.not.isnloqcdparton(i))cycle
        do j=1,nexternal
          if(.not.isnloqcdparton(j))cycle
          if(j.eq.i)cycle
          if(pmass(i).eq.0d0.and.pmass(j).eq.0d0)then
            call epem_guux_me_accessor_hook(p,hel,alphas,ans)
            ccrnlo = epem_guux_get_ccblo(i,j)
            lij = log(snlo(i,j)/mur_current**2)
            i1nnlo(-1) = i1nnlo(-1) + lij*ccrnlo
          endif
        enddo
      enddo

      I1NNLO=pref*I1NNLO

c
c     recoilers
c      do i=1,len_iref
c         iref1(iref(1,i)) = iref(2,i)
c      enddo
c
c     Hard-collinear contribution
c      do k=1,nexternal
c         if(pmass(k).ne.0d0)cycle
c         if(iref1(k).eq.0)cycle
c         if(iref1(k).eq.isec.or.iref1(k).eq.jsec.or.iref1(k).eq.k)cycle
c         Lkr = log(sNLO(k,iref1(k))/muR_current**2)
c         if(abs(leg_pdgs_%(NLO_process)s(k)).le.6) then
c            I1NNLO(-2) = 0d0
c            I1NNLO(-1) = I1NNLO(-1) + gamma_hc_q
c            I1NNLO( 0) = I1NNLO( 0) - gamma_hc_q*Lkr + phi_hc_q
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
c            I1NNLO( 0) = I1NNLO( 0) + pi**2/12d0 * CF
c         elseif(leg_pdgs_%(NLO_process)s(k).eq.21) then
c            I1NNLO(-2) = 0d0
c            I1NNLO(-1) = I1NNLO(-1) + gamma_hc_g
c            I1NNLO( 0) = I1NNLO( 0) - gamma_hc_g*Lkr + phi_hc_g
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
c            I1NNLO( 0) = I1NNLO( 0) + pi**2/12d0 * CF
c         endif
c      enddo
c
c     TODO: check (gamma/gamma_hc)
c     Include damping factors
c      A20a=A20(alpha)
c      A21a=A21(alpha)
c      A20b=A20(beta_FF)
c      do i=1,nexternal
c         if(pmass(i).ne.0d0)cycle
c         if(leg_pdgs_%(NLO_process)s(i).eq.21)I1NNLO(0) = I1NNLO(0) + CA*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_g-2d0*CA)*A20b
c         if(leg_pdgs_%(NLO_process)s(i).ne.0.and.abs(leg_pdgs_%(NLO_process)s(i)).le.6)I1NNLO(0) = I1NNLO(0) + CF*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_q-2d0*CF)*A20b
c      enddo
c
c     real
c      call %(NLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
c      RNLO = ANS(0)
c      if(RNLO.lt.0d0.or.abs(RNLO).ge.huge(1d0).or.isnan(RNLO))then
c         write(77,*) 'int_real: '
c         write(77,*) 'Wrong RNLO', RNLO
c         goto 999
c      endif
c      I1NNLO=I1NNLO*RNLO
c
c     Soft contribution
c      do i=1,nexternal
c         if(.not.ISNLOQCDPARTON(i))cycle
c         do j=1,nexternal
c            if(.not.ISNLOQCDPARTON(j))cycle
c            if(j.eq.i)cycle
c            if(pmass(i).eq.0d0.and.pmass(j).eq.0d0)then
c               Lij = log(sNLO(i,j)/muR_current**2)
c               I1S(-2) = 1d0
c               I1S(-1) = 2d0 - Lij
c               I1S( 0) = 6d0-7d0/2d0*zeta2 - 2d0*Lij + Lij**2d0/2d0
c               call %(NLO_proc_str)sME_ACCESSOR_HOOK(p,hel,alphas,ans)
c               CCRNLO = %(NLO_proc_str)sGET_CCBLO(I,J)
c               I1S = -I1S*CCRNLO
c               I1NNLOS = I1NNLOS + I1S
c            elseif(pmass(i).eq.0d0.and.pmass(j).ne.0d0)then
c               continue
c            elseif(pmass(i).ne.0d0.and.pmass(j).eq.0d0)then
c               continue
c            elseif(pmass(i).ne.0d0.and.pmass(j).ne.0d0)then
c               ML2=PMASS(I)**2
c               MK2=PMASS(J)**2
c               SS=SLO(I,J)
c               VV=DSQRT(SS**2-4D0*ML2*MK2)/SS
c               Q2=SS+ML2+MK2
c               YPL=1D0+(DSQRT(ML2)-DSQRT(Q2))*2D0*DSQRT(ML2)/SS
c               call nnlo_rv_sub(ss,vv,mk2,ml2,muR_current,ccRNLO,res)
c               I1NNLO(0) = res
c               I1NNLO(-1) = I1NNLO(-1) + CCRNLO*(-1D0/2D0)*(2D0 - 1D0/VV*DLOG((1D0+VV)/(1D0-VV)) )
c            endif
c         enddo
c      enddo
c
c      I1NNLO=I1NNLO+I1NNLOS
c      I1NNLO=pref*I1NNLO
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
