      subroutine int_counter_NLO(p,sLO,INLO,ierr)
c     MSbar integrated counterterm
c     FINITE_PART = INLO(0)
c     SINGLE_POLE = INLO(-1)
c     DOUBLE_POLE = INLO(-2)
      implicit none
      INCLUDE 'nexternal.inc'
      INCLUDE 'damping_factors.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'virtual_recoilers.inc'
      INCLUDE 'leg_PDGs_%(proc_prefix)s.inc'
      INCLUDE 'colored_partons.inc'
      integer i,j,r
      integer ierr
      double precision p(0:3,nexternal)
      double precision sLO(nexternal,nexternal)
      double precision INLO(-2:0),pref
      double precision BLO,ccBLO
      double precision A20a,A21a,A20b,A20,A21
      DOUBLE PRECISION ALPHAS,ANS(0:NSQSO_BORN)
      DOUBLE PRECISION ALPHA_QCD
      INTEGER, PARAMETER :: HEL = - 1
      DOUBLE PRECISION  GET_CCBLO
      integer iref1(nexternal)
      double precision vv,ypl,Q2,ddilog
      double precision pmass(nexternal)
      DOUBLE PRECISION SS,MK2,ML2
      DOUBLE PRECISION FF1,FF2,FF3
      PARAMETER(FF1=1D0,FF2=1D0,FF3=0D0)
      double precision res
c     massive-massless (2503.14629, final-state radiation)
      double precision LKL,LSM,rhom
      double precision IS_M,IHC_0G,IHC_1G
      double precision IHC_2G
      logical massive_coloured
      include 'pmass.inc'
c
c     initialise
      ierr = 0
c     renormalisation scale for the n-body kinematics
      call set_mur_from_momenta(p,nexternal,leg_pdgs_%(proc_prefix)s)
      alphas=alphas_current
      pref=alphas/(2d0*pi)
      INLO = 0d0
      iref1 = 0
      CCBLO = 0d0
      BLO = 0d0
      res = 0d0

      CALL ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BLO = ANS(0)
c
c     TODO: add check 
      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo
c
c     Born contribution
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(proc_prefix)s(i).eq.21) then
            INLO(0) = INLO(0) + (CA/6d0+2*TR*Nf/3d0)*(log(sLO(i,iref1(i))/muR_current**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])      
            INLO(0) = INLO(0) + pi**2/12d0 * CA
            INLO(-1) = INLO(-1) + gamma_g
            INLO(-2) = INLO(-2) + CA
         elseif(leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6) then
            INLO(0) = INLO(0) + (CF/2d0)*(10d0-7d0*zeta2+log(sLO(i,iref1(i))/muR_current**2))
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            INLO(0) = INLO(0) + pi**2/12d0 * CF
            INLO(-1) = INLO(-1) + gamma_q
            INLO(-2) = INLO(-2) + CF
         endif
c     hard-collinear counterterm with a massive recoiler r=iref1(i):
c     extra finite term -sum_n gamma^hc_(i,n) I^(ng)_hc,M(s_ir,m_r)
c     (2503.14629, Sec. 2.5 and eq. IhcM)
         if(pmass(iref1(i)).ne.0d0)then
            rhom=pmass(iref1(i))**2/sLO(i,iref1(i))
            if(leg_pdgs_%(proc_prefix)s(i).eq.21)then
               INLO(0) = INLO(0) + 2d0/3d0*TR*Nf*IHC_0G(rhom) + CA/6d0*IHC_2G(rhom)
            elseif(abs(leg_pdgs_%(proc_prefix)s(i)).le.6)then
               INLO(0) = INLO(0) + CF/2d0*IHC_1G(rhom)
            endif
         endif
      enddo
c
c     The massive integrated counterterms are computed without damping
c     factors: stop if damping is switched on with massive partons
      massive_coloured=.false.
      do i=1,nexternal
         if(ISLOQCDPARTON(i).and.pmass(i).ne.0d0)massive_coloured=.true.
      enddo
      if(massive_coloured.and.(alpha.ne.0d0.or.beta_FF.ne.0d0))then
         write(*,*)'int_counter_NLO: damping factors (alpha, beta_FF)'
         write(*,*)'are not available with massive coloured partons.'
         write(*,*)'Set them to 0 in Cards/damping_factors.inc'
         stop
      endif
c
c     Include damping factors
      A20a=A20(alpha)
      A21a=A21(alpha)
      A20b=A20(beta_FF)
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(proc_prefix)s(i).eq.21)INLO(0) = INLO(0) + CA*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_g-2d0*CA)*A20b
         if(leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6)INLO(0) = INLO(0) + CF*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_q-2d0*CF)*A20b
      enddo
c
      INLO=INLO*BLO
c
c     Colour-linked-Born contribution
      do i=1,nexternal
         if(.not.ISLOQCDPARTON(i))cycle
         do j=1,nexternal
            if(.not.ISLOQCDPARTON(j))cycle
            if(j.eq.i)cycle
            CCBLO = GET_CCBLO(i,j)
            if(pmass(i).eq.0d0.and.pmass(j).eq.0d0)then
               INLO(0) = INLO(0) + ccBLO*log(sLO(i,j)/muR_current**2)*(2d0-log(sLO(i,j)/muR_current**2)/2d0)
               INLO(-1) = INLO(-1) + ccBLO*log(sLO(i,j)/muR_current**2)
            elseif(pmass(i).eq.0d0.and.pmass(j).ne.0d0)then
c     massless i, massive j: integrated massive-massless soft
c     counterterm (2503.14629, eq. IsMF), per ordered pair. The sign of
c     the ln(s_ij/m_j^2) terms follows the paper's integrals (the
c     assembled formula for I in the paper has them with opposite sign)
               LKL=log(sLO(i,j)/muR_current**2)
               LSM=log(sLO(i,j)/pmass(j)**2)
               INLO(0) = INLO(0) + ccBLO*(2d0*LKL-LKL**2/2d0-LKL*LSM -2d0*IS_M(pmass(j)**2/sLO(i,j)))
               INLO(-1) = INLO(-1) + ccBLO*(LKL+LSM)
            elseif(pmass(i).ne.0d0.and.pmass(j).eq.0d0)then
c     massive i, massless j: L_ij term, plus the share
c     -(1/eps+4(1-zeta2)) B_ij of the massive-leg term
c     C_i (1/eps+4(1-zeta2)) B, written through colour conservation
c     as in the massive-massive case below
               LKL=log(sLO(i,j)/muR_current**2)
               INLO(0) = INLO(0) + ccBLO*(LKL-4d0*(1d0-zeta2))
               INLO(-1) = INLO(-1) - ccBLO
            elseif(pmass(i).ne.0d0.and.pmass(j).ne.0d0)then
               ML2=PMASS(I)**2
               MK2=PMASS(J)**2
               SS=SLO(I,J)
               VV=DSQRT(SS**2-4D0*ML2*MK2)/SS
               Q2=SS+ML2+MK2
               YPL=1D0+(DSQRT(ML2)-DSQRT(Q2))*2D0*DSQRT(ML2)/SS
               CALL NLO_I_MASS(ss,vv,mk2,ml2,muR_current,ccBLO,res)
               INLO(0) = INLO(0) + res
               INLO(-1) = INLO(-1) + CCBLO*(-1D0/2D0)*(2D0 - 1D0/VV*DLOG((1D0+VV)/(1D0-VV)) )
            endif
         enddo
      enddo
      INLO = INLO*pref
c
      if(abs(INLO(0)).ge.huge(1d0).or.isnan(INLO(0)))then
         write(77,*)'Exception caught in int_counter_NLO',INLO(0)
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



      subroutine NLO_I_MASS(s,v,mk2,ml2,mu,ccBLO,INLO_MASS)
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      double precision s,v,yp
      double precision mk2,ml2,mu
      double precision INLO_MASS,ddilog,ccBLO,Q2
      double precision lam, eta, alpha, beta
      double precision Q, mk, ml
c
      Q2=s+ml2+mk2
      lam=(s*v)**2d0
      Q = dsqrt(Q2)
      mk = dsqrt(mk2)
      ml = dsqrt(ml2)
      eta = (1d0-v)/(1d0+v)
      alpha = (Q-ml)/(Q+ml)
      beta = dsqrt( ((Q-ml)**2-mk2)/((Q+ml)**2-mk2) )
c
      if(alpha.lt.beta.or.alpha.gt.1d0.or.beta.lt.0d0)then
         write(*,*)'Wrong alpha, beta, values in NLO_I_MASS'
         write(*,*)alpha,beta
         stop
      endif
c
      INLO_MASS =-0.5d0*CCBLO*(8d0+(4d0*mk2*dlog((1d0-beta)/(1d0+beta)))/dsqrt(lam)-(4d0*ml2*dlog((1d0-beta)/(1d0+beta)))/dsqrt(lam) - (2d0*s*dlog((1d0-beta)/(1d0+beta))*dlog((1d0-beta**2)/(2d0*(1d0+alpha))))/dsqrt(lam)+(2d0*ml2*dlog(eta))/dsqrt(lam) - (2d0*Q2*dlog(eta))/dsqrt(lam)+(s*dlog(eta))/dsqrt(lam)-(s*dlog(eta)**2)/(2d0*dsqrt(lam)) - 2d0*dlog(lam**2/(mk2*ml2*Q2*s))-(s*dlog(eta)*dlog(lam**2/(mk2*ml2*Q2*s)))/dsqrt(lam) - 2d0*dlog(s/mu**2)-(s*dlog(eta)*dlog(s/mu**2))/dsqrt(lam)-(2d0*s*ddilog((1d0-alpha)/(1d0-beta)))/dsqrt(lam) - (2d0*s*ddilog((1d0-beta)/2d0))/dsqrt(lam) -(2d0*s*ddilog((1d0-beta)/(1d0+alpha)))/dsqrt(lam) + (2d0*s*ddilog((-2d0*beta)/(1d0-beta)))/dsqrt(lam) +(2d0*s*ddilog((-2d0*beta)/(alpha-beta)))/dsqrt(lam) + (2d0*s*ddilog((1d0-alpha)/(1d0+beta)))/dsqrt(lam) -(2d0*s*ddilog((2d0*beta)/(1d0+beta)))/dsqrt(lam) + (2d0*s*ddilog((1d0+beta)/2d0))/dsqrt(lam) +(2d0*s*ddilog((1d0+beta)/(1d0+alpha)))/dsqrt(lam) - (2d0*s*ddilog((2d0*beta)/(alpha+beta)))/dsqrt(lam) -(2d0*s*ddilog(1d0-eta))/dsqrt(lam))
c
      end


      double precision function IS_M(rho)
c     finite part of the massive-massless integrated soft counterterm,
c     rho=m^2/s (2503.14629, eq. IsMF)
      implicit none
      double precision rho,ddilog
      external ddilog
      IS_M = -ddilog(-rho) - log(rho)**2/4d0 + (1d0+rho)*log(1d0+rho) + (0.5d0-rho)*log(rho)
      return
      end


      double precision function IHC_0G(rho)
c     massive-recoiler hard-collinear finite part, g -> q qbar,
c     rho=m_r^2/s (2503.14629, eq. IhcM)
      implicit none
      double precision rho
      IHC_0G = rho + 1.5d0*rho*log(rho) - (1d0+1.5d0*rho)*log(1d0+rho) - 2d0*rho*sqrt(rho)*atan(sqrt(1d0+rho)-sqrt(rho))
      return
      end


      double precision function IHC_1G(rho)
c     massive-recoiler hard-collinear finite part, q -> q g
      implicit none
      double precision rho
      IHC_1G = rho*log(rho) - (1d0+rho)*log(1d0+rho)
      return
      end


      double precision function IHC_2G(rho)
c     massive-recoiler hard-collinear finite part, g -> g g
      implicit none
      double precision rho
      IHC_2G = -2d0*rho - log(1d0+rho) + 4d0*rho*sqrt(rho)*atan(sqrt(1d0+rho)-sqrt(rho))
      return
      end
