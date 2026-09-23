      subroutine int_counter_NLO(p,sLO,INLO,ierr)
c     MSbar integrated counterterm
c     FINITE_PART = INLO(1)
c     SINGLE_POLE = INLO(2)
c     DOUBLE_POLE = INLO(3)
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
      double precision INLO(3),pref
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
      include 'pmass.inc'
c
c     initialise
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
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
            INLO(1) = INLO(1) + (CA/6d0+2*TR*Nf/3d0)*(log(sLO(i,iref1(i))/MU_R**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])      
            INLO(1) = INLO(1) + pi**2/12d0 * CA
            INLO(2) = INLO(2) + gamma_g
            INLO(3) = INLO(3) + CA
         elseif(leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6) then
            INLO(1) = INLO(1) + (CF/2d0)*(10d0-7d0*zeta2+log(sLO(i,iref1(i))/MU_R**2))
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            INLO(1) = INLO(1) + pi**2/12d0 * CF
            INLO(2) = INLO(2) + gamma_q
            INLO(3) = INLO(3) + CF
         endif
      enddo
c
c     Include damping factors
      A20a=A20(alpha)
      A21a=A21(alpha)
      A20b=A20(beta_FF)
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(proc_prefix)s(i).eq.21)INLO(1) = INLO(1) + CA*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_g-2d0*CA)*A20b
         if(leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6)INLO(1) = INLO(1) + CF*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_q-2d0*CF)*A20b
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
               INLO(1) = INLO(1) + ccBLO*log(sLO(i,j)/MU_R**2)*(2d0-log(sLO(i,j)/MU_R**2)/2d0)
               INLO(2) = INLO(2) + ccBLO*log(sLO(i,j)/MU_R**2)
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
               CALL NLO_I_MASS(ss,vv,mk2,ml2,mu_r,ccBLO,res)
               INLO(1) = INLO(1) + res
               INLO(2) = INLO(2) + CCBLO*(-1D0/2D0)*(2D0 - 1D0/VV*DLOG((1D0+VV)/(1D0-VV)) )
            endif
         enddo
      enddo
      INLO = INLO*pref
c
      if(abs(INLO(1)).ge.huge(1d0).or.isnan(INLO(1)))then
         write(77,*)'Exception caught in int_counter_NLO',INLO(1)
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
