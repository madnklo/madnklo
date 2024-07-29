      subroutine int_counter2_NNLO(p,sLO,INLO,ierr)
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
      double precision INLO(5),pref
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
c      INLO(i) = coeff*eps^-(i-1)


c     Born contribution
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(proc_prefix)s(i).eq.21) then
            INLO(5) = INLO(5) - 1d0/2d0*CA
            INLO(4) = INLO(4) + CA*(3d0/8d0*beta0-gamma_g)
            INLO(3) = INLO(3) + 1d0/4d0*(2d0*beta0*gamma_g-1d0/4d0*CA*((8d0/3d0-4d0*zeta2)*CA+10d0/3d0*beta0)-2d0*gamma_g**2)
            INLO(2) = INLO(2) - 1d0/8d0*4d0*gamma2_g
            

            
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
c     EQ. (39) (+ FF-DEPENDENT PIECE)
c               INLO(1) = INLO(1) - CCBLO*1D0/VV*(DLOG((1D0+VV)/(1D0-VV))*(-YPL+1D0/2D0*DLOG(SS*YPL**2/MK2)+1D0/2D0*DLOG(SS/MU_R**2))+1D0/4D0*DLOG((1D0+VV)/(1D0-VV))**2+DDILOG(-2D0*VV/(1D0-VV)) )
c               INLO(1) = INLO(1) - CCBLO*1D0/VV*DLOG((1D0+VV)/(1D0-VV))*(-(SS+MK2)/SS*DLOG((SS*YPL+MK2)/MK2)+YPL)*(1D0-FF1)
c     EQ. (40) (MADE FF-DEPENDENT)
c               INLO(1) = INLO(1) - CCBLO*(-(SS+MK2)/SS*DLOG((SS*YPL+MK2)/MK2)+YPL) * FF2
c     EQ. (43)
c               INLO(1) = INLO(1) - CCBLO*(-1D0/2D0*DLOG(SS/MU_R**2) -1D0/2D0/SS*((2D0*MK2+SS)*DLOG(MK2)-2D0*SS+SS*DLOG(SS)+2D0*SS*DLOG(YPL)-2D0*(MK2+SS)*DLOG(MK2+SS*YPL)) )
c     EQ. (44) (+ FF-DEPENDENT PIECE)
c               INLO(1) = INLO(1) - CCBLO*(YPL-SS*YPL/MK2+SS*YPL**2/2D0/MK2+1D0/2D0*DLOG(MK2)-1D0/2D0*DLOG(SS)-DLOG(YPL)+1D0/2D0/VV*DLOG((1D0+VV)/(1D0-VV))-1D0/2D0*DLOG(SS/MU_R**2))
c
c     EQ. (44) (+ FF-DEPENDENT PIECE)
c               INLO(1) = INLO(1) - CCBLO*((((8d0*ml2**2*dsqrt(ss**2)*VV-2d0*ml2*ss**2*(-1d0+VV**2)-(ss**3-(ss**2)**1.5d0*VV)*(-1d0+VV**2))*dlog(1d0-(ss*(2d0*ml2+ss-dsqrt(ss**2)*VV)*ypl)/(2d0*ml2*dsqrt(ss**2)*VV)))/((2d0*ml2+ss-dsqrt(ss**2)*VV)*(ss+dsqrt(ss**2)*VV))+((8d0*ml2**2*dsqrt(ss**2)*VV+2d0*ml2*ss**2*(-1d0+VV**2)+(ss**3+(ss**2)**1.5d0*VV)*(-1d0+VV**2))*dlog(1d0+(ss*(2d0*ml2+ss+dsqrt(ss**2)*VV)*ypl)/(2d0*ml2*dsqrt(ss**2)*VV)))/((-ss+dsqrt(ss**2)*VV)*(2d0*ml2+ss+dsqrt(ss**2)*VV)))/(2d0*ss)+(2d0*mk2*ml2*(-(ss*dlog(1d0+(2d0*ss*VV)/(dsqrt(ss**2)-ss*VV)))+2d0*dsqrt(ss**2)*VV*dlog((ss*ypl)/(dsqrt(mk2)*mu_r))))/((ss**2)**1.5d0*VV*(-1d0+VV**2)))
c               
c               INLO(1) = INLO(1) - CCBLO*(FF3-1D0)/2D0/SS/MK2*(SS*YPL*(2D0*(1D0-FF3)*MK2-(FF3+1D0)*SS*(2D0-YPL))+2D0*(FF3-1D0)*MK2*(MK2+SS)*DLOG((SS*YPL+MK2)/MK2))
C     
               
            call nlo_v_sub(ss,vv,mk2,ml2,mu_r,ccBLO,res)

            INLO(1) = res
               
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

