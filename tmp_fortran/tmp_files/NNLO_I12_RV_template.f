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
      INCLUDE 'leg_PDGs_%(proc_prefix_real)s.inc'
      INCLUDE 'colored_partons.inc'
      integer i,j,r
      integer l,m,n,q
      integer lb,mb,nb,qb
      integer ierr
      double precision p(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sNLO(nexternal,nexternal),sLO(nexternal-1,nexternal-1)
      double precision I12NNLO(-2:0),pref
      double precision BLO,ccBLO,triBLO,quadBLO
      double precision A20a,A21a,A20b,A20,A21
      DOUBLE PRECISION ALPHAS,ANS(0:NSQSO_BORN)
      DOUBLE PRECISION ALPHA_QCD
      INTEGER, PARAMETER :: HEL = - 1
      double precision  %(proc_prefix_S_RV_g)s_GET_CCBLO
      double precision  %(proc_prefix_S_RV_g)s_GET_TRIBLO
      double precision  %(proc_prefix_S_RV_g)s_GET_QUADBLO
      integer iref1(nexternal)
      double precision vv,ypl,Q2,ddilog
      DOUBLE PRECISION SS,MK2,ML2
      DOUBLE PRECISION FF1,FF2,FF3
      PARAMETER(FF1=1D0,FF2=1D0,FF3=0D0)
      double precision res
c      integer isec,jsec,ksec,lsec,iref
c      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
      pref=alphas/(2d0*pi)
      I12NNLO = 0d0
      iref1 = 0
      BLO = 0d0
      CCBLO = 0d0
      TRIBLO = 0d0
      QUADBLO = 0d0
      res = 0d0  
c
c     TODO: add check
      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo
c     
c     TODO: Born and colour correlations
      xpb=0d0
      lb=mapped_labels(l)
      mb=mapped_labels(m)
      nb=mapped_labels(n)
      qb=mapped_labels(q)
c      call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = 0d0 !ANS(0)
      ccBLO = 0d0 !%(proc_prefix_S_RV_g)s_GET_CCBLO(lb,mb)
c      call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      triBLO = 0d0 !%(proc_prefix_S_RV_g)s_GET_TRIBLO(lb,mb,nb)
c      call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      quadBLO = 0d0 !%(proc_prefix_S_RV_g)s_GET_QUADBLO(lb,mb,qb,nb)
c
c     Soft contribution
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(proc_prefix_real)s(i).eq.21) then
            I12NNLO(0) = I12NNLO(0) + (CA/6d0+2*TR*Nf/3d0)*(log(sNLO(i,iref1(i))/MU_R**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CA
            I12NNLO(-1) = I12NNLO(-1) + gamma_g
            I12NNLO(-2) = I12NNLO(-2) + CA
         elseif(leg_pdgs_%(proc_prefix_real)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix_real)s(i)).le.6) then
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
         if(leg_pdgs_%(proc_prefix_real)s(i).eq.21) then
            I12NNLO(0) = I12NNLO(0) + (CA/6d0+2*TR*Nf/3d0)*(log(sNLO(i,iref1(i))/MU_R**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CA
            I12NNLO(-1) = I12NNLO(-1) + gamma_g
            I12NNLO(-2) = I12NNLO(-2) + CA
         elseif(leg_pdgs_%(proc_prefix_real)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix_real)s(i)).le.6) then
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
         if(leg_pdgs_%(proc_prefix_real)s(i).eq.21) then
            I12NNLO(0) = I12NNLO(0) + (CA/6d0+2*TR*Nf/3d0)*(log(sNLO(i,iref1(i))/MU_R**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CA
            I12NNLO(-1) = I12NNLO(-1) + gamma_g
            I12NNLO(-2) = I12NNLO(-2) + CA
         elseif(leg_pdgs_%(proc_prefix_real)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix_real)s(i)).le.6) then
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
         if(leg_pdgs_%(proc_prefix_real)s(i).eq.21)I12NNLO(0) = I12NNLO(0) + CA*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_g-2d0*CA)*A20b
         if(leg_pdgs_%(proc_prefix_real)s(i).ne.0.and.abs(leg_pdgs_%(proc_prefix_real)s(i)).le.6)I12NNLO(0) = I12NNLO(0) + CF*(A20a*(A20a-2d0*A20b)-A21a)+(gamma_q-2d0*CF)*A20b
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
            call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
            CCBLO = %(proc_prefix_S_RV_g)s_GET_CCBLO(i,j)
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
c            call nnlo_irv_sub(ss,vv,mk2,ml2,mu_r,ccBLO,res)

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
