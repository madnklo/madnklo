      subroutine test_RR(iunit,x0,alpha2,beta2)
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      integer i,isec,jsec,ksec,lsec,iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref
      common/cNNLOsecindices/isec,jsec,ksec,lsec
      common/cNNLOmaplabels/iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref
      integer iunit,ievnt
      double precision x0(mxdim)
      double precision alpha2,beta2
      double precision e1,e2,e4,e5
      character*10 dash10
      save ievnt
c
      dash10='----------'
      ievnt=ievnt+1
c
      write(iunit,*)dash10//dash10//dash10//dash10
      write(iunit,*)dash10//dash10//dash10//dash10
      write(iunit,*)' EVENT NUMBER ',ievnt
      write(iunit,*)dash10//dash10//dash10//dash10
      write(iunit,*)dash10//dash10//dash10//dash10
c
c     double-soft limit
      if(is_NNLO_SS_singular_sec(isec,jsec,ksec,lsec))then
         e1=0d0
         e2=1d0
         e4=1d0
         e5=1d0
         call do_limit_RR(iunit,'SS      ',x0,isec,jsec,ksec,lsec,alpha2,beta2,e1,e2,e4,e5)
      endif
c
c     double-collinear limit
      if(is_NNLO_CC1_singular_sec(isec,jsec,ksec,lsec).and.
     &   topology(isec,jsec,ksec,lsec).ne.3)then
         e1=0d0
         e2=1d0
         e4=0d0
         e5=1d0
         call do_limit_RR(iunit,'CC      ',x0,isec,jsec,ksec,lsec,alpha2,beta2,e1,e2,e4,e5)
      endif
c
c     single-collinear limit
      if(is_NNLO_C_singular_sec(isec,jsec,ksec,lsec))then
         e1=0d0
         e2=1d0
         e4=0d0
         e5=0d0
         call do_limit_RR(iunit,'C       ',x0,isec,jsec,ksec,lsec,alpha2,beta2,e1,e2,e4,e5)
      endif
c
c     double-soft double-collinear limit
      if(is_NNLO_SS_singular_sec(isec,jsec,ksec,lsec).and.
     &   is_NNLO_CC1_singular_sec(isec,jsec,ksec,lsec).and.
     &   topology(isec,jsec,ksec,lsec).ne.3)then
         e1=0d0
         e2=2d0
         e4=1d0
         e5=2d0
         call do_limit_RR(iunit,'SS CC   ',x0,isec,jsec,ksec,lsec,alpha2,beta2,e1,e2,e4,e5)
      endif
c
c     double-soft single-collinear limit
      if(is_NNLO_SS_singular_sec(isec,jsec,ksec,lsec).and.
     &   is_NNLO_C_singular_sec(isec,jsec,ksec,lsec).and.
     &   topology(isec,jsec,ksec,lsec).eq.1)then
         e1=0d0
         e2=2d0
         e4=1d0
         e5=1d0
         call do_limit_RR(iunit,'SS C    ',x0,isec,jsec,ksec,lsec,alpha2,beta2,e1,e2,e4,e5)
      endif
c
c     double-collinear single-collinear limit
      if(is_NNLO_CC1_singular_sec(isec,jsec,ksec,lsec).and.
     &   is_NNLO_C_singular_sec(isec,jsec,ksec,lsec).and.
     &   topology(isec,jsec,ksec,lsec).ne.3)then
         e1=0d0
         e2=2d0
         e4=0d0
         e5=1d0
         call do_limit_RR(iunit,'CC C    ',x0,isec,jsec,ksec,lsec,alpha2,beta2,e1,e2,e4,e5)
      endif
c
c     double-soft double-collinear single-collinear limit
      if(is_NNLO_SS_singular_sec(isec,jsec,ksec,lsec).and.
     &   is_NNLO_CC1_singular_sec(isec,jsec,ksec,lsec).and.
     &   is_NNLO_C_singular_sec(isec,jsec,ksec,lsec).and.
     &   topology(isec,jsec,ksec,lsec).eq.1)then
         e1=0d0
         e2=3d0
         e4=1d0
         e5=2d0
         call do_limit_RR(iunit,'SS CC C  ',x0,isec,jsec,ksec,lsec,alpha2,beta2,e1,e2,e4,e5)
      endif
c
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)
c
      return
      end


      subroutine do_limit_RR(iunit,limstr,x0,isec,jsec,ksec,lsec,alpha2,beta2,e1,e2,e4,e5)
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      integer iitn,i,j,maxitn,iunit,ierr
      integer isec,jsec,ksec,lsec
      integer iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref,itopo
      parameter(maxitn=30)
      double precision x0(mxdim),x(mxdim)
      double precision xsNNLO(-2:maxdim,-2:maxdim)
      double precision xsNLO(-2:maxdim,-2:maxdim)
      double precision xsLO(-2:maxdim,-2:maxdim)
      double precision KSS,KHHCC,KSSC,KHHCCC,KS,KHC,K1NNLO,K2NNLO,K12NNLO
      double precision W_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision WS_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision WC_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision WSC_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision WSS_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision WCC_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision WSSCC_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision WCCC_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision WSSCCC_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision WSSC_NNLO(maxdim,maxdim,maxdim,maxdim)
      double precision W_NLO(maxdim,maxdim)
      double precision Wsum,WsumSS,WsumCC,WsumSSCC,WsumCCC,WsumSSCCC,WsumSSC,WsumS,WsumC
      double precision WS_NLO(maxdim,maxdim),WC_NLO(maxdim,maxdim)
      double precision W_NLO_2(maxdim,maxdim)
      double precision WS_NLO_2(maxdim,maxdim),WC_NLO_2(maxdim,maxdim)
      double precision lam,lim,RRNNLO,double_real
      double precision alpha2,beta2
      double precision e1,e2,e4,e5
      character*5 str5
      character*8 limstr
      character*10 str10
      double precision p(0:3,-2:npartNNLO)
      double precision pb(0:3,-2:npartNLO)
      double precision pbb(0:3,-2:npartLO)
      double precision xjac,xjacB
      double precision BLO
c
c     initialise
      x=x0
      str5 ='     '
      str10='          '
      xjac=0d0
      xsNNLO=0d0
      xsNLO=0d0
      xsLO=0d0
      itopo=topology(isec,jsec,ksec,lsec)
c
c     assign iU1, iS1, iB1, iA1, iA2, iref
      call assign_phsp_labels_npt_test(limstr,isec,jsec,ksec,lsec,
     &iU1,iS1,iB1,iA1,iA2,iref)
c
c     start testing
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)'LIM = '//trim(limstr)
      write(iunit,*)str10//'lambda'//str10//str10//'RR'//
     &str10//str10//str5//'LIM'//
     &str10//str10//'|RR-LIM|/|LIM|'
c$$$c
c$$$c     possibility to set by hand the starting point
c$$$c     for the limiting procedure
c$$$      x0(1)=0.5d0
c$$$      x0(2)=0.5d0
c$$$      x0(4)=0.5d-7
c$$$      x0(5)=0.5d0
c
c     loop to get closer and closer to the limit
      do iitn=1,maxitn
         lam=10d0**(1d0-iitn)
c
c     initialise counterterms
         KS=0d0
         KHC=0d0
         KSS=0d0
         KHHCC=0d0
         KSSC=0d0
         KHHCCC=0d0
         K1NNLO=0d0
         K2NNLO=0d0
         K12NNLO=0d0
c
c     rescale relevant x random numbers
c     x(1) is zCS', while x(2) is yCS'
c     x(4) is zCS , while x(5) is yCS
         x(1)=x0(1)*lam**e1
         x(2)=x0(2)*lam**e2
         x(4)=x0(4)*lam**e4
         x(5)=x0(5)*lam**e5
c
c     recompute momenta after rescaling
         call phase_space_npt(x,sCM,iU1,iS1,iB1,iA1,iA2,p,pb,pbb,xjac,xjacB,iU2,iS2,iB2)
         if(xjac.eq.0d0)cycle
         call invariants_from_p(p,npartNNLO,xsNNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(pb,npartNLO,xsNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(pbb,npartLO,xsLO,ierr)
         if(ierr.eq.1)cycle
c
c     double real
         call double_real_NNLO(xsNNLO,RRNNLO,iref,ierr)
         if(ierr.eq.1)cycle
         call get_W_NNLO(xsNNLO,p,xsNLO,pb,alpha_spec,beta_spec,alpha2,beta2,
     &   isec,jsec,ksec,lsec,iref,is_NNLO_1unres_singular_sec(isec,jsec,ksec,lsec),
     &   Wsum,WsumS,WsumC,WsumSS,WsumSSC,ierr)
         if(ierr.eq.1)cycle
c
c     double-unresolved counterterm
         if(is_NNLO_2unres_singular_sec(isec,jsec,ksec,lsec))then
            call local_pure_2unres_counter_NNLO(isec,jsec,ksec,lsec,iref,
     &      xsNNLO,p,xsNLO,pb,xsLO,pbb,iA2,1d0,WsumSS,1d0,1d0,
     &      1d0,KSS,KHHCC,K2NNLO,ierr)
            if(ierr.eq.1)cycle
            call local_mixed_2unres_counter_NNLO(isec,jsec,ksec,lsec,iref,
     &      xsNNLO,p,xsNLO,pb,xsLO,pbb,iU1,iS1,iB1,iU2,iS2,iB2,1d0,1d0,
     &      1d0,WsumSSC,1d0,KSSC,KHHCCC,K12NNLO,ierr)
            if(ierr.eq.1)cycle
         endif
c
c     single-unresolved counterterm
         if(is_NNLO_1unres_singular_sec(isec,jsec,ksec,lsec))then
            call local_1unres_counter_NNLO(isec,jsec,ksec,lsec,iref,
     &      xsNNLO,p,xsNLO,pb,xsLO,pbb,1d0,WsumS,WsumC,1d0,KS,KHC,K1NNLO,ierr)
            if(ierr.eq.1)cycle
         endif
c
         lim=K1NNLO+K2NNLO+K12NNLO
         double_real=RRNNLO*Wsum
         if(abs(lim).gt.0d0)then
            write(iunit,*)lam,double_real,lim,abs(double_real-lim)/abs(lim)
         else
            write(iunit,*)lam,double_real,lim,double_real,' *** '
         endif
      enddo
      x=x0
c
      return
      end
