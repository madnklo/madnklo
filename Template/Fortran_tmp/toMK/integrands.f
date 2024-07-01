      function int_double_real(x,wgt)
c     (n+2)-body NNLO integrand for vegas
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      integer ndim,ierr,itopo
      integer i,j,k,l,ievt,nthres,ntest
      integer iunit
      common/ciunitNNLO/iunit
      integer ntested(-2:maxdim,-2:maxdim,-2:maxdim,-2:maxdim)
      parameter(ntest=30)
      save ievt,nthres,ntested
      double precision int_double_real,int_double_real_no_cnt
      double precision sNNLO(-2:maxdim,-2:maxdim),sminNNLO
      double precision sNLO(-2:maxdim,-2:maxdim)
      double precision sLO(-2:maxdim,-2:maxdim)
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
      double precision alpha2,beta2
      parameter(alpha2=1d0,beta2=1d0)
      double precision RRNNLO,K1NNLO,K2NNLO,K12NNLO
      double precision KSS,KHHCC,KSSC,KHHCCC,KS,KHC
      double precision x(mxdim)
      double precision wgt,wgtpl
      common/cdim/ndim
      logical dotechcut
      double precision tinycut
      logical doplot
      common/cdoplot/doplot
      integer isec,jsec,ksec,lsec,iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref
      common/cNNLOsecindices/isec,jsec,ksec,lsec
      common/cNNLOmaplabels/iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref
      double precision p(0:3,-2:npartNNLO)
      double precision pb(0:3,-2:npartNLO)
      double precision pbb(0:3,-2:npartLO)
      double precision xjac,xjacB
c
c     initialise
      xjac=0d0
      int_double_real=0d0
      int_double_real_no_cnt=0d0
      sNNLO=0d0
      sNLO=0d0
      sLO=0d0
      K1NNLO=0d0
      K2NNLO=0d0
      K12NNLO=0d0
      KSS=0d0
      KHHCC=0d0
      KSSC=0d0
      KHHCCC=0d0
      KS=0d0
      KHC=0d0
      itopo=topology(isec,jsec,ksec,lsec)
c
c     phase space and invariants
      call phase_space_npt(x,sCM,iU1,iS1,iB1,iA1,iA2,p,pb,pbb,xjac,xjacB,iU2,iS2,iB2)
      if(xjac.eq.0d0)goto 999
      call invariants_from_p(p,npartNNLO,sNNLO,ierr)
      if(ierr.eq.1)goto 999
      call invariants_from_p(pb,npartNLO,sNLO,ierr)
      if(ierr.eq.1)goto 999
      call invariants_from_p(pbb,npartLO,sLO,ierr)
      if(ierr.eq.1)goto 999
c
c     check that phase-space labels and x variables are as expected
      call check_phsp_consistency(x,npartNNLO,sNNLO,sNLO,iU1,iS1,iB1,iA1,ierr)
      if(ierr.eq.1)goto 999
c
c     tiny technical phase-space cut to avoid fluctuations
      if(.not.is_NNLO_singular_sec(isec,jsec,ksec,lsec))then
         tinycut=0d0
      else
         tinycut=1d-7
         tinycut=1d-8
      endif
      if(dotechcut(x,npartNNLO,tinycut))goto 999
c
c     possible cuts
      if(docut(p,npartNNLO))goto 555
c
c     test matrix elements
      if(is_NNLO_singular_sec(isec,jsec,ksec,lsec).and.
     &   ntested(isec,jsec,ksec,lsec).le.ntest)then
         call test_RR(iunit,x,alpha2,beta2)
         ntested(isec,jsec,ksec,lsec)=ntested(isec,jsec,ksec,lsec)+1
      endif
c
c     double real
      call double_real_NNLO(sNNLO,RRNNLO,iref,ierr)
      if(ierr.eq.1)goto 999
      call get_W_NNLO(sNNLO,p,sNLO,pb,alpha_spec,beta_spec,alpha2,beta2,
     &isec,jsec,ksec,lsec,iref,is_NNLO_1unres_singular_sec(isec,jsec,ksec,lsec),
     &Wsum,WsumS,WsumC,WsumSS,WsumSSC,ierr)
      if(ierr.eq.1)goto 999
c
c     full double real in the combination of sectors
      int_double_real_no_cnt=RRNNLO*Wsum*xjac
c
c     plot double real
      wgtpl=int_double_real_no_cnt*wgt/nitRR
      if(doplot)call histo_fill(p,sNNLO,npartNNLO,wgtpl)
 555  continue
c
c     double-unresolved counterterm
      if(is_NNLO_2unres_singular_sec(isec,jsec,ksec,lsec))then
         call local_pure_2unres_counter_NNLO(isec,jsec,ksec,lsec,iref,
     &   sNNLO,p,sNLO,pb,sLO,pbb,iA2,wgt,WsumSS,1d0,1d0,
     &   xjac,KSS,KHHCC,K2NNLO,ierr)
         if(ierr.eq.1)goto 999
         call local_mixed_2unres_counter_NNLO(isec,jsec,ksec,lsec,iref,
     &   sNNLO,p,sNLO,pb,sLO,pbb,iU1,iS1,iB1,iU2,iS2,iB2,wgt,1d0,
     &   1d0,WsumSSC,xjac,KSSC,KHHCCC,K12NNLO,ierr)
         if(ierr.eq.1)goto 999
      endif
c
c     single-unresolved counterterm
      if(is_NNLO_1unres_singular_sec(isec,jsec,ksec,lsec))then
         call local_1unres_counter_NNLO(isec,jsec,ksec,lsec,iref,
     &   sNNLO,p,sNLO,pb,sLO,pbb,wgt,WsumS,WsumC,xjac,KS,KHC,K1NNLO,ierr)
         if(ierr.eq.1)goto 999
      endif
c
c     subtraction
c     all counterterms are multiplied by the same jacobian factor
c     xjac since they all start from the same NNLO phase space point
c     and are mapped down differently. This jacobian is the one for
c     the change of variable passing from the Vegas variables x to the
c     NNLO phase-space ones, so it has to be common to all contributions
      int_double_real=int_double_real_no_cnt-(K1NNLO+K2NNLO+K12NNLO)*xjac
ccccccccccccccccccc
      if(abs(int_double_real).gt.1d10)then
         write(31,*)int_double_real
         write(31,*)sNNLO(1,2),sNNLO(1,3),sNNLO(1,4),sNNLO(2,3),sNNLO(2,4),sNNLO(3,4)
         write(31,*)int_double_real_no_cnt,(K1NNLO+K2NNLO+K12NNLO)*xjac,
     &              K1NNLO*xjac,K2NNLO*xjac,K12NNLO*xjac
         write(31,*)Wsum,WsumC,WsumSS,WsumSSC
         write(31,*)
      endif
ccccccccccccccccccc
c
c     print out current run progress
 999  ievt=ievt+1
      if(ievt.gt.nthres)then
         write(*,111)char(13),int(1d2*nthres/(nprodRR*1d0)),'% done'
         nthres=nthres+int(nprodRR/rfactRR)
      endif
 111  format(a1,i3,a6,$)
c
      return
      end
