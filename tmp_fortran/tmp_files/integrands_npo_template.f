      function int_real_%(isec)d_%(jsec)d(x,wgt)
c     (n+1)-body NLO integrand for vegas
      implicit none
      include 'math.inc'
      include 'nexternal.inc'
      integer ierr
      integer i,j,l,m,ievt,nthres,ntest
      integer iunit
      common/ciunitNLO/iunit
      integer ntested
      parameter(ntest=30)
      save ievt,nthres,ntested
      double precision int_real,int_real_no_cnt
      double precision sNLO(nexternal,nexternal),sminNLO
      double precision sLO(nexternal-1,nexternal-1)
      double precision W_NLO(nexternal,nexternal)
      double precision Wsum,WsumSi,WsumSj
      double precision WS_NLO(nexternal,nexternal),WC_NLO(nexternal,nexternal)
      double precision alpha,beta
      parameter(alpha=1d0,beta=1d0)
      double precision RNLO,KNLO,KS,KHC,RNLO_MG
c     TODO: understand x(mxdim) definition by Vegas
      double precision x(mxdim)
      double precision wgt,wgtpl
      logical dotechcut
      double precision tinycut
      logical doplot
      common/cdoplot/doplot
      integer isec,jsec,iU,iS,iB,iA,iref
      common/cNLOsecindices/isec,jsec
      common/cNLOmaplabels/iU,iS,iB,iA,iref
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision xjac,xjacB
      double precision xsave(3)
      common/cxsave/xsave
      integer counter
      save counter
c
c     initialise
      isec = %(isec)d
      jsec = %(jsec)d
      xjac=0d0
      int_real=0d0
      int_real_no_cnt=0d0
      sNLO=0d0
      sLO=0d0
      KNLO=0d0
      KS=0d0
      KHC=0d0
      do i=1,3
         xsave(i)=x(i)
      enddo
c     
c     phase space and invariants
      call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)
      if(xjac.eq.0d0)goto 999
      call invariants_from_p(p,nexternal,sNLO,ierr)
      if(ierr.eq.1)goto 999
      call invariants_from_p(pb,nexternal-1,sLO,ierr)
      if(ierr.eq.1)goto 999
c
c     check that phase-space labels and x variables are as expected
c      call check_phsp_consistency(x,npartNLO,sNLO,sLO,iU,iS,iB,iA,ierr)
c      if(ierr.eq.1)goto 999
c
c     tiny technical phase-space cut to avoid fluctuations
      tinycut=tiny1
c     TODO: look at dotechcut
      if(dotechcut(x,nexternal,tinycut))goto 999
c
c     possible cuts
c      if(docut(p,nexternal))goto 555
c
c     test matrix elements
      if(ntested.le.ntest)then
         call test_R(iunit,x,alpha,beta)
         ntested=ntested+1
      endif
c
c     real
c     TODO: look at real_NLO
      call real_NLO(sNLO,RNLO,ierr)
      if(ierr.eq.1)goto 999
c
c     TODO: look at get_W_NLO()
      call get_W_NLO(sNLO,alpha,beta,isec,jsec,W_NLO,WS_NLO,WC_NLO,ierr)
      if(ierr.eq.1)goto 999
c
c     combinations of W functions
      Wsum=W_NLO(isec,jsec)+W_NLO(jsec,isec)
      WsumSi=WS_NLO(isec,jsec)
      WsumSj=WS_NLO(jsec,isec)
c
c     full real in the combination of sectors
      int_real_no_cnt=RNLO*Wsum*xjac
c
c     plot real
      wgtpl=int_real_no_cnt*wgt/nitR
c      if(doplot)call histo_fill(p,sNLO,npartNLO,wgtpl)
c 555  continue
c
c     counterterm
c     TODO: adapt to new convention
      call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,iA,
     &wgt,WsumSi,WsumSj,xjac,KS,KHC,KNLO,ierr)
c
c     subtraction
      int_real=int_real_no_cnt-KNLO*xjac
c
c     print out current run progress
c     TODO: adapt progress bar
 999  ievt=ievt+1
      if(ievt.gt.nthres)then
         write(*,111)char(13),int(1d2*nthres/(nprodR*1d0)),'% done'
         nthres=nthres+int(nprodR/rfactR)
      endif
 111  format(a1,i3,a6,$)
c
      return
      end
