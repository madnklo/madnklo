      function int_real_%(isec)d_%(jsec)d(x,wgt)
c     (n+1)-body NLO integrand for vegas
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      include 'nexternal.inc'
      integer ndim,ierr
      integer i,j,l,m,ievt,nthres,ntest
      integer iunit
      common/ciunitNLO/iunit
      integer ntested(-2:maxdim,-2:maxdim)
      parameter(ntest=30)
      save ievt,nthres,ntested
      double precision int_real,int_real_no_cnt
      double precision sNLO(npartNLO,npartNLO),sminNLO
      double precision sLO(npartLO,npartLO)
      double precision W_NLO(maxdim,maxdim)
      double precision Wsum,WsumSi,WsumSj
      double precision WS_NLO(maxdim,maxdim),WC_NLO(maxdim,maxdim)
      double precision alpha,beta
      parameter(alpha=1d0,beta=1d0)
      double precision RNLO,KNLO,KS,KHC,RNLO_MG
      double precision x(mxdim)
      double precision wgt,wgtpl
      common/cdim/ndim
      logical dotechcut
      double precision tinycut
      logical doplot
      common/cdoplot/doplot
      integer isec,jsec,iU,iS,iB,iA,iref
      common/cNLOsecindices/isec,jsec
      common/cNLOmaplabels/iU,iS,iB,iA,iref
      double precision p(0:3,npartNLO)
      double precision pMG(0:3,npartNLO)
      double precision pb(0:3,npartLO)
      double precision xjac,xjacB
      double precision xsave(3)
      common/cxsave/xsave
      logical R_from_MG
      parameter(R_from_MG=.false.)
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
c     TODO: changed momenta convention!!!
c     phase space and invariants
      call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)
      if(xjac.eq.0d0)goto 999
      call invariants_from_p(p,npartNLO,sNLO,ierr)
      if(ierr.eq.1)goto 999
      call invariants_from_p(pb,npartLO,sLO,ierr)
      if(ierr.eq.1)goto 999
c
c     check that phase-space labels and x variables are as expected
      call check_phsp_consistency(x,npartNLO,sNLO,sLO,iU,iS,iB,iA,ierr)
      if(ierr.eq.1)goto 999
c
c     tiny technical phase-space cut to avoid fluctuations
      if(.not.is_NLO_singular_sec(isec,jsec))then
         tinycut=0d0
      else
         tinycut=tiny1
      endif
      if(dotechcut(x,npartNLO,tinycut))goto 999
c
c     possible cuts
      if(docut(p,npartNLO))goto 555
c
c     test matrix elements
      if(is_NLO_singular_sec(isec,jsec).and.ntested(isec,jsec).le.ntest)then
         call test_R(iunit,x,alpha,beta)
         ntested(isec,jsec)=ntested(isec,jsec)+1
      endif
c
c     real
      if(.not.R_from_MG)then
         call real_NLO(sNLO,RNLO,ierr)
         if(ierr.eq.1)goto 999
      else
c=========================================
c     CALLS TO MG5 MATRIX ELEMENTS
c=========================================

         call smatrix(pMG,RNLO_MG)
c$$$      write(*,*)
c$$$      write(*,*)RNLO,' in-house'
c$$$      write(*,*)RNLO_MG*Gevtopb/(2d0*sCM)*(alphas/0.118d0),' MG5'
c$$$      write(*,*)
c$$$      counter=counter+1
c$$$      if(counter.eq.11)stop
         RNLO=RNLO_MG*Gevtopb/(2d0*sCM)*(alphas/0.118d0)
      endif
c
c     add function from python
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
      if(doplot)call histo_fill(p,sNLO,npartNLO,wgtpl)
 555  continue
c
c     counterterm
      call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,iA,
     &wgt,WsumSi,WsumSj,xjac,KS,KHC,KNLO,ierr)
c
c     subtraction
      int_real=int_real_no_cnt-KNLO*xjac
c
c     print out current run progress
 999  ievt=ievt+1
      if(ievt.gt.nthres)then
         write(*,111)char(13),int(1d2*nthres/(nprodR*1d0)),'% done'
         nthres=nthres+int(nprodR/rfactR)
      endif
 111  format(a1,i3,a6,$)
c
      return
      end
