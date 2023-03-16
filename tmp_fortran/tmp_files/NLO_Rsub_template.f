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
      double precision Z_NLO,ZsumSi,ZsumSj
      double precision alpha
      parameter(alpha=1d0)
      double precision RNLO,KNLO,KS,KHC
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
      double precision sCM
      common/cxsave/xsave
      integer counter
      save counter
c     
c     initialise
      isec = %(isec)d
      jsec = %(jsec)d
      ZsumSi=0d0
      ZsumSj=0d0
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
c     TODO: pass sCM information
      if(sCM.le.0d0)then
         write(*,*) 'Wrong sCM', sCM
         stop
      endif
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
c      if(dotechcut(x,nexternal,tinycut))goto 999
c
c     possible cuts
c      if(docut(p,nexternal))goto 555
c
c     TODO: reinstate test routine
c     test matrix elements
c      if(ntested.le.ntest)then
c         call test_R(iunit,x,alpha,beta)
c         ntested=ntested+1
c      endif
c
c     real
c     TODO: look at real_NLO
      call real_NLO(sNLO,RNLO,ierr)
      if(ierr.eq.1)goto 999
c
c     TODO: introduce symmetrised sectors 
      call get_Z_NLO(sNLO,sCM,alpha,isec,jsec,Z_NLO,'F',ierr)
      if(ierr.eq.1)goto 999
c
c     full real in the combination of sectors
      int_real_no_cnt=RNLO*Z_NLO*xjac
c
c     plot real
      wgtpl=int_real_no_cnt*wgt/nitR
c      if(doplot)call histo_fill(p,sNLO,npartNLO,wgtpl)
c 555  continue
c
c     call sector functions
      %(str_int_real)s
c
c     counterterm
      call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,iA,wgt,ZsumSi,ZsumSj,xjac,KS,KHC,KNLO,ierr)
c
c     subtraction
      int_real=int_real_no_cnt-KNLO*xjac
c
c     print out current run progress
c     TODO: adapt progress bar
c 999  ievt=ievt+1
c      if(ievt.gt.nthres)then
c         write(*,111)char(13),int(1d2*nthres/(nprodR*1d0)),' done'
c         nthres=nthres+int(nprodR/rfactR)
c      endif
c 111  format(a1,i3,a6,$)
c
      return
      end
