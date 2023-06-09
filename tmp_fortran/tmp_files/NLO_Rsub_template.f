      double precision function int_real_%(isec)d_%(jsec)d(x,wgt)
c     (n+1)-body NLO integrand for vegas
      implicit none
      include 'math.inc'
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer i
      integer ierr
      integer ievt,nthres,ntest
      integer iunit
      common/ciunitNLO/iunit
      integer ntested
      parameter(ntest=20)
      save ievt,nthres,ntested
      double precision int_real_no_cnt
      double precision sNLO(nexternal,nexternal),sminNLO
      double precision sLO(nexternal-1,nexternal-1)
      double precision Z_NLO,ZsumSi,ZsumSj
      double precision alpha
      parameter(alpha=1d0)
      double precision RNLO,KNLO,KS,KHC
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgtpl
      logical dotechcut
      double precision tinycut
      logical doplot
      common/cdoplot/doplot
      logical docut
      integer isec,jsec,iU,iS,iB,iA,iref
      common/cNLOsecindices/isec,jsec
      common/cNLOmaplabels/iU,iS,iB,iA,iref
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision xjac,xjacB
      double precision xsave(3)
      double precision sCM
      common/cscm/sCM
      common/cxsave/xsave
      integer counter
      save counter
      integer nitr
      common/iterations/nitr
      double precision ans(0:1) !TODO SET CORRECTLY RANGE OF ANS 
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
c     TODO: convert to partonic sCM 
      sCM = (2d0*EBEAM(1))**2
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
c     
c     initialise
      xjac = 0d0
      xjacB = 0d0
      isec = %(isec)d
      jsec = %(jsec)d
      iref = %(iref)d
      int_real_%(isec)d_%(jsec)d=0d0
      int_real_no_cnt=0d0
      ZsumSi=0d0
      ZsumSj=0d0
      Z_NLO=0d0
      do i=1,3
         xsave(i)=x(i)
      enddo
c
c     specify phase-space mapping
      %(mapping_str)s

      if(isec.le.2.or.jsec.le.2)then
         write(*,*)'update sCM in int_real'
         stop
      endif
c
c     phase space and invariants
      if(sCM.le.0d0)then
         write(*,*) 'Wrong sCM', sCM
         stop
      endif
      call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)
      if(xjac.eq.0d0.or.xjacB.eq.0d0)goto 999
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
      if(dotechcut(snlo,nexternal,tinycut)) goto 999
c     TODO: look at dotechcut
c      if(dotechcut(x,nexternal,tinycut))goto 999
c
c     possible cuts
      if(docut(p,nexternal))goto 999
c
c     test matrix elements
      if(ntested.le.ntest)then
         ntested=ntested+1
         call test_R_%(isec)d_%(jsec)d(iunit,x)
      endif
c     TODO: implement flag 'test_only' to stop here
c
c     real
      call %(NLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      RNLO = ANS(0)
      if(RNLO.lt.0d0.or.abs(RNLO).ge.huge(1d0).or.isnan(RNLO))goto 999
c
c     real sector function
      call get_Z_NLO(sNLO,sCM,alpha,isec,jsec,Z_NLO,'F',ierr)
      if(ierr.eq.1)goto 999
c
c     full real in the combination of sectors
      int_real_no_cnt=RNLO*Z_NLO*xjac

c     plot real
      wgtpl=int_real_no_cnt*wgt/nitR/2D0/SCM
      if(doplot)call histo_fill(p,sNLO,nexternal,wgtpl)
c 555  continue
c
      %(str_int_real)s
c
c     counterterm
      call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,ZsumSi,ZsumSj,xjac,KS,KHC,KNLO,ierr)
      if(ierr.eq.1)goto 999
c
c     subtraction
      int_real_%(isec)d_%(jsec)d=int_real_no_cnt-KNLO*xjac
c     add flux factor
c     TODO: add the general case for the flux factor
      int_real_%(isec)d_%(jsec)d = int_real_%(isec)d_%(jsec)d/2d0/sCM
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
 999  return
      end
