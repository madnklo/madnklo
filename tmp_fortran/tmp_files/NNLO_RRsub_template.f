      double precision function int_double_real_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(x,wgt)
c     (n+2)-body NNLO integrand for vegas
      USE SECTORS4_MODULE
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'leg_PDGs.inc'
      INCLUDE 'ngraphs_%(UBgraphs)s.inc'
      integer i
      integer ierr
      integer ievt,nthres,ntest
      integer iunit
      common/ciunitNLO/iunit
      integer ntested
      parameter(ntest=20)
      save ievt,nthres,ntested
      double precision int_double_real_no_cnt
      double precision sNNLO(nexternal,nexternal)
      double precision sNLO(nexternal-1,nexternal-1),sminNLO
      double precision sLO(nexternal-2,nexternal-2)
c      double precision Z_NNLO
      double precision alphaZ
      parameter(alphaZ=2d0)
      double precision RNNLO,KNNLO
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgtpl,wgt_chan
      logical dotechcut
      double precision tinycut
      logical doplot
      common/cdoplot/doplot
      logical docut
      integer isec,jsec,ksec,lsec,iref
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer asec,bsec,csec,dsec
      common/csecindices/asec,bsec,csec,dsec
      integer iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
      common/cNNLOmaplabels/iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision ptilde(0:3,nexternal-2)
      double precision xjac,xjacB,xjacCS1
      double precision xsave(3)
      double precision sCM
      common/cscm/sCM
      common/cxsave/xsave
      integer counter
      save counter
      integer nitr
      common/iterations/nitr
      integer %(proc_prefix_rr)sfl_factor
      common/%(proc_prefix_rr)sflavour_factor/%(proc_prefix_rr)sfl_factor
      double precision dummy_ans(0:1),ans(0:1) !TODO SET CORRECTLY RANGE OF ANS
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
      integer ich
      common/comich/ich
      double precision  amp2(n_max_cg)
      common/to_amp2/amp2
      double precision K1,K2,K12
c     TODO: convert to partonic sCM
      sCM = (2d0*EBEAM(1))**2
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
c
c     initialise
      xjac = 0d0
      xjacB = 0d0
c     cpartindices: 
c     each sector-relevant particle -> one index
      isec = %(isec)d
      jsec = %(jsec)d
      ksec = %(ksec)d
      lsec = %(lsec)d
      iref = %(iref)d
c     csecindices: 
c     ordered list of particle indices identifying the sector
      asec = %(isec)d
      bsec = %(jsec)d
      csec = %(c3p)d
      dsec = %(d3p)d
      int_double_real_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d=0d0
      int_double_real_no_cnt=0d0
      Z_NNLO=0d0
      RNNLO=0d0
      do i=1,3
         xsave(i)=x(i)
      enddo
c
c     specify phase-space mapping 
      %(mapping_str)s

      if(asec.le.2.or.bsec.le.2.or.csec.le.2.or.dsec.le.2)then
         write(*,*)'ISR: update sCM in int_real'
         stop
      endif
c
c     phase space and invariants
      if(sCM.le.0d0)then
         write(*,*) 'Wrong sCM', sCM
         stop
      endif

      call configs_%(str_UBorn)s
      call props_%(str_UBorn)s
      call decaybw_%(str_UBorn)s
      call getleshouche_%(str_UBorn)s


c     call to phase space
      call phase_space_npt(x,sCM,iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,p,pb,ptilde,xjac,xjacB,xjacCS1)
      if(xjac.eq.0d0.or.xjacB.eq.0d0 .or. xjacCS1 .eq. 0d0) then
         write(77,*) 'int_double_real: '
         write(77,*) 'Jacobians = 0 in phase space ', xjac, xjacB, xjacCS1
         goto 999
      endif
      call invariants_from_p(p,nexternal,sNNLO,ierr)
      if(ierr.eq.1) then
         write(77,*) 'int_double_real: '
         write(77,*) 'Wrong NNLO invariants ', sNNLO
         goto 999
      endif
      call invariants_from_p(pb,nexternal-1,sNLO,ierr)
      if(ierr.eq.1) then
         write(77,*) 'int_double_real: '
         write(77,*) 'Wrong NLO invariants ', sNLO
         goto 999
      endif
      call invariants_from_p(ptilde,nexternal-2,sLO,ierr)
      if(ierr.eq.1) then
         write(77,*) 'int_double_real: '
         write(77,*) 'Wrong LO invariants ', sLO
         goto 999
      endif
c
c     tiny technical phase-space cut to avoid fluctuations
      tinycut=tiny1
      if(dotechcut(snnlo,nexternal,tinycut)) goto 999
C
c     Call the Underlying Born matrix element to fill the amp2 array,
c     in order to implement the multi channel
      call %(str_UBorn)s_ME_ACCESSOR_HOOK(PB,HEL,ALPHAS,dummy_ANS)
      WGT_CHAN=AMP2(ICH)
c
c     possible cuts
      IF(DOCUT(P,NEXTERNAL,leg_pdgs,2))GOTO 555
c
c     test matrix elements
      if(ntested.lt.ntest)then
         ntested=ntested+1
         call test_RR_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(iunit,x)
      endif
c     TODO: implement flag 'test_only' to stop here
c
c     double real
      call %(proc_prefix_rr)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      RNNLO = ANS(0) * %(proc_prefix_rr)sfl_factor
      if(RNNLO.lt.0d0.or.abs(RNNLO).ge.huge(1d0).or.isnan(RNNLO))then
         write(77,*) 'int_double_real: '
         write(77,*) 'Wrong RNNLO', RNNLO
         goto 999
      endif
c
c     double real sector function
c      call  get_Z_NNLO(sNNLO,sCM,alphaZ,asec,bsec,csec,dsec,Z_NNLO,ierr)
      call get_sigNNLO(SNNLO,alphaz,nexternal)
c      call get_Z_NNLO(asec,bsec,csec,dsec).  !!! GB: move to get_W_NNLO

      if(ierr.eq.1)then
         write(77,*) 'int_double_real: '
         write(77,*) 'Wrong Z_NNLO', Z_NNLO
         goto 999
      endif
c
c     full real in the combination of sectors
      int_double_real_no_cnt=RNNLO*Z_NNLO*xjac
c
c     plot real
      wgtpl=int_double_real_no_cnt*wgt/nitR*wgt_chan
      if(doplot)call histo_fill(p,sNNLO,nexternal,leg_pdgs,wgtpl)
 555  continue
c
c     counterterm
      call local_counter_NNLO_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(sNNLO,p,sNLO,pb,sLO,ptilde,wgt,xjac,xjacB,x,K1,K2,K12,KNNLO,wgt_chan,ierr)
      if(ierr.eq.1)then
         write(77,*) 'int_double_real: '
         write(77,*) 'Something wrong in the counterterm', KNNLO
         goto 999
      endif
c
c     subtraction (phase-space jacobian included in counterterm definition)
      int_double_real_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d = int_double_real_no_cnt-KNNLO
      int_double_real_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d = int_double_real_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d*wgt_chan
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
