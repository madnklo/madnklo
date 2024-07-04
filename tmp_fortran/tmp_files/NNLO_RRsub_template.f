      double precision function int_double_real_%(isec)d_%(jsec)d(x,wgt)
c     (n+1)-body NLO integrand for vegas
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
      double precision Z_NLO
      double precision alphaZ
      parameter(alphaZ=1d0)
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
      integer iS1,iB1,iA1,iU2,iS2,iB2,iA2
      integer isec,jsec,iref
      common/cNLOsecindices/isec,jsec
      common/cNLOmaplabels/iU,iS,iB,iA,iref
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision ptilde(0:3,nexternal-2)
      double precision xjac,xjacB
      double precision xsave(3)
      double precision sCM
      common/cscm/sCM
      common/cxsave/xsave
      integer counter
      save counter
      integer nitr
      common/iterations/nitr
      integer %(NLO_proc_str)sfl_factor 
      common/%(NLO_proc_str)sflavour_factor/%(NLO_proc_str)sfl_factor
      double precision dummy_ans(0:1),ans(0:1) !TODO SET CORRECTLY RANGE OF ANS 
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
      integer ich
      common/comich/ich
      double precision  amp2(n_max_cg)
      common/to_amp2/amp2
c     TODO: convert to partonic sCM 
      sCM = (2d0*EBEAM(1))**2
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
c     
c     initialise
      xjac = 0d0
      xjacB = 0d0
      isec = %(isec)d
      jsec = %(jsec)d
      iref = %(iref)d
      int_double_real_%(isec)d_%(jsec)d=0d0
      int_double_real_no_cnt=0d0
      Z_NLO=0d0
      RNNLO=0d0
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

      call configs_%(strUB)s
      call props_%(strUB)s
      call decaybw_%(strUB)s
      call getleshouche_%(strUB)s



      call phase_space_npt(x,sCM,iU1,iS1,iB1,iA1,iA2,p,pbar,ptilde,xjac,xjacB,iU2,iS2,iB2)
      if(xjac.eq.0d0.or.xjacB.eq.0d0) then
         write(77,*) 'int_double_real: '
         write(77,*) 'Jacobians = 0 in phase space ', xjac, xjacB
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
      call %(strUB)s_ME_ACCESSOR_HOOK(PB,HEL,ALPHAS,dummy_ANS)
      WGT_CHAN=AMP2(ICH)
c
c     possible cuts
      IF(DOCUT(P,NEXTERNAL,leg_pdgs,1))GOTO 555
c
c     test matrix elements
      if(ntested.lt.ntest)then
         ntested=ntested+1
         call test_R_%(isec)d_%(jsec)d(iunit,x)
      endif
c     TODO: implement flag 'test_only' to stop here
c
c     real
      call %(NLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      RNNLO = ANS(0) * %(NLO_proc_str)sfl_factor
      if(RNNLO.lt.0d0.or.abs(RNNLO).ge.huge(1d0).or.isnan(RNNLO))then
         write(77,*) 'int_double_real: '
         write(77,*) 'Wrong RNNLO', RNNLO
         goto 999
      endif
c
c     real sector function
      call get_Z_NLO(sNLO,sCM,alphaZ,isec,jsec,Z_NLO,ierr)
      if(ierr.eq.1)then
         write(77,*) 'int_real: '
         write(77,*) 'Wrong Z_NLO', Z_NLO
         goto 999
      endif
c
c     full real in the combination of sectors
      int_double_real_no_cnt=RNNLO*Z_NLO*xjac
c
c     plot real
      wgtpl=int_double_real_no_cnt*wgt/nitR*wgt_chan
      if(doplot)call histo_fill(p,sNNLO,nexternal,wgtpl)
 555  continue
c
      %(str_int_real)s
c
c     counterterm
      call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,xjac,xjacB,x,KNLO,wgt_chan,ierr)
      if(ierr.eq.1)then
         write(77,*) 'int_real: '
         write(77,*) 'Something wrong in the counterterm', KNLO
         goto 999
      endif
c
c     subtraction (phase-space jacobian included in counterterm definition)
      int_double_real_%(isec)d_%(jsec)d=int_double_real_no_cnt-KNNLO
      int_double_real_%(isec)d_%(jsec)d = int_double_real_%(isec)d_%(jsec)d*wgt_chan
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
