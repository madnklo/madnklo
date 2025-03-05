      double precision function int_real_%(isec)d_%(jsec)d(x,wgt)
c     (n+1)-body NLO integrand for vegas
      use sectors2_module
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
      integer ievt,nthres,ntest,ntested
      save ievt,nthres,ntested
      parameter(ntest=20)
      integer iunit
      common/ciunitNLO/iunit
      double precision int_real_no_cnt
      double precision sNLO(nexternal,nexternal),sminNLO
      double precision sLO(nexternal-1,nexternal-1)
      double precision alphaZ
      parameter(alphaZ=1d0)
      double precision RNLO,KNLO,KS,KHC
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgtpl,wgt_chan
      logical dotechcut
      double precision tinycut
      logical doplot
      common/cdoplot/doplot
      logical docut
      integer iU,iS,iB,iA,iref
      integer isec,jsec,ksec,lsec
      common/csecindices/isec,jsec,ksec,lsec
      common/cNLOmaplabels/iU,iS,iB,iA,iref
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision xjac,xjacB
      double precision sCM
      common/cscm/sCM
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
      ksec = 0d0
      lsec = 0d0
      iref = %(iref)d
      int_real_%(isec)d_%(jsec)d=0d0
      int_real_no_cnt=0d0
      RNLO=0d0
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
c
c     test IR limits matrix elements
c     TODO: test_only option
      if(ntested.lt.ntest)then
         ntested=ntested+1
         call test_R_%(isec)d_%(jsec)d(iunit,ntested)
      endif
c      
c     call phase space
      call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)
      if(xjac.eq.0d0.or.xjacB.eq.0d0) then
         write(77,*) 'int_real: '
         write(77,*) 'Jacobians = 0 in phase space ', xjac, xjacB
         goto 999
      endif
      call invariants_from_p(p,nexternal,sNLO,ierr)
      if(ierr.eq.1) then
         write(77,*) 'int_real: '
         write(77,*) 'Wrong NLO invariants ', sNLO
         goto 999
      endif
      call invariants_from_p(pb,nexternal-1,sLO,ierr)  
      if(ierr.eq.1) then
         write(77,*) 'int_real: '
         write(77,*) 'Wrong LO invariants ', sLO
         goto 999
      endif
      call get_sig2(SNLO,alphaZ,nexternal)
c
c     tiny technical phase-space cut to avoid fluctuations
      tinycut=tiny1
      if(dotechcut(snlo,nexternal,tinycut)) goto 999
C
c     Call the Underlying Born matrix element to fill the amp2 array,
c     in order to implement the multi channel
      call %(strUB)s_ME_ACCESSOR_HOOK(PB,HEL,ALPHAS,dummy_ANS)
      WGT_CHAN=AMP2(ICH)
c
c     possible cuts
      IF(DOCUT(P,NEXTERNAL,leg_pdgs,1))GOTO 555
c
c     real
      call %(NLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      RNLO = ANS(0) * %(NLO_proc_str)sfl_factor
      if(RNLO.lt.0d0.or.abs(RNLO).ge.huge(1d0).or.isnan(RNLO))then
         write(77,*) 'int_real: '
         write(77,*) 'Wrong RNLO', RNLO
         goto 999
      endif
c
c     real sector function
      call get_Z_NLO(isec,jsec)
      if(ierr.eq.1)then
         write(77,*) 'int_real: '
         write(77,*) 'Wrong Z_NLO', Z_NLO
         goto 999
      endif
c
c     full real in the combination of sectors
      int_real_no_cnt=RNLO*Z_NLO*xjac
c
c     plot real
      wgtpl=int_real_no_cnt*wgt/nitR*wgt_chan
      if(doplot)call histo_fill(p,sNLO,nexternal,wgtpl)
 555  continue
c
      %(str_int_real)s
c
c     counterterm
      call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,xjac,xjacB,KNLO,wgt_chan,ierr)
      if(ierr.eq.1)then
         write(77,*) 'int_real: '
         write(77,*) 'Something wrong in the counterterm', KNLO
         goto 999
      endif
c
c     subtraction (phase-space jacobian included in counterterm definition)
      int_real_%(isec)d_%(jsec)d=int_real_no_cnt-KNLO
      int_real_%(isec)d_%(jsec)d = int_real_%(isec)d_%(jsec)d*wgt_chan
c
 999  return
      end
