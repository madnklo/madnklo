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
      integer ievt,nthres,ntest
      integer iunit
      common/ciunitNLO/iunit
      integer ntested
      parameter(ntest=20)
      save ievt,nthres,ntested
      double precision int_real_no_cnt
      double precision sNLO(nexternal,nexternal)
      double precision sLO(nexternal-1,nexternal-1)
      double precision alphaZ
      parameter(alphaZ=1d0)
      double precision RNLO,KNLO
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgtpl,wgts(1),wgt_chan
      logical dotechcut
      logical doplot
      common/cdoplot/doplot
      logical docut
      integer iU,iS,iB,iA
      common/cNLOmaplabels/iU,iS,iB,iA
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision xjac,xjacB
      double precision xsave(3)
      double precision sCM
      common/cscm/sCM
      common/cxsave/xsave
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
      logical firsttime
      data firsttime/.true./
      save firsttime
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
c
c     call initialisation function
      if(firsttime)then
         call initialise_sector()
         firsttime=.false.
      endif
c
c     TODO: convert to partonic sCM 
      sCM = (2d0*EBEAM(1))**2
      if(sCM.le.0d0)then
         write(*,*) 'Wrong sCM', sCM
         stop
      endif
c
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
c     
c     initialise
      xjac = 0d0
      xjacB = 0d0
      int_real_%(isec)d_%(jsec)d=0d0
      int_real_no_cnt=0d0
      RNLO=0d0
      xsave(1:3)=x(1:3)
      WGT_CHAN=1d0
c
c     phase space and invariants
      call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)
      if(xjac*xjacB.eq.0d0) then
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
c
c     tiny technical phase-space cut to avoid fluctuations
      if(dotechcut(snlo,nexternal,tiny1)) goto 999
c
c     possible cuts
      IF(DOCUT(P,NEXTERNAL,leg_pdgs,1))GOTO 555
c
c     test phase-space singularities of matrix elements
      if(ntested.lt.ntest)then
         ntested=ntested+1
         call test_R_%(isec)d_%(jsec)d(iunit,x)
      endif
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
      call get_sig2(SNLO,alphaZ,nexternal)
      call get_W_NLO(isec,jsec)
      if(ierr.eq.1)then
         write(77,*) 'int_real: '
         write(77,*) 'Wrong W_NLO', W_NLO
         goto 999
      endif
c
c     full real in sector Wij
      int_real_no_cnt=RNLO*W_NLO*xjac
c
c     plot real
      wgtpl=int_real_no_cnt*wgt/nitR*wgt_chan
c      if(doplot)call histo_fill(p,sNLO,nexternal,leg_pdgs,wgtpl)
      wgts=wgtpl
      if(doplot)call analysis_fill(p,sNLO,nexternal,leg_pdgs,wgts)
 555  continue
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
      int_real_%(isec)d_%(jsec)d=int_real_no_cnt-KNLO
      int_real_%(isec)d_%(jsec)d = int_real_%(isec)d_%(jsec)d*wgt_chan
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



      subroutine initialise_sector()
      implicit none
      include 'nexternal.inc'
      include 'leg_PDGs.inc'
      integer iU,iS,iB,iA
      common/cNLOmaplabels/iU,iS,iB,iA
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
c
      isec = %(isec)d
      jsec = %(jsec)d
      ksec = 0
      lsec = 0
      iref = %(iref)d
c
c     check we are not in the ISR case
      if(isec.le.2.or.jsec.le.2)then
         write(*,*)'ISR indices',isec,jsec
         stop
      endif
c
c     specify phase-space mapping
      %(mapping_str)s
c
c     configuration files
      call configs_%(strUB)s
      call props_%(strUB)s
      call decaybw_%(strUB)s
      call getleshouche_%(strUB)s
c
c     fill underlying pdgs, labels and flavours
      call get_underlying_pdgs(isec,jsec,ksec,lsec,nexternal-1,underlying_leg_pdgs)
      call get_mapped_labels(nexternal,isec,jsec,leg_pdgs,underlying_leg_pdgs,mapped_labels)
      return
      end
