      double precision function int_double_real_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(x,wgt)
c     (n+2)-body NNLO integrand for vegas
      USE SECTORS4_MODULE
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'nexternal.inc'
      include 'input.inc'
      include 'run.inc'
      include 'cuts.inc'
      include 'leg_PDGs.inc'
      include 'ngraphs_%(UBgraphs)s.inc'
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
c      double precision W_NNLO
      double precision RNNLO,KNNLO
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgts(1),wgtpl,wgt_chan
      logical dotechcut
      logical doplot
      common/cdoplot/doplot
      logical docut
      integer iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
      common/cNNLOmaplabels/iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
      integer isec,jsec,ksec,lsec,iref
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer asec,bsec,csec,dsec
      common/csecindices/asec,bsec,csec,dsec
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision ptilde(0:3,nexternal-2)
      double precision xjac,xjacB,xjacCS1
      double precision xsave(3)
      double precision sCM
      common/cscm/sCM
      common/cxsave/xsave
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
      logical firsttime
      data firsttime/.true./
      save firsttime
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
      logical consistency_check
      common/cconscheck/consistency_check
c
C     call initialisation function
      if(firsttime)then
        call initialise_sector()
        firsttime=.false.
      endif
c
c     TODO: convert to partonic sCM
      sCM = (2d0*EBEAM(1))**2
      IF(SCM.LE.0D0)THEN
        WRITE(*,*) 'Wrong sCM', SCM
        STOP
      ENDIF
c
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
c
c     initialise
      xjac = 0d0
      xjacB = 0d0
      int_double_real_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d=0d0
      int_double_real_no_cnt=0d0
      RNNLO=0d0
      xsave(1:3)=x(1:3)
      WGT_CHAN=1D0
c
c     phase space and invariants
      call phase_space_npt(x,sCM,iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,p,pb,ptilde,xjac,xjacB,xjacCS1)
      if(xjac*xjacB*xjacCS1.eq.0d0) then
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
      if(dotechcut(snnlo,nexternal,tiny1)) goto 999
c
c     possible cuts
      IF(DOCUT(P,NEXTERNAL,leg_pdgs,2))GOTO 555
c
c     test phase-space singularities of matrix elements
      test_sector_function = .false.
      consistency_check = .false.
      if(ntested.lt.ntest)then
         ntested=ntested+1
         call test_RR_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(iunit,x)
      endif
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
      call get_sigNNLO(snnlo,nexternal)
      call get_W_NNLO(asec,bsec,csec,dsec)
c
c     full double real in sector Wijkl
      int_double_real_no_cnt=RNNLO*W_NNLO*xjac
c
c     plot double real
      wgtpl=int_double_real_no_cnt*wgt/nitR*wgt_chan
c      if(doplot)call histo_fill(p,sNNLO,nexternal,leg_pdgs,wgtpl)
      wgts=wgtpl
      if(doplot) call analysis_fill(p,sNNLO,nexternal,leg_pdgs,wgts)
 555  continue
c
c     counterterm
      call local_counter_NNLO_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(sNNLO,p,sNLO,pb,sLO,ptilde,wgt,xjac,xjacB,x,KNNLO,wgt_chan,ierr)
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


      subroutine initialise_sector()
      implicit none
      include 'nexternal.inc'
      include 'leg_PDGs.inc'
      integer iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
      integer i,j,i1,i2,ib3,ib4
      common/cNNLOmaplabels/iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
      integer isec,jsec,ksec,lsec,iref
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer asec,bsec,csec,dsec
      common/csecindices/asec,bsec,csec,dsec
      integer map1,map2
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
      integer real_ss_mapped_labels(nexternal),Born_ss_mapped_labels(nexternal-1)
      common/c_NNLO_ss_mapped_labels/real_ss_mapped_labels,Born_ss_mapped_labels
      integer real_sc_mapped_labels(nexternal),Born_sc1_mapped_labels(nexternal-1),Born_sc2_mapped_labels(nexternal-1)
      common/c_NNLO_sc_mapped_labels/real_sc_mapped_labels,Born_sc1_mapped_labels,Born_sc2_mapped_labels
      integer real_s_sc_1_mapped_labels(nexternal),real_s_sc_2_mapped_labels(nexternal),Born_s_sc_1_mapped_labels(nexternal-1),Born_s_sc_2_mapped_labels(nexternal-1)
      common/c_nnlo_s_sc_mapped_labels/real_s_sc_1_mapped_labels,real_s_sc_2_mapped_labels,Born_s_sc_1_mapped_labels,Born_s_sc_2_mapped_labels
      integer Born_2_mapped_labels(nexternal-1)
      common/c_NNLO_2_mapped_labels/Born_2_mapped_labels
      integer real_hcc_ia_mapped_labels(nexternal),Born_hcc_ia_mapped_labels(nexternal-1)
      common/c_NNLO_hcc_mapped_labels/real_hcc_ia_mapped_labels,Born_hcc_ia_mapped_labels
      integer real_hcc_ib_mapped_labels(nexternal),Born_hcc_ib_mapped_labels(nexternal-1)
      common/c_NNLO_hcc_mapped_labels/real_hcc_ib_mapped_labels,Born_hcc_ib_mapped_labels

      include 'all_sector_list.inc'
C
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
c
c     check we are not in the ISR case
      if(asec.le.2.or.bsec.le.2.or.csec.le.2.or.dsec.le.2)then
         write(*,*)'ISR indices',asec,bsec,csec,dsec
      endif
c
c     specify phase-space mapping
      %(mapping_str)s
c
c     configuration files
      call configs_%(str_UBorn)s
      call props_%(str_UBorn)s
      call decaybw_%(str_UBorn)s
      call getleshouche_%(str_UBorn)s
c
c     fill underlying pdgs, labels and flavours
      call get_underlying_pdgs(asec,bsec,csec,dsec,nexternal-1,real_leg_pdgs)
      call get_mapped_labels(nexternal,asec,bsec,leg_pdgs,real_leg_pdgs,real_mapped_labels)
c     for mapped n+1 -> n mapped labels:
c     if lsec =0 the unresolved pair is jsec, ksec,
c     if lsec!=0 the unresolved pair is ksec, lsec.
      call get_underlying_pdgs(asec,bsec,csec,dsec,nexternal-2,Born_leg_pdgs)
      map1=real_mapped_labels(ksec)
      map2=real_mapped_labels(jsec)
      if(lsec.ne.0)map2=real_mapped_labels(lsec)
      call get_mapped_labels(nexternal-1,map1,map2,real_leg_pdgs,Born_leg_pdgs,Born_mapped_labels)
c
c     ad-hoc Born mapped labels (first appeared in S_SS)
c     TODO: checked for IJJK only
c     TODO: generalise construction of mapped labels
      map1=real_mapped_labels(iref)
      map2=real_mapped_labels(jsec)
      call get_mapped_labels(nexternal-1,map1,map2,real_leg_pdgs,Born_leg_pdgs,born_2_mapped_labels)
c
c     fill mapped labels for double-soft mapping
c     for ijkj and ijkl if pdg(i)+pdg(k)=0
      if((bsec.ne.csec).and.(leg_pdgs(asec)+leg_pdgs(ksec).eq.0)) then
         call get_mapped_labels(nexternal,asec,csec,leg_pdgs,real_leg_pdgs,real_ss_mapped_labels)
         map1=real_ss_mapped_labels(ksec)
         map2=real_ss_mapped_labels(jsec)
         if(lsec.ne.0)map2=real_ss_mapped_labels(lsec)
         call get_mapped_labels(nexternal-1,map1,map2,real_leg_pdgs,Born_leg_pdgs,Born_ss_mapped_labels)
      endif
c
c     fill mapped labels for soft-collinear mapping
      call get_mapped_labels(nexternal,csec,dsec,leg_pdgs,real_leg_pdgs,real_sc_mapped_labels)
      map1=real_sc_mapped_labels(isec)
      map2=real_sc_mapped_labels(ksec)
      call get_mapped_labels(nexternal-1,map1,map2,real_leg_pdgs,Born_leg_pdgs,Born_sc1_mapped_labels)
      map1=real_sc_mapped_labels(isec)
      map2=real_sc_mapped_labels(jsec)
      call get_mapped_labels(nexternal-1,map1,map2,real_leg_pdgs,Born_leg_pdgs,Born_sc2_mapped_labels)
c
c     fill mapped labels for soft soft-collinear mapping
      call get_mapped_labels(nexternal,asec,csec,leg_pdgs,real_leg_pdgs,real_s_sc_1_mapped_labels)
      map1=real_s_sc_1_mapped_labels(dsec)
      map2=real_s_sc_1_mapped_labels(csec)
      call get_mapped_labels(nexternal-1,map1,map2,real_leg_pdgs,Born_leg_pdgs,Born_s_sc_1_mapped_labels)

      call get_mapped_labels(nexternal,asec,dsec,leg_pdgs,real_leg_pdgs,real_s_sc_2_mapped_labels)
      map1=real_s_sc_2_mapped_labels(csec)
      map2=real_s_sc_2_mapped_labels(dsec)
      call get_mapped_labels(nexternal-1,map1,map2,real_leg_pdgs,Born_leg_pdgs,Born_s_sc_2_mapped_labels)
c
c     fill mapped labels for soft/double-soft hard double-collinear mapping
      call get_mapped_labels(nexternal,asec,csec,leg_pdgs,real_leg_pdgs,real_hcc_ia_mapped_labels)
      map1=real_hcc_ia_mapped_labels(dsec)
      map2=real_hcc_ia_mapped_labels(csec)
      call get_mapped_labels(nexternal-1,map1,map2,real_leg_pdgs,Born_leg_pdgs,Born_hcc_ia_mapped_labels)

      call get_mapped_labels(nexternal,asec,dsec,leg_pdgs,real_leg_pdgs,real_hcc_ib_mapped_labels)
      map1=real_hcc_ib_mapped_labels(csec)
      map2=real_hcc_ib_mapped_labels(dsec)
      call get_mapped_labels(nexternal-1,map1,map2,real_leg_pdgs,Born_leg_pdgs,Born_hcc_ib_mapped_labels)
c
c     fill bar_indices for barred sector functions
      j=1
      do i=1,lensectors
         i1=all_sector_list(1,i)
         i2=all_sector_list(2,i)
         if(i1.eq.asec.and.i2.eq.bsec) then
            ib3=real_mapped_labels(all_sector_list(3,i))
            ib4=real_mapped_labels(all_sector_list(4,i))
            bar_indices(j) = ib3
            bar_indices(j+1) = ib4
         endif
         j=j+2
      enddo

      return
      end
