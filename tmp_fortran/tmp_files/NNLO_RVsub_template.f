      double precision function int_real_virtual_%(isec)d_%(jsec)d(x,wgt)
c     (n+1)-body NNLO integrand for vegas
      use sectors2_module
      implicit none
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'leg_PDGs.inc'
      INCLUDE 'ngraphs_%(UBgraphs)s.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'nsquaredSO.inc'
      integer i,ierr
      integer ntest
      parameter(ntest=20)
      integer ntested
      save ntested
      integer iunit
      common/ciunitNLO/iunit
      double precision sNLO(nexternal,nexternal)
      double precision sLO(nexternal-1,nexternal-1)
      double precision RVNNLO(-2:0),KRVNNLO(-2:0)
      double precision I1NNLO(-2:0),I12NNLO(-2:0)
      double precision int_rv_i1, int_krv_i12
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgts(1),wgtpl,wgt_chan
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
      double precision sCM
      common/cscm/sCM
      integer nitRV
      common/iterations/nitRV
      integer %(NNLO_RV_proc_str)sfl_factor
      common/%(NNLO_RV_proc_str)sflavour_factor/%(NNLO_RV_proc_str)sfl_factor
      double precision dummy_ans(0:1),ans(0:1) !TODO SET CORRECTLY RANGE OF ANS
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
      logical firsttime
      data firsttime/.true./
      save firsttime
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      logical init
      data init/.true./
      integer matelem_array_dim
      real*8 , allocatable :: matelem(:,:)
      integer returncode
      integer nsquaredso_loop
      real*8 , allocatable :: prec_found(:)
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     call initialisation function
      if(firsttime)then
         call initialise_sector()
         call %(NNLO_RV_proc_str)sme_accessor_hook(p,hel,alphas,ans)
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
      int_real_virtual_%(isec)d_%(jsec)d=0d0
      WGT_CHAN=1d0
c
c     phase-space tests (for double pole, single pole, finite part)
      if(ntested.lt.ntest)then
         ntested=ntested+1
         call test_RV_%(isec)d_%(jsec)d(ntested)
      endif
c
c     phase space and invariants
      call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB,mapped_labels)
      if(xjac.eq.0d0.or.xjacB.eq.0d0) then
         write(77,*) 'int_real_virtual: '
         write(77,*) 'Jacobians = 0 in phase space ', xjac, xjacB
         goto 999
      endif
      call invariants_from_p(p,nexternal,sNLO,ierr)
      if(ierr.eq.1) then
         write(77,*) 'int_real_virtual: '
         write(77,*) 'Wrong (n+1)-body invariants ', sNLO
         goto 999
      endif
      call invariants_from_p(pb,nexternal-1,sLO,ierr)
      if(ierr.eq.1) then
         write(77,*) 'int_real_virtual: '
         write(77,*) 'Wrong n-body invariants ', sLO
         goto 999
      endif
c
c     tiny technical phase-space cut to avoid fluctuations
      if(dotechcut(snlo,nexternal,tiny1)) goto 999
c
c     possible cuts
      if(docut(p,nexternal,leg_pdgs,1))goto 555
c
c     real virtual
      if (init) then
         init=.false.
         call %(long_proc_prefix)sget_answer_dimension(matelem_array_dim)
         allocate(matelem(0:3,0:matelem_array_dim))
         call %(long_proc_prefix)sget_nsqso_loop(nsquaredso_loop)
         allocate(prec_found(0:nsquaredso_loop))
      endif
c
      CALL %(long_proc_prefix)sSLOOPMATRIX_THRES(P,MATELEM,-1.0D0,PREC_FOUND,RETURNCODE)
      RVNNLO(-2:0) = [(MATELEM(i,0), i=3,1,-1)]
      do i=-2,0
         if(abs(RVNNLO(i)).ge.huge(1d0).or.isnan(RVNNLO(i)))then
            write(77,*) 'int_real_virtual: '
            write(77,*) 'Wrong RVNNLO at eps^',i,' : ', RVNNLO(i)
            goto 999
         endif
      enddo
c
c     real sector function
      call get_sig2(snlo,nexternal)
      call get_W_NLO(isec,jsec)
      if(ierr.eq.1)then
         write(77,*) 'int_real_virtual: '
         write(77,*) 'Wrong W_NLO', W_NLO
         goto 999
      endif
c
c     1-unresolved integrated counterterm
      call int_counter_I1_NNLO_%(isec)d_%(jsec)d(p,sNLO,sLO,I1NNLO,ierr)
      if(ierr.eq.1)goto 999
c
c     pole-free combination RV+I1
      int_rv_i1 = (RVNNLO(0)+I1NNLO(0))*w_nlo*xjac
      int_rv_i1 = int_rv_i1 * %(NNLO_RV_proc_str)sfl_factor
c
c     plot real virtual
      wgtpl=int_rv_i1*wgt/nitRV*wgt_chan
c     if(doplot)call histo_fill(p,sNLO,nexternal,leg_pdgs,wgtpl)
      wgts=wgtpl
      if(doplot)call analysis_fill(p,sNLO,nexternal,leg_pdgs,wgts)
 555  continue
c
c     real-virtual counterterm
      call local_RV_counter_NNLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,xjac,xjacB,x,KRVNNLO,wgt_chan,ierr)
      if(ierr.eq.1)then
         write(77,*) 'int_real_virtual: '
         write(77,*) 'Something wrong in the RV counterterm', KRVNNLO
         goto 999
      endif
c
c     12-unresolved integrated counterterm
      call int_counter_I12_NNLO_%(isec)d_%(jsec)d(p,sNLO,sLO,I12NNLO,ierr)
      if(ierr.eq.1)goto 999
c
c     test coefficients of epsilon poles
      if(ntested.lt.ntest)then
         ntested=ntested+1
         write(50,*)
         write(50,*)'Testing point # ', ntested
         write(50,*)'Double pole RV, I1, sum', RVNNLO(-2), I1NNLO(-2), RVNNLO(-2)+I1NNLO(-2)
         write(50,*)'Single pole RV, I1, sum', RVNNLO(-1), I1NNLO(-1), RVNNLO(-1)+I1NNLO(-1)
         write(50,*)
         write(50,*)'Double pole KRV, I12, sum', KRVNNLO(-2), I12NNLO(-2), KRVNNLO(-2)+I12NNLO(-2)
         write(50,*)'Single pole KRV, I12, sum', KRVNNLO(-1), I12NNLO(-1), KRVNNLO(-1)+I12NNLO(-1)
         write(50,*)
         write(50,*)
      endif
c
c     subtraction (phase-space jacobian included in counterterm definition)
      int_krv_i12 = KRVNNLO(0)+I12NNLO(0)
      int_real_virtual_%(isec)d_%(jsec)d = int_rv_i1 - int_krv_i12
      int_real_virtual_%(isec)d_%(jsec)d = int_real_virtual_%(isec)d_%(jsec)d*wgt_chan
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
         write(*,*)'update sCM in int_real'
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
