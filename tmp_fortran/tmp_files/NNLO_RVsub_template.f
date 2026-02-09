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
      integer i
      integer ierr
      integer ievt,nthres,ntest
      integer iunit
      common/ciunitNLO/iunit
      integer ntested
      parameter(ntest=20)
      save ievt,nthres,ntested
      double precision int_real_virtual_no_cnt
      double precision sNLO(nexternal,nexternal)
      double precision sLO(nexternal-1,nexternal-1)
      double precision alphaZ
      parameter(alphaZ=1d0)
      double precision RVNNLO(-2:0),KRVNNLO(-2:0)
      double precision I1NNLO(-2:0),I12NNLO(-2:0)
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgtpl,wgt_chan
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
      integer nitRV
      common/niterationsrv/nitRV
      integer %(NNLO_RV_proc_str)sfl_factor
      common/%(NNLO_RV_proc_str)sflavour_factor/%(NNLO_RV_proc_str)sfl_factor
      double precision dummy_ans(0:1),ans(0:1) !TODO SET CORRECTLY RANGE OF ANS
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
      integer ich
      common/comich/ich
      integer ngraphs2
      double precision amp2(n_max_cg)
      common/to_amp2/amp2,ngraphs2
      logical firsttime
      data firsttime/.true./
      save firsttime
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
C
C     block for virtual matel
C
      LOGICAL INIT
      DATA INIT/.TRUE./
C      COMMON/INITCHECKSA/INIT
      INTEGER MATELEM_ARRAY_DIM
      REAL*8 , ALLOCATABLE :: MATELEM(:,:)
      REAL*8 SQRTS,AO2PI,TOTMASS
C     sqrt(s)= center of mass energy
      REAL*8 PIN(0:3), POUT(0:3)
      CHARACTER*120 BUFF(NEXTERNAL)
      INTEGER RETURNCODE, UNITS, TENS, HUNDREDS
      INTEGER NSQUAREDSO_LOOP
      REAL*8 , ALLOCATABLE :: PREC_FOUND(:)
      REAL*8 BLO
C
C     GLOBAL VARIABLES
C
C     This is from ML code for the list of split orders selected by
C     the process definition
C
      INTEGER NLOOPCHOSEN
      CHARACTER*20 CHOSEN_LOOP_SO_INDICES(NSQUAREDSO)
      LOGICAL CHOSEN_LOOP_SO_CONFIGS(NSQUAREDSO)
      COMMON/%(long_proc_prefix)sCHOSEN_LOOP_SQSO/CHOSEN_LOOP_SO_CONFIGS
      INTEGER NBORNCHOSEN
      CHARACTER*20 CHOSEN_BORN_SO_INDICES(NSQSO_BORN)
      LOGICAL CHOSEN_BORN_SO_CONFIGS(NSQSO_BORN)
      COMMON/%(long_proc_prefix)sCHOSEN_BORN_SQSO/CHOSEN_BORN_SO_CONFIGS
      integer iconfig,mincfig,maxcfig,invar
C
      double precision pmass(nexternal)
      INCLUDE 'pmass.inc'

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
      iref = %(iref)d
      int_real_virtual_%(isec)d_%(jsec)d=0d0
      int_real_virtual_no_cnt=0d0
      RVNNLO  = 0d0
      KRVNNLO = 0d0
      I1NNLO  = 0d0
      I12NNLO = 0d0
c      pb_to_ME=0d0
      xsave(1:3)=x(1:3)
      WGT_CHAN=1d0
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
      if(docut(p,nexternal,leg_pdgs,0))goto 555
c
c     Call the Underlying Born matrix element to fill the amp2 array,
c     in order to implement the multi channel
c      call %(strUB)s_ME_ACCESSOR_HOOK(PB,HEL,ALPHAS,dummy_ANS)
c      WGT_CHAN=AMP2(ICH)
c
c     test phase-space singularities of matrix elements
      if(ntested.lt.ntest)then
         ntested=ntested+1
         call test_RV_%(isec)d_%(jsec)d(iunit,x)
      endif
c     TODO: implement flag 'test_only' to stop here
c
c     real virtual
      IF (INIT) THEN
         INIT=.FALSE.
         CALL %(long_proc_prefix)sGET_ANSWER_DIMENSION(MATELEM_ARRAY_DIM)
         ALLOCATE(MATELEM(0:3,0:MATELEM_ARRAY_DIM))
         CALL %(long_proc_prefix)sGET_NSQSO_LOOP(NSQUAREDSO_LOOP)
         ALLOCATE(PREC_FOUND(0:NSQUAREDSO_LOOP))
      ENDIF
c
      CALL ML5_1_1_SLOOPMATRIX_THRES(P,MATELEM,-1.0D0,PREC_FOUND,RETURNCODE)
      RVNNLO(-2:0) = MATELEM(1:3,0) * %(NNLO_RV_proc_str)sfl_factor
      do i=-2,0
         if(abs(RVNNLO(i)).ge.huge(1d0).or.isnan(RVNNLO(i)))then
            write(77,*) 'int_real_virtual: '
            write(77,*) 'Wrong RVNNLO at eps^',i,' : ', RVNNLO(i)
            goto 999
         endif
      enddo
c
c     real sector function
      call get_sig2(SNLO,alphaZ,nexternal)
      call get_W_NLO(isec,jsec)
      if(ierr.eq.1)then
         write(77,*) 'int_real_virtual: '
         write(77,*) 'Wrong W_NLO', W_NLO
         goto 999
      endif
c
c     full real virtual in sector Wij
      int_real_virtual_no_cnt=RVNNLO(0)*W_NLO*xjac
c
c     plot real virtual
      wgtpl=int_real_virtual_no_cnt*wgt/nitRV*wgt_chan
      if(doplot)call histo_fill(p,sNLO,nexternal,leg_pdgs,wgtpl)
 555  continue
c
c
c     real-virtual counterterm
      call local_RV_counter_NNLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,xjac,xjacB,x,KRVNNLO,wgt_chan,ierr)
      if(ierr.eq.1)then
         write(77,*) 'int_real_virtual: '
         write(77,*) 'Something wrong in the RV counterterm', KRVNNLO
         goto 999
      endif
c
c     1-unresolved integrated counterterm
      call int_counter_I1_NNLO_%(isec)d_%(jsec)d(p,sNLO,sLO,I1NNLO,ierr)
      if(ierr.eq.1)goto 999
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
      int_real_virtual_%(isec)d_%(jsec)d=(int_real_virtual_no_cnt+I1NNLO(0)) - (KRVNNLO(0)+I12NNLO(0))
      int_real_virtual_%(isec)d_%(jsec)d = int_real_virtual_%(isec)d_%(jsec)d*wgt_chan
c
c     print out current run progress
c     TODO: adapt progress bar
c     999  ievt=ievt+1
c      if(ievt.gt.nthres)then
c         write(*,111)char(13),int(1d2*nthres/(nprodRV*1d0)),' done'
c         nthres=nthres+int(nprodRV/rfactRV)
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
c      integer mapped_labels_s(nexternal)
c      integer mapped_flavours_s(nexternal-1),mapped_indices_shuff_s(nexternal-1)
c      common/c_mapped_quantities_s/mapped_labels_s,mapped_flavours_s,mapped_indices_shuff_s
c      integer mapped_labels_c(nexternal)
c      integer mapped_flavours_c(nexternal-1),mapped_indices_shuff_c(nexternal-1)
c      common/c_mapped_quantities_c/mapped_labels_c,mapped_flavours_c,mapped_indices_shuff_c
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
