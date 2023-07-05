      function int_virtual(x,wgt)
c     n-body NLO integrand for vegas
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'ngraphs.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'nsquaredSO.inc'
      INCLUDE 'leg_PDGs.inc'
      integer ndim,ierr,ievt,nthres,i
      save ievt,nthres
      integer nitV
      common/niterationsv/nitV
      double precision int_virtual,VNLO(3),INLO(3)
      double precision sLO(nexternal,nexternal)
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgtpl
      logical doplot, docut
      common/cdoplot/doplot
      double precision p(0:3,nexternal)
      double precision xjac
      double precision sCM
      integer fl_factor 
      common/flavour_factor/fl_factor
      double precision ans(0:1) !TODO SET CORRECTLY RANGE OF ANS 
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
      LOGICAL INIT
      DATA INIT/.TRUE./
      COMMON/INITCHECKSA/INIT
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
      integer ntested,ntest
      parameter(ntest=20)
      save ntested
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
C     EXTERNAL
C
c     TODO: convert to partonic sCM 
      sCM = (2d0*EBEAM(1))**2
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
c
c     initialise
      xjac = Gevtopb
      int_virtual = 0d0
      sLO = 0d0
      VNLO = 0d0
      INLO = 0d0
C
C     BEGIN CODE
C
      IF (INIT) THEN
        INIT=.FALSE.
        CALL %(long_proc_prefix)sGET_ANSWER_DIMENSION(MATELEM_ARRAY_DIM)
        ALLOCATE(MATELEM(0:3,0:MATELEM_ARRAY_DIM))
        CALL %(long_proc_prefix)sGET_NSQSO_LOOP(NSQUAREDSO_LOOP)
        ALLOCATE(PREC_FOUND(0:NSQUAREDSO_LOOP))
!        INCLUDE 'pmass.inc'
      ENDIF
c
c     phase space and invariants
      if(sCM.le.0d0)then
         write(*,*) 'Wrong sCM', sCM
         stop
      endif
C     Hard coded settings for gen_mom
      iconfig = 1
      mincfig = 1
      maxcfig = 1
      invar = 2
      call gen_mom(iconfig,mincfig,maxcfig,invar,xjac,x,p,nexternal)

!      call phase_space_n(x,sCM,p,nexternal,xjac)
      if(xjac.eq.0d0) then
         write(*,*)'Wrong jacobian in NLO_V'
         goto 999
      endif
      call invariants_from_p(p,nexternal,sLO,ierr)
      if(ierr.eq.1) then
         write(*,*)'Wrong invariants in NLO_V', sLO
         goto 999
      endif
c
c     possible cuts
      if(docut(p,nexternal,leg_pdgs,0))goto 999
c
c     call virtual
      CALL %(long_proc_prefix)sSLOOPMATRIX_THRES(p,MATELEM,-1.0D0,PREC_FOUND,RETURNCODE)
      VNLO(1:3) = MATELEM(1:3,0)
c
c     call Born
      CALL ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BLO = ANS(0)
c
c     call counterterm
      call int_counter_NLO(p,sLO,INLO,ierr)
      if(ierr.eq.1)goto 999
c
c     test coefficients of epsilon poles
      if(ntested.lt.ntest)then
         ntested=ntested+1
         write(50,*)
         write(50,*)'Testing point # ', ntested
         write(50,*)'Double pole V, I, sum', VNLO(3), INLO(3), VNLO(3)+INLO(3)
         write(50,*)'Single pole V, I, sum', VNLO(2), INLO(2), VNLO(2)+INLO(2)
         write(50,*)
      endif
c
c     subtracted vrtual
      int_virtual=(VNLO(1)+INLO(1))*xjac
!      int_virtual = int_virtual/2d0/sCM
c
c     apply flavour multiplicity factor
      int_virtual=int_virtual*fl_factor
c
c     plot
      wgtpl=int_virtual*wgt/nitV
      if(doplot)call histo_fill(p,sLO,nexternal,wgtpl)
c
c     print out current run progress
c     999  ievt=ievt+1
c      if(ievt.gt.nthres)then
c         write(*,111)char(13),int(1d2*nthres/(nprodV*1d0)),' done'
c         nthres=nthres+int(nprodV/rfactV)
c      endif
c 111  format(a1,i3,a6,$)
c
 999  return
      end
