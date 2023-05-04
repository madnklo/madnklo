      function int_virtual(x,wgt)
c     n-body NLO integrand for vegas
      implicit none
      include 'math.inc'
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'ngraphs.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'nsquaredSO.inc'
      integer ndim,ierr,ievt,nthres,i
      save ievt,nthres
      double precision int_virtual,VNLO,INLO,VNLO_MG
      double precision sLO(nexternal,nexternal)
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgtpl
      logical doplot
      common/cdoplot/doplot
      double precision p(0:3,nexternal)
      double precision xjac
      double precision sCM
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
      REAL * 8 finite_part, single_pole, double_pole, diffeps0, diffeps1, diffeps2
      REAL * 8 finite_part_torino, single_pole_torino, double_pole_torino
      REAL * 8 Q2
C     
C     GLOBAL VARIABLES
C     
C     This is from ML code for the list of split orders selected by
C     the process definition
C     
      INTEGER NLOOPCHOSEN
      CHARACTER*20 CHOSEN_LOOP_SO_INDICES(NSQUAREDSO)
      LOGICAL CHOSEN_LOOP_SO_CONFIGS(NSQUAREDSO)
      COMMON/ML5_1_2_CHOSEN_LOOP_SQSO/CHOSEN_LOOP_SO_CONFIGS
      INTEGER NBORNCHOSEN
      CHARACTER*20 CHOSEN_BORN_SO_INDICES(NSQSO_BORN)
      LOGICAL CHOSEN_BORN_SO_CONFIGS(NSQSO_BORN)
      COMMON/ML5_1_2_CHOSEN_BORN_SQSO/CHOSEN_BORN_SO_CONFIGS

C     
C     EXTERNAL
C     
c     TODO: convert to partonic sCM 
      sCM = (2d0*EBEAM(1))**2
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)

c
c     initialise
      xjac=0d0
      int_virtual=0d0
      sLO=0d0
      VNLO = 0d0
      INLO = 0d0
c

C     
C     BEGIN CODE
C     
C     


      IF (INIT) THEN
        INIT=.FALSE.
        CALL ML5_1_2_GET_ANSWER_DIMENSION(MATELEM_ARRAY_DIM)
        ALLOCATE(MATELEM(0:3,0:MATELEM_ARRAY_DIM))
        CALL ML5_1_2_GET_NSQSO_LOOP(NSQUAREDSO_LOOP)
        ALLOCATE(PREC_FOUND(0:NSQUAREDSO_LOOP))
!        INCLUDE 'pmass.inc'
      ENDIF


c     phase space and invariants
      if(sCM.le.0d0)then
         write(*,*) 'Wrong sCM', sCM
         stop
      endif
      call phase_space_n(x,sCM,p,nexternal,xjac)
      if(xjac.eq.0d0)goto 999
      call invariants_from_p(p,nexternal,sLO,ierr)
      if(ierr.eq.1)goto 999
c
c     possible cuts
c      if(docut(p,npartLO))goto 999
c
c     virtual
C      call virtual_NLO(sLO,VNLO,0,ierr)
      CALL ML5_1_2_SLOOPMATRIX_THRES(p,MATELEM,-1.0D0,PREC_FOUND
     $   ,RETURNCODE)

      
      finite_part = MATELEM(1,0)
      single_pole = MATELEM(2,0)
      double_pole = MATELEM(3,0)

      Q2 = sCM ! Ellis-Exton scale, to be found!

!     Attempt to convert from ML convention for the virtual to the Torino paper one

      
      diffeps1 = -(double_pole*EulerGamma) + 
     -  double_pole*(dlog((4d0*MU_R**2*Pi)/Q2) - dlog(MU_R**2/sCM))
      

      single_pole_torino = single_pole + diffeps1

      diffeps0 = (6d0*double_pole*EulerGamma**2 - 
     -    12d0*EulerGamma*single_pole - 
     -    12d0*(double_pole*EulerGamma - single_pole)*
     -     dlog((4d0*MU_R**2*Pi)/Q2) + 
     -    6d0*double_pole*
     -     dlog((4d0*MU_R**2*Pi)/Q2)**2 - 
     -    12d0*single_pole_torino*dlog(MU_R**2/sCM) - 
     -     6d0*double_pole*dlog(MU_R**2/sCM)**2)/12d0

      finite_part_torino = finite_part + diffeps0
      
c$$$      diffeps1 = -eulergamma+dlog((4d0*MU_R**2*Pi)/Q2)-dlog(MU_R**2/sCM)
c$$$      
c$$$      diffeps2 =  (eulergamma**2+dlog(4d0*Pi)**2 - 
c$$$     -    2d0*eulergamma*dlog((4d0*MU_R**2*Pi)/Q2) + 
c$$$     -    dlog(MU_R**2/Q2)*
c$$$     -     dlog((16d0*MU_R**2*Pi**2)/Q2) - 
c$$$     -    dlog(MU_R**2/sCM)**2)/2d0
      
C      VNLO = VNLO - (single_pole * diffeps1 + double_pole * diffeps2)

      VNLO = finite_part_torino
      
      if(ierr.eq.1)goto 999
c
c     counterterm
      call int_counter_NLO(p,sLO,INLO,ierr)
      if(ierr.eq.1)goto 999
c
c     subtraction
      int_virtual=(VNLO+INLO)*xjac
      int_virtual = int_virtual/2d0/sCM
c
c     plot
c$$$      wgtpl=int_virtual*wgt/nitV
c$$$      if(doplot)call histo_fill(p,sLO,npartLO,wgtpl)
c$$$c
c$$$c     print out current run progress
c$$$ 999  ievt=ievt+1
c$$$      if(ievt.gt.nthres)then
c$$$         write(*,111)char(13),int(1d2*nthres/(nprodV*1d0)),'% done'
c$$$         nthres=nthres+int(nprodV/rfactV)
c$$$      endif
c$$$ 111  format(a1,i3,a6,$)
c$$$c
 999  return
      end
