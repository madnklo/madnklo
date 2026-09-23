      FUNCTION INT_DOUBLE_VIRTUAL(X,WGT)
C     n-body NNLO integrand for vegas
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'ngraphs.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'nsquaredSO.inc'
      INCLUDE 'leg_PDGs.inc'
      INTEGER NDIM,IERR,IEVT,NTHRES,I
      SAVE IEVT,NTHRES
      INTEGER NITVV
      COMMON/NITERATIONSV/NITVV
      DOUBLE PRECISION INT_VV,VVNNLO(-4:0)
      DOUBLE PRECISION I2NNLO(-4:0),IRVNNLO(-4:0)
      DOUBLE PRECISION SLO(NEXTERNAL,NEXTERNAL)
C     TODO: understand x(mxdim) definition by Vegas
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      DOUBLE PRECISION WGT,WGTS(1),WGTPL
      LOGICAL DOPLOT, DOCUT
      COMMON/CDOPLOT/DOPLOT
      DOUBLE PRECISION P(0:3,NEXTERNAL)
      DOUBLE PRECISION XJAC
      DOUBLE PRECISION SCM
      INTEGER FL_FACTOR
      COMMON/FLAVOUR_FACTOR/FL_FACTOR
      DOUBLE PRECISION ANS(0:1)  !TODO SET CORRECTLY RANGE OF ANS 
      DOUBLE PRECISION ALPHAS, ALPHA_QCD
      INTEGER, PARAMETER :: HEL=-1
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


      integer NGRAPHS2
      double precision amp2(N_MAX_CG)
      COMMON/TO_AMP2/AMP2,NGRAPHS2
      INTEGER ICH
      COMMON/COMICH/ICH
      DOUBLE PRECISION PMASS(NEXTERNAL)
      INTEGER NCOLORCORRELATORS
      PARAMETER (NCOLORCORRELATORS=4)
C     
C     Index 0 is the number of correlators to consider and the next
C     indices are which one to consider
      INTEGER COLOR_CORRELATORS_TO_CONSIDER(0:NCOLORCORRELATORS)
      REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS, 0:3
     $ ,0:NSQUAREDSO)
      COMMON/%(long_proc_prefix)sCOLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER
     $ ,COLOR_CORRELATED_EVALS

      INCLUDE 'pmass.inc'
C     
C     EXTERNAL
C
c     TODO: convert to partonic sCM 
      sCM = (2d0*EBEAM(1))**2
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
C     
C     initialise
      XJAC = GEVTOPB
      INT_VV = 0D0
      SLO = 0D0
      VVNNLO = 0D0
      I2NNLO = 0D0
      IRVNNLO = 0D0
C     
C     BEGIN CODE
C     
      IF (INIT) THEN
        INIT=.FALSE.
c$$$        CALL %(long_proc_prefix)sGET_ANSWER_DIMENSION(MATELEM_ARRAY_DIM)
c$$$        ALLOCATE(MATELEM(0:3,0:MATELEM_ARRAY_DIM))
c$$$        CALL %(long_proc_prefix)sGET_NSQSO_LOOP(NSQUAREDSO_LOOP)
c$$$        ALLOCATE(PREC_FOUND(0:NSQUAREDSO_LOOP))
      ENDIF
C     
C     phase space and invariants
      IF(SCM.LE.0D0)THEN
        WRITE(*,*) 'Wrong sCM', SCM
        STOP
      ENDIF
C     Hard coded settings for gen_mom
      iconfig = ich
      mincfig = 1
      maxcfig = 1
      invar = 2
      call configs_born
      call props_born
      call decaybw_born
      call getleshouche_born
      call gen_mom(iconfig,mincfig,maxcfig,invar,xjac,x,p,nexternal)
      IF(XJAC.EQ.0D0) THEN
        WRITE(77,*)'Wrong jacobian in NNLO_VV'
        GOTO 999
      ENDIF

      CALL INVARIANTS_FROM_P(P,NEXTERNAL,SLO,IERR)
      IF(IERR.EQ.1) THEN
        WRITE(77,*)'Wrong invariants in NNLO_VV', SLO
        GOTO 999
      ENDIF
C     
C     possible cuts
      IF(DOCUT(P,NEXTERNAL,LEG_PDGS,0))GOTO 999
C     
C     call virtual
c$$$      COLOR_CORRELATED_EVALS = 0D0
c$$$      CALL V_ML5_1_1_SLOOPMATRIX_THRES(P,MATELEM,-1.0D0,PREC_FOUND
c$$$     $ ,RETURNCODE)
c$$$      VVNNLO(-4:0) = [(MATELEM(I,0), I=5,1,-1)]
      DO I=-4,0
         IF(ABS(VVNNLO(I)).GE.HUGE(1D0).OR.ISNAN(VVNNLO(I)))THEN
            WRITE(77,*) 'int_VV: '
            WRITE(77,*) 'Wrong VVNNLO at eps^',I,' : ', VVNNLO(I)
            GOTO 999
         ENDIF
      ENDDO
C     
C     call counterterm
      CALL INT_COUNTER_I2_NNLO(P,SLO,I2NNLO,IERR)
      IF(IERR.EQ.1)GOTO 999
      CALL INT_COUNTER_IRV_NNLO(P,SLO,IRVNNLO,IERR)
      IF(IERR.EQ.1)GOTO 999
C     
C     test coefficients of epsilon poles
      IF(NTESTED.LT.NTEST)THEN
        NTESTED=NTESTED+1
        WRITE(50,*)
        WRITE(50,*)'Testing point # ', NTESTED
        WRITE(50,*)'Quadruple pole VV, I2, IRV, sum', VVNLO(-4),
     $    I2NNLO(-4), IRVNNLO(-4), VVNLO(-4) + I2NNLO(-4) + IRVNNLO(-4)
        WRITE(50,*)'Triple    pole VV, I2, IRV, sum', VVNLO(-3),
     $    I2NNLO(-3), IRVNNLO(-3), VVNLO(-3) + I2NNLO(-3) + IRVNNLO(-3)
        WRITE(50,*)'Double    pole VV, I2, IRV, sum', VVNLO(-2),
     $    I2NNLO(-2), IRVNNLO(-2), VVNLO(-2) + I2NNLO(-2) + IRVNNLO(-2)
        WRITE(50,*)'Single    pole VV, I2, IRV, sum', VVNLO(-1),
     $    I2NNLO(-1), IRVNNLO(-1), VVNLO(-1) + I2NNLO(-1) + IRVNNLO(-1)
        WRITE(50,*)
      ENDIF
C     
C     subtracted VV
      INT_VV=(VVNNLO(0)+I2NNLO(0)+IRVNNLO(0))*XJAC
C     
C     apply flavour multiplicity factor
      INT_VV=INT_VV*FL_FACTOR
C     Multi channeling
      INT_VV = INT_VV * AMP2(ICH)
C     
C     plot
      WGTPL=INT_VV*WGT
      WGTS=WGTPL
      IF(DOPLOT)CALL ANALYSIS_FILL(P,SLO,NEXTERNAL,LEG_PDGS,WGTS)
C     
C     print out current run progress
C     999  ievt=ievt+1
C     if(ievt.gt.nthres)then
C     write(*,111)char(13),int(1d2*nthres/(nprodVV*1d0)),' done'
C     nthres=nthres+int(nprodVV/rfactVV)
C     endif
C     111  format(a1,i3,a6,$)
C     
 999  RETURN
      END
