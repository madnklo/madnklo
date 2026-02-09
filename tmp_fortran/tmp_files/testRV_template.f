      subroutine test_RV_%(isec)d_%(jsec)d(iunit,x0)
      implicit none
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer isec,jsec,ksec,lsec
      common/csecindices/isec,jsec,ksec,lsec
      integer i,iU,iS,iB,iA,iref
      common/cNLOmaplabels/iU,iS,iB,iA,iref
      integer iUtmp,iStmp
      integer iunit,ievnt
      INTEGER, PARAMETER :: MXDIM = 30
      double precision x0(mxdim)
      character*10 dash10
      save ievnt
      double precision xsave(3)
      common/cxsave/xsave
      double precision e(2), l(2)
c
      dash10='----------'
      ievnt=ievnt+1
c
      write(iunit,*)dash10//dash10//dash10//dash10
      write(iunit,*)dash10//dash10//dash10//dash10
      write(iunit,*)' EVENT NUMBER ',ievnt
      write(iunit,*)dash10//dash10//dash10//dash10
      write(iunit,*)dash10//dash10//dash10//dash10
%(limit_str)s
c
c     reinstate original xsave after testing
      do i=1,3
         xsave(i)=x0(i)
      enddo
c
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)
c
      return
      end


      subroutine do_limit_RV_%(isec)d_%(jsec)d(iunit,limstr,x0,e,l)
      use sectors2_module
      implicit none
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'leg_PDGs.inc'
      INCLUDE 'ngraphs.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'nsquaredSO.inc'
      integer iitn,i,j,maxitn,iunit,ierr
      integer isec,jsec,ksec,lsec
      common/csecindices/isec,jsec,ksec,lsec
      integer iU,iS,iB,iA,iref
      common/cnlomaplabels/iU,iS,iB,iA,iref
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      integer, parameter :: mxdim=30
      parameter(maxitn=12)
      double precision x0(mxdim),x(mxdim)
      double precision sNLO(nexternal,nexternal)
      double precision sLO(nexternal-1,nexternal-1)
      double precision KRVNNLO(-2:0)
      double precision lam,lim,RVNNLO(-2:0),single_real
      character*5 str5
      character*8 limstr
      character*10 str10
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision xjac,xjacB
      double precision xsave(3)
      DOUBLE PRECISION ANS(0:1) !TODO SET CORRECTLY RANGE OF ANS
      DOUBLE PRECISION ALPHAS, ALPHA_QCD
      DOUBLE PRECISION WGT,WGTPL,wgt_chan
      DOUBLE PRECISION SCM
      INTEGER, PARAMETER :: HEL=-1
      integer %(NNLO_RV_proc_str)sfl_factor
      common/%(NNLO_RV_proc_str)sflavour_factor/%(NNLO_RV_proc_str)sfl_factor
      DOUBLE PRECISION ALPHAZ
      PARAMETER(ALPHAZ=1D0)
      common/cxsave/xsave
      double precision e(2),l(2)
C
C     block for virtual matel
C
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

      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
      SCM = (2D0*EBEAM(1))**2
c
c     initialise
      x=x0
      str5 ='     '
      str10='          '
      xjac=0d0
      sNLO=0d0
      sLO=0d0
      wgt_chan=1d0
c
c     TODO: MAP SOFT LIMIT AS (ilm), I.E. ONE MAPPING PER DIPOLE
c
c     start testing
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)'LIM = '//trim(limstr)
      write(iunit,*)str10//'lambda'//str10//str10//'R'//str10//str10//str5//'LIM'//str10//str10//'|R-LIM|/|LIM|'
c
c     possibility to set by hand the starting point
c     for the limiting procedure
c      x0(1)=0.5d0
c      x0(2)=0.5d0
c
c     loop to get closer and closer to the limit
      do iitn=1,maxitn
         lam=10d0**(1d0-iitn)
c
c     initialise
         KRVNNLO=0d0
c
c     rescale relevant x random numbers
c     x(1) is zCS, while x(2) is yCS
c     TODO: this rescaling is specific for (ijr) mapping; generalise
         x(1)=abs(l(1)-x0(1)*lam**e(1))
         x(2)=abs(l(2)-x0(2)*lam**e(2))
c
c     set xsave so that the counterterms will be called with
c     more and more singular kinematics
         do i=1,3
            xsave(i)=x(i)
         enddo
c
c     recompute momenta after rescaling
         call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB,mapped_labels)
         if(xjac.eq.0d0.or.xjacb.eq.0d0)cycle
         call invariants_from_p(p,nexternal,sNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(pb,nexternal-1,sLO,ierr)
         if(ierr.eq.1)cycle
c
c     real virtual
      IF (INIT) THEN
         INIT=.FALSE.
         CALL %(long_proc_prefix)sGET_ANSWER_DIMENSION(MATELEM_ARRAY_DIM)
         ALLOCATE(MATELEM(0:3,0:MATELEM_ARRAY_DIM))
         CALL %(long_proc_prefix)sGET_NSQSO_LOOP(NSQUAREDSO_LOOP)
         ALLOCATE(PREC_FOUND(0:NSQUAREDSO_LOOP))
      ENDIF
      RVNNLO(-2:0) = MATELEM(1:3,0) * %(NNLO_RV_proc_str)sfl_factor
      if(abs(RVNNLO(0)).ge.huge(1d0).or.isnan(RVNNLO(0)))cycle
c
      call get_sig2(SNLO,alphaZ,nexternal)
      call get_W_NLO(isec,jsec)
c
c     counterterm
         call local_RV_counter_NNLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,xjac,xjacB,x,KRVNNLO,wgt_chan,ierr)
         if(ierr.eq.1)cycle

         lim=KRVNNLO(0)
         single_real=RVNNLO(0)*W_NLO*xjac

         if(abs(lim).gt.0d0)then
            write(iunit,*)lam,single_real,lim,abs(single_real-lim)/abs(lim)
         else
            write(iunit,*)lam,single_real,lim,single_real,' *** '
         endif
      enddo
      x=x0
c
      return
      end
