c      double precision function int_Born_multich(x,ich,wgt)
      double precision function int_Born(x,wgt)
c     n-body LO integrand for vegas
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'leg_PDGs.inc'
      
      integer ich
      integer ierr
      integer ievt,nthres
      save ievt,nthres
      double precision sLO(nexternal,nexternal),sminLO
      double precision BLO
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgtpl
      logical doplot
      common/cdoplot/doplot
      logical docut
      integer nitB
      common/iterations/nitB
      integer fl_factor 
      common/flavour_factor/fl_factor
      double precision p(0:3,nexternal)
      double precision xjac
      double precision sCM
      double precision ans(0:1) !TODO SET CORRECTLY RANGE OF ANS 
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
      integer iconfig,mincfig,maxcfig,invar
      double precision dot
      integer NGRAPHS2
      double precision amp2(8)
      COMMON/TO_AMP2/AMP2,NGRAPHS2
      common/comich/ich

      double precision pmass(nexternal)
      include 'pmass.inc'


      
c     TODO: convert to partonic sCM 
      sCM = (2d0*EBEAM(1))**2
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
c
c     initialise
      xjac=Gevtopb
      int_Born=0d0
c
c     phase space and invariants
      if(sCM.le.0d0)then
         write(*,*) 'Wrong sCM', sCM
         stop
      endif
C     Hard coded settings for gen_mom
c     TODO: At the moment the variables mincfig,maxcfig,invar seem no to be used
c      Check if we actually need them!
      iconfig = ich
c      iconfig = 1
      mincfig = 1
      maxcfig = 1
      invar = 1
      call gen_mom(iconfig,mincfig,maxcfig,invar,xjac,x,p,nexternal)
c      write(*,*) p(:,1),dot(p(:,1),p(:,1)),pmass(1)**2
c      write(*,*) p(:,2),dot(p(:,2),p(:,2)),pmass(2)**2
c      write(*,*) p(:,3),dot(p(:,3),p(:,3)),pmass(3)**2
c      write(*,*) p(:,4),dot(p(:,4),p(:,4)),pmass(4)**2
c      write(*,*) p(:,5),dot(p(:,5),p(:,5)),pmass(5)**2
c      pause
!      call phase_space_n(x,sCM,p,nexternal,xjac)
      if(xjac.eq.0d0)goto 999
c      call invariants_from_p(p,nexternal,sLO,ierr)
c      if(ierr.eq.1)goto 999
c
c     possible cuts
      if(docut(p,nexternal,leg_pdgs,0))goto 999
c
c     Born
      call ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
c      BLO = ANS(0)*AMP2(ich)
c      BLO = ANS(0)

      if(BLO.lt.0d0.or.abs(BLO).ge.huge(1d0).or.isnan(BLO))goto 999
c      int_Born=BLO*xjac
c     add flux factor
C      int_Born = int_Born/2d0/sCM
c     apply flavour factor
      int_Born=int_Born * fl_factor
      int_Born = xjac
c
c     plot
      wgtpl=int_Born*wgt/nitB
      if(doplot)call histo_fill(p,sLO,nexternal,wgtpl)
c
c     print out current run progress
c 999  ievt=ievt+1
c      if(ievt.gt.nthres)then
c         write(*,111)char(13),int(1d2*nthres/(nprodB*1d0)),' done'
c         nthres=nthres+int(nprodB/rfactB)
c      endif
c 111  format(a1,i3,a6,$)
c
c
 999  return
      end