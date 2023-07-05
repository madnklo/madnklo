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
      iconfig = 1
      mincfig = 1
      maxcfig = 1
      invar = 2
      call gen_mom(iconfig,mincfig,maxcfig,invar,xjac,x,p,nexternal)
!      call phase_space_n(x,sCM,p,nexternal,xjac)
      if(xjac.eq.0d0)goto 999
      call invariants_from_p(p,nexternal,sLO,ierr)
      if(ierr.eq.1)goto 999
c
c     possible cuts
      if(docut(p,nexternal,leg_pdgs,0))goto 999
c
c     Born
      call ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BLO = ANS(0)
      if(BLO.lt.0d0.or.abs(BLO).ge.huge(1d0).or.isnan(BLO))goto 999
      int_Born=BLO*xjac
c     add flux factor
C      int_Born = int_Born/2d0/sCM
c     apply flavour factor
      int_Born = int_Born * fl_factor
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
