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
      INCLUDE 'ngraphs.inc'
      
      integer ich
      integer ierr
      integer ievt,nthres
      save ievt,nthres
      double precision sLO(nexternal,nexternal),sminLO
      double precision BLO
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgts(1),wgtpl
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
      double precision amp2(N_MAX_CG)
      COMMON/TO_AMP2/AMP2,NGRAPHS2
      common/comich/ich
      double precision mass2
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
      call configs_born
      call props_born
      call decaybw_born
      call getleshouche_born
      call gen_mom(iconfig,mincfig,maxcfig,invar,xjac,x,p,nexternal)
      if(xjac.eq.0d0)goto 999
      call invariants_from_p(p,nexternal,sLO,ierr)
      if(ierr.eq.1)goto 999
c
c     possible cuts
      if(docut(p,nexternal,leg_pdgs,0))goto 999
c
      
c     Born
      call ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BLO = ANS(0)*AMP2(ich)

      if(BLO.lt.0d0.or.abs(BLO).ge.huge(1d0).or.isnan(BLO))goto 999
      int_Born=BLO*xjac
c     apply flavour factor
      int_Born=int_Born * fl_factor
c     plot
      wgtpl=int_Born*wgt/nitB
      wgts(1)=wgtpl
c      if(doplot)call histo_fill(p,sLO,nexternal,wgtpl)
      if(doplot)call analysis_fill(p,sLO,nexternal,wgts)
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
