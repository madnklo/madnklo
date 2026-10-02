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
      double precision ans(0:1) !TODO SET CORRECTLY RANGE OF ANS 
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
      integer ich
      common/comich/ich
      integer iconfig,mincfig,maxcfig,invar
      common/cfig/iconfig,mincfig,maxcfig,invar
      double precision dot
      integer NGRAPHS2
      double precision amp2(N_MAX_CG)
      COMMON/TO_AMP2/AMP2,NGRAPHS2
      integer last_ich
      data last_ich /-1/
      save last_ich
c
c     call initialisation function, once per channel
      if (ich.ne.last_ich) then
         call initialise_born_channel()
         last_ich=ich
      endif
c
c     TODO: muR from card
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
c
c     initialise
      xjac=Gevtopb
      int_Born=0d0
c
      call gen_mom(iconfig,mincfig,maxcfig,invar,xjac,x,p,nexternal)
      if(xjac.eq.0d0)goto 999
c
c     possible cuts
      if(docut(p,nexternal,leg_pdgs,0))goto 999
c
c     Born
      call ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BLO = ANS(0)*AMP2(ich)
      if(BLO.lt.0d0.or.abs(BLO).ge.huge(1d0).or.isnan(BLO))goto 999
      int_Born=BLO*xjac*fl_factor
c
c     plot Born
      wgtpl=int_Born*wgt
      wgts=wgtpl
      if(doplot) then
         call invariants_from_p(p,nexternal,sLO,ierr)
         if(ierr.eq.1)goto 999
         call analysis_fill(p,slo,nexternal,leg_pdgs,wgts)
      endif
c
 999  return
      end


      subroutine initialise_born_channel()
      implicit none
      integer ich
      common/comich/ich
      integer iconfig,mincfig,maxcfig,invar
      common/cfig/iconfig,mincfig,maxcfig,invar
c
      iconfig = ich
      mincfig = 1
      maxcfig = 1
      invar = 1
      call configs_born
      call props_born
      call decaybw_born
      call getleshouche_born
c
      return
      end
