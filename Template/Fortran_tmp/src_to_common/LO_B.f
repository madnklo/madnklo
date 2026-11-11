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
      INCLUDE 'qcd_power.inc'
      INCLUDE 'initial_states.inc'
      integer ierr
      integer ievt,nthres
      save ievt,nthres
      double precision sLO(nexternal,nexternal),sminLO
      double precision BLO
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgts(sv_maxwgt),wgtpl
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
      integer nfail
      character*8 failure_stage
      common/cfailure_count/nfail
      common/cfailure_stage/failure_stage
      integer iconfig,mincfig,maxcfig,invar
      common/cfig/iconfig,mincfig,maxcfig,invar
      double precision dot
      integer NGRAPHS2
      double precision amp2(N_MAX_CG)
      COMMON/TO_AMP2/AMP2,NGRAPHS2
      integer last_ich
      data last_ich /-1/
      save last_ich
c     hadron beams (lpp != 0): PDFs and mirror processes (pdf_lumi.f)
      logical hadrons,firsttime
      data firsttime/.true./
      save hadrons,firsttime
      double precision pm(0:3,nexternal),BLOm,muf,lum,lumm,den,num
      double precision rf(10),facf(10)
      integer nvf,ii
      double precision get_muf,pdf_lumi,pdf_swap_initial
      external get_muf,pdf_lumi,pdf_swap_initial
      double precision born_chan_weight
      external born_chan_weight
c
      if(firsttime)then
         hadrons=(lpp(1).ne.0.or.lpp(2).ne.0)
         if(hadrons)call set_tau_min(nexternal,leg_pdgs)
         firsttime=.false.
      endif
c
c     call initialisation function, once per channel
      if (ich.ne.last_ich) then
         call initialise_born_channel()
         last_ich=ich
      endif
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
c     renormalisation scale for the n-body kinematics
      call set_mur_from_momenta(p,nexternal,leg_pdgs)
      alphas=alphas_current
c
c     Born
      call ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BLO = ANS(0)*born_chan_weight(ich)
      if(BLO.lt.0d0.or.abs(BLO).ge.huge(1d0).or.isnan(BLO))goto 999
      if(.not.hadrons)then
         int_Born=BLO*xjac*fl_factor
      else
c        hadron beams: sum over the processes of this directory of
c        f_a(x1) f_b(x2) B_ab(p), plus the mirror processes, with the
c        matrix element evaluated with the initial momenta exchanged
         muf=get_muf(p,nexternal,leg_pdgs)
         q2fact(1)=muf**2
         q2fact(2)=muf**2
         lum=pdf_lumi(xbk(1),xbk(2),muf,1)
         BLOm=0d0
         lumm=0d0
         if(has_mirror)then
            lumm=pdf_lumi(xbk(1),xbk(2),muf,2)*
     &           pdf_swap_initial(p,pm,nexternal)
            call ME_ACCESSOR_HOOK(PM,HEL,ALPHAS,ANS)
            BLOm = ANS(0)*born_chan_weight(ich)
            if(BLOm.lt.0d0.or.abs(BLOm).ge.huge(1d0).or.isnan(BLOm))
     &           then
               int_Born=0d0
               goto 999
            endif
         endif
         den=BLO*lum+BLOm*lumm
         int_Born=den*xjac
c        factorisation-scale variations: PDF ratios at k*muF
         call sv_get_muf_factors(nvf,facf)
         do ii=1,nvf
            num=BLO*pdf_lumi(xbk(1),xbk(2),facf(ii)*muf,1)
            if(has_mirror)num=num+
     &           BLOm*pdf_lumi(xbk(1),xbk(2),facf(ii)*muf,2)
            rf(ii)=1d0
            if(den.ne.0d0)rf(ii)=num/den
         enddo
         call sv_set_muf_ratios(rf)
      endif
c     renormalisation-scale variations: B ~ alpha_s^b
      call sv_add(int_Born,0d0,born_qcd_power)
c
c     plot Born
      wgtpl=int_Born*wgt
      call sv_wgts(wgtpl,0d0,born_qcd_power,wgts)
      if(doplot) then
         call invariants_from_p(p,nexternal,sLO,ierr)
         if(ierr.ne.0)then
            nfail=nfail+1
            write(77,*)'ierr failure: B',ierr,ich,failure_stage
            goto 999
         endif
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

