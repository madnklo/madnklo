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
      INCLUDE 'qcd_power.inc'
      INCLUDE 'initial_states.inc'
      integer ndim,ierr,ievt,nthres,i
      save ievt,nthres
      integer nitV
      common/niterationsv/nitV
      double precision int_virtual,VNLO(-2:0),INLO(-2:0)
      double precision sLO(nexternal,nexternal)
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgts(sv_maxwgt),wgtpl
      double precision BLOI,dvsv
      logical doplot, docut
      common/cdoplot/doplot
      double precision p(0:3,nexternal)
      double precision xjac
      integer fl_factor 
      common/flavour_factor/fl_factor
      double precision ans(0:1) !TODO SET CORRECTLY RANGE OF ANS 
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
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
c     number of points whose poles were tested (also read by the driver
c     in test-only mode)
      common/c_test_poles_v/ntested
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
      integer ich
      common/comich/ich
      integer nfail
      character*8 failure_stage
      common/cfailure_count/nfail
      common/cfailure_stage/failure_stage
      integer iconfig,mincfig,maxcfig,invar
      common/cfig/iconfig,mincfig,maxcfig,invar
      integer NGRAPHS2
      double precision amp2(N_MAX_CG)
      COMMON/TO_AMP2/AMP2,NGRAPHS2
c     multi-channel weight of channel ich (born_channels.f)
      double precision chw,born_chan_weight
      external born_chan_weight
      INTEGER NCOLORCORRELATORS
      PARAMETER (NCOLORCORRELATORS=4)     
C     Index 0 is the number of correlators to consider and the next
C      indices are which one to consider
      INTEGER COLOR_CORRELATORS_TO_CONSIDER(0:NCOLORCORRELATORS)
      REAL*8 COLOR_CORRELATED_EVALS(NCOLORCORRELATORS, 0:3,0:NSQUAREDSO)
      COMMON/V_ML5_1_1_COLOR_CORRELATIONS/COLOR_CORRELATORS_TO_CONSIDER,COLOR_CORRELATED_EVALS
      integer last_ich
      logical loop_arrays_initialised
      data last_ich /-1/
      data loop_arrays_initialised /.false./
      save last_ich,loop_arrays_initialised
c     hadron beams (lpp != 0): PDFs at the Born momentum fractions and
c     mirror processes, as in LO_B.f (pdf_lumi.f)
      logical hadrons,firsttime
      data firsttime/.true./
      save hadrons,firsttime
      double precision pm(0:3,nexternal),sLOm(nexternal,nexternal)
      double precision VNLOm(-2:0),INLOm(-2:0),BLOIm
      double precision muf,lum,lumm
c     factorisation-scale variations (rw_fscale): nvf factors facf (0:
c     none), luminosities, J + C (index 0: central muF), V+I+J+C and
c     the coefficient of the muR log at each muF (at most sv_max = 10
c     factors, ren_scale.f)
      integer nfmax
      parameter (nfmax=10)
      integer nvf,ii
      double precision facf(nfmax),lumf(nfmax),lummf(nfmax)
      double precision conv(0:nfmax),convm(0:nfmax)
      double precision vf(nfmax),dvf(nfmax),VIp,VIm
c     number of vegas variables: the last one is x of J(x) + C(x)
      integer ndimv
      common/cdim/ndimv
c     integrand (driver_v.f): vmode = 0 loop + I + J + C, 1 loop only,
c     2 I + J + C only; floop, fct = 1 or 0 select the two parts
      integer vmode
      common/c_vmode/vmode
      double precision floop,fct
      double precision get_muf,pdf_lumi,pdf_swap_initial
      external get_muf,pdf_lumi,pdf_swap_initial
c
      if(firsttime)then
         hadrons=(lpp(1).ne.0.or.lpp(2).ne.0)
         if(hadrons)call set_tau_min(nexternal,leg_pdgs)
         firsttime=.false.
      endif
c
c     Allocate the integrand-owned loop work arrays once.
      if (.not.loop_arrays_initialised) then
         CALL %(long_proc_prefix)sGET_ANSWER_DIMENSION(MATELEM_ARRAY_DIM)
         ALLOCATE(MATELEM(0:3,0:MATELEM_ARRAY_DIM))
         CALL %(long_proc_prefix)sGET_NSQSO_LOOP(NSQUAREDSO_LOOP)
         ALLOCATE(PREC_FOUND(0:NSQUAREDSO_LOOP))
         loop_arrays_initialised=.true.
      endif
c
c     Refresh channel-dependent phase-space data only on an ich change.
      if (ich.ne.last_ich) then
         call initialise_born_channel()
         last_ich=ich
      endif
C
C     EXTERNAL
c     initialise
      floop=1d0
      if(vmode.eq.2)floop=0d0
      fct=1d0
      if(vmode.eq.1)fct=0d0
      xjac = Gevtopb
      int_virtual = 0d0
      VNLO = 0d0
      INLO = 0d0
c
      call gen_mom(iconfig,mincfig,maxcfig,invar,xjac,x,p,nexternal)
      if(xjac.eq.0d0)goto 999
      call invariants_from_p(p,nexternal,sLO,ierr)
      if(ierr.ne.0)then
         nfail=nfail+1
         write(77,*)'ierr failure: V',ierr,ich,failure_stage
         goto 999
      endif
c
c     possible cuts
      if(docut(p,nexternal,leg_pdgs,0))goto 999
c
c     renormalisation scale for the n-body kinematics
c     (set_mur also updates MU_R and the couplings used by MadLoop)
      call set_mur_from_momenta(p,nexternal,leg_pdgs)
      alphas=alphas_current
c
c     call virtual (not for I + J + C only)
      if(vmode.ne.2)then
         COLOR_CORRELATED_EVALS = 0d0
         CALL %(long_proc_prefix)sSLOOPMATRIX_THRES(p,MATELEM,-1.0D0,PREC_FOUND,RETURNCODE)
         VNLO(-2:0) = [(MATELEM(i,0), i=3,1,-1)]
      endif
c
c     call counterterm
      call int_counter_NLO(p,sLO,INLO,BLOI,ierr)
      if(ierr.ne.0)then
         nfail=nfail+1
         write(77,*)'ierr failure: V',ierr,ich,failure_stage
         goto 999
      endif
c
c     test coefficients of epsilon poles (needs the loop)
      if(ntested.lt.ntest.and.vmode.ne.2)then
         ntested=ntested+1
         write(50,*)
         write(50,*)'Testing point # ', ntested
         write(50,*)'Double pole V, I, sum', VNLO(-2), INLO(-2), VNLO(-2)+INLO(-2)
         write(50,*)'Single pole V, I, sum', VNLO(-1), INLO(-1), VNLO(-1)+INLO(-1)
         write(50,*)
      endif
c
c     subtracted virtual
      if(.not.hadrons)then
         int_virtual=(floop*VNLO(0)+fct*INLO(0))*xjac
c
c     apply flavour multiplicity factor
         int_virtual=int_virtual*fl_factor
      else
c     hadron beams: sum over the processes of this directory of
c     f_a(x1) f_b(x2) (V+I)_ab(p), plus the mirror processes, with V+I
c     evaluated with the initial momenta exchanged
         muf=get_muf(p,nexternal,leg_pdgs)
         lum=pdf_lumi(xbk(1),xbk(2),muf,1)
         lumm=0d0
         BLOIm=0d0
         VNLOm=0d0
         INLOm=0d0
         if(has_mirror)then
            lumm=pdf_lumi(xbk(1),xbk(2),muf,2)*pdf_swap_initial(p,pm,nexternal)
            call invariants_from_p(pm,nexternal,sLOm,ierr)
            if(ierr.ne.0)goto 999
            if(vmode.ne.2)then
               COLOR_CORRELATED_EVALS = 0d0
               CALL %(long_proc_prefix)sSLOOPMATRIX_THRES(pm,MATELEM,-1.0D0,PREC_FOUND,RETURNCODE)
               VNLOm(-2:0) = [(MATELEM(i,0), i=3,1,-1)]
            endif
            call int_counter_NLO(pm,sLOm,INLOm,BLOIm,ierr)
            if(ierr.ne.0)then
               nfail=nfail+1
               write(77,*)'ierr failure: V (mirror)',ierr,ich,failure_stage
               goto 999
            endif
         endif
         VIp=floop*VNLO(0)+fct*INLO(0)
         VIm=floop*VNLOm(0)+fct*INLOm(0)
         int_virtual=(VIp*lum+VIm*lumm)*xjac
c     factorisation-scale variations: luminosities at facf(ii)*muF
         call sv_get_muf_factors(nvf,facf)
         do ii=1,nvf
            lumf(ii)=pdf_lumi(xbk(1),xbk(2),facf(ii)*muf,1)
            lummf(ii)=0d0
            if(has_mirror)lummf(ii)=pdf_lumi(xbk(1),xbk(2),facf(ii)*muf,2)
            vf(ii)=(VIp*lumf(ii)+VIm*lummf(ii))*xjac
         enddo
c     x-dependent contribution J(x) + C(x) (2209.09123, eq. sigNLOsub),
c     with the PDFs inside (int_convolution_NLO, NLO_I_template.f);
c     ~ alpha_s^(b+1), without explicit log(muR) (they cancel between
c     its soft terms by colour conservation). It depends on muF also
c     through log(muF) (PDF counterterm): computed at each facf(ii)*muF
         if(vmode.ne.1)then
            call int_convolution_NLO(p,sLO,1,x(ndimv),nvf,facf,conv,ierr)
            if(ierr.ne.0)goto 999
            convm=0d0
            if(has_mirror)then
               call int_convolution_NLO(pm,sLOm,2,x(ndimv),nvf,facf,convm,ierr)
               if(ierr.ne.0)goto 999
            endif
            int_virtual=int_virtual+(conv(0)+convm(0))*xjac
            do ii=1,nvf
               vf(ii)=vf(ii)+(conv(ii)+convm(ii))*xjac
            enddo
         endif
c        renormalisation-scale variations below: the Born summed in the
c        same way (at each muF)
         do ii=1,nvf
            dvf(ii)=BLOI*lumf(ii)+BLOIm*lummf(ii)
         enddo
         BLOI=BLOI*lum+BLOIm*lumm
      endif
c
c     Multi channeling
      chw = born_chan_weight(ich)
      int_virtual = int_virtual * chw
c
c     renormalisation-scale variations, with b = born_qcd_power:
c     (V+I)(k muR) = r_k^(b+1) [ (V+I)(muR) + b beta0/2 alphas/(2 pi) B log(k^2) ]
c     (the explicit infrared logarithms of V and I cancel, since their
c     poles do; the remaining log comes from the UV renormalisation, so
c     it goes with the loop: floop = 0 for I + J + C only)
c     With hadron beams, the same at each factorisation scale facf(ii)*muF
c     (vf, dvf).
      if(.not.hadrons)then
         dvsv=floop*born_qcd_power*beta0/2d0*alphas/(2d0*pi)*BLOI*xjac*fl_factor*chw
         call sv_add(int_virtual,dvsv,born_qcd_power+1)
      else
         dvsv=floop*born_qcd_power*beta0/2d0*alphas/(2d0*pi)*BLOI*xjac*chw
         do ii=1,nvf
            vf(ii)=vf(ii)*chw
            dvf(ii)=floop*born_qcd_power*beta0/2d0*alphas/(2d0*pi)*dvf(ii)*xjac*chw
         enddo
         call sv_add_muf(int_virtual,dvsv,vf,dvf,born_qcd_power+1)
      endif
c
c     plot
      wgtpl=int_virtual*wgt
      if(.not.hadrons)then
         call sv_wgts(wgtpl,dvsv*wgt,born_qcd_power+1,wgts)
      else
         do ii=1,nvf
            vf(ii)=vf(ii)*wgt
            dvf(ii)=dvf(ii)*wgt
         enddo
         call sv_wgts_muf(wgtpl,dvsv*wgt,vf,dvf,born_qcd_power+1,wgts)
      endif
      if(doplot)call analysis_fill(p,sLO,nexternal,leg_pdgs,wgts)
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
c      iconfig = 1
      mincfig = 1
      maxcfig = 1
      invar = 2
      call configs_born
      call props_born
      call decaybw_born
      call getleshouche_born
c
      return
      end
