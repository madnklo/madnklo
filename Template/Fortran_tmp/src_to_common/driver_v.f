      program driver_n_body
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'colored_partons.inc'
      INCLUDE 'ngraphs.inc'
      integer nchan,n_born_channels
      external n_born_channels
      integer mxdim
      parameter(mxdim=30)
c     test-only mode (command-line argument --test_only)
      logical test_only
      character*40 arg
      integer iarg
      double precision xtest(mxdim),wtest,ftest,ran2
      external ran2
      integer ntested_v,ntest_v
      parameter(ntest_v=20)
      common/c_test_poles_v/ntested_v
      integer ndim,i,j,idum
      double precision s_had
      integer iu,iu1,iu2,iu7,iu8
      common/cdim/ndim
      double precision int_virtual
      double precision res_V,err_V
      external int_virtual
      integer order
      logical doplot
      common/cdoplot/doplot
      double precision rescale_plot_V
      character*100 line
      integer nitVth,nclVth,nitV,nclV
      integer nitVth0,nclVth0,nclV0,nclVth1,nclV1
      integer iseed
      common /to_seed/iseed
      common/niterationsv/nitV
c
c     vegas declarations
      integer ndmx,nprn,ndo,init,it
      parameter (ndmx=50,nprn=0)
      double precision chi2a,acc,si,swgt,schi
      double precision region(2*mxdim),xi(ndmx,mxdim)
      parameter(acc=1d-10)
      common/rand/idum
c
      integer ich
      common/comich/ich
      integer nfail
      character*8 failure_stage
      common/cfailure_count/nfail
      common/cfailure_stage/failure_stage
      double precision sum_v, sum_err_v
      double precision sum_err_v_a,err_v_a(N_MAX_CG)
      double precision xi_ch(ndmx,mxdim,N_MAX_CG)
      integer ndo_ch(N_MAX_CG),ncl_ch(N_MAX_CG),nclw_ch(N_MAX_CG)
      integer k,itw,nflw
      double precision s_ch(N_MAX_CG),s_prev(N_MAX_CG),s_main,sum_s
      double precision rcalls
      common/abchia/rcalls
c     phases (--warmup, --main: see below)
      integer iuw,ncl_tot,nch_w,ndim_w,nit_w
      double precision s_mn(N_MAX_CG),ncalls_w,t_w0,t_w1,cpu_call
      integer rstate(35)
      logical do_warm,do_main
      character*7 fpos
      integer nwgt
      character*80 weights_info(100)
c     integrand (common with int_virtual, NLO_V_*.f):
c       vmode = 0: loop + I + J + C;
c       vmode = 1: loop only (I + J + C integrated by './virtual
c                  --counterterms', n_vc_folds > 1 in run_card.dat);
c       vmode = 2: I + J + C only ('--counterterms'), with n_vc_folds
c                  times the calls of V; tag of the output files: V or VC
      integer vmode
      common/c_vmode/vmode
      logical counterterms
      character*2 tag
c     virtual trick (loop computed only at a fraction of the points
c     from the second warm-up iteration on, see below and
c     NLO_V_template.f): vtrick = .true. switches it on by default; the
c     environment variable MADNKLO_VTRICK = 0 (1) switches it off (on)
c     at run time. At least vt_nlmin loops per channel and main
c     iteration (reliable errors of the iterations, combined with
c     weights 1/sigma^2), vt_nlminw per warm-up iteration
      logical vtrick
      parameter (vtrick=.true.)
      integer vt_nlmin,vt_nlminw
      parameter (vt_nlmin=1000,vt_nlminw=200)
      double precision vt_c(4,N_MAX_CG),vt_p(N_MAX_CG),vt_n(N_MAX_CG)
      double precision vt_g(4,4,N_MAX_CG),vt_b(4,N_MAX_CG)
      double precision vt_yy(N_MAX_CG)
      double precision vt_tl(N_MAX_CG),vt_tt(N_MAX_CG),vt_nl(N_MAX_CG)
      integer vt_mode,vt_s1,vt_s2
      common/c_vtrick/vt_c,vt_p,vt_n,vt_g,vt_b,vt_yy,vt_tl,vt_tt,vt_nl,
     &     vt_mode,vt_s1,vt_s2
      logical vt_on,vt_main
      integer iversion,vt_ncl
      character*40 envv
      double precision vt_popt(N_MAX_CG),vt_d(N_MAX_CG)
      double precision vt_sraw(N_MAX_CG)
      double precision vt_tlc(N_MAX_CG),vt_trc(N_MAX_CG),vt_gain
      double precision vt_cf,vt_sw,vt_pe,vt_pmain,vt_popt_f
      external vt_pmain,vt_popt_f
      double precision vt_pw(N_MAX_CG),vt_sr1(N_MAX_CG),vt_sr0(N_MAX_CG)
      double precision vt_d1(N_MAX_CG),vt_d0(N_MAX_CG)
      double precision vt_di(N_MAX_CG),vt_ni(N_MAX_CG),vt_nlw(N_MAX_CG)
      common/c_vtrick_w/vt_di,vt_ni,vt_nlw
c
      nfail=0
      sum_v=0d0
      sum_err_v=0d0
      res_v=0d0
      err_v=0d0
c
      call SETPARA('param_card.dat')
      call SETRUN('run_card.dat')
c
c     read inputs
      region=0d0
      order=0
      idum = -max(abs(iseed),1)
      s_had = (EBEAM(1)+EBEAM(2))**2
c
c     command-line arguments: '--test_only' (tests only, see below),
c     '--counterterms' (I + J + C only), and the phases '--warmup' and
c     '--main' (see below)
      test_only=.false.
      counterterms=.false.
      do_warm=.true.
      do_main=.true.
      do iarg=1,command_argument_count()
         call get_command_argument(iarg,arg)
         if(trim(arg).eq.'--test_only')test_only=.true.
         if(trim(arg).eq.'--counterterms')counterterms=.true.
         if(trim(arg).eq.'--warmup')do_main=.false.
         if(trim(arg).eq.'--main')do_warm=.false.
      enddo
      iuw=97
      fpos='rewind'
      if(.not.do_warm)fpos='append'
      vmode=0
      if(n_vc_folds.gt.1)vmode=1
      if(counterterms)then
         if(n_vc_folds.le.1)then
            write(*,*)'virtual --counterterms: n_VC_folds <= 1 in the',
     &           ' run_card (I + J + C integrated together with V)'
            stop 1
         endif
         vmode=2
      endif
      if(test_only)vmode=0
c     virtual trick: only with the loop (not for I + J + C only, nor in
c     the tests)
      vt_on=vtrick
      call get_environment_variable('MADNKLO_VTRICK',envv)
      if(trim(envv).eq.'0')vt_on=.false.
      if(trim(envv).eq.'1')vt_on=.true.
      if(vmode.eq.2.or.test_only)vt_on=.false.
      vt_mode=0
      vt_main=.false.
      vt_n=0d0
      vt_g=0d0
      vt_b=0d0
      vt_yy=0d0
      vt_tl=0d0
      vt_tt=0d0
      vt_nl=0d0
      vt_c=0d0
      vt_p=1d0
      vt_popt=1d0
      vt_d=0d0
      vt_sraw=0d0
      vt_tlc=0d0
      vt_trc=0d0
      vt_pw=1d0
      vt_di=0d0
      vt_ni=0d0
      vt_nlw=0d0
      vt_sr1=0d0
      vt_sr0=0d0
      vt_d1=0d0
      vt_d0=0d0
      if(vt_on)call vt_seed(iseed,1)
      tag='V'
      if(vmode.eq.2)tag='VC'
      NITVTH = NITERS_V_GRID
      NCLVTH = NPOINTS_V_GRID
      NITV = NITERS_V
      NCLV = NPOINTS_V
      if(vmode.eq.2)then
         NCLVTH = N_VC_FOLDS*NPOINTS_V_GRID
         NCLV = N_VC_FOLDS*NPOINTS_V
      endif
c     TODO: understand muR input fixed/dyn scale
c
c     initialise physics parameters
      iu1=44
      iu2=50
      iu=55
      iu7=77
      iu8=88
c
c     phase-space dimension, same for all contributions to this folder:
c     one more variable per hadron beam (momentum fractions, as in
c     gen_mom and driver_n.f)
      ndim=3*(nexternal-2)-4
      do i=1,2
         if(lpp(i).ne.0) ndim = ndim + 1
      enddo
c     hadron beams: one more variable, the momentum fraction x of the
c     x-dependent contribution J(x) + C(x) (int_convolution_NLO); it is
c     the last one (gen_mom uses its own number of Born variables)
      if(lpp(1).ne.0.or.lpp(2).ne.0) ndim = ndim + 1
      do i=1,ndim
         region(i)=0d0
         region(i+ndim)=1d0
      enddo
c
c     test-only mode: './virtual --test_only' only runs the tests
      if(test_only)then
c        pole cancellation V+I on the first ntest_v points passing the
c        cuts (test_poles_V.log), without integration
         open(unit=iu2,file='test_poles_V.log')
         open(unit=iu7,file='failures_V_test.log',status='replace')
         ich=1
         doplot=.false.
         failure_stage='test'
         ntested_v=0
         do i=1,100000
            do j=1,ndim
               xtest(j)=ran2(idum)
            enddo
            wtest=1d0
            ftest=int_virtual(xtest,wtest)
            if(ntested_v.ge.ntest_v)exit
         enddo
         close(iu2)
         close(iu7)
         write(*,*)'Test only:',ntested_v,
     &        ' points with V+I poles written to test_poles_V.log'
         stop
      endif
c
c     initialise histograms and open output files
c     scale variations (rw_rscale, rw_fscale in run_card.dat): one
c     histogram weight per (muF,muR) pair; with hadron beams, muF is
c     varied by PDF reweighting:
c     the integrand computes V+I and J+C at each muF (NLO_V_template.f)
      call sv_allow_muf()
      call sv_init()
      call sv_weights_info(nwgt,weights_info,100)
      call analysis_begin(nwgt,weights_info)
c     phases: './virtual --warmup' (also with --counterterms) only
c     runs the warm-up and writes warmup_<tag>.dat (grids, standard
c     deviation per point of each channel, cost per call, state of the
c     random numbers); './virtual --main' reads it and runs the main
c     integration with the calls per iteration of main_calls_<tag>.dat,
c     written by run_all.py, which shares the calls among the P1_*
c     directories of the contribution (optimal for the total error at
c     fixed CPU time); without these arguments, both, with npoints_V
c     calls per iteration (times n_VC_folds for VC). --warmup followed
c     by --main with these calls gives the same results (up to the
c     initialisation of MadLoop, done in the first process).
      open(unit=iu1,file='integration_'//trim(tag)//'.log',
     &     position=fpos)
c     pole test (V + I) on the first points: needs the loop
      if(vmode.ne.2)open(unit=iu2,file='test_poles_V.log',position=fpos)
      if(do_warm)then
         open(unit=iu7,file='failures_'//trim(tag)//'.log',
     &        status='replace',action='write')
      else
         open(unit=iu7,file='failures_'//trim(tag)//'.log',
     &        position='append',action='write')
      endif
      open(unit=iu8,file=trim(tag)//'_chan.log',position=fpos)
      if(do_main)then
         open(unit=iu,file='results_'//trim(tag)//'.log')
         if(vmode.eq.0)write(iu,*)' Virtual contribution '
         if(vmode.eq.1)write(iu,*)' Virtual contribution, loop only '
         if(vmode.eq.2)write(iu,*)
     &        ' Virtual contribution, I + J + C only '
      endif
      line='=================================================='
c
c     channels: the phase-space configurations (configs.f), each
c     integrating the virtual times its multi-channel weight
c     (born_chan_weight). Same scheme as the Born (driver_n.f).
c     Warm-up, iteration by iteration over all the channels: each
c     channel runs one vegas iteration on its own grid (flat at the
c     first one), which refines the grid; its standard deviation per
c     point on that grid, s = (error of the iteration)*sqrt(calls), is
c     measured. The calls of an iteration (npoints_V_grid) are shared
c     equally at the first one, then in proportion to s of the
c     previous one (optimal for the total error), at least
c     max(1000,npoints_V_grid/(10 nchan)) each, so that every channel
c     stays measured (at least 1000 calls per channel and iteration:
c     with fewer, the combination of the iterations by vegas, weighted
c     by 1/sigma^2, is biased for integrands with large fluctuations,
c     e.g. by 20% for the virtual of e+e- > 3 jets).
c     Main runs: the calls (npoints_V) are shared in proportion to s
c     averaged over the last two warm-up iterations, at least
c     max(1000,npoints_V/(10 nchan)) each; each main run continues on
c     the grid of its warm-up.
      nchan=n_born_channels()
      if(do_warm)then
      call cpu_time(t_w0)
      ncalls_w=0d0
      nflw=max(1000,nclVth/(10*nchan))
      do ich=1,nchan
         nclw_ch(ich)=max(1000,nclVth/nchan)
         s_ch(ich)=0d0
         s_prev(ich)=0d0
      enddo
      do itw=1,max(nitVth,1)
         do ich=1,nchan
            if(itw.eq.1)then
               init=0
            else
               do k=1,ndim
                  do j=1,ndmx
                     xi(j,k)=xi_ch(j,k,ich)
                  enddo
               enddo
               ndo=ndo_ch(ich)
               init=1
            endif
            if(itw.eq.1)then
               write(*,*)'Virtual warmup for channel',ich
               write(iu1,*)
               write(iu1,*)'============================='
               write(iu1,*)' VIRTUAL WARMUP, CHANNEL',ich
               write(iu1,*)'============================='
            endif
            doplot=.false.
            failure_stage='warmup'
c           virtual trick: in the first iteration the loop is always
c           computed; then with probability vt_pw(ich) (and the
c           reference fitted on the previous iterations); sums for the
c           fit in all the iterations
            vt_mode=0
            if(vt_on)then
               vt_mode=1
               vt_p(ich)=vt_pw(ich)
               call vt_iter_start(ich)
            endif
            call vegas(region,ndim,int_virtual,init,nclw_ch(ich),1,
     &      nprn,res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
            vt_mode=0
            do k=1,ndim
               do j=1,ndmx
                  xi_ch(j,k,ich)=xi(j,k)
               enddo
            enddo
            ndo_ch(ich)=ndo
            ncalls_w=ncalls_w+rcalls
            s_prev(ich)=s_ch(ich)
            s_ch(ich)=err_v*dsqrt(rcalls)
            if(.not.(s_ch(ich).ge.0d0.and.s_ch(ich).lt.huge(1d0)))
     &           s_ch(ich)=0d0
c           virtual trick: s without the trick (vt_sr1, squared), the
c           residual D (vt_d1), new fit of the reference, probability
c           of the loop and CPU times per loop and per call for the rest
            if(vt_on)then
               vt_sr0(ich)=vt_sr1(ich)
               vt_d0(ich)=vt_d1(ich)
               call vt_iter_end(ich,itw,s_ch(ich),vt_sr1(ich),
     &              vt_d1(ich),vt_popt(ich),vt_tlc(ich),vt_trc(ich))
               write(iu8,'(a,i4,i4,a,f8.5,a,es10.3,a,es10.3,a,f8.5)')
     &              ' '//trim(tag)//' trick warmup: channel, itn',ich,
     &              itw,', p =',vt_p(ich),', sd/point without trick =',
     &              sqrt(vt_sr1(ich)),', D/S =',vt_d1(ich)
     &              /max(vt_sr1(ich),1d-300),', p_opt =',vt_popt(ich)
            endif
            write(iu8,*)trim(tag)//
     &           ' warmup: channel, itn, calls, sd/point = ',
     &           ich,itw,nclw_ch(ich),s_ch(ich)
         enddo
c        calls of the next warm-up iteration
         sum_s=0d0
         do ich=1,nchan
            sum_s=sum_s+s_ch(ich)
         enddo
         if(nchan.gt.1.and.sum_s.gt.0d0)then
            do ich=1,nchan
               nclw_ch(ich)=max(nflw,nint(nclVth*s_ch(ich)/sum_s))
            enddo
         endif
c        virtual trick: probability of the loop in the next iteration
         if(vt_on)then
            do ich=1,nchan
               vt_pw(ich)=vt_pmain(vt_popt(ich),nclw_ch(ich),vt_nlminw)
            enddo
         endif
      enddo
      call cpu_time(t_w1)
      cpu_call=0d0
      if(ncalls_w.gt.0d0)cpu_call=(t_w1-t_w0)/ncalls_w
      do ich=1,nchan
         s_mn(ich)=s_ch(ich)
         if(nitVth.ge.2)s_mn(ich)=0.5d0*(s_ch(ich)+s_prev(ich))
      enddo
c     virtual trick: s without the trick (vt_sraw) and residual D
c     (vt_d) of each channel, averaged over the last two warm-up
c     iterations as s_mn, and probability of the loop for the main runs;
c     s_mn and the CPU time per call become those expected with the
c     trick, so that the calls are shared accordingly (here and by
c     run_all.py), with the probabilities expected for the calls per
c     channel of a main run with the default calls (see vt_pmain)
      vt_cf=0d0
      vt_sw=0d0
      sum_s=0d0
      do ich=1,nchan
         vt_sraw(ich)=s_mn(ich)
         if(vt_on)then
            vt_sraw(ich)=sqrt(vt_sr1(ich))
            vt_d(ich)=vt_d1(ich)
            if(nitVth.ge.2)then
               vt_sraw(ich)=0.5d0*(sqrt(vt_sr1(ich))+sqrt(vt_sr0(ich)))
               vt_d(ich)=0.5d0*(vt_d1(ich)+vt_d0(ich))
            endif
            vt_popt(ich)=vt_popt_f(vt_sraw(ich)**2,vt_d(ich),
     &           vt_tlc(ich),vt_trc(ich),vt_n(ich))
         endif
         sum_s=sum_s+vt_sraw(ich)
      enddo
      do ich=1,nchan
         if(vt_on)then
            vt_ncl=nclV
            if(nchan.gt.1.and.sum_s.gt.0d0)vt_ncl=max(1000,
     &           nclV/(10*nchan),nint(nclV*vt_sraw(ich)/sum_s))
            vt_pe=vt_pmain(vt_popt(ich),vt_ncl,vt_nlmin)
            s_mn(ich)=sqrt(vt_sraw(ich)**2
     &           +(1d0/vt_pe-1d0)*vt_d(ich))
            if(vt_tlc(ich)+vt_trc(ich).gt.0d0)then
               vt_cf=vt_cf+s_mn(ich)*(vt_trc(ich)+vt_pe*vt_tlc(ich))
               vt_sw=vt_sw+s_mn(ich)
            endif
         endif
      enddo
c     (the CPU time per call measured in the warm-up mixes iterations
c     with and without the trick: the one expected in the main runs)
      if(vt_sw.gt.0d0.and.vt_cf.gt.0d0)cpu_call=vt_cf/vt_sw
      if(.not.do_main)then
c        warm-up only: save the state and stop
         call ran2_get_state(rstate)
         open(unit=iuw,file='warmup_'//trim(tag)//'.dat',
     &        status='replace')
         write(iuw,*)'MADNKLO_WARMUP 2'
         write(iuw,*)nchan,ndim,ndmx,nclV,nitV,ncalls_w,cpu_call
         do ich=1,nchan
            write(iuw,'(es26.17e3,i6)')s_mn(ich),ndo_ch(ich)
            write(iuw,'(4es26.17e3)')
     &           ((xi_ch(j,k,ich),j=1,ndmx),k=1,ndim)
         enddo
         write(iuw,'(5i14)')rstate
c        virtual trick (version 2): per channel the probability of the
c        loop, the coefficients of the reference, the residual and the
c        standard deviation per point without the trick, the CPU times
c        per call of the loops and of the rest
         do ich=1,nchan
            write(iuw,'(10es26.17e3)')vt_popt(ich),(vt_c(k,ich),k=1,4),
     &           vt_d(ich),vt_sraw(ich),vt_tlc(ich),vt_trc(ich)
         enddo
         close(iuw)
         close(iu1)
         if(vmode.ne.2)close(iu2)
         close(iu7)
         close(iu8)
         write(*,*)'Warm-up written to warmup_'//trim(tag)//'.dat'
         stop
      endif
      else
c     main only: state of the warm-up
      open(unit=iuw,file='warmup_'//trim(tag)//'.dat',status='old')
      read(iuw,'(a)')line
      iversion=1
      if(index(line,'MADNKLO_WARMUP 2').gt.0)iversion=2
      read(iuw,*)nch_w,ndim_w,j,k,nit_w,ncalls_w,cpu_call
      if(nch_w.ne.nchan.or.ndim_w.ne.ndim.or.j.ne.ndmx)then
         write(*,*)'warmup_'//trim(tag)//'.dat does not match'
         stop 1
      endif
      do ich=1,nchan
         read(iuw,*)s_mn(ich),ndo_ch(ich)
         read(iuw,*)((xi_ch(j,k,ich),j=1,ndmx),k=1,ndim)
      enddo
      read(iuw,*)rstate
      if(iversion.ge.2)then
         do ich=1,nchan
            read(iuw,*)vt_popt(ich),(vt_c(k,ich),k=1,4),vt_d(ich),
     &           vt_sraw(ich),vt_tlc(ich),vt_trc(ich)
         enddo
      else
         vt_on=.false.
      endif
      close(iuw)
      call ran2_set_state(rstate)
c     virtual trick switched off for the main runs only: the calls are
c     shared with the standard deviations without the trick
      if(.not.vt_on.and.iversion.ge.2)then
         do ich=1,nchan
            if(vt_sraw(ich).gt.0d0)s_mn(ich)=vt_sraw(ich)
         enddo
      endif
      endif
c     calls per iteration of the main runs: npoints_V (times n_VC_folds
c     for VC), or those given by run_all.py (main_calls_<tag>.dat)
      ncl_tot=nclV
      if(.not.do_warm)then
         open(unit=iuw,file='main_calls_'//trim(tag)//'.dat',
     &        status='old',err=11)
         read(iuw,*)ncl_tot
         close(iuw)
 11      continue
      endif
      sum_s=0d0
      do ich=1,nchan
         sum_s=sum_s+s_mn(ich)
      enddo
      do ich=1,nchan
         if(nchan.eq.1)then
            ncl_ch(ich)=ncl_tot
         elseif(sum_s.gt.0d0)then
            ncl_ch(ich)=max(1000,ncl_tot/(10*nchan),
     &           nint(ncl_tot*s_mn(ich)/sum_s))
         else
            ncl_ch(ich)=max(1000,ncl_tot/nchan)
         endif
      enddo
c     virtual trick: probability of the loop of each channel, at least
c     vt_nlmin loops per iteration (1: no trick in this channel); the
c     choices are drawn from their own random numbers (vt_rand), so that
c     those of vegas are not changed
      do ich=1,nchan
         vt_p(ich)=1d0
         if(vt_on)vt_p(ich)=vt_pmain(vt_popt(ich),ncl_ch(ich),vt_nlmin)
         if(vt_p(ich).lt.1d0)vt_main=.true.
      enddo
      if(vt_main)call vt_seed(iseed,2)
c
c     main runs, on the grids of the warm-ups
      do ich=1,nchan
         do k=1,ndim
            do j=1,ndmx
               xi(j,k)=xi_ch(j,k,ich)
            enddo
         enddo
         ndo=ndo_ch(ich)
         write(*,*)'Virtual for channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' VIRTUAL, CHANNEL',ich
         write(iu1,*)'============================='
         init=1
         doplot=.true.
         failure_stage='main'
         vt_mode=0
         if(vt_p(ich).lt.1d0)vt_mode=2
         call vegas(region,ndim,int_virtual,init,ncl_ch(ich),nitV,nprn,
     &   res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
         vt_mode=0
         call sv_add_channel()
         write(iu8,*)trim(tag)//': channel, itns, calls = ',ich,nitV,
     &        ncl_ch(ich)
         if(vt_on)then
c           virtual trick: summary of the channel; expected gain in
c           error^2 x CPU time
            vt_gain=1d0
            if(vt_p(ich).ge.1d0)vt_nl(ich)=dble(nitV)*ncl_ch(ich)
            if(vt_p(ich).lt.1d0.and.vt_sraw(ich).gt.0d0)
     &           vt_gain=vt_sraw(ich)**2*(vt_trc(ich)+vt_tlc(ich))
     &           /((vt_sraw(ich)**2+(1d0/vt_p(ich)-1d0)*vt_d(ich))
     &           *(vt_trc(ich)+vt_p(ich)*vt_tlc(ich)))
            write(iu8,'(a,i4,a,f8.5,a,4es12.4,a,es10.3,a,f8.4,a,f7.2)')
     &           ' '//trim(tag)//' trick: channel',ich,', p =',
     &           vt_p(ich),', c =',(vt_c(k,ich),k=1,4),
     &           ', D/S =',vt_d(ich)/max(vt_sraw(ich)**2,1d-300),
     &           ', loops/calls =',vt_nl(ich)/dble(max(1,nitV
     &           *ncl_ch(ich))),', expected gain =',vt_gain
         endif
         rescale_plot_V=dble(nitV)/min(dble(nitV),dble(it))
         sum_v = sum_v + res_v
         sum_err_v = sum_err_v + err_v**2
         write(iu8,*)' sigma '//trim(tag)//' [pb], channel',ich,' = ',
     &   res_v,' +-',err_v
         write(*,*)'...done'
      enddo
c
c     finalise histograms and output files
      call analysis_end(1d0,'plot_'//trim(tag)//'.dat')
      sum_err_v = dsqrt(sum_err_v)
      write(iu,*)' sigma '//trim(tag)//' [pb]  = ',sum_v,' +-',sum_err_v
      call sv_write_results(iu,trim(tag),sum_v,sum_err_v)
      close(iu)
      close(iu1)
      if(vmode.ne.2)close(iu2)
      if(nfail.gt.0)write(iu7,*)'Total ierr failures: ',nfail
      close(iu7)
      close(iu8)
c
      end


c     ------------------------------------------------------------------
c     Virtual trick (see also NLO_V_template.f). The loop L of a point
c     is computed only with probability p (vt_p, per channel); the
c     integrand uses instead the reference
c       R = sum_k c_k X_k,  X = (Born, finite part, single and double
c                                pole of the integrated counterterm I),
c     available at every point, plus (L - R)/p when the loop is
c     computed. The average is unchanged; per call the CPU time drops
c     from t_r + t_l (rest, loop) to t_r + p t_l, and the variance grows
c     from S to S + (1/p - 1) D, D = <(J (L - R))^2> (J: weight of the
c     point; the added variance is not reduced by the stratification).
c     The poles of I are minus those of the loop: they carry the
c     logarithms of the scales of V, so the reference follows the loop
c     at any choice of the scales. The scale variations are linear in
c     the loop, so they stay unbiased.
c     Warm-up: in the first iteration the loop is always computed; after
c     each iteration the c_k are fitted (least squares, minimum D) on
c     all the warm-up points so far (those without the loop weighted
c     as in the integrand: the loop terms of the sums by 1/p), and used
c     in the next iteration, where the loop is computed with probability
c     p = sqrt(D/(S - D) t_r/t_l) (minimum of variance x CPU time, at
c     least vt_pmin and vt_nlminw loops per iteration). S (without the
c     trick) is obtained from the measured s^2 = S + (1/p - 1) D, D from
c     the points of the iteration (with the c_k fixed before it). The c_k
c     and p are fixed within each iteration, so that every iteration is
c     unbiased. Main runs: c_k of the end of the warm-up, p from S and D
c     of its last two iterations (at least vt_nlmin loops per
c     iteration).
c     ------------------------------------------------------------------
      subroutine vt_accum(fac,wgt,v0,v0m,xp,xm,lum,lumm,loop)
c     warm-up: sums for the fit of the reference of channel ich; the
c     loop part of the integrand is fac*(v0*lum + v0m*lumm) (process
c     and mirror; if computed, loop), its weight J = wgt*rcalls (vegas),
c     the reference quantities xp, xm (process, mirror) enter in the
c     same way. Also the residual with the present reference (vt_di)
      implicit none
      include 'ngraphs.inc'
      double precision fac,wgt,v0,v0m,xp(4),xm(4),lum,lumm
      logical loop
      double precision vt_c(4,N_MAX_CG),vt_p(N_MAX_CG),vt_n(N_MAX_CG)
      double precision vt_g(4,4,N_MAX_CG),vt_b(4,N_MAX_CG)
      double precision vt_yy(N_MAX_CG)
      double precision vt_tl(N_MAX_CG),vt_tt(N_MAX_CG),vt_nl(N_MAX_CG)
      integer vt_mode,vt_s1,vt_s2
      common/c_vtrick/vt_c,vt_p,vt_n,vt_g,vt_b,vt_yy,vt_tl,vt_tt,vt_nl,
     &     vt_mode,vt_s1,vt_s2
c     sums of the current warm-up iteration: residual, calls; loops
c     computed in the warm-up
      double precision vt_di(N_MAX_CG),vt_ni(N_MAX_CG),vt_nlw(N_MAX_CG)
      common/c_vtrick_w/vt_di,vt_ni,vt_nlw
      integer ich
      common/comich/ich
      double precision rcalls
      common/abchia/rcalls
      double precision f,y,x(4),w,r
      integer k,l
      f=fac*wgt*rcalls
      do k=1,4
         x(k)=f*(xp(k)*lum+xm(k)*lumm)
      enddo
      do k=1,4
         do l=1,4
            vt_g(k,l,ich)=vt_g(k,l,ich)+x(k)*x(l)
         enddo
      enddo
      if(.not.loop)return
      w=1d0/vt_p(ich)
      y=f*(v0*lum+v0m*lumm)
      r=0d0
      do k=1,4
         vt_b(k,ich)=vt_b(k,ich)+w*x(k)*y
         r=r+vt_c(k,ich)*x(k)
      enddo
      vt_yy(ich)=vt_yy(ich)+w*y*y
      vt_di(ich)=vt_di(ich)+w*(y-r)**2
      return
      end


      subroutine vt_time(dt,dtl,did)
c     warm-up: one more call of channel ich, CPU time dt, of which dtl
c     in the loops (if computed, did)
      implicit none
      include 'ngraphs.inc'
      double precision dt,dtl
      logical did
      double precision vt_c(4,N_MAX_CG),vt_p(N_MAX_CG),vt_n(N_MAX_CG)
      double precision vt_g(4,4,N_MAX_CG),vt_b(4,N_MAX_CG)
      double precision vt_yy(N_MAX_CG)
      double precision vt_tl(N_MAX_CG),vt_tt(N_MAX_CG),vt_nl(N_MAX_CG)
      integer vt_mode,vt_s1,vt_s2
      common/c_vtrick/vt_c,vt_p,vt_n,vt_g,vt_b,vt_yy,vt_tl,vt_tt,vt_nl,
     &     vt_mode,vt_s1,vt_s2
c     sums of the current warm-up iteration: residual, calls; loops
c     computed in the warm-up
      double precision vt_di(N_MAX_CG),vt_ni(N_MAX_CG),vt_nlw(N_MAX_CG)
      common/c_vtrick_w/vt_di,vt_ni,vt_nlw
      integer ich
      common/comich/ich
      vt_n(ich)=vt_n(ich)+1d0
      vt_ni(ich)=vt_ni(ich)+1d0
      vt_tt(ich)=vt_tt(ich)+dt
      if(did)then
         vt_tl(ich)=vt_tl(ich)+dtl
         vt_nlw(ich)=vt_nlw(ich)+1d0
      endif
      return
      end


      subroutine vt_iter_start(ich)
c     warm-up: start of an iteration of channel ich
      implicit none
      include 'ngraphs.inc'
      integer ich
      double precision vt_c(4,N_MAX_CG),vt_p(N_MAX_CG),vt_n(N_MAX_CG)
      double precision vt_g(4,4,N_MAX_CG),vt_b(4,N_MAX_CG)
      double precision vt_yy(N_MAX_CG)
      double precision vt_tl(N_MAX_CG),vt_tt(N_MAX_CG),vt_nl(N_MAX_CG)
      integer vt_mode,vt_s1,vt_s2
      common/c_vtrick/vt_c,vt_p,vt_n,vt_g,vt_b,vt_yy,vt_tl,vt_tt,vt_nl,
     &     vt_mode,vt_s1,vt_s2
c     sums of the current warm-up iteration: residual, calls; loops
c     computed in the warm-up
      double precision vt_di(N_MAX_CG),vt_ni(N_MAX_CG),vt_nlw(N_MAX_CG)
      common/c_vtrick_w/vt_di,vt_ni,vt_nlw
      vt_di(ich)=0d0
      vt_ni(ich)=0d0
      return
      end


      subroutine vt_iter_end(ich,itw,s,sr,d,popt,tl,tr)
c     warm-up: end of iteration itw of channel ich, with standard
c     deviation per call s (with the trick, probability vt_p(ich)).
c     Returns sr = S (s^2 without the trick), d = D (with the reference
c     used in the iteration; after the first, with the new fit), popt
c     (vt_popt_f), tl, tr: CPU time per loop and per call for the
c     rest. The reference (vt_c) is fitted again on all the points so
c     far
      implicit none
      include 'ngraphs.inc'
      integer ich,itw
      double precision s,sr,d,popt,tl,tr
      double precision vt_c(4,N_MAX_CG),vt_p(N_MAX_CG),vt_n(N_MAX_CG)
      double precision vt_g(4,4,N_MAX_CG),vt_b(4,N_MAX_CG)
      double precision vt_yy(N_MAX_CG)
      double precision vt_tl(N_MAX_CG),vt_tt(N_MAX_CG),vt_nl(N_MAX_CG)
      integer vt_mode,vt_s1,vt_s2
      common/c_vtrick/vt_c,vt_p,vt_n,vt_g,vt_b,vt_yy,vt_tl,vt_tt,vt_nl,
     &     vt_mode,vt_s1,vt_s2
c     sums of the current warm-up iteration: residual, calls; loops
c     computed in the warm-up
      double precision vt_di(N_MAX_CG),vt_ni(N_MAX_CG),vt_nlw(N_MAX_CG)
      common/c_vtrick_w/vt_di,vt_ni,vt_nlw
      double precision dit,dfit,p,vt_popt_f
      external vt_popt_f
      dit=-1d0
      if(itw.gt.1.and.vt_ni(ich).gt.0d0)dit=vt_di(ich)/vt_ni(ich)
      p=vt_p(ich)
      sr=s*s
      if(dit.ge.0d0.and.p.lt.1d0)
     &     sr=max(s*s-(1d0/p-1d0)*dit,0.1d0*s*s)
      call vt_fit(vt_g(1,1,ich),vt_b(1,ich),vt_yy(ich),vt_c(1,ich),
     &     dfit)
      if(vt_n(ich).gt.0d0)dfit=dfit/vt_n(ich)
      d=dfit
      if(dit.ge.0d0)d=dit
      tl=0d0
      tr=0d0
      if(vt_nlw(ich).gt.0d0)tl=vt_tl(ich)/vt_nlw(ich)
      if(vt_n(ich).gt.0d0)tr=max((vt_tt(ich)-vt_tl(ich))/vt_n(ich),0d0)
      popt=vt_popt_f(sr,d,tl,tr,vt_n(ich))
      return
      end


      double precision function vt_popt_f(sr,d,tl,tr,n)
c     probability of the loop that minimises (variance x CPU time):
c     sr = S, d = D, tl, tr: CPU times per loop and per call for the
c     rest, n: calls of the sums; 1 (no trick) with fewer than 1000
c     calls, or if the trick does not pay off (above 0.7)
      implicit none
      double precision sr,d,tl,tr,n
      double precision vt_pmin
      parameter (vt_pmin=0.02d0)
      vt_popt_f=1d0
      if(n.lt.1000d0.or.d.ge.sr.or.tl.le.0d0)return
      vt_popt_f=min(1d0,max(vt_pmin,sqrt(d/(sr-d)*tr/tl)))
      if(vt_popt_f.gt.0.7d0)vt_popt_f=1d0
      return
      end


      subroutine vt_apply(v0,x,loop)
c     the loop v0 of a point of channel ich replaced by R + (v0 - R)/p
c     if it was computed (loop), by R otherwise (p = vt_p(ich) < 1)
      implicit none
      include 'ngraphs.inc'
      double precision v0,x(4)
      logical loop
      double precision vt_c(4,N_MAX_CG),vt_p(N_MAX_CG),vt_n(N_MAX_CG)
      double precision vt_g(4,4,N_MAX_CG),vt_b(4,N_MAX_CG)
      double precision vt_yy(N_MAX_CG)
      double precision vt_tl(N_MAX_CG),vt_tt(N_MAX_CG),vt_nl(N_MAX_CG)
      integer vt_mode,vt_s1,vt_s2
      common/c_vtrick/vt_c,vt_p,vt_n,vt_g,vt_b,vt_yy,vt_tl,vt_tt,vt_nl,
     &     vt_mode,vt_s1,vt_s2
      integer ich
      common/comich/ich
      double precision r
      integer k
      if(vt_p(ich).ge.1d0)return
      r=0d0
      do k=1,4
         r=r+vt_c(k,ich)*x(k)
      enddo
      if(loop)then
         v0=r+(v0-r)/vt_p(ich)
      else
         v0=r
      endif
      return
      end


      double precision function vt_pmain(popt,ncl,nlmin)
c     probability of the loop in a main run with ncl calls per iteration
c     of the channel: popt (vt_popt_f), but at least nlmin loops per
c     iteration; 1 (no trick) above 0.7
      implicit none
      double precision popt
      integer ncl,nlmin
      vt_pmain=1d0
      if(popt.ge.1d0)return
      vt_pmain=max(popt,min(1d0,dble(nlmin)/dble(max(ncl,1))))
      if(vt_pmain.gt.0.7d0)vt_pmain=1d0
      return
      end


      subroutine vt_fit(g,b,yy,c,d)
c     least squares: c minimises sum (y - sum_k c_k x_k)^2, from the
c     sums g = sum x x^T, b = sum x y, yy = sum y^2; d is the minimum.
c     Cholesky decomposition of g (columns normalised), the columns
c     that depend on the previous ones (pivot < 1d-10) dropped, e.g.
c     the double pole proportional to the Born for a single Born colour
c     structure
      implicit none
      double precision g(4,4),b(4),yy,c(4),d
      double precision a(4,4),sc(4),z(4),cs(4),t
      logical keep(4)
      integer j,k,l
      do k=1,4
         sc(k)=sqrt(max(g(k,k),0d0))
         keep(k)=sc(k).gt.0d0
         z(k)=0d0
         cs(k)=0d0
         c(k)=0d0
      enddo
      a=0d0
      do k=1,4
         if(.not.keep(k))cycle
         t=1d0
         do j=1,k-1
            if(keep(j))t=t-a(k,j)**2
         enddo
         if(t.le.1d-10)then
            keep(k)=.false.
            cycle
         endif
         a(k,k)=sqrt(t)
         do l=k+1,4
            if(.not.keep(l))cycle
            t=g(l,k)/(sc(l)*sc(k))
            do j=1,k-1
               if(keep(j))t=t-a(l,j)*a(k,j)
            enddo
            a(l,k)=t/a(k,k)
         enddo
      enddo
c     L z = b (normalised), then L^T c = z
      d=yy
      do k=1,4
         if(.not.keep(k))cycle
         t=b(k)/sc(k)
         do j=1,k-1
            if(keep(j))t=t-a(k,j)*z(j)
         enddo
         z(k)=t/a(k,k)
         d=d-z(k)**2
      enddo
      do k=4,1,-1
         if(.not.keep(k))cycle
         t=z(k)
         do l=k+1,4
            if(keep(l))t=t-a(l,k)*cs(l)
         enddo
         cs(k)=t/a(k,k)
         c(k)=cs(k)/sc(k)
      enddo
      d=max(d,0d0)
      return
      end


      subroutine vt_seed(iseed,k)
c     random numbers of the virtual trick (vt_rand), from the seed;
c     k = 1: warm-up, 2: main runs (the same with --warmup and --main)
      implicit none
      include 'ngraphs.inc'
      integer iseed,k
      double precision vt_c(4,N_MAX_CG),vt_p(N_MAX_CG),vt_n(N_MAX_CG)
      double precision vt_g(4,4,N_MAX_CG),vt_b(4,N_MAX_CG)
      double precision vt_yy(N_MAX_CG)
      double precision vt_tl(N_MAX_CG),vt_tt(N_MAX_CG),vt_nl(N_MAX_CG)
      integer vt_mode,vt_s1,vt_s2
      common/c_vtrick/vt_c,vt_p,vt_n,vt_g,vt_b,vt_yy,vt_tl,vt_tt,vt_nl,
     &     vt_mode,vt_s1,vt_s2
c     sums of the current warm-up iteration: residual, calls; loops
c     computed in the warm-up
      double precision vt_di(N_MAX_CG),vt_ni(N_MAX_CG),vt_nlw(N_MAX_CG)
      common/c_vtrick_w/vt_di,vt_ni,vt_nlw
      integer*8 is
      is=abs(iseed)+1000003_8*k
      vt_s1=int(1+mod(7919_8*is+12345_8,2147483562_8))
      vt_s2=int(1+mod(104729_8*is+67890_8,2147483398_8))
      return
      end


      double precision function vt_rand()
c     uniform random number in (0,1) for the virtual trick: L'Ecuyer's
c     combination of two multiplicative congruential generators
c     (Schrage's method; ran2 without the shuffle), separate from the
c     random numbers of vegas
      implicit none
      include 'ngraphs.inc'
      double precision vt_c(4,N_MAX_CG),vt_p(N_MAX_CG),vt_n(N_MAX_CG)
      double precision vt_g(4,4,N_MAX_CG),vt_b(4,N_MAX_CG)
      double precision vt_yy(N_MAX_CG)
      double precision vt_tl(N_MAX_CG),vt_tt(N_MAX_CG),vt_nl(N_MAX_CG)
      integer vt_mode,vt_s1,vt_s2
      common/c_vtrick/vt_c,vt_p,vt_n,vt_g,vt_b,vt_yy,vt_tl,vt_tt,vt_nl,
     &     vt_mode,vt_s1,vt_s2
      integer k,z
      k=vt_s1/53668
      vt_s1=40014*(vt_s1-k*53668)-k*12211
      if(vt_s1.lt.0)vt_s1=vt_s1+2147483563
      k=vt_s2/52774
      vt_s2=40692*(vt_s2-k*52774)-k*3791
      if(vt_s2.lt.0)vt_s2=vt_s2+2147483399
      z=vt_s1-vt_s2
      if(z.lt.1)z=z+2147483562
      vt_rand=z*4.656613057d-10
      return
      end
