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
            call vegas(region,ndim,int_virtual,init,nclw_ch(ich),1,
     &      nprn,res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
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
      enddo
      call cpu_time(t_w1)
      cpu_call=0d0
      if(ncalls_w.gt.0d0)cpu_call=(t_w1-t_w0)/ncalls_w
      do ich=1,nchan
         s_mn(ich)=s_ch(ich)
         if(nitVth.ge.2)s_mn(ich)=0.5d0*(s_ch(ich)+s_prev(ich))
      enddo
      if(.not.do_main)then
c        warm-up only: save the state and stop
         call ran2_get_state(rstate)
         open(unit=iuw,file='warmup_'//trim(tag)//'.dat',
     &        status='replace')
         write(iuw,*)'MADNKLO_WARMUP 1'
         write(iuw,*)nchan,ndim,ndmx,nclV,nitV,ncalls_w,cpu_call
         do ich=1,nchan
            write(iuw,'(es26.17e3,i6)')s_mn(ich),ndo_ch(ich)
            write(iuw,'(4es26.17e3)')
     &           ((xi_ch(j,k,ich),j=1,ndmx),k=1,ndim)
         enddo
         write(iuw,'(5i14)')rstate
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
      read(iuw,*)
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
      close(iuw)
      call ran2_set_state(rstate)
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
         call vegas(region,ndim,int_virtual,init,ncl_ch(ich),nitV,nprn,
     &   res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
         call sv_add_channel()
         write(iu8,*)trim(tag)//': channel, itns, calls = ',ich,nitV,
     &        ncl_ch(ich)
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
