      program driver_%(isec)d_%(jsec)d
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      include 'colored_partons.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'ngraphs_%(UBgraphs)s.inc'
      integer mxdim
      parameter(mxdim=30)
c     test-only mode (command-line argument --test_only)
      logical test_only
      character*40 arg
      integer iarg
      double precision xtest(mxdim),wtest,ftest,ran2
      external ran2
      integer ntest_r
      parameter(ntest_r=20)
      integer ndim,i,j,idum
      double precision s_had
      integer iu,iu1,iu7,iu8,iu9,iuw
      common/cdim/ndim
      double precision int_real_%(isec)d_%(jsec)d
      double precision err_r,res_r
      external int_real_%(isec)d_%(jsec)d
      common/ciunitNLO/iu8
      integer order
      logical doplot
      common/cdoplot/doplot
      double precision rescale_plot_R
      character*100 line
      integer nitRth,nclRth,nitR,nclR
      integer nitRth0,nclRth0,nclR0,nclRth1,nclR1
      integer iseed
      common /to_seed/iseed
      COMMON/iterations/NITR
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
      double precision sum_r,sum_err_r
c     double precision sum_err_r_a,err_r_a(N_MAX_CG)
      integer nwgt
      character*80 weights_info(100)
c     channel decomposition (run_card channel_decomposition)
      integer maxchan
      parameter(maxchan=100)
      integer ichan,nch,k
      common/c_nlo_chan_cur/ichan
      double precision err_ch(maxchan),sum_err_ch
      double precision xi_ch(ndmx,mxdim,maxchan)
      integer ndo_ch(maxchan),ncl_ch(maxchan)
c     warm-up iteration by iteration over the channels, and phases
c     (see below)
      integer nclw_ch(maxchan),itw,nflw,ncl_tot,nit_w,nch_w,ndim_w
      double precision s_ch(maxchan),s_prev(maxchan),s_main(maxchan)
      double precision sum_s,ncalls_w,t_w0,t_w1,cpu_call
      double precision rcalls
      common/abchia/rcalls
      integer rstate(35)
      logical do_warm,do_main
      character*7 fpos
      integer ntested_c,nchecked_c
      logical protect_rng_c
      common/c_tests_R/ntested_c,nchecked_c,protect_rng_c
      INCLUDE 'all_sector_list.inc'
c
      nfail=0
      sum_r=0d0
      sum_err_r=0d0
      res_r=0d0
      err_r=0d0
c
      call SETPARA('param_card.dat')
      call SETRUN('run_card.dat')
c     hadron beams: lower bound on x1 x2 of the Born generation
      if(lpp(1).ne.0.or.lpp(2).ne.0)call set_tau_min_npo(nexternal,leg_pdgs)
c
c     read inputs
      region=0d0
      order=1
      idum = -max(abs(iseed),1)
      s_had = (EBEAM(1)+EBEAM(2))**2
      NITRTH = NITERS_R_GRID
      NCLRTH = max(1,int(NPOINTS_R_GRID/dble(lensectors)))
      NITR = NITERS_R
      NCLR = max(1,int(NPOINTS_R/dble(lensectors)))
c     TODO: understand muR input fixed/dyn scale
c
c     initialise physics parameters and set sector parametrisation
      iu1=44
      iu=55
      iu7=77
      iu8=88
      iu9=99
      iuw=97
c
c     phase-space dimension, same for all contributions to this folder:
c     one more variable per hadron beam (momentum fractions, as in
c     gen_mom and driver_n.f)
      ndim=3*(nexternal-2)-4
      do i=1,2
         if(lpp(i).ne.0) ndim = ndim + 1
      enddo
      do i=1,ndim
         region(i)=0d0
         region(i+ndim)=1d0
      enddo
c
c     channel decomposition: list of channels (momentum mappings)
      nch=1
      ichan=1
      if(channel_decomposition)call setup_channels_%(isec)d_%(jsec)d(nch)
c
c     test-only mode: './sector_i_j --test_only' only runs the tests
      test_only=.false.
      do iarg=1,command_argument_count()
         call get_command_argument(iarg,arg)
         if(trim(arg).eq.'--test_only')test_only=.true.
      enddo
      if(test_only)then
c        IR-limit tests of the first ntest_r calls of the integrand
c        (testR_i_j.log), without integration
         open(unit=iu8,file='testR_%(isec)d_%(jsec)d.log')
         open(unit=iu7,file='failures_R_%(isec)d_%(jsec)d_test.log',status='replace')
         ich=1
         doplot=.false.
         failure_stage='test'
         do i=1,ntest_r
            do j=1,ndim
               xtest(j)=ran2(idum)
            enddo
            wtest=1d0
            ftest=int_real_%(isec)d_%(jsec)d(xtest,wtest)
         enddo
         close(iu8)
         close(iu7)
         write(*,*)'Test only: IR limits written to testR_%(isec)d_%(jsec)d.log'
         stop
      endif
c
c     phases: './sector_%(isec)d_%(jsec)d --warmup' only runs the
c     warm-up and writes warmup_R_%(isec)d_%(jsec)d.dat (grids,
c     standard deviation per point of each channel, cost per call,
c     state of the random numbers); './sector_%(isec)d_%(jsec)d --main'
c     reads it and runs the main integration with the calls per
c     iteration of main_calls_R_%(isec)d_%(jsec)d.dat, written by
c     run_all.py, which shares the calls among all the sectors of the
c     contribution (optimal for the total error at fixed CPU time);
c     without these arguments, both, with npoints_R/(number of sectors)
c     calls per iteration. --warmup followed by --main with these calls
c     gives exactly the same results.
      do_warm=.true.
      do_main=.true.
      do iarg=1,command_argument_count()
         call get_command_argument(iarg,arg)
         if(trim(arg).eq.'--warmup')do_main=.false.
         if(trim(arg).eq.'--main')do_warm=.false.
      enddo
      fpos='rewind'
      if(.not.do_warm)fpos='append'
c
c     initialise histograms and open output files
c      call histo_init
c     scale variations (rw_rscale, rw_fscale in run_card.dat): one
c     histogram weight per (muF,muR) pair; with hadron beams, muF is
c     varied by PDF reweighting:
c     R and the counterterms set their PDF ratios in lumi_npo
      call sv_allow_muf()
      call sv_init()
      call sv_weights_info(nwgt,weights_info,100)
      call analysis_begin(nwgt,weights_info)
      open(unit=iu1,file='integration_R_%(isec)d_%(jsec)d.log',position=fpos)
      if(do_warm)then
         open(unit=iu7,file='failures_R_%(isec)d_%(jsec)d.log',status='replace',action='write')
      else
         open(unit=iu7,file='failures_R_%(isec)d_%(jsec)d.log',position='append',action='write')
      endif
      open(unit=iu8,file='testR_%(isec)d_%(jsec)d.log',position=fpos)
      open(unit=iu9,file='chan_R_%(isec)d_%(jsec)d.log',position=fpos)
      if(do_main)open(unit=iu,file='results_R_%(isec)d_%(jsec)d.log')
      line='=================================================='
      if(do_warm)write(iu9,*)' Real contribution '
c
c     single phase-space channel for the Born (ich=1). With the
c     channel decomposition (run_card channel_decomposition), one
c     integration per channel ichan (momentum mapping of a
c     counterterm piece), with its own vegas grid.
c     Warm-up, iteration by iteration over the channels: each channel
c     runs one vegas iteration on its own grid (flat at the first one),
c     which refines the grid; its standard deviation per point on that
c     grid, s = (error of the iteration)*sqrt(calls), is measured. With
c     more than one channel, the nclRth calls of an iteration are
c     shared equally at the first one, then in proportion to s of the
c     previous one, at least max(1000,nclRth/(10 nch)) each (with
c     fewer calls per iteration, the combination of the iterations by
c     vegas is biased for integrands with large fluctuations).
c     Main runs: the calls are shared among the channels in proportion
c     to s averaged over the last two warm-up iterations, at least
c     max(1000,calls/(10 nch)) each; each main run continues on the
c     grid of its warm-up. The results, histograms and scale variations
c     are summed over the channels.
      ich=1
      if(do_warm)then
         call cpu_time(t_w0)
         ncalls_w=0d0
         nflw=max(1000,nclRth/(10*nch))
         do ichan=1,nch
            nclw_ch(ichan)=nclRth
            if(nch.gt.1)nclw_ch(ichan)=max(1000,nclRth/nch)
            s_ch(ichan)=0d0
            s_prev(ichan)=0d0
         enddo
         do itw=1,max(nitRth,1)
            do ichan=1,nch
               if(itw.eq.1)then
                  init=0
                  write(*,*)'Real %(isec)d%(jsec)d warmup for channel',ichan
                  write(iu1,*)
                  write(iu1,*)'============================='
                  write(iu1,*)' REAL_%(isec)d_%(jsec)d WARMUP, CHANNEL',ichan
                  write(iu1,*)'============================='
               else
                  do k=1,ndim
                     do j=1,ndmx
                        xi(j,k)=xi_ch(j,k,ichan)
                     enddo
                  enddo
                  ndo=ndo_ch(ichan)
                  init=1
               endif
               doplot=.false.
               failure_stage='warmup'
               call vegas(region,ndim,int_real_%(isec)d_%(jsec)d,init,nclw_ch(ichan),1,nprn,res_r,err_r,chi2a,acc,xi,it,ndo,si,swgt,schi)
               do k=1,ndim
                  do j=1,ndmx
                     xi_ch(j,k,ichan)=xi(j,k)
                  enddo
               enddo
               ndo_ch(ichan)=ndo
               ncalls_w=ncalls_w+rcalls
               s_prev(ichan)=s_ch(ichan)
               s_ch(ichan)=err_r*dsqrt(rcalls)
               if(.not.(s_ch(ichan).ge.0d0.and.s_ch(ichan).lt.huge(1d0)))s_ch(ichan)=0d0
               write(iu9,*)'R%(isec)d%(jsec)d warmup: channel, itn, calls, sd/point = ',ichan,itw,nclw_ch(ichan),s_ch(ichan)
            enddo
c           calls of the next warm-up iteration
            if(nch.gt.1)then
               sum_s=0d0
               do ichan=1,nch
                  sum_s=sum_s+s_ch(ichan)
               enddo
               if(sum_s.gt.0d0)then
                  do ichan=1,nch
                     nclw_ch(ichan)=max(nflw,nint(nclRth*s_ch(ichan)/sum_s))
                  enddo
               endif
            endif
         enddo
         call cpu_time(t_w1)
         cpu_call=0d0
         if(ncalls_w.gt.0d0)cpu_call=(t_w1-t_w0)/ncalls_w
         do ichan=1,nch
            s_main(ichan)=s_ch(ichan)
            if(nitRth.ge.2)s_main(ichan)=0.5d0*(s_ch(ichan)+s_prev(ichan))
         enddo
         if(.not.do_main)then
c           warm-up only: save the state and stop
            call ran2_get_state(rstate)
            open(unit=iuw,file='warmup_R_%(isec)d_%(jsec)d.dat',status='replace')
            write(iuw,*)'MADNKLO_WARMUP 1'
            write(iuw,*)nch,ndim,ndmx,nclR,nitR,ncalls_w,cpu_call
            do ichan=1,nch
               write(iuw,'(es26.17e3,i6)')s_main(ichan),ndo_ch(ichan)
               write(iuw,'(4es26.17e3)')((xi_ch(j,k,ichan),j=1,ndmx),k=1,ndim)
            enddo
            write(iuw,'(5i14)')rstate
            close(iuw)
            close(iu1)
            close(iu7)
            close(iu8)
            close(iu9)
            write(*,*)'Warm-up written to warmup_R_%(isec)d_%(jsec)d.dat'
            stop
         endif
      else
c        main only: state of the warm-up
         open(unit=iuw,file='warmup_R_%(isec)d_%(jsec)d.dat',status='old')
         read(iuw,*)
         read(iuw,*)nch_w,ndim_w,j,k,nit_w,ncalls_w,cpu_call
         if(nch_w.ne.nch.or.ndim_w.ne.ndim.or.j.ne.ndmx)then
            write(*,*)'warmup_R_%(isec)d_%(jsec)d.dat does not match this executable'
            stop 1
         endif
         do ichan=1,nch
            read(iuw,*)s_main(ichan),ndo_ch(ichan)
            read(iuw,*)((xi_ch(j,k,ichan),j=1,ndmx),k=1,ndim)
         enddo
         read(iuw,*)rstate
         close(iuw)
         call ran2_set_state(rstate)
c        the channel checks of the first points have been done in the
c        warm-up phase; the limit tests are run again (they initialise
c        the matrix elements) without changing the random numbers
         nchecked_c=ntest_r
         protect_rng_c=.true.
      endif
c
c     calls per iteration of the main runs: npoints_R/(number of
c     sectors), or those given by run_all.py (main_calls_R_i_j.dat)
      ncl_tot=nclR
      if(.not.do_warm)then
         open(unit=iuw,file='main_calls_R_%(isec)d_%(jsec)d.dat',status='old',err=11)
         read(iuw,*)ncl_tot
         close(iuw)
 11      continue
      endif
      if(nch.eq.1)then
         ncl_ch(1)=ncl_tot
      else
         sum_s=0d0
         do ichan=1,nch
            sum_s=sum_s+s_main(ichan)
         enddo
         do ichan=1,nch
            if(sum_s.gt.0d0)then
               ncl_ch(ichan)=max(1000,ncl_tot/(10*nch),nint(ncl_tot*s_main(ichan)/sum_s))
            else
               ncl_ch(ichan)=max(1000,ncl_tot/nch)
            endif
         enddo
      endif
      do ichan=1,nch
         do k=1,ndim
            do j=1,ndmx
               xi(j,k)=xi_ch(j,k,ichan)
            enddo
         enddo
         ndo=ndo_ch(ichan)
         write(*,*)'Real %(isec)d%(jsec)d for channel',ichan
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' REAL_%(isec)d_%(jsec)d, CHANNEL',ichan
         write(iu1,*)'============================='
         init=1
         doplot=.true.
         failure_stage='main'
         call vegas(region,ndim,int_real_%(isec)d_%(jsec)d,init,ncl_ch(ichan),nitR,nprn,res_r,err_r,chi2a,acc,xi,it,ndo,si,swgt,schi)
         call sv_add_channel()
         rescale_plot_R=dble(nitR)/min(dble(nitR),dble(it))
         sum_r = sum_r + res_r
         sum_err_r = sum_err_r + err_r**2
         write(iu9,*)'R%(isec)d%(jsec)d: channel, itns, calls = ',ichan,nitR,ncl_ch(ichan)
         write(iu9,*)' sigma R%(isec)d_%(jsec)d [pb], channel',ichan,' = ',res_r,' +-',err_r
         write(iu9,*)
      enddo
c
      write(*,*)'...done'
c
c     finalise histograms and output files
      call analysis_end(1d0,'plot_R_%(isec)d_%(jsec)d.dat')
      sum_err_r = dsqrt(sum_err_r)
      write(iu,*)' sigma R%(isec)d%(jsec)d [pb]  = ',sum_r,' +-',sum_err_r
      call sv_write_results(iu,'R%(isec)d%(jsec)d',sum_r,sum_err_r)
      close(iu)
      close(iu1)
      if(nfail.gt.0)write(iu7,*)'Total ierr failures: ',nfail
      close(iu7)
      close(iu8)
      close(iu9)
c
      end
