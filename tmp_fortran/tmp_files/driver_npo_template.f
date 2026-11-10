      program driver_%(isec)d_%(jsec)d
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      include 'colored_partons.inc'
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
      integer iu,iu1,iu7,iu8,iu9
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
c
c     phase-space dimension, same for all contributions to this folder
      ndim=3*(nexternal-2)-4
      do i=1,2
         if(ISNLOQCDPARTON(i)) ndim = ndim + 1
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
c     initialise histograms and open output files
c      call histo_init
c     renormalisation-scale variations (rw_rscale in run_card.dat):
c     one histogram weight per scale factor
      call sv_init()
      call sv_weights_info(nwgt,weights_info,100)
      call analysis_begin(nwgt,weights_info)
      open(unit=iu1,file='integration_R_%(isec)d_%(jsec)d.log')
      open(unit=iu7,file='failures_R_%(isec)d_%(jsec)d.log',status='replace',action='write')
      open(unit=iu8,file='testR_%(isec)d_%(jsec)d.log')
      open(unit=iu9,file='chan_R_%(isec)d_%(jsec)d.log')
      open(unit=iu,file='results_R_%(isec)d_%(jsec)d.log')
      line='=================================================='
      write(iu9,*)' Real contribution '
c
c     quickly get integration error per channel so to modulate
c     number of points thrown per channel in the main loop
c      nclRth0=max(10000,int(nclRth/5d0))
c      nitRth0=max(5,int(nitRth/2d0))
c      sum_err_r_a=0d0
c      do i=1,N_MAX_CG
c         ich=i
c         init=0
c         doplot=.false.
c         call vegas(region,ndim,int_real_%(isec)d_%(jsec)d,init,nclRth0,nitRth0,nprn,res_r,err_r,chi2a,acc,xi,it,ndo,si,swgt,schi)
c         err_r_a(ich) = err_r
c         sum_err_r_a = sum_err_r_a + err_r_a(ich)
c      enddo
c
c     single phase-space channel for the Born (ich=1). With the
c     channel decomposition (run_card channel_decomposition), one
c     integration per channel ichan (momentum mapping of a
c     counterterm piece), with its own vegas grid: warm-ups of all
c     channels first (nclRth/nch calls each), then the main runs,
c     the nclR calls being shared among the channels in proportion
c     to their warm-up errors (optimal for the total error). The
c     results, histograms and scale variations are summed over the
c     channels.
      ich=1
      do ichan=1,nch
         write(*,*)'Real %(isec)d%(jsec)d warmup for channel',ichan
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' REAL_%(isec)d_%(jsec)d WARMUP, CHANNEL',ichan
         write(iu1,*)'============================='
         init=0
         doplot=.false.
         failure_stage='warmup'
         call vegas(region,ndim,int_real_%(isec)d_%(jsec)d,init,max(1,nclRth/nch),nitRth,nprn,res_r,err_r,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu9,*)'R%(isec)d%(jsec)d warmup: channel, itns, calls = ',ichan,nitRth,max(1,nclRth/nch)
         err_ch(ichan)=err_r
         if(nch.gt.1)then
            do k=1,ndim
               do j=1,ndmx
                  xi_ch(j,k,ichan)=xi(j,k)
               enddo
            enddo
            ndo_ch(ichan)=ndo
         endif
      enddo
c     calls of the main runs
      if(nch.eq.1)then
         ncl_ch(1)=nclR
      else
         sum_err_ch=0d0
         do ichan=1,nch
            sum_err_ch=sum_err_ch+err_ch(ichan)
         enddo
         do ichan=1,nch
            if(sum_err_ch.gt.0d0)then
               ncl_ch(ichan)=max(nclR/(10*nch),nint(nclR*err_ch(ichan)/sum_err_ch))
            else
               ncl_ch(ichan)=max(1,nclR/nch)
            endif
         enddo
      endif
      do ichan=1,nch
         if(nch.gt.1)then
            do k=1,ndim
               do j=1,ndmx
                  xi(j,k)=xi_ch(j,k,ichan)
               enddo
            enddo
            ndo=ndo_ch(ichan)
         endif
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
c      call histo_final('plot_R_%(isec)d_%(jsec)d.dat',rescale_plot_R)
c      write(iu,*)
c      write(iu,*)' '//line
c      write(iu,*)
      write(iu,*)' sigma R%(isec)d%(jsec)d [pb]  = ',sum_r,' +-',sum_err_r
      call sv_write_results(iu,'R%(isec)d%(jsec)d',sum_r,sum_err_r)
c      write(iu,*)
c      write(iu,*)' '//line
c      write(iu,*)
      close(iu)
      close(iu1)
      if(nfail.gt.0)write(iu7,*)'Total ierr failures: ',nfail
      close(iu7)
      close(iu8)
      close(iu9)
c
      end
