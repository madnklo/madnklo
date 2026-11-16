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
c     command-line arguments: '--test_only' (tests only, see below) and
c     '--counterterms' (I + J + C only)
      test_only=.false.
      counterterms=.false.
      do iarg=1,command_argument_count()
         call get_command_argument(iarg,arg)
         if(trim(arg).eq.'--test_only')test_only=.true.
         if(trim(arg).eq.'--counterterms')counterterms=.true.
      enddo
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
      open(unit=iu1,file='integration_'//trim(tag)//'.log')
c     pole test (V + I) on the first points: needs the loop
      if(vmode.ne.2)open(unit=iu2,file='test_poles_V.log')
      open(unit=iu7,file='failures_'//trim(tag)//'.log',status='replace',
     &     action='write')
      open(unit=iu8,file=trim(tag)//'_chan.log')
      open(unit=iu,file='results_'//trim(tag)//'.log')
      line='=================================================='
      if(vmode.eq.0)write(iu,*)' Virtual contribution '
      if(vmode.eq.1)write(iu,*)' Virtual contribution, loop only '
      if(vmode.eq.2)write(iu,*)' Virtual contribution, I + J + C only '
c
c     quickly get integration error per channel so to modulate
c     number of points thrown per channel in the main loop
c     (probe: npoints_V_grid/5 calls per channel, at least 500; with
c     the former floor of 10000 the probe dominated the run time of
c     slow loops at low statistics)
      nclVth0=max(500,int(nclVth/5d0))
      nitVth0=max(5,int(nitVth/2d0))
c     channels: the phase-space configurations (configs.f)
      nchan=n_born_channels()
      sum_err_v_a=0d0
      do ich=1,nchan
         init=0
         doplot=.false.
         failure_stage='probe'
         call vegas(region,ndim,int_virtual,init,nclVth0,nitVth0,nprn,
     &   res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
         err_v_a(ich) = err_v
         sum_err_v_a = sum_err_v_a + err_v_a(ich)
      enddo
c
c     main loop over channels
      do ich=1,nchan
         write(*,*)'Virtual warmup for channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' VIRTUAL WARMUP, CHANNEL',ich
         write(iu1,*)'============================='
         init=0
         doplot=.false.
         nclVth1=max(1000,int(nclVth*err_v_a(ich)/sum_err_v_a))
         failure_stage='warmup'
         call vegas(region,ndim,int_virtual,init,nclVth1,nitVth,nprn,
     &   res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu8,*)trim(tag)//' warmup: channel, itns, calls = ',ich,
     &        nitVth,nclVth1
c
         write(*,*)'Virtual for channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' VIRTUAL, CHANNEL',ich
         write(iu1,*)'============================='
         init=1
         doplot=.true.
         nclV1=max(1000,int(nclV*err_v_a(ich)/sum_err_v_a))
         failure_stage='main'
         call vegas(region,ndim,int_virtual,init,nclV1,nitV,nprn,
     &   res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
         call sv_add_channel()
         write(iu8,*)trim(tag)//': channel, itns, calls = ',ich,nitV,
     &        nclV1
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
