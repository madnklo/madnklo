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
      integer ndim,i,j,idum
      double precision s_had
      integer iu,iu1,iu7,iu8
      common/cdim/ndim
      double precision int_Born
      double precision res_B,err_B
      external int_Born
      integer order
      logical doplot
      common/cdoplot/doplot
      double precision rescale_plot_B
      character*100 line
      integer nitBth,nclBth,nitB,nclB
      integer nitBth0,nclBth0,nclB0,nclBth1,nclB1
      integer iseed
      common /to_seed/iseed
      common/iterations/nitB
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
      double precision sum_b,sum_err_b
      double precision sum_err_b_a,err_b_a(N_MAX_CG)
      double precision xi_ch(ndmx,mxdim,N_MAX_CG)
      integer ndo_ch(N_MAX_CG),ncl_ch(N_MAX_CG),nclw_ch(N_MAX_CG)
      integer k,itw,nflw
      double precision s_ch(N_MAX_CG),s_prev(N_MAX_CG),s_main,sum_s
      double precision rcalls
      common/abchia/rcalls
      integer nwgt
      character*80 weights_info(100)
c
      nfail=0
      sum_b=0d0
      sum_err_b=0d0
      res_b=0d0
      err_b=0d0
c
      call SETPARA('param_card.dat')
      call SETRUN('run_card.dat')
c
c     read inputs
      region=0d0
      order=0
      idum = -max(abs(iseed),1)
      s_had = (EBEAM(1)+EBEAM(2))**2
      NITBTH = NITERS_B_GRID
      NCLBTH = NPOINTS_B_GRID
      NITB = NITERS_B
      NCLB = NPOINTS_B
c     TODO: understand muR input fixed/dyn scale
c
c     initialise physics parameters
      iu1=44
      iu=55
      iu7=77
      iu8=88
c
c     phase-space dimension, same for all contributions to this folder
      ndim=max(3*(nexternal-2)-4,0)
c     plus one Bjorken x per hadron beam (same rule as gen_mom, genps.f)
      do i=1,2
         if(lpp(i).ne.0) ndim = ndim + 1
      enddo
      do i=1,ndim
         region(i)=0d0
         region(i+ndim)=1d0
      enddo
c
c     initialise histograms and open output files
c     renormalisation-scale variations (rw_rscale in run_card.dat):
c     one histogram weight per scale factor
c     LO: the PDF ratios of the muF variations are computed (LO_B.f)
      call sv_allow_muf()
      call sv_init()
      call sv_weights_info(nwgt,weights_info,100)
      call analysis_begin(nwgt,weights_info)
      open(unit=iu1,file='integration_B.log')
      open(unit=iu7,file='failures_B.log',status='replace',action='write')
      open(unit=iu8 ,file='B_chan.log')
      open(unit=iu ,file='results_B.log')
      line='=================================================='
      write(iu,*)' Born contribution '
c
c     channels: the phase-space configurations (configs.f), each
c     integrating the Born times its multi-channel weight
c     (born_chan_weight).
c     Warm-up, iteration by iteration over all the channels: each
c     channel runs one vegas iteration on its own grid (flat at the
c     first one), which refines the grid; its standard deviation per
c     point on that grid, s = (error of the iteration)*sqrt(calls), is
c     measured. The npoints_B_grid calls of an iteration are shared
c     equally at the first one, then in proportion to s of the
c     previous one (optimal for the total error), at least
c     max(1000,npoints_B_grid/(10 nchan)) each, so that every channel
c     stays measured (at least 1000 calls per channel and iteration:
c     with fewer, the combination of the iterations by vegas, weighted
c     by 1/sigma^2, is biased for integrands with large fluctuations,
c     e.g. by 20% for the virtual of e+e- > 3 jets).
c     Main runs: the npoints_B calls are shared in proportion to s
c     averaged over the last two warm-up iterations, at least
c     max(1000,npoints_B/(10 nchan)) each; each main run continues on
c     the grid of its warm-up.
      nchan=n_born_channels()
      nflw=max(1000,nclBth/(10*nchan))
      do ich=1,nchan
         nclw_ch(ich)=max(1000,nclBth/nchan)
         s_ch(ich)=0d0
         s_prev(ich)=0d0
      enddo
      do itw=1,max(nitBth,1)
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
               write(iu1,*)
               write(iu1,*)'============================='
               write(iu1,*)' BORN WARMUP, CHANNEL',ich
               write(iu1,*)'============================='
            endif
            doplot=.false.
            failure_stage='warmup'
            call vegas(region,ndim,int_Born,init,nclw_ch(ich),1,nprn,
     &      res_b,err_b,chi2a,acc,xi,it,ndo,si,swgt,schi)
            do k=1,ndim
               do j=1,ndmx
                  xi_ch(j,k,ich)=xi(j,k)
               enddo
            enddo
            ndo_ch(ich)=ndo
            s_prev(ich)=s_ch(ich)
            s_ch(ich)=err_b*dsqrt(rcalls)
            if(.not.(s_ch(ich).ge.0d0.and.s_ch(ich).lt.huge(1d0)))
     &           s_ch(ich)=0d0
            write(iu8,*)'B warmup: channel, itn, calls, sd/point = ',
     &           ich,itw,nclw_ch(ich),s_ch(ich)
         enddo
c        calls of the next warm-up iteration
         sum_s=0d0
         do ich=1,nchan
            sum_s=sum_s+s_ch(ich)
         enddo
         if(nchan.gt.1.and.sum_s.gt.0d0)then
            do ich=1,nchan
               nclw_ch(ich)=max(nflw,nint(nclBth*s_ch(ich)/sum_s))
            enddo
         endif
      enddo
c     calls of the main runs
      sum_s=0d0
      do ich=1,nchan
         s_main=s_ch(ich)
         if(nitBth.ge.2)s_main=0.5d0*(s_ch(ich)+s_prev(ich))
         s_prev(ich)=s_main
         sum_s=sum_s+s_main
      enddo
      do ich=1,nchan
         if(nchan.eq.1)then
            ncl_ch(ich)=nclB
         elseif(sum_s.gt.0d0)then
            ncl_ch(ich)=max(1000,nclB/(10*nchan),
     &           nint(nclB*s_prev(ich)/sum_s))
         else
            ncl_ch(ich)=max(1000,nclB/nchan)
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
         write(*,*)'Born for channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' BORN, CHANNEL',ich
         write(iu1,*)'============================='
         init=1
         doplot=.true.
         failure_stage='main'
         call vegas(region,ndim,int_Born,init,ncl_ch(ich),nitB,nprn,
     &        res_b,err_b,chi2a,acc,xi,it,ndo,si,swgt,schi)
         call sv_add_channel()
         write(iu8,*)'B: channel, itns, calls = ',ich,nitB,ncl_ch(ich)
         rescale_plot_B=dble(nitB)/min(dble(nitB),dble(it))
         sum_b = sum_b + res_b
         sum_err_b = sum_err_b + err_b**2
         write(iu8,*)' sigma B [pb], channel',ich,' = ',
     &   res_b,' +-',err_b
         write(iu8,*)
         write(*,*)'...done'
      enddo
c
c     finalise histograms and output files
      call analysis_end(1d0,'plot_B.dat')
      sum_err_b = dsqrt(sum_err_b)
      write(iu,*)' sigma B [pb]  = ',sum_b,' +-',sum_err_b
      call sv_write_results(iu,'B',sum_b,sum_err_b)
      close(iu)
      close(iu1)
      if(nfail.gt.0)write(iu7,*)'Total ierr failures: ',nfail
      close(iu7)
c
      end
