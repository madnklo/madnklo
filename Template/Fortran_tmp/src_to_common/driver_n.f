      program driver_n_body
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'colored_partons.inc'
      INCLUDE 'ngraphs.inc'
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
      double precision sum_b,sum_err_b
      double precision sum_err_b_a,err_b_a(N_MAX_CG)
c
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
      s_had = (EBEAM(1)+EBEAM(2))**2
      NITBTH =  NITERS_FO_GRID
      NCLBTH = NPOINTS_FO_GRID
      NITB = NITERS_FO
      NCLB = NPOINTS_FO

c     TODO: understand muR input fixed/dyn scale
c
c     initialise physics parameters
      iu1=44
      iu=55
      iu7=77
      iu8=88
c
c     phase-space dimension, same for all contributions to this folder
      ndim=3*(nexternal-2)-4
      do i=1,2
         if(ISLOQCDPARTON(i)) ndim = ndim + 1
      enddo
      do i=1,ndim
         region(i)=0d0
         region(i+ndim)=1d0
      enddo
c
c     initialise histograms and open output files
c      call histo_init
      call analysis_begin(1,'central')
      open(unit=iu1,file='integration_B.log')
      open(unit=iu7,file='failures_B.log')
      open(unit=iu8 ,file='B_chan.log')
      open(unit=iu ,file='results_B.log')
      line='=================================================='
c      write(iu,*)' Born contribution '
c      write(iu,*)
c
c     quickly get integration error per channel so to modulate
c     number of points thrown per channel in the main loop
      nclBth0=max(10000,int(nclBth/5d0))
      nitBth0=max(5,int(nitBth/2d0))
      sum_err_b_a=0d0
      do i=1,N_MAX_CG
         ich=i
         init=0
         doplot=.false.
         call vegas(region,ndim,int_Born,init,nclBth0,nitBth0,nprn,
     &   res_b,err_b,chi2a,acc,xi,it,ndo,si,swgt,schi)
         err_b_a(ich) = err_b
         sum_err_b_a = sum_err_b_a + err_b_a(ich)
      enddo
c
c     main loop over channels
      do i=1,N_MAX_CG
         ich=i
         write(*,*)'Born warmup for channel',ich
         write(iu7,*)'Failures for Born warmup, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' BORN WARMUP, CHANNEL',ich
         write(iu1,*)'============================='
         init=0
         doplot=.false.
         nclBth1=max(1000,int(nclBth*err_b_a(ich)/sum_err_b_a))
         call vegas(region,ndim,int_Born,init,nclBth1,nitBth,nprn,
     &   res_b,err_b,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu8,*)'B warmup: channel, itns, calls = ',ich,nitBth,nclBth1
c
         write(*,*)'Born for channel',ich
         write(iu7,*)'Failures for Born, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' BORN, CHANNEL',ich
         write(iu1,*)'============================='
         init=1
         doplot=.true.
         nclB1=max(1000,int(nclB*err_b_a(ich)/sum_err_b_a))
         call vegas(region,ndim,int_Born,init,nclB1,nitB,nprn,
     &   res_b,err_b,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu8,*)'B: channel, itns, calls = ',ich,nitB,nclB1
         rescale_plot_B=dble(nitB)/min(dble(nitB),dble(it))
         sum_b = sum_b + res_b
         sum_err_b = sum_err_b + err_b**2
         write(iu8,*)' sigma B [pb], channel',ich,' = ',
     &   res_b,' +-',err_b
         write(iu8,*)
c     
         write(*,*)'...done'
      enddo
c
c     finalise histograms and output files
      sum_err_b = dsqrt(sum_err_b)
      call analysis_end('plot_B.dat',rescale_plot_B)
c      write(iu,*)
c      write(iu,*)' '//line
c      write(iu,*)
      write(iu,*)' sigma B [pb]  = ',sum_b,' +-',sum_err_b
c      write(iu,*)
c      write(iu,*)' '//line
c      write(iu,*)
      close(iu)
      close(iu1)
      close(iu7)
c
      stop
      end
