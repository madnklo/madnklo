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
      integer iu,iu1,iu2,iu7,iu8
      common/cdim/ndim
      double precision int_VV
      double precision res_VV,err_VV
      external int_VV
      integer order
      logical doplot
      common/cdoplot/doplot
      double precision rescale_plot_VV
      character*100 line
      integer nitVVth,nclVVth,nitVV,nclVV
      integer nitVVth0,nclVVth0,nclVV0,nclVVth1,nclVV1
      common/iterations/nitVV
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
      double precision sum_vv,sum_err_vv
      double precision sum_err_vv_a,err_vv_a(N_MAX_CG)
c
      sum_vv=0d0
      sum_err_vv=0d0
      res_vv=0d0
      err_vv=0d0
c
      call SETPARA('param_card.dat')
      call SETRUN('run_card.dat')
c
c     read inputs
      region=0d0
      order=0
      s_had = (EBEAM(1)+EBEAM(2))**2
      NITVVTH =  NITERS_FO_GRID
      NCLVVTH = NPOINTS_FO_GRID
      NITVV = NITERS_FO
      NCLVV = NPOINTS_FO

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
      call histo_init
c      call analysis_begin(1,'central')
      open(unit=iu1,file='integration_VV.log')
      open(unit=iu7,file='failures_VV.log')
      open(unit=iu8 ,file='VV_chan.log')
      open(unit=iu ,file='results_VV.log')
      line='=================================================='
c      write(iu,*)' VV contribution '
c      write(iu,*)
c
c     quickly get integration error per channel so to modulate
c     number of points thrown per channel in the main loop
      nclVVth0=max(10000,int(nclVVth/5d0))
      nitVVth0=max(5,int(nitVVth/2d0))
      sum_err_VV_a=0d0
      do i=1,N_MAX_CG
         ich=i
         init=0
         doplot=.false.
         call vegas(region,ndim,int_VV,init,nclVVth0,nitVVth0,nprn,
     &   res_vv,err_vv,chi2a,acc,xi,it,ndo,si,swgt,schi)
         err_vv_a(ich) = err_vv
         sum_err_vv_a = sum_err_vv_a + err_vv_a(ich)
      enddo
c
c     main loop over channels
      do i=1,N_MAX_CG
         ich=i
         write(*,*)'VV warmup for channel',ich
         write(iu7,*)'Failures for VV warmup, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' VV WARMUP, CHANNEL',ich
         write(iu1,*)'============================='
         init=0
         doplot=.false.
         nclVVth1=max(1000,int(nclVVth*err_vv_a(ich)/sum_err_vv_a))
         call vegas(region,ndim,int_VV,init,nclVVth1,nitVVth,nprn,
     &   res_vv,err_vv,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu8,*)'VV warmup: channel, itns, calls = ',ich,nitVVth,nclVVth1
c
         write(*,*)'VV for channel',ich
         write(iu7,*)'Failures for VV, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' VV, CHANNEL',ich
         write(iu1,*)'============================='
         init=1
         doplot=.true.
         nclVV1=max(1000,int(nclVV*err_vv_a(ich)/sum_err_vv_a))
         call vegas(region,ndim,int_VV,init,nclVV1,nitVV,nprn,
     &   res_vv,err_vv,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu8,*)'VV: channel, itns, calls = ',ich,nitVV,nclVV1
         rescale_plot_VV=dble(nitVV)/min(dble(nitVV),dble(it))
         sum_vv = sum_vv + res_vv
         sum_err_vv = sum_err_vv + err_vv**2
         write(iu8,*)' sigma VV [pb], channel',ich,' = ',
     &   res_vv,' +-',err_vv
         write(iu8,*)
c     
         write(*,*)'...done'
      enddo
c
c     finalise histograms and output files
      sum_err_vv = dsqrt(sum_err_vv)
      call histo_final('plot_VV.dat',rescale_plot_VV)
c      write(iu,*)
c      write(iu,*)' '//line
c      write(iu,*)
      write(iu,*)' sigma VV [pb]  = ',sum_vv,' +-',sum_err_vv
c      write(iu,*)
c      write(iu,*)' '//line
c      write(iu,*)
      close(iu)
      close(iu1)
      close(iu7)
c
      stop
      end
