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
      integer iu,iu1,iu2,iu7
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
      double precision sum_v, sum_err_v
      double precision sum_err_v_a,err_v_a(N_MAX_CG)
c
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
      s_had = (EBEAM(1)+EBEAM(2))**2
C     TODO: read from run_card.inc
      NITVTH = 10
      NCLVTH = 10000
      NITV = 10 
      NCLV = 10000
c     TODO: understand muR input fixed/dyn scale
c
c     initialise physics parameters
      iu1=44
      iu2=50
      iu=55
      iu7=77
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
      open(unit=iu1,file='integration_V.log')
      open(unit=iu2,file='test_poles_V.log')
      open(unit=iu7,file='failures_V.log')
      open(unit=iu,file='results_V.log')
      line='=================================================='
      write(iu,*)' Virtual contribution '
      write(iu,*)
c
c     quickly get integration error per channel so to modulate
c     number of points thrown per channel in the main loop
      nclVth0=max(10000,int(nclVth/5d0))
      nitVth0=max(5,int(nitVth/2d0))
      sum_err_v_a=0d0
      do i=1,N_MAX_CG
         ich=i
         init=0
         doplot=.false.
         call vegas(region,ndim,int_virtual,init,nclVth0,nitVth0,nprn,
     &   res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
         err_v_a(ich) = err_v
         sum_err_v_a = sum_err_v_a + err_v_a(ich)
      enddo
c
c     main loop over channels
      do i=1,N_MAX_CG
         ich=i
         write(*,*)'Virtual warmup for channel',ich
         write(iu7,*)'Failures for virtual warmup, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' VIRTUAL WARMUP, CHANNEL',ich
         write(iu1,*)'============================='
         init=0
         doplot=.false.
         nclVth1=max(1000,int(nclVth*err_v_a(ich)/sum_err_v_a))
         call vegas(region,ndim,int_virtual,init,nclVth1,nitVth,nprn,
     &   res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu,*)'V warmup: channel, itns, calls = ',ich,nitVth,nclVth1
c
         write(*,*)'Virtual for channel',ich
         write(iu7,*)'Failures for virtual, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' VIRTUAL, CHANNEL',ich
         write(iu1,*)'============================='
         init=1
         doplot=.true.
         nclV1=max(1000,int(nclV*err_v_a(ich)/sum_err_v_a))
         call vegas(region,ndim,int_virtual,init,nclV1,nitV,nprn,
     &   res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu,*)'V: channel, itns, calls = ',ich,nitV,nclV1
         rescale_plot_V=dble(nitV)/min(dble(nitV),dble(it))
         sum_v = sum_v + res_v
         sum_err_v = sum_err_v + err_v**2
         write(iu,*)' sigma V [pb], channel',ich,' = ',
     &   res_v,' +-',err_v
         write(iu,*)
c     
         write(*,*)'...done'
      enddo
c
c     finalise histograms and output files
      sum_err_v = dsqrt(sum_err_v)      
      call histo_final('plot_V.dat',rescale_plot_V)
      write(iu,*)
      write(iu,*)' '//line
      write(iu,*)
      write(iu,*)' sigma V [pb]  = ',sum_v,' +-',sum_err_v
      write(iu,*)
      write(iu,*)' '//line
      write(iu,*)
      close(iu)
      close(iu1)
      close(iu2)
      close(iu7)
c
      stop
      end
