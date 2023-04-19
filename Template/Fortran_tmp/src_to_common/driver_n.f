      program driver_n_body
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'colored_partons.inc'
      integer mxdim
      parameter(mxdim=30)
      integer ndim,i,j,idum
      double precision s_had
      integer iu,iu1,iu7
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
c
c     vegas declarations
      integer ndmx,nprn,ndo,init,it
      parameter (ndmx=50,nprn=0)
      double precision chi2a,acc,si,swgt,schi
      double precision region(2*mxdim),xi(ndmx,mxdim)
      parameter(acc=1d-10)
      common/rand/idum
c
      call SETPARA('param_card.dat')
      call SETRUN('run_card.dat')
c
c     read inputs
      region=0d0
      order=0
      s_had = (EBEAM(1)+EBEAM(2))**2
C     TODO: read from run_card.inc
      NITBTH = 10
      NCLBTH = 10000
      NITB = 10 
      NCLB = 10000
c     TODO: understand muR input fixed/dyn scale
c
c     initialise physics parameters
      iu1=44
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

      open(unit=iu1,file='integration_B.log')
      open(unit=iu7,file='failures_B.log')

      write(*,'(a)')'Warm up B'
      write(iu1,'(a)')'============================='
      write(iu1,'(a)')' BORN WARMUP                 '
      write(iu1,'(a)')'============================='
      init=0
      doplot=.false.
      call vegas(region,ndim,int_Born,init,nclBth,nitBth,nprn,res_b,err_b,chi2a,acc,xi,it,ndo,si,swgt,schi)
c
      write(iu1,'(a)')
      write(iu1,'(a)')'============================='
      write(iu1,'(a)')' BORN                        '
      write(iu1,'(a)')'============================='
      init=1
c      doplot=.true.
c      call histo_init
      call vegas(region,ndim,int_Born,init,nclB,nitB,nprn,res_b,err_b,chi2a,acc,xi,it,ndo,si,swgt,schi)
c      rescale_plot_B=dble(nitB)/min(dble(nitB),dble(it))
      write(*,110)char(13),'...done     '
      write(*,*)
c      call histo_final('plot_B.dat',rescale_plot_B)
c
      open(unit=iu,file='results_B.log')
      line='=================================================='
      write(iu,*)' B '
      write(iu,*)
      write(iu,*)' itns and calls for B warmup  = ',nitBth,nclBth
      write(iu,*)' itns and calls for B integration = ',nitB,nclB
      write(iu,*)
c      write(iu,*)' '//line//line
      write(iu,*)' '//line
      write(iu,*)' '//line
      write(iu,*)' sigma B [pb]  = ',res_b,' +-',err_b
      write(iu,*)' '//line
      write(iu,*)' '//line
c      write(iu,*)' '//line//line
      write(iu,*)
      close(iu)
      close(iu1)
      close(iu7)
c
c
 110  format(a1,a12,$)
c
      stop
      end
