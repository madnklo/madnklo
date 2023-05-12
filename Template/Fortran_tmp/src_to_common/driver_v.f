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
      double precision int_virtual
      double precision res_V,err_V
      external int_virtual
      integer order
      logical doplot
      common/cdoplot/doplot
      double precision rescale_plot_V
      character*100 line
      integer nitVth,nclVth,nitV,nclV
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

      open(unit=iu1,file='integration_V.log')
      open(unit=iu7,file='failures_V.log')

      write(*,'(a)')'Warm up Virtual'
      write(iu1,'(a)')'============================='
      write(iu1,'(a)')' VIRTUAL WARMUP                 '
      write(iu1,'(a)')'============================='
      init=0
      doplot=.false.
      call vegas(region,ndim,int_virtual,init,nclVth,nitVth,nprn,res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
c
      write(*,'(a)')'Integrating Virtual'
      write(iu1,'(a)')
      write(iu1,'(a)')'============================='
      write(iu1,'(a)')' VIRTUAL                        '
      write(iu1,'(a)')'============================='
      init=1
      doplot=.true.
      call histo_init
      call vegas(region,ndim,int_virtual,init,nclV,nitV,nprn,res_v,err_v,chi2a,acc,xi,it,ndo,si,swgt,schi)
      rescale_plot_V=dble(nitV)/min(dble(nitV),dble(it))
      write(*,110)char(13),'...done     '
      write(*,*)
      call histo_final('plot_V.dat',rescale_plot_V)
c
      open(unit=iu,file='results_V.log')
      line='=================================================='
      write(iu,*)' V '
      write(iu,*)
      write(iu,*)' itns and calls for V warmup  = ',nitVth,nclVth
      write(iu,*)' itns and calls for V integration = ',nitV,nclV
      write(iu,*)
c      write(iu,*)' '//line//line
      write(iu,*)' '//line
      write(iu,*)' '//line
      write(iu,*)' sigma V [pb]  = ',res_v,' +-',err_v
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
