      program driver_%(isec)d_%(jsec)d
      implicit none
      include 'nexternal.inc'
      include 'colored_partons.inc'
      integer mxdim
      parameter(mxdim=30)
      integer ndim,i,j,idum
      integer isec,jsec
      double precision s_had
      double precision err_r,errnlo,res_r,totnlo
      integer iu,iu1,iu4,iu6,iu7,iu8
      common/ciunitNLO/iu8
      double precision int_real_%(isec)d_%(jsec)d
      external int_real_%(isec)d_%(jsec)d
      integer order
      logical doplot
      common/cdoplot/doplot
      double precision rescale_plot_R
      character*100 line
      integer nitRth,nclRth,nitR,nclR
      double precision ebeam(2)
c
c     vegas declarations
      integer ndmx,nprn,ndo,init,it
      parameter (ndmx=50,nprn=0)
      double precision chi2a,acc,si,swgt,schi
      double precision region(2*mxdim),xi(ndmx,mxdim)
      parameter(acc=1d-10)
      common/rand/idum
c
      call setpara('param_card.dat')
c     TODO: include run.inc, setrun.f in Source directory, write run_card.inc
c      call setpara('run_card.dat')
c
c     read inputs
      region=0d0
      order=1
      s_had = (EBEAM(1)+EBEAM(2))**2
c     TODO: understand muR input fixed/dyn scale
c      muR
c     TODO: include nitR in run_card as in MG5
c      nitRth,nclRth
c      nitR,nclR

      isec=%(isec)d
      jsec=%(jsec)d
c
c     initialise physics parameters and set sector parametrisation
      iu1=44
      iu4=54
      iu=55
      iu6=56
      iu7=77
      iu8=88
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
      open(unit=iu1,file='integration_R_%(isec)d_%(jsec)d.log')
      open(unit=iu7,file='failures_R_%(isec)d_%(jsec)d.log')
      open(unit=iu8,file='test_R_%(isec)d_%(jsec)d.log')
c
      write(*,'(a)')'Warm up Rsub_%(isec)d_%(jsec)d'
      write(iu1,'(a)')'============================='
      write(iu1,'(a)')' REAL_%(isec)d_%(jsec)d WARMUP '
      write(iu1,'(a)')'============================='
      init=0
      doplot=.false.
      call vegas(region,ndim,int_real_%(isec)d_%(jsec)d,init,nclRth,nitRth,nprn,res_r,err_r,chi2a,acc,xi,it,ndo,si,swgt,schi)
c
      write(*,'(a)')'Integrating Rsub_%(isec)d_%(jsec)d'
      write(iu1,'(a)')
      write(iu1,'(a)')'============================='
      write(iu1,'(a)')' REAL_%(isec)d_%(jsec)d '
      write(iu1,'(a)')'============================='
      init=1
c      doplot=.true.
c      call histo_init
      call vegas(region,ndim,int_real_%(isec)d_%(jsec)d,init,nclR,nitR,nprn,res_r,err_r,chi2a,acc,xi,it,ndo,si,swgt,schi)
c      rescale_plot_R=dble(nitR)/min(dble(nitR),dble(it))
      write(*,110)char(13),'...done     '
      write(*,*)
c      call histo_final('plot_R_%(isec)d_%(jsec)d'.dat',rescale_plot_R)
c
      totNLO=res_r
      errNLO=err_r
      open(unit=iu,file='results_R_%(isec)d_%(jsec)d.log')
c
      line='=================================================='
      write(iu,*)' Rsub_%(isec)d_%(jsec)d '
      write(iu,*)
      write(iu,*)' itns and calls for Rsub_%(isec)d_%(jsec)d warmup  = ',nitRth,nclRth
      write(iu,*)' itns and calls for Rsub_%(isec)d_%(jsec)d integration = ',nitR,nclR
      write(iu,*)
      write(iu,*)' '//line//line
      write(iu,*)' sigma Rsub_%(isec)d_%(jsec)d [pb]  = ',res_r,' +-',err_r
      write(iu,*)' '//line//line
      write(iu,*)
      close(iu)
      close(iu1)
      close(iu7)
      close(iu8)
c
c
 110  format(a1,a12,$)
c
      stop
      end

