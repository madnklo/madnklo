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
      integer ndim,i,j,idum
      integer isec,jsec
      double precision s_had
      integer iu,iu1,iu7,iu8,iu9,iu0
      common/cdim/ndim
      double precision int_real_virtual_%(isec)d_%(jsec)d
      double precision err_rv,res_rv
      external int_real_virtual_%(isec)d_%(jsec)d
      common/ciunitNLO/iu8
      integer order
      logical doplot
      common/cdoplot/doplot
      double precision rescale_plot_RV
      character*100 line
      integer nitRVth,nclRVth,nitRV,nclRV
      integer nitRVth0,nclRVth0,nclRV0,nclRVth1,nclRV1
      COMMON/iterations/NITRV
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
      double precision sum_rv,sum_err_rv
c      double precision sum_err_rv_a,err_rv_a(N_MAX_CG)
      integer nwgt
      character*20 weights_info(1)
c
      sum_rv=0d0
      sum_err_rv=0d0
      res_rv=0d0
      err_rv=0d0
c
      call SETPARA('param_card.dat')
      call SETRUN('run_card.dat')
c
c     read inputs
      region=0d0
      order=1
      s_had = (EBEAM(1)+EBEAM(2))**2
      NITRVTH = NITERS_FO_GRID
      NCLRVTH = NPOINTS_FO_GRID
      NITRV = NITERS_FO
      NCLRV = NPOINTS_FO
c     TODO: understand muR input fixed/dyn scale
c
c     initialise physics parameters and set sector parametrisation
      iu1=44
      iu=55
      iu7=77
      iu8=88
      iu9=99
      iu0=50
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
c     initialise histograms and open output files
      isec=%(isec)d
      jsec=%(jsec)d
c     call histo_init
      nwgt=1
      weights_info(1)='central'
      call analysis_begin(nwgt,weights_info)
      open(unit=iu1,file='integration_RV_%(isec)d_%(jsec)d.log')
      open(unit=iu7,file='failures_RV_%(isec)d_%(jsec)d.log')
      open(unit=iu8,file='testRV_%(isec)d_%(jsec)d.log')
      open(unit=iu9,file='chan_RV_%(isec)d_%(jsec)d.log')
      open(unit=iu,file='results_RV_%(isec)d_%(jsec)d.log')
      open(unit=iu0,file='poles_RV_%(isec)d_%(jsec)d.log')
      line='=================================================='
      write(iu9,*)' Real-Virtual contribution '
c
c
c     quickly get integration error per channel so to modulate
c     number of points thrown per channel in the main loop
c      nclRVth0=max(10000,int(nclRVth/5d0))
c      nitRVth0=max(5,int(nitRVth/2d0))
c      sum_err_rv_a=0d0
c      do i=1,N_MAX_CG
c         ich=i
c         init=0
c         doplot=.false.
c         call vegas(region,ndim,int_real_virtual_%(isec)d_%(jsec)d,init,nclRVth0,nitRVth0,nprn,res_rv,err_rv,chi2a,acc,xi,it,ndo,si,swgt,schi)
c         err_rv_a(ich) = err_rv
c         sum_err_rv_a = sum_err_rv_a + err_rv_a(ich)
c      enddo
c
c     main loop over channels
c     do i=1,N_MAX_CG
      do i=1,1
         ich=i
         write(*,*)'Real Virtual %(isec)d%(jsec)d warmup for channel',ich
         write(iu7,*)'Failures for RV%(isec)d%(jsec)d warmup, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' REAL_VIRTUAL_%(isec)d_%(jsec)d WARMUP, CHANNEL',ich
         write(iu1,*)'============================='
         init=0
         doplot=.false.
c        nclRVth1=max(1000,int(nclRVth*err_rv_a(ich)/sum_err_rv_a))
c         call vegas(region,ndim,int_real_virtual_%(isec)d_%(jsec)d,init,nclRVth1,nitRVth,nprn,res_rv,err_rv,chi2a,acc,xi,it,ndo,si,swgt,schi)
         call vegas(region,ndim,int_real_virtual_%(isec)d_%(jsec)d,init,nclRVth,nitRVth,nprn,res_rv,err_rv,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu9,*)'RV%(isec)d%(jsec)d warmup: channel, itns, calls = ',ich,nitRVth,nclRVth1
c
         write(*,*)'Real Virtual %(isec)d%(jsec)d for channel',ich
         write(iu7,*)'Failures for RV%(isec)d%(jsec)d, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' REAL_VIRTUAL_%(isec)d_%(jsec)d, CHANNEL',ich
         write(iu1,*)'============================='
         init=1
         doplot=.true.
c         nclRV1=max(1000,int(nclRV*err_rv_a(ich)/sum_err_rv_a))
c         call vegas(region,ndim,int_real_virtual_%(isec)d_%(jsec)d,init,nclRV1,nitRV,nprn,res_rv,err_rv,chi2a,acc,xi,it,ndo,si,swgt,schi)
         call vegas(region,ndim,int_real_virtual_%(isec)d_%(jsec)d,init,nclRV,nitRV,nprn,res_rv,err_rv,chi2a,acc,xi,it,ndo,si,swgt,schi)
         rescale_plot_RV=dble(nitRV)/min(dble(nitRV),dble(it))
         sum_rv = sum_rv + res_rv
         sum_err_rv = sum_err_rv + err_rv**2
         write(iu9,*)'RV%(isec)d%(jsec)d: channel, itns, calls = ',ich,nitRV,nclRV1
         write(iu9,*)' sigma RV%(isec)d_%(jsec)d [pb], channel',ich,' = ',res_rv,' +-',err_rv
         write(iu9,*)
c
         write(*,*)'...done'
      enddo
c
c     finalise histograms and output files
      call analysis_end(1d0)
      sum_err_rv = dsqrt(sum_err_rv)
c      call histo_final('plot_RV_%(isec)d_%(jsec)d.dat',rescale_plot_RV)
c      write(iu,*)
c      write(iu,*)' '//line
c      write(iu,*)
      write(iu,*)' sigma RV%(isec)d%(jsec)d [pb]  = ',sum_rv,' +-',sum_err_rv
c      write(iu,*)
c      write(iu,*)' '//line
c      write(iu,*)
      close(iu)
      close(iu1)
      close(iu7)
      close(iu8)
      close(iu9)
c
      stop
      end
