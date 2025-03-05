      program driver_%(isec)d_%(jsec)d
      use init_R_module
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
      common/cdim/ndim
      double precision int_real_%(isec)d_%(jsec)d
      double precision err_r,res_r
      external int_real_%(isec)d_%(jsec)d
      integer order
      double precision rescale_plot_R
      character*100 line
      integer nitRth,nclRth,nclR
      integer nitRth0,nclRth0,nclR0,nclRth1,nclR1
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
      double precision sum_r,sum_err_r
      double precision sum_err_r_a,err_r_a(N_MAX_CG)
c
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
      NITRTH = NITERS_FO_GRID
      NCLRTH = NPOINTS_FO_GRID
      NITR = NITERS_FO
      NCLR = NPOINTS_FO
c     TODO: understand muR input fixed/dyn scale
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
      call histo_init
      open(unit=44,file='integration_R_%(isec)d_%(jsec)d.log')
      open(unit=77,file='failures_R_%(isec)d_%(jsec)d.log')
      open(unit=88,file='testR_%(isec)d_%(jsec)d.log')
      open(unit=99,file='chan_R_%(isec)d_%(jsec)d.log')
      open(unit=55,file='results_R_%(isec)d_%(jsec)d.log')
      line='=================================================='
      write(99,*)' Real contribution '
c
c     quickly get integration error per channel so to modulate
c     number of points thrown per channel in the main loop
      nclRth0=max(10000,int(nclRth/5d0))
      nitRth0=max(5,int(nitRth/1d0))
      sum_err_r_a=0d0
      do i=1,N_MAX_CG
         ich=i
         init=0
         doplot=.false.
         call vegas(region,ndim,int_real_%(isec)d_%(jsec)d,init,nclRth0,nitRth0,nprn,res_r,err_r,chi2a,acc,xi,it,ndo,si,swgt,schi)
         err_r_a(ich) = err_r**2d0
         sum_err_r_a = sum_err_r_a + err_r_a(ich)
      enddo
c
c     main loop over channels
      do i=1,N_MAX_CG
         ich=i
         write(*,*)'Real %(isec)d%(jsec)d warmup for channel',ich
         write(77,*)'Failures for R%(isec)d%(jsec)d warmup, channel',ich
         write(44,*)
         write(44,*)'============================='
         write(44,*)' REAL_%(isec)d_%(jsec)d WARMUP, CHANNEL',ich
         write(44,*)'============================='
         init=0
         doplot=.false.
         nclRth1=max(2000,int(nclRth*err_r_a(ich)/sum_err_r_a))
         call vegas(region,ndim,int_real_%(isec)d_%(jsec)d,init,nclRth1,nitRth,nprn,res_r,err_r,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(99,*)'R%(isec)d%(jsec)d warmup: channel, itns, calls = ',ich,nitRth,nclRth1
c
         write(*,*)'Real %(isec)d%(jsec)d for channel',ich
         write(77,*)'Failures for R%(isec)d%(jsec)d, channel',ich
         write(44,*)
         write(44,*)'============================='
         write(44,*)' REAL_%(isec)d_%(jsec)d, CHANNEL',ich
         write(44,*)'============================='
         init=1
         doplot=.true.
         nclR1=max(4000,int(nclR*err_r_a(ich)/sum_err_r_a))
         call vegas(region,ndim,int_real_%(isec)d_%(jsec)d,init,nclR1,nitR,nprn,res_r,err_r,chi2a,acc,xi,it,ndo,si,swgt,schi)
         rescale_plot_R=dble(nitR)/min(dble(nitR),dble(it))
         sum_r = sum_r + res_r
         sum_err_r = sum_err_r + err_r**2
         write(99,*)' sigma R%(isec)d_%(jsec)d [pb], channel',ich,' = ',res_r,' +-',err_r
         write(99,*)
c     
         write(*,*)'...done'
      enddo
c
c     finalise histograms and output files
      sum_err_r = dsqrt(sum_err_r)      
      call histo_final('plot_R_%(isec)d_%(jsec)d.dat',rescale_plot_R)
      write(55,*)' sigma R%(isec)d%(jsec)d [pb]  = ',sum_r,' +-',sum_err_r
      close(55)
      close(44)
      close(77)
      close(88)
      close(99)
c
      stop
      end


