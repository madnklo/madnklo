      program driver_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d
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
      integer isec,jsec,ksec,lsec
      double precision s_had
      integer iu,iu1,iu7,iu8,iu9
      common/cdim/ndim
      double precision int_double_real_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d
      double precision err_rr,res_rr
      external int_double_real_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d
      common/ciunitNLO/iu8
      integer order
      logical doplot
      common/cdoplot/doplot
      double precision rescale_plot_RR
      character*100 line
      integer nitRRth,nclRRth,nitRR,nclRR
      integer nitRRth0,nclRRth0,nclRR0,nclRRth1,nclRR1
      COMMON/iterations/NITRR
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
      double precision sum_rr,sum_err_rr
      double precision sum_err_rr_a,err_rr_a(N_MAX_CG)
c
      sum_rr=0d0
      sum_err_rr=0d0
      res_rr=0d0
      err_rr=0d0
c
      call SETPARA('param_card.dat')
      call SETRUN('run_card.dat')
c
c     read inputs
      region=0d0
      order=1
      s_had = (EBEAM(1)+EBEAM(2))**2
      NITRRTH = NITERS_FO_GRID
      NCLRRTH = NPOINTS_FO_GRID
      NITRR = NITERS_FO
      NCLRR = NPOINTS_FO
c     TODO: understand muR input fixed/dyn scale
c
c     initialise physics parameters and set sector parametrisation
      iu1=44
      iu=55
      iu7=77
      iu8=88
      iu9=99
c
c     phase-space dimension, same for all contributions to this folder
      ndim=3*(nexternal-2)-4
      do i=1,2
         if(ISNNLOQCDPARTON(i)) ndim = ndim + 1
      enddo
      do i=1,ndim
         region(i)=0d0
         region(i+ndim)=1d0
      enddo
c
c     initialise histograms and open output files
      isec=%(isec)d
      jsec=%(jsec)d
      ksec=%(ksec)d
      lsec=%(lsec)d
      call histo_init
      open(unit=iu1,file='integration_RR_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d.log')
      open(unit=iu7,file='failures_RR_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d.log')
      open(unit=iu8,file='testRR_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d.log')
      open(unit=iu9,file='chan_RR_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d.log')
      open(unit=iu,file='results_RR_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d.log')
      line='=================================================='
      write(iu9,*)' Double-real contribution '
c      write(iu,*)
c
c     quickly get integration error per channel so to modulate
c     number of points thrown per channel in the main loop
      nclRRth0=max(10000,int(nclRRth/5d0))
      nitRRth0=max(5,int(nitRRth/2d0))
      sum_err_rr_a=0d0
      do i=1,N_MAX_CG
         ich=i
         init=0
         doplot=.false.
         call vegas(region,ndim,int_double_real_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d,init,nclRRth0,nitRRth0,nprn,res_rr,err_rr,chi2a,acc,xi,it,ndo,si,swgt,schi)
         err_rr_a(ich) = err_rr
         sum_err_rr_a = sum_err_rr_a + err_rr_a(ich)
      enddo
c
c     main loop over channels
      do i=1,N_MAX_CG
         ich=i
         write(*,*)'Double-real %(isec)d%(jsec)d%(ksec)d%(lsec)d warmup for channel',ich
         write(iu7,*)'Failures for RR%(isec)d%(jsec)d%(ksec)d%(lsec)d  warmup, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' DOUBLE-REAL_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d WARMUP, CHANNEL',ich
         write(iu1,*)'============================='
         init=0
         doplot=.false.
         nclRRth1=max(1000,int(nclRRth*err_rr_a(ich)/sum_err_rr_a))
         call vegas(region,ndim,int_double_real_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d,init,nclRRth1,nitRRth,nprn,res_rr,err_rr,chi2a,acc,xi,it,ndo,si,swgt,schi)
         write(iu9,*)'RR%(isec)d%(jsec)d%(ksec)d%(lsec)d warmup:channel, itns, calls = ',ich,nitRRth,nclRRth1
c
         write(*,*)'Double-real %(isec)d%(jsec)d%(ksec)d%(lsec)d for channel',ich
         write(iu7,*)'Failures for RR%(isec)d%(jsec)d%(ksec)d%(lsec)d, channel',ich
         write(iu1,*)
         write(iu1,*)'============================='
         write(iu1,*)' DOUBLE-REAL_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d, CHANNEL',ich
         write(iu1,*)'============================='
         init=1
         doplot=.true.
         nclRR1=max(1000,int(nclRR*err_rr_a(ich)/sum_err_rr_a))
         call vegas(region,ndim,int_double_real_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d,init,nclRR1,nitRR,nprn,res_rr,err_rr,chi2a,acc,xi,it,ndo,si,swgt,schi)
         rescale_plot_RR=dble(nitRR)/min(dble(nitRR),dble(it))
         sum_rr = sum_rr + res_rr
         sum_err_rr = sum_err_rr + err_rr**2
         write(iu9,*)' sigma RR%(isec)d_%(jsec)d%(ksec)d_%(lsec)d [pb], channel',ich,' = ',res_rr,' +-',err_rr
         write(iu9,*)
c     
         write(*,*)'...done'
      enddo
c
c     finalise histograms and output files
      sum_err_rr = dsqrt(sum_err_rr)      
      call histo_final('plot_RR_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d.dat', rescale_plot_RR)
c      write(iu,*)
c      write(iu,*)' '//line
c      write(iu,*)
      write(iu,*)' sigma RR%(isec)d%(jsec)d%(ksec)d%(lsec)d [pb]  = ', sum_rr,' +-',sum_err_rr
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


