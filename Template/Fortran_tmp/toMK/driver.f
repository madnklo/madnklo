      program driver
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      include 'output.inc'
      integer ndim,i,j,idum
      integer isec,jsec,ksec,lsec,iUU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref
      common/cNNLOsecindices/isec,jsec,ksec,lsec
      common/cNNLOmaplabels/iUU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref
      integer iu0,iu1,iu4,iu5,iu6,iu7,iu8
      common/cdim/ndim
      common/ciunitNNLO/iu8
      double precision int_double_real
      external int_double_real
      character*80 tmpstr,runname,infilename
      character*10 order
      character*1 strisec,strjsec,strksec,strlsec
      common/cNNLOstrsecindices/strisec,strjsec,strksec,strlsec
      logical doplot
      common/cdoplot/doplot
      double precision rescale_plot_RR
c
c     vegas declarations
      integer ndmx,nprn,ndo,init,it
      parameter (ndmx=50,nprn=0)
      double precision chi2a,acc,si,swgt,schi
      double precision region(2*mxdim),xi(ndmx,mxdim)
      parameter(acc=1d-10)
      common/rand/idum
c
c     read inputs
      iu0=33
      region=0d0
      read(*,*)infilename
      open(unit=iu0,file=infilename)
      read(iu0,*)order,tmpstr
      read(iu0,*)sCM,tmpstr
      read(iu0,*)muR,tmpstr
      read(iu0,*)nitBth,nclBth,tmpstr
      read(iu0,*)nitB,nclB,tmpstr
      read(iu0,*)nitVth,nclVth,tmpstr
      read(iu0,*)nitV,nclV,tmpstr
      read(iu0,*)nitVVth,nclVVth,tmpstr
      read(iu0,*)nitVV,nclVV,tmpstr
      read(iu0,*)nitRth,nclRth,tmpstr
      read(iu0,*)nitR,nclR,tmpstr
      read(iu0,*)nitRVth,nclRVth,tmpstr
      read(iu0,*)nitRV,nclRV,tmpstr
      read(iu0,*)nitRRth,nclRRth,tmpstr
      read(iu0,*)nitRR,nclRR,tmpstr
      read(iu0,*)idum,tmpstr
      read(iu0,*)runname,tmpstr
      read(iu0,*)isec,jsec,ksec,lsec
      close(iu0)
c
c     input checks
      if(order.ne.'NLO'.and.order.ne.'NNLO'.and.
     &   order.ne.'NLOonly'.and.order.ne.'NNLOonly')then
         write(*,*)'Wrong order ',order
         stop
      endif
      if(sCM.le.0d0)then
         write(*,*)'Wrong input sCM ',sCM
         stop
      endif
      if(muR.le.0d0)then
         write(*,*)'Wrong input muR ',muR
         stop
      endif
c
c     initialise physics parameters and set sector parametrisation
      call setup
      call assign_phsp_labels_npt(isec,jsec,ksec,lsec,
     &iUU1,iS1,iB1,iA1,iA2,iref)
c
      iu1=44
      iu4=54
      iu5=55
      iu6=56
      iu7=77
      iu8=88
c
c     phase-space dimension, same for all contributions to this folder
      ndim=ndimLO+6
      do i=1,ndim
         region(i)=0d0
         region(i+ndim)=1d0
      enddo
c
c     sector indices, same for all contributions to this folder
      write(strisec,'(i1)')isec
      write(strjsec,'(i1)')jsec
      write(strksec,'(i1)')ksec
      write(strlsec,'(i1)')lsec
c
c     NNLO contribution
      if(order.eq.'NNLO'.or.order.eq.'NNLOonly')then
         open(unit=iu1,file='integration_RR'//strisec//strjsec//strksec//strlsec//'.log')
         open(unit=iu7,file='failures_RR'//strisec//strjsec//strksec//strlsec//'.log')
c
c     test if there is a counterterm
         if(is_NNLO_singular_sec(isec,jsec,ksec,lsec))then
            open(unit=iu8,file='test_RR'//strisec//strjsec//strksec//strlsec//'.log')
         endif
c
         write(*,'(a)')'Integrating RR '//strisec//' '//strjsec//' '//strksec//' '//strlsec
         write(iu1,'(a)')'============================='
         write(iu1,'(a)')' DOUBLE REAL '//strisec//' '//strjsec//' '//strksec//' '//strlsec//' WARMUP'
         write(iu1,'(a)')'============================='
         init=0
         doplot=.false.
         call vegas(region,ndim,int_double_real,init,nclRRth,nitRRth,
     &        nprn,res_rr,err_rr,chi2a,acc,xi,it,ndo,si,swgt,schi)
c
         write(iu1,'(a)')
         write(iu1,'(a)')'============================='
         write(iu1,'(a)')' DOUBLE REAL '//strisec//' '//strjsec//' '//strksec//' '//strlsec
         write(iu1,'(a)')'============================='
         init=1
         doplot=.true.
         call histo_init
         call vegas(region,ndim,int_double_real,init,nclRR,nitRR,
     &        nprn,res_rr,err_rr,chi2a,acc,xi,it,ndo,si,swgt,schi)
         rescale_plot_RR=dble(nitRR)/min(dble(nitRR),dble(it))
         write(*,110)char(13),'...done     '
         write(*,*)
         call histo_final('plot_RR'//strisec//strjsec//strksec//strlsec//'.dat',rescale_plot_RR)
c
         totNNLO=res_rr
         errNNLO=err_rr
         open(unit=iu6,
     &   file='results_RR'//strisec//strjsec//strksec//strlsec//'.log')
         call outputNNLO(iu6)
         close(iu6)
         close(iu1)
         close(iu7)
         if(is_NNLO_singular_sec(isec,jsec,ksec,lsec))close(iu8)
      endif
c
 110  format(a1,a12,$)
c
c$$$c     combine plots
c$$$      call system('gfortran -o combine_plots combine_plots.f')
c$$$      call system('cp combine_plots.in '//resdir)
c$$$      call system('mv combine_plots '//resdir)
c
      stop
      end
