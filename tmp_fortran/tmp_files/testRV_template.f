      subroutine test_RV_%(isec)d_%(jsec)d(ievnt)
      implicit none
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
c      integer isec,jsec,ksec,lsec
c      common/csecindices/isec,jsec,ksec,lsec
c      integer i,iU,iS,iB,iA,iref
c      common/cNLOmaplabels/iU,iS,iB,iA,iref
c      integer iUtmp,iStmp
c      integer iunit,ievnt
      integer ievnt
c      INTEGER, PARAMETER :: MXDIM = 30
c      double precision x0(mxdim)
      character*10 dash10
c      save ievnt
c      double precision xsave(3)
c      common/cxsave/xsave
      double precision e(2), l(2)
c
      dash10='----------'
c      ievnt=ievnt+1
c
      write(88,*)dash10//dash10//dash10//dash10
      write(88,*)dash10//dash10//dash10//dash10
      write(88,*)' EVENT NUMBER ',ievnt
      write(88,*)dash10//dash10//dash10//dash10
      write(88,*)dash10//dash10//dash10//dash10
%(limit_str)s
c
c     reinstate original xsave after testing
c      do i=1,3
c         xsave(i)=x0(i)
c      enddo
c
      write(88,*)
      write(88,*)
      write(88,*)
c
      return
      end


c      subroutine do_limit_RV_%(isec)d_%(jsec)d(iunit,limstr,x0,e,l)
      subroutine do_limit_RV_%(isec)d_%(jsec)d(limstr,e,l)
      use sectors2_module
      implicit none
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'leg_PDGs.inc'
      INCLUDE 'ngraphs.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'nsquaredSO.inc'
      integer iitn,i,j,maxitn,iunit,ierr
      integer isec,jsec,ksec,lsec
      common/csecindices/isec,jsec,ksec,lsec
      integer iU,iS,iB,iA,iref
      common/cnlomaplabels/iU,iS,iB,iA,iref
      integer, parameter :: mxdim=30
      parameter(maxitn=12)
c      double precision x0(mxdim),x(mxdim)
      double precision x0(3*nexternal-10),xr(3*nexternal-10)
      double precision sNLO(nexternal,nexternal)
      double precision sLO(nexternal-1,nexternal-1)
      double precision KRVNNLO(-2:0)
      double precision lam,lim,RVNNLO(-2:0),single_real
      character*5 str5
      character*8 limstr
      character*10 str10
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision xjac,xjacB
      double precision xsave(3)
      DOUBLE PRECISION ANS(0:1) !TODO SET CORRECTLY RANGE OF ANS
      DOUBLE PRECISION ALPHAS, ALPHA_QCD
      DOUBLE PRECISION WGT,WGTPL,wgt_chan
      DOUBLE PRECISION SCM
      INTEGER, PARAMETER :: HEL=-1
      integer %(NNLO_RV_proc_str)sfl_factor
      common/%(NNLO_RV_proc_str)sflavour_factor/%(NNLO_RV_proc_str)sfl_factor
      common/cxsave/xsave
      double precision e(2),l(2)
      logical inittest
      data inittest/.true./
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      double precision ran2
      double precision CSpow(2)
      common /cCSpow/CSpow
      integer matelem_array_dim
      real*8 , allocatable :: matelem(:,:)
      integer returncode
      integer nsquaredso_loop
      real*8 , allocatable :: prec_found(:)

      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
      SCM = (2D0*EBEAM(1))**2
c
c     initialise
      str5 ='     '
      str10='          '
      xjac=0d0
      sNLO=0d0
      sLO=0d0
      wgt_chan=1d0
      do i=1,3*nexternal-10
         x0(i)=ran2(33+10*i)
      enddo
      xr=x0

c
c
c     start testing
      write(88,*)
      write(88,*)
      write(88,*)'LIM = '//trim(limstr)
      write(88,*)str10//'lambda'//str10//str10//'R'//str10//str10//str5//'LIM'//str10//str10//'|R-LIM|/|LIM|'
c
c     loop to get closer and closer to the limit
      do iitn=1,maxitn
         lam=10d0**(1d0-iitn)
c
c     initialise
         KRVNNLO=0d0
c
c     rescale relevant x random numbers
c     x(1) is zCS, while x(2) is yCS
c     TODO: this rescaling is specific for (ijr) mapping; generalise
         xr(1:2)=abs(l(1:2)-x0(1:2)*lam**(e(1:2)/CSpow(1:2)))
                  
c
c
c     recompute momenta after rescaling
         call phase_space_npo(xr,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)
         if(xjac.eq.0d0.or.xjacb.eq.0d0)cycle
         call invariants_from_p(p,nexternal,sNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(pb,nexternal-1,sLO,ierr)
         if(ierr.eq.1)cycle
c
c     real virtual
      if (inittest) then
         inittest=.false.
         call %(long_proc_prefix)sget_answer_dimension(matelem_array_dim)
         allocate(matelem(0:3,0:matelem_array_dim))
         call %(long_proc_prefix)sget_nsqso_loop(nsquaredso_loop)
         allocate(prec_found(0:nsquaredso_loop))
      endif
c
      call %(long_proc_prefix)ssloopmatrix_thres(p,matelem,-1.0d0,prec_found,returncode)


      do j=-2,0
         RVNNLO(j) = MATELEM(1-j,0) 
      enddo
       RVNNLO = RVNNLO* %(NNLO_RV_proc_str)sfl_factor

      if(abs(RVNNLO(0)).ge.huge(1d0).or.isnan(RVNNLO(0)))cycle
c
      call get_sig2(snlo,nexternal)
      call get_W_NLO(isec,jsec)
c
c     counterterm
      call local_RV_counter_NNLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,xjac,xjacB,xr,KRVNNLO,wgt_chan,ierr)
         if(ierr.eq.1)cycle

         lim=KRVNNLO(0)
         single_real=RVNNLO(0)*W_NLO*xjac

         if(abs(lim).gt.0d0)then
            write(88,*)lam,single_real,lim,abs(single_real-lim)/abs(lim)
         else
            write(88,*)lam,single_real,lim,single_real,' *** '
         endif
      enddo
c
      return
      end
