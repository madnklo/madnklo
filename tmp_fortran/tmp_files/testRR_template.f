      subroutine test_RR_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(iunit,x0)
      implicit none
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
      integer isec,jsec,ksec,lsec,iref
      common/cNNLOsecindices/isec,jsec,ksec,lsec
      common/cNNLOmaplabels/iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref
      integer iunit,ievnt
      INTEGER, PARAMETER :: MXDIM = 30
      double precision x0(mxdim)
      character*10 dash10
      save ievnt
      double precision xsave(3)
      common/cxsave/xsave
      double precision e(5),l(5)
c
      dash10='----------'
      ievnt=ievnt+1
c
      write(iunit,*)dash10//dash10//dash10//dash10
      write(iunit,*)dash10//dash10//dash10//dash10
      write(iunit,*)' EVENT NUMBER ',ievnt
      write(iunit,*)dash10//dash10//dash10//dash10
      write(iunit,*)dash10//dash10//dash10//dash10
%(limit_str)s
c     
c     reinstate original xsave after testing
      do i=1,5
         xsave(i)=x0(i)
      enddo
c
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)
c
      return
      end


      subroutine do_limit_RR_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(iunit,limstr,x0,e,l)
      implicit none
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer iitn,i,j,maxitn,iunit,ierr
      integer isec,jsec,ksec,lsec
      common/cnlosecindices/isec,jsec,ksec,lsec
      integer iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref
      common/cNNLOmaplabels/iS1,iB1,iA1,iU2,iS2,iB2,iA2,iref
C      common/cnlomaplabels/iU,iS,iB,iA,iref
      integer, parameter :: mxdim=30
      parameter(maxitn=12)
      double precision x0(mxdim),x(mxdim)
      double precision sNNLO(nexternal,nexternal)
      double precision sNLO(nexternal-1,nexternal-1)
      double precision sLO(nexternal-2,nexternal-2)
      double precision KNNLO
      double precision lam,lim,RNNLO,double_real
      character*5 str5
      character*8 limstr
      character*10 str10
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision ptilde(0:3,nexternal-2)
      double precision xjac,xjacB
      double precision xsave(5)
      DOUBLE PRECISION ANS(0:1) !TODO SET CORRECTLY RANGE OF ANS
      DOUBLE PRECISION ALPHAS, ALPHA_QCD
      DOUBLE PRECISION Z_NNLO
      DOUBLE PRECISION WGT,WGTPL,wgt_chan
      DOUBLE PRECISION SCM
      INTEGER, PARAMETER :: HEL=-1
      integer %(NNLO_proc_str)sfl_factor 
      common/%(NNLO_proc_str)sflavour_factor/%(NNLO_proc_str)sfl_factor
      DOUBLE PRECISION ALPHAZ
      PARAMETER(ALPHAZ=1D0)
      common/cxsave/xsave
      double precision e(5),l(5)
      integer i 
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
      SCM = (2D0*EBEAM(1))**2
c     
c     initialise
      x=x0
      str5 ='     '
      str10='          '
      xjac=0d0
      sNNLO=0d0
      sNLO=0d0
      sLO=0d0
      Z_NNLO=0d0
      wgt_chan=1d0
c
c     TODO: MAP SOFT LIMIT AS (ilm), I.E. ONE MAPPING PER DIPOLE
c
c     start testing
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)'LIM = '//trim(limstr)
      write(iunit,*)str10//'lambda'//str10//str10//'RR'//str10//str10//str5//'LIM'//str10//str10//'|RR-LIM|/|LIM|'
c
c     possibility to set by hand the starting point
c     for the limiting procedure
c      x0(1)=0.5d0
c      x0(2)=0.5d0
c
c     loop to get closer and closer to the limit
      do iitn=1,maxitn
         lam=10d0**(1d0-iitn)
c
c     initialise
         KNNLO=0d0
c
c     rescale relevant x random numbers
c     x(1) is zCS, while x(2) is yCS
c     TODO: this rescaling is specific for (ijr) mapping; generalise
         x(1)=abs(l(1)-x0(1))*lam**e(1)
         x(2)=abs(l(2)-x0(2))*lam**e(2)
         x(4)=abs(l(4)-x0(4))*lam**e(4)
         x(5)=abs(l(5)-x0(5))*lam**e(5)

         
c
c     set xsave so that the counterterms will be called with
c     more and more singular kinematics
         do i=1,5
            xsave(i)=x(i)
         enddo
c
c     recompute momenta after rescaling
         call phase_space_npt(x,sCM,iU1,iS1,iB1,iA1,iA2,p,pb,ptilde,xjac,xjacB,iU2,iS2,iB2)
         if(xjac.eq.0d0.or.xjacB.eq.0d0) cycle
         call invariants_from_p(p,nexternal,sNNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(pb,nexternal-1,sNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(ptilde,nexternal-2,sLO,ierr)
         if(ierr.eq.1)cycle
c
c     double real
         call %(NNLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
         RNNLO = ANS(0) * %(NNLO_proc_str)sfl_factor
         if(RNNLO.lt.0d0.or.abs(RNNLO).ge.huge(1d0).or.isnan(RNNLO))cycle
         call  get_Z_NNLO(sNNLO,sCM,alphaZ,isec,jsec,ksec,lsec,Z_NNLO,ierr)
         if(ierr.eq.1)cycle
c
c     counterterm
         call local_counter_NNLO_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(sNNLO,p,sNLO,pb,sLO,ptilde,wgt,xjac,xjacB,x,KNNLO,wgt_chan,ierr)
         if(ierr.eq.1)cycle
         
         lim=KNNLO
         double_real=RNNLO*Z_NNLO*xjac
         
         if(abs(lim).gt.0d0)then
            write(iunit,*)lam,double_real,lim,abs(double_real-lim)/abs(lim)
         else
            write(iunit,*)lam,double_real,lim,double_real,' *** '
         endif
      enddo
      x=x0
c
      return
      end
