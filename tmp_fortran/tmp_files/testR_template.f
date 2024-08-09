      subroutine test_R_%(isec)d_%(jsec)d(iunit,x0)
      implicit none
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer isec,jsec,ksec,lsec
      common/csecindices/isec,jsec,ksec,lsec
      integer i,iU,iS,iB,iA,iref
      common/cNLOmaplabels/iU,iS,iB,iA,iref
      integer iUtmp,iStmp
      integer iunit,ievnt
      INTEGER, PARAMETER :: MXDIM = 30
      double precision x0(mxdim)
      character*10 dash10
      save ievnt
      double precision xsave(3)
      common/cxsave/xsave
      double precision e(2), l(2)
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
      do i=1,3
         xsave(i)=x0(i)
      enddo
c
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)
c
      return
      end


      subroutine do_limit_R_%(isec)d_%(jsec)d(iunit,limstr,x0,e,l)
      implicit none
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer iitn,i,j,maxitn,iunit,ierr
      integer isec,jsec,ksec,lsec
      common/csecindices/isec,jsec,ksec,lsec
      integer iU,iS,iB,iA,iref
      common/cnlomaplabels/iU,iS,iB,iA,iref
      integer, parameter :: mxdim=30
      parameter(maxitn=12)
      double precision x0(mxdim),x(mxdim)
      double precision sNLO(nexternal,nexternal)
      double precision sLO(nexternal-1,nexternal-1)
      double precision KS,KHC,KNLO
      double precision lam,lim,RNLO,single_real
      character*5 str5
      character*8 limstr
      character*10 str10
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision xjac,xjacB
      double precision xsave(3)
      DOUBLE PRECISION ANS(0:1) !TODO SET CORRECTLY RANGE OF ANS
      DOUBLE PRECISION ALPHAS, ALPHA_QCD
      DOUBLE PRECISION Z_NLO
      DOUBLE PRECISION WGT,WGTPL,wgt_chan
      DOUBLE PRECISION SCM
      INTEGER, PARAMETER :: HEL=-1
      integer %(NLO_proc_str)sfl_factor 
      common/%(NLO_proc_str)sflavour_factor/%(NLO_proc_str)sfl_factor
      DOUBLE PRECISION ALPHAZ
      PARAMETER(ALPHAZ=1D0)
      common/cxsave/xsave
      double precision e(2),l(2)
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
      SCM = (2D0*EBEAM(1))**2
c     
c     initialise
      x=x0
      str5 ='     '
      str10='          '
      xjac=0d0
      sNLO=0d0
      sLO=0d0
      Z_NLO=0d0
      wgt_chan=1d0
c
c     TODO: MAP SOFT LIMIT AS (ilm), I.E. ONE MAPPING PER DIPOLE
c
c     start testing
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)'LIM = '//trim(limstr)
      write(iunit,*)str10//'lambda'//str10//str10//'R'//str10//str10//str5//'LIM'//str10//str10//'|R-LIM|/|LIM|'
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
         KNLO=0d0
c
c     rescale relevant x random numbers
c     x(1) is zCS, while x(2) is yCS
c     TODO: this rescaling is specific for (ijr) mapping; generalise 
         x(1)=abs(l(1)-x0(1)*lam**e(1))
         x(2)=abs(l(2)-x0(2)*lam**e(2))
c
c     set xsave so that the counterterms will be called with
c     more and more singular kinematics
         do i=1,3
            xsave(i)=x(i)
         enddo
c
c     recompute momenta after rescaling
         call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)
         if(xjac.eq.0d0.or.xjacb.eq.0d0)cycle
         call invariants_from_p(p,nexternal,sNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(pb,nexternal-1,sLO,ierr)
         if(ierr.eq.1)cycle
c
c     real
         call %(NLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
         RNLO = ANS(0) * %(NLO_proc_str)sfl_factor
         if(RNLO.lt.0d0.or.abs(RNLO).ge.huge(1d0).or.isnan(RNLO))cycle
         CALL GET_Z_NLO(SNLO,SCM,ALPHAZ,%(isec)d,%(jsec)d,Z_NLO,IERR)
         if(ierr.eq.1)cycle
c
c     counterterm
         call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,xjac,xjacB,x,KNLO,wgt_chan,ierr)
         if(ierr.eq.1)cycle
         
         lim=KNLO
         single_real=RNLO*Z_NLO*xjac
         
         if(abs(lim).gt.0d0)then
            write(iunit,*)lam,single_real,lim,abs(single_real-lim)/abs(lim)
         else
            write(iunit,*)lam,single_real,lim,single_real,' *** '
         endif
      enddo
      x=x0
c
      return
      end
