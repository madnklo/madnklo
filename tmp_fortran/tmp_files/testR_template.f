      subroutine test_R_%(isec)d_%(jsec)d(ievnt)
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
      integer ievnt
      character*10 dash10
      double precision e(2), l(2)
c
      dash10='----------'
      write(88,*)dash10//dash10//dash10//dash10
      write(88,*)dash10//dash10//dash10//dash10
      write(88,*)' EVENT NUMBER ',ievnt
      write(88,*)dash10//dash10//dash10//dash10
      write(88,*)dash10//dash10//dash10//dash10
%(limit_str)s
c
      write(88,*)
      write(88,*)
      write(88,*)
c
      return
      end


      subroutine do_limit_R_%(isec)d_%(jsec)d(limstr,e,l)
      use sectors2_module
      implicit none
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer iitn,i,j,maxitn,ierr
      integer isec,jsec,ksec,lsec
      common/csecindices/isec,jsec,ksec,lsec
      integer iU,iS,iB,iA,iref
      common/cnlomaplabels/iU,iS,iB,iA,iref
      parameter(maxitn=12)
      DOUBLE PRECISION x0(3*nexternal-10),xr(3*nexternal-10)
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
      DOUBLE PRECISION ANS(0:1) !TODO SET CORRECTLY RANGE OF ANS
      DOUBLE PRECISION ALPHAS, ALPHA_QCD
      DOUBLE PRECISION WGT,WGTPL,wgt_chan
      DOUBLE PRECISION SCM
      INTEGER, PARAMETER :: HEL=-1
      integer %(NLO_proc_str)sfl_factor 
      common/%(NLO_proc_str)sflavour_factor/%(NLO_proc_str)sfl_factor
      DOUBLE PRECISION ALPHAZ
      PARAMETER(ALPHAZ=1D0)
      double precision e(2),l(2)
      double precision ran2
c
      SCM = (2D0*EBEAM(1))**2
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
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
c     rescale relevant xr random numbers
c     xr(1) is zCS, while xr(2) is yCS
         xr(1)=abs(l(1)-x0(1)*lam**e(1))
         xr(2)=abs(l(2)-x0(2)*lam**e(2))
c
c     recompute momenta after rescaling
         call phase_space_npo(xr,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)
         if(xjac.eq.0d0.or.xjacb.eq.0d0)cycle
         call invariants_from_p(p,nexternal,sNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(pb,nexternal-1,sLO,ierr)
         if(ierr.eq.1)cycle
         call get_sig2(SNLO,alphaZ,nexternal)
c
c     real
         call %(NLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
         RNLO = ANS(0) * %(NLO_proc_str)sfl_factor
         if(RNLO.lt.0d0.or.abs(RNLO).ge.huge(1d0).or.isnan(RNLO))cycle
         CALL GET_Z_NLO(%(isec)d,%(jsec)d)
c
c     counterterm
         call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,xjac,xjacB,KNLO,wgt_chan,ierr)
         if(ierr.eq.1)cycle
         
         lim=KNLO
         single_real=RNLO*Z_NLO*xjac
         
         if(abs(lim).gt.0d0)then
            write(88,*)lam,single_real,lim,abs(single_real-lim)/abs(lim)
         else
            write(88,*)lam,single_real,lim,single_real,' *** '
         endif
      enddo
c
      return
      end
