      subroutine test_R_%(isec)d_%(jsec)d(iunit,x0)
      implicit none
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer i,iU,iS,iB,iA,iref
      common/cNLOmaplabels/iU,iS,iB,iA,iref
      integer iunit,ievnt
      INTEGER, PARAMETER :: MXDIM = 30
      double precision x0(mxdim)
      double precision e1,e2
      character*10 dash10
      save ievnt
      double precision xsave(3)
      common/cxsave/xsave
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


      subroutine do_limit_R_%(isec)d_%(jsec)d(iunit,limstr,x0,iU,iS,e1,e2)
      implicit none
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer iitn,i,j,maxitn,iunit,ierr
      integer isec,jsec
      common/cnlosecindices/isec,jsec
      integer iU,iS,iB,iA,iref
      integer, parameter :: mxdim=30
      parameter(maxitn=12)
      double precision x0(mxdim),x(mxdim)
      double precision sNLO(nexternal,nexternal)
      double precision sLO(nexternal-1,nexternal-1)
      double precision KS,KHC,KNLO
      double precision Wsum,WsumSi,WsumSj
      double precision lam,lim,RNLO,single_real
      double precision e1,e2
      character*5 str5
      character*8 limstr
      character*10 str10
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision xjac,xjacB
      double precision xsave(3)
      DOUBLE PRECISION ANS(0:1) !TODO SET CORRECTLY RANGE OF ANS
      DOUBLE PRECISION ALPHAS, ALPHA_QCD
      DOUBLE PRECISION Z_NLO,ZSUMSI,ZSUMSJ
      DOUBLE PRECISION WGT,WGTPL
      DOUBLE PRECISION SCM
      INTEGER, PARAMETER :: HEL=-1
      DOUBLE PRECISION W34,W35,W45
      DOUBLE PRECISION e3,e4,e5
      DOUBLE PRECISION ZSUMSIUSR,Z_NLOUSR
      DOUBLE PRECISION ALPHA
      PARAMETER(ALPHA=1D0)
      common/cxsave/xsave
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
      ZSUMSI=0d0
      ZSUMSJ=0d0
      Z_NLO=0d0
c
c     TODO: MAP SOFT LIMIT AS (ilm), I.E. ONE MAPPING PER DIPOLE
      IB = %(iref)d
      IA = 1
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
         KS=0d0
         KHC=0d0
         KNLO=0d0
c
c     rescale relevant x random numbers
c     x(1) is zCS, while x(2) is yCS 
         x(1)=x0(1)*lam**e1
         x(2)=x0(2)*lam**e2
c
c     set xsave so that the counterterms will be called with
c     more and more singular kinematics
         do i=1,3
            xsave(i)=x(i)
         enddo
c
c     recompute momenta after rescaling
         call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)
         if(xjac.eq.0d0)cycle
         call invariants_from_p(p,nexternal,sNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(pb,nexternal-1,sLO,ierr)
         if(ierr.eq.1)cycle
c
c     real
c         CALL EPEM_GDDX_ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
         call %(NLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
         RNLO = ANS(0)
         if(RNLO.lt.0d0.or.abs(RNLO).ge.huge(1d0).or.isnan(RNLO))cycle
         CALL GET_Z_NLO(SNLO,SCM,1D0,IU,IS,Z_NLO,'F',IERR)
         if(ierr.eq.1)cycle
c
c     counterterm
%(str_Zsum)s

         call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,ZsumSi,ZsumSj,xjac,KS,KHC,KNLO,ierr)
         if(ierr.eq.1)cycle
c
c         write(iunit,*) 'KS= ', KS
c         write(iunit,*) 'KHC= ', KHC
c         
         lim=KNLO
         single_real=RNLO*Z_NLO
c         e3=(sNLO(1,3)+sNLO(2,3))/sNLO(1,2)
c         e4=(sNLO(1,4)+sNLO(2,4))/sNLO(1,2)
c         e5=(sNLO(1,5)+sNLO(2,5))/sNLO(1,2)
c
c         W34 = sNLO(3,4)/e3/e4/sNLO(1,2)
c         W35 = sNLO(3,5)/e3/e5/sNLO(1,2)
c         W45 = sNLO(4,5)/e4/e5/sNLO(1,2)
c
c         Z_NLOusr=(1d0/e3/W34+1d0/e4/W34)/
c     $        (1d0/e3/W34+1d0/e3/W35+  
c     $        1d0/e4/W34+1d0/e4/W45+
c     $        1d0/e5/W35+1d0/e5/W45)
c
c         ZSUMSIusr=(1d0/e3/W34)/
c     $        (1d0/e3/W34+1d0/e3/W35) 
         
c         lim=ZSUMSI
c         single_real=Z_NLO

c         lim=ZSUMSIusr
c         single_real=Z_NLOusr


c         write(iunit, *)  'RATIOS=', ZSUMSI/ZSUMSIusr, Z_NLO/Z_NLOusr
         
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
