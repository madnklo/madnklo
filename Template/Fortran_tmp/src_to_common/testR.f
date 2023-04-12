      subroutine test_R(iunit,x0,isec,jsec)
      implicit none
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer i,isec,jsec,iU,iS,iB,iA,iref
!      common/cNLOsecindices/isec,jsec
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
c
c     soft limit
c      if(is_NLO_S_singular_sec(isec,jsec))then
      e1=1d0
      e2=1d0
      call do_limit_R(iunit,'S       ',x0,isec,jsec,e1,e2)
c      endif
c
c     collinear limit
c      if(is_NLO_C_singular_sec(isec,jsec))then
c         e1=0d0
c         e2=1d0
c         call do_limit_R(iunit,'C       ',x0,isec,jsec,alpha,beta,e1,e2)
c      endif
c
c     soft-collinear limit
c      if(is_NLO_S_singular_sec(isec,jsec).and.
c     &   is_NLO_C_singular_sec(isec,jsec))then
c         e1=1d0
c         e2=2d0
c         call do_limit_R(iunit,'SC      ',x0,isec,jsec,alpha,beta,e1,e2)
c      endif
c
c     reinstate original xsave after testing
      do i=1,3
         xsave(i)=x0(i)
      enddo
c
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)
c
      return
      end


      subroutine do_limit_R(iunit,limstr,x0,isec,jsec,e1,e2)
      implicit none
      INCLUDE 'math.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      integer iitn,i,j,maxitn,iunit,ierr
      integer isec,jsec
      integer iU,iS,iB,iA,iref
      integer, parameter :: mxdim=30
      parameter(maxitn=12)
      double precision x0(mxdim),x(mxdim)
      double precision sNLO(nexternal,nexternal)
      double precision sLO(nexternal-1,nexternal-1)
      double precision KS,KHC,KNLO
!      double precision W_NLO(maxdim,maxdim)
!      double precision WS_NLO(maxdim,maxdim),WC_NLO(maxdim,maxdim)
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
c
c     assign iU, iS, iB, iA, iref
c      call assign_phsp_labels_npo(isec,jsec,iU,iS,iB,iA,iref)

      IREF = 5
      IU = ISEC
      IS = JSEC
      IB = IREF
      IA=1
c
c     start testing
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)'LIM = '//trim(limstr)
      write(iunit,*)str10//'lambda'//str10//str10//'R'//
     &str10//str10//str5//'LIM'//
     &str10//str10//'|R-LIM|/|LIM|'
c
c     possibility to set by hand the starting point
c     for the limiting procedure
c$$$      x0(1)=0.5d0
c$$$      x0(2)=0.5d0
c
c     loop to get closer and closer to the limit
      do iitn=1,maxitn
         lam=10d0**(1d0-iitn)
c
c     initialise counterterms
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


c     recompute momenta after rescaling
         call phase_space_npo(x,sCM,iU,iS,iB,iA,p,pb,xjac,xjacB)


         if(xjac.eq.0d0)cycle
         call invariants_from_p(p,nexternal,sNLO,ierr)
         if(ierr.eq.1)cycle
         call invariants_from_p(pb,nexternal-1,sLO,ierr)
         if(ierr.eq.1)cycle
c
c     real
         CALL EPEM_GDDX_ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
         RNLO = ANS(0)
!         if(ierr.eq.1)cycle
!     call get_W_NLO(xsNLO,alpha,beta,isec,jsec,W_NLO,WS_NLO,WC_NLO,ierr)
         CALL GET_Z_NLO(SNLO,SCM,1D0,ISEC,JSEC,Z_NLO,'F',IERR)
         if(ierr.eq.1)cycle
c
c     combinations of W functions
c         Wsum=W_NLO(isec,jsec)+W_NLO(jsec,isec)
c         WsumSi=WS_NLO(isec,jsec)
c         WsumSj=WS_NLO(jsec,isec)
c
c     counterterm
!         call local_counter_NLO(isec,jsec,iref,xsNLO,p,xsLO,pb,iA,
!     &        1d0,WsumSi,WsumSj,1d0,KS,KHC,KNLO,ierr)

         
         CALL GET_Z_NLO(SNLO,SCM,1d0,ISEC,JSEC,ZSUMSI,'S',IERR)
 
         
         CALL LOCAL_COUNTER_NLO_3_4(SNLO,P,SLO,PB,WGT,ZSUMSI,ZSUMSJ,XJAC
     $        ,KS,KHC,KNLO,IERR)


!        write(iunit,*) 'ZSUMSI=', ZSUMSI
         
         if(ierr.eq.1)cycle
c
         write(iunit,*) 'KS= ', KS
         write(iunit,*) 'KHC= ', KHC
         
         lim=KNLO
         single_real=RNLO*Z_NLO


c$$$         e3=(sNLO(1,3)+sNLO(2,3))/sNLO(1,2)
c$$$         e4=(sNLO(1,4)+sNLO(2,4))/sNLO(1,2)
c$$$         e5=(sNLO(1,5)+sNLO(2,5))/sNLO(1,2)
c$$$
c$$$         W34 = sNLO(3,4)/e3/e4/sNLO(1,2)
c$$$         W35 = sNLO(3,5)/e3/e5/sNLO(1,2)
c$$$         W45 = sNLO(4,5)/e4/e5/sNLO(1,2)
c$$$
c$$$         Z_NLOusr=(1d0/e3/W34+1d0/e4/W34)/
c$$$     $        (1d0/e3/W34+1d0/e3/W35+  
c$$$     $        1d0/e4/W34+1d0/e4/W45+
c$$$     $        1d0/e5/W35+1d0/e5/W45)
c$$$
c$$$
c$$$
c$$$         ZSUMSIusr=(1d0/e3/W34)/
c$$$     $        (1d0/e3/W34+1d0/e3/W35) 
         
!         lim=ZSUMSI
!         single_real=Z_NLO

!         lim=ZSUMSIusr
!         single_real=Z_NLOusr


!         write(iunit, *)  'RATIOS=', ZSUMSI/ZSUMSIusr, Z_NLO/Z_NLOusr
         
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
