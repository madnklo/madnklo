      logical function docut(p,npart,part_pdgs,nUnres)
      implicit none
      include 'jets.inc'
      include 'nexternal.inc'
      integer npart,i,j,nQCD,nUnres
      double precision p(0:3,npart)
      double precision rfj,sycut,palg,pQCD(0:3,maxdim),etamax
      logical dojets
      parameter(dojets=.true.)
      double precision dot
      double precision, parameter :: tiny=1d-8
      logical alliso
c     Sort array of results: ismode>0 for real, isway=0 for ascending order
      integer ismode,isway,izero,isorted(npart)
      parameter (ismode=1)
      parameter (isway=0)
      parameter (izero=0)
      integer k,nem,nin,nph
      double precision ptg
      logical is_a_lp(npart),is_a_lm(npart),is_a_ph(npart)
      double precision pgamma(0:3,npart),pem(0:3,npart),pt,eta
      double precision drlist(npart)
      double precision Etsum(0:npart)
      real * 8 chi_gamma_iso
      double precision r2,invm2,iso_getdrv
      integer part_pdgs(npart)
c
c     local variables
      double precision x1,x2
c     
c     initialise
      docut=.false.
      pjet=0d0
      jet=0
      njet=0
      ptjet=0d0
      ptg=0d0
      pgamma=0d0



c     


C***************************************************************
C***************************************************************
C     Cuts from the run_card.dat
C***************************************************************
C***************************************************************
c     

c     TODO:
c     leptons cuts must be read from the run_card


      ptl=10d0                  ! minimum pt for the charged leptons 
      etal=2.5d0                ! max rap for the charged leptons 
      drll=0.4d0                ! min distance between leptons 
!      drlal=0.4d0                ! min distance between gamma and lepton 
      mll=30d0
      drll_sf=0.4d0
      mll_sf=30d0
      maxjetflavor=4
c     CHARGED LEPTON CUTS
c     
c     find the charged leptons (also used in the photon isolation cuts below)
      do i=nincoming+1,npart
         if(part_pdgs(i).eq.11 .or. part_pdgs(i).eq.13 .or. 
     $        part_pdgs(i).eq.15) then
            is_a_lm(i)=.true.
         else
            is_a_lm(i)=.false.
         endif
         if(part_pdgs(i).eq.-11 .or. part_pdgs(i).eq.-13 .or. 
     $        part_pdgs(i).eq.-15) then
            is_a_lp(i)=.true.
         else
            is_a_lp(i)=.false.
         endif
      enddo
c     apply the charged lepton cuts
      do i=nincoming+1,npart
         if (is_a_lp(i).or.is_a_lm(i)) then
c     transverse momentum
            if (ptl.gt.0d0) then
               if (pt(p(0,i)).lt.ptl) then
                  docut=.true.
                  return
               endif
            endif
c     pseudo-rapidity
            if (etal.gt.0d0) then
               if (abs(eta(p(0,i))).gt.etal) then
                  docut=.true.
                  return
               endif
            endif
c     DeltaR and invariant mass cuts
            if (is_a_lp(i)) then
               do j=nincoming+1,npart
                  if (is_a_lm(j)) then
                     if (drll.gt.0d0) then
                        if (R2(p(0,i),p(0,j)).lt.drll**2) then
                           docut=.true.
                           return
                        endif
                     endif
                     if (mll.gt.0d0) then
                        if (invm2(p(0,i),p(0,j),1d0).lt.mll**2) then
                           docut=.true.
                           return
                        endif
                     endif
                     if (part_pdgs(i).eq.-part_pdgs(j)) then
                        if (drll_sf.gt.0d0) then
                           if (R2(p(0,i),p(0,j)).lt.drll_sf**2) then
                              docut=.true.
                              return
                           endif
                        endif
                        if (mll_sf.gt.0d0) then
                           if (invm2(p(0,i),p(0,j),1d0).lt.mll_sf**2)
     $                          then
                              docut=.true.
                              return
                           endif
                        endif
                     endif
                  endif
               enddo
            endif
         endif
      enddo

c
c     cluster partons into jets
      if(dojets)then
         nQCD=0
         do j=nincoming+1,npart
            if(abs(part_pdgs(j)) .le. maxjetflavor .or. part_pdgs(j)
     $           .eq. 21) then
               nQCD=nQCD+1
               do i=0,3
                  pQCD(i,nQCD)=p(i,j)
               enddo
            endif
         enddo
c     
c     clustering parameters
         palg=-1d0
         rfj=0.4d0
         sycut=20d0
         etamax = 5d0

c******************************************************************************
c     call FASTJET to get all the jets
c     
c     INPUT:
c     input momenta:               pQCD(0:3,n), energy is 0th component
c     number of input momenta:     nQCD
c     radius parameter:            rfj
c     minumum jet pt:              sycut
c     jet algorithm:               palg, 1.0=kt, 0.0=C/A, -1.0 = anti-kt
c     maximum jet rapidity:        etamax    
c     
c     OUTPUT:
c     jet momenta:                           pjet(0:3,n), E is 0th cmpnt
c     the number of jets (with pt > SYCUT):  njet
c     the jet for a given particle 'i':      jet(i),   note that this is the
c     particle in pQCD, which doesn't
c     necessarily correspond to the particle
c     label in the process
c
         call fastjetppgenkt_etamax(pQCD,nQCD,rfj,sycut,etamax,palg,pjet,njet,jet)
         
c         call fastjetppgenkt(pQCD,nQCD,rfj,sycut,palg,pjet,njet,jet)

c     
c******************************************************************************
         do i=1,njet
            ptjet(i)=sqrt(pjet(1,i)**2+pjet(2,i)**2)
            etajet(i)=eta(pjet(0,i))
            if(i.gt.1)then
               if (ptjet(i)-ptjet(i-1).gt.tiny) then
                  write (*,*) 'Error 1 in docut: jets unordered in pt'
                  stop
               endif
            endif
         enddo
      endif
c
c
c check
      if(nQCD.gt.nexternal-nincoming) then
         write(*,*)'Wrong nQCD in cuts.f ',nQCD,nexternal,nincoming
         stop
      endif
      if (nUnres .ne. 0 .and. nUnres .ne. 1) then
         write(*,*)'Wrong nUnres in cuts.f ',nUnres
         stop
      endif

      if (njet.lt.nQCD-nUnres) then
         docut = .true.
         return
      endif


c$$$
c$$$
c$$$      if (nUnres.eq.0) then
c$$$         if(njet.ne.nQCD) then
c$$$            docut = .true.
c$$$            return
c$$$         endif
c$$$      elseif(nUnres.eq.1) then
c$$$         if(njet.ne.nQCD .and. njet.ne.nQCD-1) then
c$$$            docut = .true.
c$$$            return
c$$$         endif
c$$$      else
c$$$         write(*,*)'Unknown number of unresolved particles',nUnres
c$$$         stop
c$$$      endif


c PHOTON (ISOLATION) CUTS
c


!     TODO: These parameters must be read in the run_card.dat

c$$$   ptgmin    ! Min photon transverse momentum
c$$$   etagamma  ! Max photon abs(pseudo-rap)
c$$$   R0gamma   ! Radius of isolation code
c$$$   xn        ! n parameter of eq.(3.4) in hep-ph/9801442
c$$$   epsgamma  ! epsilon_gamma parameter of eq.(3.4) in hep-ph/9801442
c$$$   isoEM  ! isolate photons from EM energy (photons and leptons)
      ptgmin = 20d0
      etagamma = -1d0
      R0gamma = 0.4d0
      xn = 1d0
      epsgamma = 1d0
      isoEM = .true.

c find the photons
      do i=nincoming+1,npart
         if (part_pdgs(i).eq.22) then
            is_a_ph(i)=.true.
         else
            is_a_ph(i)=.false.
         endif
      enddo
      if (ptgmin.ne.0d0) then
         nph=0
         do j=nincoming+1,npart
            if (is_a_ph(j)) then
               nph=nph+1
               do i=0,3
                  pgamma(i,nph)=p(i,j)
               enddo
            endif
         enddo
         if(nph.eq.0)goto 444
         
         if(isoEM)then
            nem=nph
            do k=1,nem
               do i=0,3
                  pem(i,k)=pgamma(i,k)
               enddo
            enddo
            do j=nincoming+1,npart
               if (is_a_lp(j).or.is_a_lm(j)) then
                  nem=nem+1
                  do i=0,3
                     pem(i,nem)=p(i,j)
                  enddo
               endif
            enddo
         endif
         
         alliso=.true.

         j=0
         do while(j.lt.nph.and.alliso)
c Loop over all photons
            j=j+1
            
            ptg=pt(pgamma(0,j))
            if(ptg.lt.ptgmin)then
               docut=.true.
               return
            endif
            if (etagamma.gt.0d0) then
               if (abs(eta(pgamma(0,j))).gt.etagamma) then
                  docut=.true.
                  return
               endif
            endif
         
c Isolate from hadronic energy
            do i=1,nQCD
               !drlist(i)=sngl(iso_getdrv(pgamma(0,j),pQCD(0,i)))
               drlist(i)=iso_getdrv(pgamma(0,j),pQCD(0,i))
            enddo
            call sortzv(drlist,isorted,nQCD,ismode,isway,izero)
            Etsum(0)=0.d0
            nin=0
            do i=1,nQCD
               if(drlist(isorted(i)).le.R0gamma)then
                  nin=nin+1
                  Etsum(nin)=Etsum(nin-1)+pt(pQCD(0,isorted(i)))
               endif
            enddo
            do i=1,nin
               alliso=alliso .and.
     $              Etsum(i).le.chi_gamma_iso(drlist(isorted(i)),
     $              R0gamma,xn,epsgamma,ptg)
            enddo

c$$$            
c$$$c Isolate from EM energy
            if(isoEM.and.nem.gt.1)then
               do i=1,nem
                  drlist(i)=iso_getdrv(pgamma(0,j),pem(0,i))
               enddo
               call sortzv(drlist,isorted,nem,ismode,isway,izero)
c First of list must be the photon: check this, and drop it
               if(isorted(1).ne.j.or.drlist(isorted(1)).gt.1.e-4)then
                  write(*,*)'Error #1 in photon isolation'
                  write(*,*)j,isorted(1),drlist(isorted(1))
                  stop
               endif
               Etsum(0)=0.d0
               nin=0
               do i=2,nem
                  if(drlist(isorted(i)).le.R0gamma)then
                     nin=nin+1
                     Etsum(nin)=Etsum(nin-1)+pt(pem(0,isorted(i)))
                  endif
               enddo
               do i=1,nin
                  alliso=alliso .and.
     $               Etsum(i).le.chi_gamma_iso(drlist(isorted(i)),
     $               R0gamma,xn,epsgamma,ptg)
               enddo
            endif
c End of loop over photons
         enddo
         if(.not.alliso)then
            docut=.true.
            return
         endif
 444     continue
c End photon isolation
      endif



      return
      end


* $Id: sortzv.F,v 1.1.1.1 1996/02/15 17:49:50 mclareni Exp $
*
* $Log: sortzv.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:50  mclareni
* Kernlib
*
*
c$$$#include "kerngen/pilot.h"
      SUBROUTINE SORTZV (A,INDEX,N1,MODE,NWAY,NSORT)
C
C CERN PROGLIB# M101    SORTZV          .VERSION KERNFOR  3.15  820113
C ORIG. 02/10/75
C
      DIMENSION A(N1),INDEX(N1)
C
C
      N = N1
      IF (N.LE.0)            RETURN
      IF (NSORT.NE.0) GO TO 2
      DO 1 I=1,N
    1 INDEX(I)=I
C
    2 IF (N.EQ.1)            RETURN
      IF (MODE)    10,20,30
   10 STOP 5 ! CALL SORTTI (A,INDEX,N)
      GO TO 40
C
   20 STOP 5 ! CALL SORTTC(A,INDEX,N)
      GO TO 40
C
   30   CALL SORTTF (A,INDEX,N)
C
   40 IF (NWAY.EQ.0) GO TO 50
      N2 = N/2
      DO 41 I=1,N2
      ISWAP = INDEX(I)
      K = N+1-I
      INDEX(I) = INDEX(K)
   41 INDEX(K) = ISWAP
   50 RETURN
      END



      double precision function pt(p)
      implicit none
      double precision p(0:3)

      pt = dsqrt(p(1)**2+p(2)**2)
      
      return
      end




      function eta(p)
      implicit none
      real*8 en,ptx,pty,pl,tiny,pt,eta,th
      real * 8 p(0:3)
      parameter (tiny=1.d-5)

      en = p(0)
      ptx = p(1)
      pty = p(2)
      pl = p(3)
c
      pt=sqrt(ptx**2+pty**2)
      if(pt.lt.tiny.and.abs(pl).lt.tiny)then
        eta=sign(1.d0,pl)*1.d8
      else
        th=atan2(pt,pl)
        eta=-log(tan(th/2.d0))
      endif
      return
      end




      function iso_getdrv(p1,p2)
      implicit none
      real*8 iso_getdrv,p1(0:3),p2(0:3)
      real*8 iso_getdr
c
      iso_getdrv=iso_getdr(p1(0),p1(1),p1(2),p1(3),
     #                     p2(0),p2(1),p2(2),p2(3))
      return
      end


      function chi_gamma_iso(dr,R0,xn,epsgamma,pTgamma)
c Eq.(3.4) of Phys.Lett. B429 (1998) 369-374 [hep-ph/9801442]
      implicit none
      real*8 chi_gamma_iso,dr,R0,xn,epsgamma,pTgamma
      real*8 tmp,axn
c
      axn=abs(xn)
      tmp=epsgamma*pTgamma
      if(axn.ne.0.d0)then
        tmp=tmp*( (1-cos(dr))/(1-cos(R0)) )**axn
      endif
      chi_gamma_iso=tmp
      return
      end


      function iso_getdr(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
      implicit none
      real*8 iso_getdr,en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2,deta,dphi,
     # iso_getpseudorap,iso_getdelphi
c
      deta=iso_getpseudorap(en1,ptx1,pty1,pl1)-
     #     iso_getpseudorap(en2,ptx2,pty2,pl2)
      dphi=iso_getdelphi(ptx1,pty1,ptx2,pty2)
      iso_getdr=sqrt(dphi**2+deta**2)
      return
      end


      function iso_getpseudorap(en,ptx,pty,pl)
      implicit none
      real*8 iso_getpseudorap,en,ptx,pty,pl,tiny,pt,eta,th
      parameter (tiny=1.d-5)
c
      pt=sqrt(ptx**2+pty**2)
      if(pt.lt.tiny.and.abs(pl).lt.tiny)then
        eta=sign(1.d0,pl)*1.d8
      else
        th=atan2(pt,pl)
        eta=-log(tan(th/2.d0))
      endif
      iso_getpseudorap=eta
      return
      end


      SUBROUTINE SORTTF (A,INDEX,N1)
C
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      SUBROUTINE SORTTI (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      SUBROUTINE SORTTC (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF(ICMPCH(AI,A(I22)))3,3,21
   21 INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (ICMPCH(A(I22),A(I222))) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (ICMPCH(AI,A(I22))) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      FUNCTION ICMPCH(IC1,IC2)
C     FUNCTION TO COMPARE TWO 4 CHARACTER EBCDIC STRINGS - IC1,IC2
C     ICMPCH=-1 IF HEX VALUE OF IC1 IS LESS THAN IC2
C     ICMPCH=0  IF HEX VALUES OF IC1 AND IC2 ARE THE SAME
C     ICMPCH=+1 IF HEX VALUES OF IC1 IS GREATER THAN IC2
      I1=IC1
      I2=IC2
      IF(I1.GE.0.AND.I2.GE.0)GOTO 40
      IF(I1.GE.0)GOTO 60
      IF(I2.GE.0)GOTO 80
      I1=-I1
      I2=-I2
      IF(I1-I2)80,70,60
 40   IF(I1-I2)60,70,80
 60   ICMPCH=-1
      RETURN
 70   ICMPCH=0
      RETURN
 80   ICMPCH=1
      RETURN
      END


      function iso_getdelphi(ptx1,pty1,ptx2,pty2)
      implicit none
      real*8 iso_getdelphi,ptx1,pty1,ptx2,pty2,tiny,pt1,pt2,tmp
      parameter (tiny=1.d-5)
c
      pt1=sqrt(ptx1**2+pty1**2)
      pt2=sqrt(ptx2**2+pty2**2)
      if(pt1.ne.0.d0.and.pt2.ne.0.d0)then
        tmp=ptx1*ptx2+pty1*pty2
        tmp=tmp/(pt1*pt2)
        if(abs(tmp).gt.1.d0+tiny)then
          write(*,*)'Cosine larger than 1'
          stop
        elseif(abs(tmp).ge.1.d0)then
          tmp=sign(1.d0,tmp)
        endif
        tmp=acos(tmp)
      else
        tmp=1.d8
      endif
      iso_getdelphi=tmp
      return
      end


      DOUBLE PRECISION FUNCTION R2(P1,P2)
c************************************************************************
c     Distance in eta,phi between two particles.
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p1(0:3),p2(0:3)
      integer i
c
c     External
c
      double precision eta,DELTA_PHI
      external eta,delta_phi
      real * 8 dphi,iso_getdelphi
      
c-----
c  Begin Code
c-----
      dphi=iso_getdelphi(p1(1),p1(2),p2(1),p2(2))
      R2 = DPHI**2+(eta(p1)-eta(p2))**2
      RETURN
      END




      DOUBLE PRECISION FUNCTION invm2(P1,P2,dsign)
c************************************************************************
c     Invarient mass of 2 particles
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p1(0:3),p2(0:3),dsign
c
c     Local
c      
      integer i
      double precision ptot(0:3)
c
c     External
c
      double precision dot
      external dot
c-----
c  Begin Code
c-----

      do i=0,3
         ptot(i)=p1(i)+dsign*p2(i)
      enddo
      invm2 = dot(ptot,ptot)
      RETURN
      END
