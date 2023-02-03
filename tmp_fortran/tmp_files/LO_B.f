      subroutine Born_LO(sLO,BLO,ierr)
c     Born matrix element
      implicit none
      include 'math.inc'
c     TODO: nexternal changes according to 
c     the directory LO/NLO_R/NLO_V in which it is
c     generated! Check
      include 'nexternal.inc'
      integer ierr
      double precision BLO,pref
      double precision sLO(nexternal-1,nexternal-1),cth
c
c     initialise
      BLO=0d0
      ierr=0
c
c     TODO: read sCM,qLO(1)
      pref=(4d0*pi*alphaEM*qLO(1))**2*Nc/(2d0*sCM)*Gevtopb
      cth=max(-1d0,min(1d0,1d0-2d0*sLO(2,3)/sCM))
      BLO=pref*(1d0+cth**2)
c
      if(abs(BLO).ge.huge(1d0).or.isnan(BLO))then
         write(77,*)'Exception caught in Born_LO',BLO
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      subroutine cc_Born_LO(sLO,ccBLO,ierr)
c     colour-correlated Born matrix element
      implicit none
      include 'math.inc'
c     TODO: nexternal changes according to 
c     the directory LO/NLO_R/NLO_V in which it is
c     generated! Check
      include 'nexternal.inc'
      integer ierr
      double precision BLO,ccBLO(nexternal-1,nexternal-1),pref
      double precision sLO(nexternal-1,nexternal-1),cth
c
c     initialise
      BLO=0d0
      ccBLO=0d0
      ierr=0
c
c     TODO: read sCM,qLO(1)     
      pref=(4d0*pi*alphaEM*qLO(1))**2*Nc/(2d0*sCM)*Gevtopb
      cth=max(-1d0,min(1d0,1d0-2d0*sLO(2,3)/sCM))
      BLO=pref*(1d0+cth**2)
c
      ccBLO(3,3)= CF*BLO
      ccBLO(3,4)=-CF*BLO
      ccBLO(4,3)=-CF*BLO
      ccBLO(4,4)= CF*BLO
c
      return
 999  ierr=1
      return
      end


      subroutine Born_LO_kp(xsLO,pkt,ktkt,BLOkp,ierr)
c     This routine returns B^munu kt_mu kt_nu
c     kt is passed as argument through pkt and ktkt
c     (i.e. p.kt and kt.kt)
      implicit none
      include 'math.inc'
c     TODO: nexternal changes according to 
c     the directory LO/NLO_R/NLO_V in which it is
c     generated! Check
      include 'nexternal.inc'
      integer ierr
      double precision xsLO(nexternal-1,nexternal-1)
      double precision pref
      double precision BLOkp,pkt(nexternal-1),ktkt
c
c     initialise
      BLOkp=0d0
      ierr=0
c
      if(abs(BLOkp).ge.huge(1d0).or.isnan(BLOkp))then
         write(77,*)'Exception caught in Born_LO_kp',BLOkp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
