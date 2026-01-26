      subroutine local_RV_counter_NNLO_%(isec)d_%(jsec)d(xs,xp,xsb,xpb,wgt,xj,xjB,x,KRVNNLO,wgt_chan,ierr)
c     local real-virtual counterterm for sector (isec,jsec)
      implicit none
      include 'nexternal.inc'
      integer isec,jsec,iref,ierr
      integer nitRV
      common/iterations/nitRV
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision wgt,xj,xjB
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      double precision KS,KC,KSC,KRVNNLO,wgt_chan
      logical default_soft
      parameter(default_soft=.true.)
      %(str_def_M2)s
c
c     initialise
      isec = %(isec)d
      jsec = %(jsec)d
      iref = %(iref)d
      KS=0d0
      KC=0d0
      KSC=0d0
      KRVNNLO=0d0
c
c     counterterms
      %(str_M2)s
c     
c     combination
      KRVNNLO=KS+KC-KSC
c
      return
 999  ierr=1
      return
      end


