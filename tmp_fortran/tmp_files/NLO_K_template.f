      subroutine local_counter_NLO_%(isec)d_%(jsec)d(xs,xp,xsb,xpb,wgt,ZSi,ZSj,xj,xjB,x,KS,KHC,KNLO,wgt_chan,ierr)
c     local NLO counterterm for sector [isec,jsec]
      implicit none
      include 'nexternal.inc'
      integer isec,jsec,iref,ierr
      integer nitR
      common/iterations/nitR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision wgt,ZSi,ZSj,xj,xjB
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      double precision M2_S_G
      double precision KS,KHC,KNLO,wgt_chan
      logical default_soft
      parameter(default_soft=.true.)
      %(str_defHC)s
c
c     initialise
      isec = %(isec)d
      jsec = %(jsec)d
      iref = %(iref)d
      KS=0d0
      KHC=0d0
      KNLO=0d0
c
c     counterterms
      %(str_M2)s
c     
c     combination
      KNLO=KS+KHC
c
      return
 999  ierr=1
      return
      end


