      subroutine local_counter_NLO_%(isec)d_%(jsec)d(xs,xp,xsb,xpb,wgt,xj,xjB,x,KNLO,wgt_chan,ierr)
c     local NLO counterterm for sector (isec,jsec)
      implicit none
      include 'nexternal.inc'
      integer ierr,nitR
      common/iterations/nitR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision wgt,xj,xjB
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      double precision KS,KC,KSC,KNLO,wgt_chan
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      logical default_soft
      parameter(default_soft=.true.)
      %(str_def_M2)s
c
c     initialise
      KS=0d0
      KC=0d0
      KSC=0d0
      KNLO=0d0
c
c     counterterms
      %(str_M2)s
c     
c     combination
      KNLO=KS+KC-KSC
c
      return
 999  ierr=1
      return
      end


