      subroutine local_counter_NLO_%(isec)d_%(jsec)d(xs,xp,xsb,xpb,wgt,WsumSi,WsumSj,xj,KS,KHC,KNLO,ierr)
c     local NLO counterterm for sector [isec,jsec]
      implicit none
c
c     TODO: pass nitR information
c
      include 'nexternal.inc'
      integer isec,jsec,iref,ierr
      integer nitR
      common/iterations/nitR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision wgt,WsumSi,WsumSj,xj
      double precision M2_S,KS,KHC,KNLO
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


