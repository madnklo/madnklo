      subroutine local_counter_NLO_%(isec)d_%(jsec)d(xs,xp,xsb,xpb,
     &wgt,WsumSi,WsumSj,xj,KS,KHC,KNLO,ierr)
c     local NLO counterterm for sector [isec,jsec]
      implicit none
c
c     TODO: pass nitR information
c
      include 'all_steps_info.inc'
c     all_steps_info.inc includes: 
c     parameters (npartNLO,npartLO)
      integer isec,jsec,ierr
      double precision xs(npartNLO,npartNLO)
      double precision xp(0:3,npartNLO)
      double precision xsb(npartLO,npartLO)
      double precision xpb(0:3,npartLO)
      double precision wgt,WsumSi,WsumSj,xj
      double precision M2_S,M2_H_C,KS,KHC,KNLO
c      double precision necessary_ct_list(5)
c     pass the list of necessary cts, in order:
c     [Si, Sj, Cij, SiCij, SjCij]
c     initialise
      isec = %(isec)d
      jsec = %(jsec)d
      KS=0d0
      KHC=0d0
      KNLO=0d0
c     counterterms
      %(str_M2)s
c     
c     soft counterterm
c      KS=KS+M2_S(isec,xs,xp,wgt,WsumSi,xj,nitR,1d0,ierr)*necessary_ct_list(1)
c      KS=KS+M2_S(jsec,xs,xp,wgt,WsumSj,xj,nitR,1d0,ierr)*necessary_ct_list(2)
c
c     hard-collinear counterterm
c      KHC=KHC+M2_H_C(isec,jsec,iref(isec,jsec),xs,xp,xsb,xpb,wgt,xj,nitR,1d0,ierr)*necessary_ct_list(3)
c
c     combination
      KNLO=KS+KHC
c
      return
 999  ierr=1
      return
      end
