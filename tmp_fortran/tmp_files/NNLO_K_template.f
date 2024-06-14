      subroutine local_counter_NNLO_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,ZSi,ZSj,xj,xjB,x,KNNLO,wgt_chan,ierr)
c     wrapper for 3/4 particle sectors; 3p sector: ijk0, 4p sector: ijkl      
      implicit none

      call local_counter_NNLO_K1_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,wgt,ZSi,ZSj,xj,xjB,x,K1,wgt_chan,ierr)
      call local_counter_NNLO_K2_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,ZSi,ZSj,xj,xjB,x,K2,wgt_chan,ierr)
      call local_counter_NNLO_K12_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,ZSi,ZSj,xj,xjB,x,K12,wgt_chan,ierr)

c     combination
      KNNLO = K1+K2-K12
      
      end subroutine


      subroutine local_counter_NNLO_K1_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,wgt,ZSi,ZSj,xj,xjB,x,K1,wgt_chan,ierr)
c     local NNLO counterterm K1 for sector [isec,jsec,ksec,lsec]
      implicit none
      include 'nexternal.inc'
      integer isec,jsec,ksec,lsec,iref,ierr
      integer nitRR
      common/iterations/nitRR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision wgt,ZSi,ZSj,xj,xjB
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      double precision M2_S
      double precision KS,KHC,K1,wgt_chan
      %(str_defHC_K1)s
c
c     initialise
      isec = %(isec)d
      jsec = %(jsec)d
      ksec = %(ksec)d
      lsec = %(lsec)d
      iref = %(iref)d
      KS=0d0
      KHC=0d0
      K1=0d0
c
c     counterterms
      %(str_M2_K1)s
c     
c     combination
      K1=KS+KHC
c
      return
 999  ierr=1
      return
      end

      subroutine local_counter_NNLO_K2_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,ZSi,ZSj,xj,xjB,x,K2,wgt_chan,ierr)
c     local NNLO counterterm K2 for sector [isec,jsec,ksec,lsec]
      implicit none
      include 'nexternal.inc'
      integer isec,jsec,ksec,lsec,iref,ierr
      integer nitRR
      common/iterations/nitRR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision xpbb(0:3,nexternal-2)
      double precision wgt,ZSi,ZSj,xj,xjB
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      double precision K2,wgt_chan
      double precision KSS,KSC,KHCC,KCSHC
      double precision KHCC
      %(str_defK2)s
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
      double precision M2_S,M2_S_ALT
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
