      subroutine local_counter_NNLO_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,x,K1,K2,K12,KNNLO,wgt_chan,ierr)
c     wrapper for 3/4 particle sectors; 3p sector: ijjk & ijkj, 4p sector: ijkl
      implicit none
      include 'nexternal.inc'
      integer nitRR, ierr
      common/iterations/nitRR
      double precision xs(nexternal,nexternal), xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1), xpb(0:3,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2), xpbb(0:3,nexternal-2)
      double precision wgt,xj,xjB
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision KNNLO, K1, K2, K12, wgt_chan

      KNNLO = 0d0
      K1 = 0d0
      K2 = 0d0
      K12 = 0d0

      call local_counter_NNLO_K1_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(xs,xp,xsb,xpb,wgt,xj,xjB,x,K1,wgt_chan,ierr)
      call local_counter_NNLO_K2_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,x,K2,wgt_chan,ierr)
      call local_counter_NNLO_K12_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,x,K12,wgt_chan,ierr)

c     combination
      KNNLO = K1+K2-K12

      end subroutine


      subroutine local_counter_NNLO_K1_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(xs,xp,xsb,xpb,wgt,xj,xjB,x,K1,wgt_chan,ierr)
c     local NNLO counterterm K1 for sector [%(isec)d,%(jsec)d,%(c3p)d,%(d3p)d]
      implicit none
      include 'nexternal.inc'
      integer isec,jsec,ksec,lsec,iref
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer nitRR
      common/iterations/nitRR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision wgt,xj,xjB
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      double precision KS,KC,KSC,K1,wgt_chan%(str_defK1)s
c
c     initialise
      KS=0d0
      KC=0d0
      KSC=0d0
      K1=0d0
c
c     counterterms
      %(str_M2_K1)s
c
c     combination
      K1=KS+KC-KSC
c
      return
 999  ierr=1
      return
      end

      subroutine local_counter_NNLO_K2_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,x,K2,wgt_chan,ierr)
c     local NNLO counterterm K2 for sector [isec,jsec,ksec,lsec]
      implicit none
      include 'nexternal.inc'
      integer isec,jsec,ksec,lsec,iref,ierr
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer nitRR
      common/iterations/nitRR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision xpbb(0:3,nexternal-2)
      double precision wgt,xj,xjB
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      double precision K2,wgt_chan
      double precision KSS,KSC,KCC%(str_defK2)s
c
c     initialise
      KSS=0d0
      KSC=0d0
      KCC=0d0
      K2=0d0
c
c     counterterms
      %(str_M2_K2)s
c
c     combination
      K2=KSS+KSC+KCC
c
      return
 999  ierr=1
      return
      end

      subroutine local_counter_NNLO_K12_%(isec)d_%(jsec)d_%(c3p)d_%(d3p)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,x,K12,wgt_chan,ierr)
c     local NNLO counterterm for sector [isec,jsec,ksec,lsec]
      implicit none
      include 'nexternal.inc'
      integer isec,jsec,ksec,lsec,iref,ierr
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer nitRR
      common/iterations/nitRR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision xpbb(0:3,nexternal-2)
      double precision wgt,xj,xjB
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      double precision K12,wgt_chan
      double precision KS_SS, KS_SC, KS_CC
      double precision KHC_SS, KHC_SC, KHC_CC%(str_defK12)s
c
c     initialise
      KS_SS=0d0
      KS_SC=0d0
      KS_CC=0d0
      KHC_SS=0d0
      KHC_SC=0d0
      KHC_CC=0d0
      K12=0d0
c
c     counterterms
      %(str_M2_K12)s
c
c     combination
      K12=KS_SS+KS_SC+KS_CC+KHC_SS+KHC_SC+KHC_CC
c
      return
 999  ierr=1
      return
      end
