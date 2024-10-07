      subroutine local_counter_NNLO_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,x,KNNLO,wgt_chan,ierr)
c     wrapper for 3/4 particle sectors; 3p sector: ijk0, 4p sector: ijkl      
      implicit none
      include 'nexternal.inc'
      integer nitRR, ierr
      common/iterations/nitRR
      double precision xs(nexternal,nexternal), xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1), xpb(0:3,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2), xpbb(0:3,nexternal-2)
      double precision ZSi, ZSj
      double precision wgt,xj,xjB
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision KNNLO, K1, K2, K12, wgt_chan

      KNNLO = 0d0
      K1 = 0d0
      K2 = 0d0
      K12 = 0d0 

      call local_counter_NNLO_K1_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,wgt,xj,xjB,x,K1,wgt_chan,ierr)
      call local_counter_NNLO_K2_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,x,K2,wgt_chan,ierr)
      call local_counter_NNLO_K12_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,ZSi,ZSj,xj,xjB,x,K12,wgt_chan,ierr)

c     combination
      KNNLO = K1+K2-K12
      
      end subroutine


      subroutine local_counter_NNLO_K1_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,wgt,xj,xjB,x,K1,wgt_chan,ierr)
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
      double precision KS,KHC,K1,wgt_chan%(str_defK1)s
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

      subroutine local_counter_NNLO_K2_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,x,K2,wgt_chan,ierr)
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
      double precision KSS,KSHC,KHCC,KCCSHC%(str_defK2)s
c
c     initialise
      isec = %(isec)d
      jsec = %(jsec)d
      ksec = %(ksec)d
      lsec = %(lsec)d
      iref = %(iref)d
      KSS=0d0
      KSHC=0d0
      KHCC=0d0
      KCCSHC=0d0
      K2=0d0
c
c     counterterms
      %(str_M2_K2)s
c     
c     combination
      K2=KSS+KSHC+KHCC-KCCSHC
c
      return
 999  ierr=1
      return
      end

      subroutine local_counter_NNLO_K12_%(isec)d_%(jsec)d_%(ksec)d_%(lsec)d(xs,xp,xsb,xpb,xsbb,xpbb,wgt,ZSi,ZSj,xj,xjB,x,K12,wgt_chan,ierr)
c     local NNLO counterterm for sector [isec,jsec,ksec,lsec]
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
      double precision K12,wgt_chan
      double precision KS_SS, KS_SHC, KS_HCC
      double precision KHC_SS, KHC_SHC, KHC_HCC%(str_defK12)s
c
c     initialise
      isec = %(isec)d
      jsec = %(jsec)d
      ksec = %(ksec)d
      lsec = %(lsec)d
      iref = %(iref)d
      KS_SS=0d0
      KS_SHC=0d0
      KS_HCC=0d0
      KHC_SS=0d0
      KHC_SHC=0d0
      KHC_HCC=0d0 
      K12=0d0
c
c     counterterms
      %(str_M2_K12)s
c     
c     combination
      K12=KS_SS+KS_SHC+KS_HCC+KHC_SS+KHC_SHC+KHC_HCC
c
      return
 999  ierr=1
      return
      end
