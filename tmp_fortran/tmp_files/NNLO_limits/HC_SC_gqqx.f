      double precision function M2_HC_SC_gqqx(k,i,j,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     Quark-pair entry point: no symmetric d_ABC term occurs in (C.33).
      implicit none
      include 'nexternal.inc'
      integer k,i,j,r,nit,ierr
      double precision xs(nexternal,nexternal),xp(0:3,nexternal)
      double precision wgt,xj,xjb,extra,wgt_chan,M2_HC_SC_gluon
      M2_HC_SC_gqqx=M2_HC_SC_gluon(k,i,j,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
      end
