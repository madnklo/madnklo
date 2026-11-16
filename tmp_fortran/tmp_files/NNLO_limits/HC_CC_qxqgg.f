      double precision function M2_HC_CC_qxqgg(i,j,k,l,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
      implicit none
      include 'nexternal.inc'
      integer i,j,k,l,nit,ierr
      double precision xs(nexternal,nexternal),xp(0:3,nexternal)
      double precision wgt,xj,xjb,extra,wgt_chan,M2_CC_gluon_parents
      M2_HC_CC_qxqgg=M2_CC_gluon_parents(i,j,k,l,2,1,.true.,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
      end
