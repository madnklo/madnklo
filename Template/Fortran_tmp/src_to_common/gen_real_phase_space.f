      subroutine phase_space_npo(x,shat,iU,iS,iB,iA,p,pbar,xjac,xjacB)
c     iU is the unresolved parton associated with the soft singularity
      implicit none
      include 'math.inc'
      include 'nexternal.inc'
      include 'leg_PDGs.inc'
      double precision x(7),shat
      double precision p(0:3,nexternal),pbar(0:3,nexternal-1)
      double precision xjac,xjacB,xjacCS
      integer i,j,iU,iS,iB,iA
c
c     initialise
      p=0d0
      pbar=0d0
      xjac=0d0
c
c     call Born phase space
      call phase_space_n(x(4),shat,pbar,nexternal-1,xjacB)
c
c     call radiation phase space
      call phase_space_CS(x(1),iU,iS,iB,iA,p,pbar,nexternal,leg_PDGs,'C',xjacCS)
c
c     total jacobian
      xjac=xjacB*xjacCS
c
      return
      end
