      subroutine phase_space_npo(x,shat,iU,iS,iB,iA,p,pbar,xjac,xjacB)
c     iU is the unresolved parton associated with the soft singularity
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'nexternal.inc'
      include 'leg_PDGs.inc'
      double precision x(3*nexternal-10),shat
      double precision p(0:3,nexternal),pbar(0:3,nexternal-1)
      double precision xjac,xjacB,xjacCS
      integer i,j,iU,iS,iB,iA
      integer iconfig,mincfig,maxcfig,invar
      integer ich
      common/comich/ich
c
c     initialise
      p=0d0
      pbar=0d0
      xjacB=Gevtopb
c
C     Hard coded settings for gen_mom
      iconfig = ich
      mincfig = 1
      maxcfig = 1
      invar = 2
      call gen_mom(iconfig,mincfig,maxcfig,invar,xjacB,x(4),pbar,nexternal-1)
c
c     call radiation phase space
      call phase_space_CS(x(1),iU,iS,iB,iA,p,pbar,nexternal,leg_PDGs,'C',xjacCS)
c
c     total jacobian
      xjac=xjacB*xjacCS
c
      return
      end
