      subroutine phase_space_npt(x,shat,iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,p,pbar,ptilde,xjac,xjacB,xjacCS1)
c     iU1 and iU2 are the unresolved partons associated with the soft singularity
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'nexternal.inc'
      include 'leg_PDGs.inc'
      double precision x(3*nexternal-10),shat
      double precision p(0:3,nexternal),pbar(0:3,nexternal-1),ptilde(0:3,nexternal-2)
      double precision xjac,xjacB,xjacCS1,xjacCS2
      integer i,j,iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
      integer iconfig,mincfig,maxcfig,invar
      integer ich
      common/comich/ich
c
c     initialise
      p=0d0
      pbar=0d0
      ptilde=0d0
      xjacB=Gevtopb
c
C     Hard coded settings for gen_mom
      iconfig = ich
      mincfig = 1
      maxcfig = 1
      invar = 2
      call gen_mom(iconfig,mincfig,maxcfig,invar,xjacB,x(7),ptilde,nexternal-2)
c
c     call radiation phase space from Born to real
      call phase_space_CS(x(4),iU2,iS2,iB2,iA2,pbar,ptilde,nexternal-1,leg_PDGs,'C',xjacCS2)
c
c     call radiation phase space from real to double real
      call phase_space_CS(x(1),iU1,iS1,iB1,iA1,p,pbar,nexternal,leg_PDGs,'C',xjacCS1)
c
c     total jacobian
      xjac=xjacB*xjacCS1*xjacCS2
c
      return
      end
