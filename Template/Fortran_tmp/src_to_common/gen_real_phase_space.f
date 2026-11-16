      subroutine phase_space_npo(x,iU,iS,iB,iA,p,pbar,xjac,xjacB)
c     iU is the unresolved parton associated with the soft singularity
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'nexternal.inc'
      include 'leg_PDGs.inc'
      double precision x(3*nexternal-10)
      double precision p(0:3,nexternal),pbar(0:3,nexternal-1)
      double precision xjac,xjacB,xjacCS
      integer i,j,iU,iS,iB,iA
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      double precision xbk_npo(2),xbk_bar(2)
      common/c_xbk_mapped/xbk_npo,xbk_bar
      include 'run.inc'
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
c     momentum fractions of the initial-state partons (labels of the
c     real emission), before the radiation: those of their Born slots
c     (phase_space_CS rescales the one of an initial-state emitter or
c     recoiler)
      do i=1,2
         if(mapped_labels(i).eq.1.or.mapped_labels(i).eq.2)then
            xbk_npo(i)=xbk(mapped_labels(i))
         else
            xbk_npo(i)=xbk(i)
         endif
         xbk_bar(i)=xbk_npo(i)
      enddo
c
c     call radiation phase space
      call phase_space_CS(x(1),iU,iS,iB,iA,p,pbar,nexternal,leg_PDGs,xjacCS,mapped_labels)
c
c     total jacobian
      xjac=xjacB*xjacCS
c
      return
      end
