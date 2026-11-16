      subroutine phase_space_npo(x,iU,iS,iB,iA,p,pbar,xjac,xjacB)
c     iU is the unresolved parton associated with the soft singularity
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'nexternal.inc'
      include 'rv_mapped_labels_common.inc'
      include 'leg_PDGs.inc'
      double precision x(3*nexternal-10)
      double precision p(0:3,nexternal),pbar(0:3,nexternal-1)
      double precision xjac,xjacB,xjacCS
      integer i,j,iU,iS,iB,iA
c     mapped labels (real-emission label -> Born slot) of the mapping
c     with unresolved iU and emitter iS: those cached for the
c     real-virtual (rv_mapped_labels_common.inc, filled by
c     fill_mapped_labels) when the pair is filled, otherwise those of
c     the NLO real-emission sector (c_mapped_labels, set by
c     initialise_sector)
      integer nlo_mapped_labels(nexternal),ml(nexternal)
      common/c_mapped_labels/nlo_mapped_labels
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
      if(mapped_pair_done(iU,iS))then
         ml(:)=mapped_labels(:,iU,iS)
      else
         ml(:)=nlo_mapped_labels(:)
      endif
c
c     momentum fractions of the initial-state partons (labels of the
c     real emission), before the radiation: those of their Born slots
c     (phase_space_CS rescales the one of an initial-state emitter or
c     recoiler)
      do i=1,2
         if(ml(i).eq.1.or.ml(i).eq.2)then
            xbk_npo(i)=xbk(ml(i))
         else
            xbk_npo(i)=xbk(i)
         endif
         xbk_bar(i)=xbk_npo(i)
      enddo
c
c     call radiation phase space
      call phase_space_CS(x(1),iU,iS,iB,iA,p,pbar,nexternal,leg_PDGs,xjacCS,ml)
c
c     total jacobian
      xjac=xjacB*xjacCS
c
      return
      end
