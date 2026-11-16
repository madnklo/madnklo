      subroutine get_mapped_labels(npart,a,b,pdgs,underlying_pdgs,mapped_labels)
c     assigns labels and flavours of particles after a mapping (a,b,y)
c     that clusters partons (a,b) in an npart-body final state
c     inputs : npart,a,b,pdgs,underlying_pdgs
c     outputs: mapped_labels,mapped_flavours,mapped_indices_shuff
c     With one of a, b in the initial state (initial-state splitting),
c     see get_mapped_labels_initial below.
      implicit none
      integer npart,a,b,i,jj,i1,i2
      integer pdgs(npart),mapped_labels(npart),target_flav(npart)
      integer underlying_pdgs(npart-1)
      logical isgg,isqq,isqg,isgq,assigned(npart-1)
c
c     initialise
      if (a.eq.b) then
         write(*,*)'mapped_labels cannot be identical'
         write(*,*)a,b,pdgs(a),pdgs(b)
         stop
      endif
      i1=min(a,b)
      i2=max(a,b)
      mapped_labels(:)=0
      target_flav(:)=pdgs(:)
c
c     initial-state splitting
      if(i1.le.2)then
         call get_mapped_labels_initial(npart,i1,i2,pdgs,
     &        underlying_pdgs,mapped_labels)
         return
      endif
c
c     preliminary checks
      if(i1.lt.3.or.i2.gt.npart)then
         write(*,*)'wrong a, b in get_collinear_mapped_labels',a,b
         stop
      endif
c     
c     possible QCD pairs
      isqq=abs(pdgs(i1)).le.6.and.pdgs(i1)+pdgs(i2).eq.0
      isgg=pdgs(i1).eq.21.and.pdgs(i2).eq.21
      isqg=abs(pdgs(i1)).le.6.and.pdgs(i2).eq.21
      isgq=pdgs(i1).eq.21.and.abs(pdgs(i2)).le.6
      if (.not.(isgg.or.isqq.or.isqg.or.isgq)) then
         write(*,*)'inconsistent pair in get_collinear_mapped_labels'
         write(*,*)i1,i2,a,b,pdgs(a),pdgs(b)
         stop
      endif
      if(isqq.or.isgg)target_flav(i2)=21
      if(isqg)target_flav(i2)=pdgs(i1)
      if(isgq)target_flav(i2)=pdgs(i2)
c
c     explicit examples of output mapped_labels
c=======================================================
c     example 1, npart = 6:
c     i                        1  (2)  3  (4)  5   6
c     pdgs(i)                  u   ub  g   g   d   db
c     jj                       1       2   3   4   5
c     underlying_pdgs(jj)      g       u   d   ub  db
c     mapped_labels(i)         2   4   1   4   3   5
c=======================================================
c=======================================================
c     example 2, npart = 6:
c     i                        1   2   3  (4)  5  (6)
c     pdgs(i)                  u   ub  g   g   d   db
c     jj                       1   2   3       4   5
c     underlying_pdgs(jj)      g   u   d       ub  db
c     mapped_labels(i)         2   4   1   5   3   5
c=======================================================
c=======================================================
c     example 3, npart = 6:
c     i                        1   2  (3) (4)  5   6
c     pdgs(i)                  u   ub  g   g   d   db
c     jj                       1   2       3   4   5
c     underlying_pdgs(jj)      g   u       d   ub  db
c     mapped_labels(i)         2   4   1   1   3   5      
c=======================================================
c=======================================================
c     example 4, npart = 6:
c     i                        1   2   3   4  (5) (6)
c     pdgs(i)                  u   ub  g   g   d   db
c     jj                       1   2   3   4       5
c     underlying_pdgs(jj)      g   u   g   ub      g
c     mapped_labels(i)         2   4   1   3   5   5      
c=======================================================
c=======================================================
c     example 5, npart = 6:
c     i                       (1) (2)  3   4   5   6
c     pdgs(i)                  u   ub  g   g   d   db
c     jj                           1   2   3   4   5
c     underlying_pdgs(jj)          g   d   g   db  g
c     mapped_labels(i)         1   1   3   5   2   4
c=======================================================
c
      assigned=.false.
      do i=1,npart
         if(i.eq.i1)cycle
         do jj=1,npart-1
            if(underlying_pdgs(jj).eq.target_flav(i).and.
     &           .not.assigned(jj))then
               assigned(jj)=.true.
               mapped_labels(i)=jj
c              stop at the first free match: without this exit, a particle
c              would take all free underlying partons of its flavour, and a
c              second identical-flavour particle would be left with label 0
               exit
            endif
         enddo
      enddo
      mapped_labels(i1)=mapped_labels(i2)
c
c     check
      if(any(.not.assigned))then
         write(*,*)'Something went wrong in mapped_labels'
         stop
      endif
c
      return
      end


      subroutine get_iref(npart,i,j,mapped_labels,iref_r)
c     recoiler iref_r of the (n+1)-body sector (i,j): the parton whose
c     Born label is the recoiler, in virtual_recoilers.inc, of the Born
c     parent of i and j (one entry (parent, recoiler) per Born parton,
c     initial-state ones included)
      implicit none
      integer npart,i,j,mapped_labels(npart)
      integer ii,iref_v,iref_r
      include 'virtual_recoilers.inc'
c
      if(mapped_labels(i).ne.mapped_labels(j))then
         write(*,*)'Wrong mapped labels in get_iref'
         write(*,*)i,j,mapped_labels(i),mapped_labels(j)
         stop
      endif
      iref_v=0
      do ii=1,len_iref
         if(iref(1,ii).eq.mapped_labels(i))iref_v=iref(2,ii)
      enddo
      iref_r=0
      do ii=1,npart
         if(ii.eq.i.or.ii.eq.j)cycle
         if(mapped_labels(ii).eq.iref_v)iref_r=ii
      enddo
      if(iref_v.eq.0.or.iref_r.eq.0)then
         write(*,*)'No recoiler found in get_iref',i,j,iref_v,iref_r
         stop
      endif
c
      return
      end


      subroutine get_mapped_labels_initial(npart,a,b,pdgs,
     &     underlying_pdgs,mapped_labels)
c     Labels of the Born (npart-1)-body kinematics for an initial-state
c     splitting: a = 1 or 2 is the initial-state parton j, b > 2 the
c     final-state parton i, removed. The parent [ij] takes the place of
c     j in the initial state, with flavour (2209.09123, Sec. 2)
c        j = q, i = g  ->  q          j = g, i = g  ->  g
c        j = g, i = q  ->  qbar       j = q, i = q  ->  g
c     (pdgs of initial-state partons as incoming particles). The other
c     partons are matched to the Born ones by flavour, in order,
c     initial-state ones to the Born initial states and final-state ones
c     to the Born final states. The Born may have the two initial states
c     in the other order (sectors.py, Born_candidates): the matching by
c     flavour then exchanges their labels.
      implicit none
      integer npart,a,b,i,jj
      integer pdgs(npart),mapped_labels(npart),target_flav(npart)
      integer underlying_pdgs(npart-1)
      logical assigned(npart-1)
c
      if(a.lt.1.or.a.gt.2.or.b.le.2.or.b.gt.npart)then
         write(*,*)'wrong a, b in get_mapped_labels_initial',a,b
         stop
      endif
      target_flav(:)=pdgs(:)
      if(abs(pdgs(a)).le.6.and.pdgs(b).eq.21)then
         target_flav(a)=pdgs(a)
      elseif(pdgs(a).eq.21.and.pdgs(b).eq.21)then
         target_flav(a)=21
      elseif(pdgs(a).eq.21.and.abs(pdgs(b)).le.6)then
         target_flav(a)=-pdgs(b)
      elseif(abs(pdgs(a)).le.6.and.pdgs(b).eq.pdgs(a))then
         target_flav(a)=21
      else
         write(*,*)'inconsistent pair in get_mapped_labels_initial'
         write(*,*)a,b,pdgs(a),pdgs(b)
         stop
      endif
c
      mapped_labels(:)=0
      assigned=.false.
      do i=1,npart
         if(i.eq.b)cycle
         do jj=1,npart-1
            if((i.le.2).neqv.(jj.le.2))cycle
            if(underlying_pdgs(jj).eq.target_flav(i).and.
     &           .not.assigned(jj))then
               assigned(jj)=.true.
               mapped_labels(i)=jj
               exit
            endif
         enddo
      enddo
      mapped_labels(b)=mapped_labels(a)
c
c     check
      if(any(.not.assigned))then
         write(*,*)'Something went wrong in mapped_labels_initial'
         write(*,*)pdgs,underlying_pdgs
         stop
      endif
c
      return
      end



      subroutine get_unp_mapped_labels(npart,a,b,mapped_labels)
      implicit none
      integer npart, a, b
      integer mapped_labels(npart)
      integer mother, daughter
      integer i
c
      mapped_labels=0
      mother = min(a,b)
      daughter = max(a,b)
c
      do i=1,npart
         if(i.eq.daughter) cycle
         if(i.lt.daughter) then
            mapped_labels(i)=i
         elseif(i.gt.daughter) then
            mapped_labels(i) = i-1
         endif
      enddo
      mapped_labels(daughter) = mapped_labels(mother)
      return
      end



c$$$      subroutine get_soft_mapped_labels(a,n,leg_pdgs,mapped_labels,mapped_flavours,ismappedQCDparton)
c$$$c     assigns labels and flavours of particles after a mapping (a,x,y) that removes gluon a from
c$$$c     an n-body final state
c$$$      implicit none
c$$$      integer a,n,i
c$$$      integer leg_pdgs(n),mapped_labels(n),mapped_flavours(n)
c$$$      logical ismappedQCDparton(n-1)
c$$$c
c$$$c     initialise
c$$$      mapped_labels=0
c$$$      mapped_flavours=0
c$$$      ismappedQCDparton=.false.
c$$$c
c$$$c     preliminary checks
c$$$      if(a.lt.3)then
c$$$         write(*,*)'get_soft_mapped_labels: wrong parton a ',a
c$$$         stop
c$$$      endif
c$$$      if(leg_pdgs(a).ne.21)then
c$$$         write(*,*)'get_soft_mapped_labels: a is not a gluon',a,leg_pdgs(a)
c$$$         stop
c$$$      endif
c$$$c
c$$$c     assign mapped labels, flavours, ismappedQCDparton
c$$$      do i=1,n
c$$$         if(i.eq.a)cycle
c$$$         mapped_flavours(i)=leg_pdgs(i)
c$$$         if(i.lt.a)then
c$$$            mapped_labels(i)=i
c$$$         else
c$$$            mapped_labels(i)=i-1
c$$$         endif
c$$$         if(abs(mapped_flavours(i)).le.6.or.mapped_flavours(i).eq.21)
c$$$     &        ismappedQCDparton(mapped_labels(i)) = .true.
c$$$      enddo
c$$$c
c$$$      return
c$$$      end
c$$$
c$$$
c$$$
c$$$      subroutine get_collinear_mapped_labels(a,b,n,leg_pdgs,mapped_labels,mapped_flavours)
c$$$c     assigns labels and flavours of particles after a mapping (a,b,y) that clusters partons
c$$$c     (a,b) in an n-body final state
c$$$      implicit none
c$$$      integer a,b,n,i
c$$$      integer leg_pdgs(n),mapped_labels(n),mapped_flavours(n)
c$$$      logical isgluon,isqqbar,isQCD
c$$$c
c$$$c     initialise
c$$$      mapped_labels=0
c$$$      mapped_flavours=0
c$$$c
c$$$c     preliminary checks
c$$$      if(a.lt.3.or.b.lt.3)then
c$$$         write(*,*)'get_collinear_mapped_labels: wrong partons a, b ',a,b
c$$$         stop
c$$$      endif
c$$$      isgluon=leg_pdgs(a).eq.21.or.leg_pdgs(b).eq.21
c$$$      isqqbar=leg_pdgs(a)+leg_pdgs(b).eq.0
c$$$      isQCD=(abs(leg_pdgs(a)).le.6.or.leg_pdgs(a).eq.21).and.
c$$$     &      (abs(leg_pdgs(b)).le.6.or.leg_pdgs(b).eq.21)
c$$$      if (.not.(isgluon.or.isqqbar.or.isQCD)) then
c$$$         write(*,*)'get_collinear_mapped_labels: inconsistent a, b '
c$$$         write(*,*)leg_pdgs(a),leg_pdgs(b)
c$$$         stop
c$$$      endif
c$$$c
c$$$      do i=1,n
c$$$         if(i.eq.a)cycle
c$$$         mapped_flavours(i)=leg_pdgs(i)
c$$$         if(i.lt.a)then
c$$$            mapped_labels(i)=i
c$$$         else
c$$$            mapped_labels(i)=i-1
c$$$         endif
c$$$       enddo
c$$$c TODO: think if a -> min(a,b), b -> max(a,b) or similar??
c$$$      if(leg_pdgs(a)+leg_pdgs(b).eq.0)mapped_flavours(b)=21
c$$$      if(leg_pdgs(b).eq.21)mapped_flavours(b)=leg_pdgs(a)
c$$$
c$$$
c$$$
c$$$      return
c$$$      end
c$$$
c$$$
c$$$
c$$$      subroutine reshuffle_momenta(n,leg_pdgs,mapped_flavours,mapped_labels,xpb)
c$$$      implicit none
c$$$      include 'nexternal.inc'
c$$$      integer i,j,n
c$$$      integer leg_pdgs(n-1), mapped_labels(nexternal),mapped_flavours(nexternal)
c$$$      double precision xpb(0:3,n-1), xpb_mapped(0:3,n-1)
c$$$      integer aux_labels(nexternal)
c$$$
c$$$      xpb_mapped(:,:) = 0d0
c$$$      aux_labels(:) = 0
c$$$
c$$$      do i=1,n-1
c$$$         do j=1,nexternal
c$$$            if(leg_pdgs(i).eq.mapped_flavours(j)) then
c$$$               if(mapped_flavours(j).eq.0.or.aux_labels(j).ne.0) cycle
c$$$               xpb_mapped(:,mapped_labels(j)) = xpb(:,i)
c$$$               aux_labels(j) = i
c$$$               exit
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$
c$$$      xpb(:,:) = xpb_mapped(:,:)
c$$$      mapped_labels(:) = aux_labels(:)
c$$$      
c$$$      return
c$$$      end


      integer function nspincol(pdg)
c     number of spin and colour states of a massless QCD parton, over
c     which the matrix elements are averaged for an initial-state one
      implicit none
      integer pdg
      if(pdg.eq.21)then
         nspincol=2*8
      elseif(abs(pdg).ge.1.and.abs(pdg).le.6)then
         nspincol=2*3
      else
         write(*,*)'nspincol: not a QCD parton',pdg
         stop 1
      endif
      return
      end
