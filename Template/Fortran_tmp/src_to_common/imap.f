      subroutine get_mapped_labels(npart,a,b,pdgs,underlying_pdgs,mapped_labels)
c     assigns labels and flavours of particles after a mapping (a,b,y)
c     that clusters partons (a,b) in an npart-body final state
c     inputs : npart,a,b,pdgs,underlying_pdgs
c     outputs: mapped_labels,mapped_flavours,mapped_indices_shuff
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
      do ii=3,npart-1
         if(iref(1,ii-2).eq.mapped_labels(i))iref_v=iref(2,ii-2)
      enddo
      do ii=3,npart
         if(mapped_labels(ii).eq.iref_v)iref_r=ii
      enddo
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
