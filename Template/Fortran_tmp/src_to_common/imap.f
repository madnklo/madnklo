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





      subroutine get_soft_mapped_labels(a,n,leg_pdgs,mapped_labels,mapped_flavours,ismappedQCDparton)
c     assigns labels and flavours of particles after a mapping (a,x,y) that removes gluon a from
c     an n-body final state
      implicit none
      integer a,n,i
      integer leg_pdgs(n),mapped_labels(n),mapped_flavours(n)
      logical ismappedQCDparton(n-1)
c
c     initialise
      mapped_labels=0
      mapped_flavours=0
      ismappedQCDparton=.false.
c
c     preliminary checks
      if(a.lt.3)then
         write(*,*)'get_soft_mapped_labels: wrong parton a ',a
         stop
      endif
      if(leg_pdgs(a).ne.21)then
         write(*,*)'get_soft_mapped_labels: a is not a gluon',a,leg_pdgs(a)
         stop
      endif
c
c     assign mapped labels, flavours, ismappedQCDparton
      do i=1,n
         if(i.eq.a)cycle
         mapped_flavours(i)=leg_pdgs(i)
         if(i.lt.a)then
            mapped_labels(i)=i
         else
            mapped_labels(i)=i-1
         endif
         if(abs(mapped_flavours(i)).le.6.or.mapped_flavours(i).eq.21)
     &        ismappedQCDparton(mapped_labels(i)) = .true.
      enddo
c
      return
      end



      subroutine get_collinear_mapped_labels(a,b,n,leg_pdgs,mapped_labels,mapped_flavours)
c     assigns labels and flavours of particles after a mapping (a,b,y) that clusters partons
c     (a,b) in an n-body final state
      implicit none
      integer a,b,n,i
      integer leg_pdgs(n),mapped_labels(n),mapped_flavours(n)
      logical isgluon,isqqbar,isQCD
c
c     initialise
      mapped_labels=0
      mapped_flavours=0
c
c     preliminary checks
      if(a.lt.3.or.b.lt.3)then
         write(*,*)'get_collinear_mapped_labels: wrong partons a, b ',a,b
         stop
      endif
      isgluon=leg_pdgs(a).eq.21.or.leg_pdgs(b).eq.21
      isqqbar=leg_pdgs(a)+leg_pdgs(b).eq.0
      isQCD=(abs(leg_pdgs(a)).le.6.or.leg_pdgs(a).eq.21).and.
     &      (abs(leg_pdgs(b)).le.6.or.leg_pdgs(b).eq.21)
      if (.not.(isgluon.or.isqqbar.or.isQCD)) then
         write(*,*)'get_collinear_mapped_labels: inconsistent a, b '
         write(*,*)leg_pdgs(a),leg_pdgs(b)
         stop
      endif
c
      do i=1,n
         if(i.eq.a)cycle
         mapped_flavours(i)=leg_pdgs(i)
         if(i.lt.a)then
            mapped_labels(i)=i
         else
            mapped_labels(i)=i-1
         endif
       enddo
c TODO: think if a -> min(a,b), b -> max(a,b) or similar??
      if(leg_pdgs(a)+leg_pdgs(b).eq.0)mapped_flavours(b)=21
      if(leg_pdgs(b).eq.21)mapped_flavours(b)=leg_pdgs(a)



      return
      end



      subroutine reshuffle_momenta(n,leg_pdgs,mapped_flavours,mapped_labels,xpb)
      implicit none
      include 'nexternal.inc'
      integer i,j,n
      integer leg_pdgs(n-1), mapped_labels(nexternal),mapped_flavours(nexternal)
      double precision xpb(0:3,n-1), xpb_mapped(0:3,n-1)
      integer aux_labels(nexternal)

      xpb_mapped(:,:) = 0d0
      aux_labels(:) = 0

      do i=1,n-1
         do j=1,nexternal
            if(leg_pdgs(i).eq.mapped_flavours(j)) then
               if(mapped_flavours(j).eq.0.or.aux_labels(j).ne.0) cycle
               xpb_mapped(:,mapped_labels(j)) = xpb(:,i)
               aux_labels(j) = i
               exit
            endif
         enddo
      enddo

      xpb(:,:) = xpb_mapped(:,:)
      mapped_labels(:) = aux_labels(:)
      
      return
      end
