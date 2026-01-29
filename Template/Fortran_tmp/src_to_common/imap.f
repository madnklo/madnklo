      subroutine get_mapped_labels(npart,a,b,pdgs,underlying_pdgs,mapped_labels)
c     assigns labels and flavours of particles after a mapping (a,b,y)
c     that clusters partons (a,b) in an npart-body final state
c     inputs : npart,a,b,pdgs,underlying_pdgs
c     outputs: mapped_labels,mapped_flavours,mapped_indices_shuff
      implicit none
      integer npart,a,b,i,jj,i1,i2
      integer pdgs(npart),mapped_labels(npart),target(npart)
      integer underlying_pdgs(npart-1)
      logical isgg,isqq,isqg,isgq,assigned(npart-1)
c
c     initialise
      i1=min(a,b)
      i2=max(a,b)
      mapped_labels(:)=0
      target(:)=pdgs(:)
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
      if(isqq.or.isgg)target(i2)=21
      if(isqg)target(i2)=pdgs(i1)
      if(isgq)target(i2)=pdgs(i2)
c
c     explicit examples of output mapped_labels
c=====================================================
c     example 1, npart = 6:
c     i                        1  (2)  3  (4)  5   6
c     pdgs(i)                  u   ub  g   g   d   db
c     jj                       1       2   3   4   5
c     underlying_pdgs(jj)      g       u   d   ub  db
c     mapped_labels(i)         2   4   1   4   3   5
c=====================================================
c=====================================================
c     example 2, npart = 6:
c     i                        1   2   3  (4)  5  (6)
c     pdgs(i)                  u   ub  g   g   d   db
c     jj                       1   2   3       4   5
c     underlying_pdgs(jj)      g   u   d       ub  db
c     mapped_labels(i)         2   4   1   5   3   5
c=====================================================
c=====================================================
c     example 3, npart = 6:
c     i                        1   2  (3) (4)  5   6
c     pdgs(i)                  u   ub  g   g   d   db
c     jj                       1   2       3   4   5
c     underlying_pdgs(jj)      g   u       d   ub  db
c     mapped_labels(i)         2   4   1   1   3   5      
c=====================================================
c=====================================================
c     example 4, npart = 6:
c     i                        1   2   3   4  (5) (6)
c     pdgs(i)                  u   ub  g   g   d   db
c     jj                       1   2   3   4       5
c     underlying_pdgs(jj)      g   u   g   ub      g
c     mapped_labels(i)         2   4   1   3   5   5      
c=====================================================
c
      assigned=.false.
      do i=1,npart
         if(i.eq.i1)cycle
         do jj=1,npart-1
            if(underlying_pdgs(jj).eq.target(i).and.
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

