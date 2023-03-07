      subroutine get_soft_mapped_labels(a,b,c,n,leg_pdgs,mapped_labels,
     $           mapped_flavours,isLOQCDparton)
      implicit none
      integer a,b,c,n
      integer leg_pdgs(n)
      integer mapped_labels(n),mapped_flavours(n)
      integer i,j
      logical isLOQCDparton(n-1)
c
c     initialise
      j = 0
      mapped_labels = 0
      mapped_flavours = 0
      isLOQCDparton = .false.
c
c     TODO: consistency check on (a,b,c) PDGs
c
c     check mapping structure
      if(a.le.2)then
         write(*,*) 'The first particle must be in the final state!'
         write(*,*) a,b,c
         write(*,*) 'Exit...'
         call exit(-1)
      endif
c
c     For NLO mapping type
      mapped_flavours = leg_pdgs
      do i=1,n
         if(i.lt.a)then
            mapped_labels(i) = i
         elseif(i.eq.a)then
            mapped_labels(i) = 0
            mapped_flavours(i) = 0
         elseif(i.gt.a)then
            mapped_labels(i) = a + j
            j = j + 1
         endif
c     write isLOQCDparton
c     exclude the mapped_flavours=0 value of the removed gluon
         if(mapped_flavours(i).ne.0) then
            if(abs(mapped_flavours(i)).lt.6.or.mapped_flavours(i).eq.21) then
               isLOQCDparton(i) = .true.
            endif
         endif 
      enddo

      end


      subroutine get_collinear_mapped_labels(a,b,c,n,leg_pdgs,
     $           mapped_labels,mapped_flavours)
      implicit none
      integer a,b,c,n
      integer leg_pdgs(n)
      integer min_ab, max_ab
      integer mapped_labels(n),mapped_flavours(n)
      integer i,j
c
c     initialise
      min_ab = 0
      max_ab = 0
      j = 0
      mapped_labels = 0
      mapped_flavours = 0
c
c     TODO: consistency check on (a,b,c) PDGs
c
c     check mapping structure
      if(a.le.2)then
         write(*,*) 'The first particle must be in the final state!'
         write(*,*) a,b,c
         write(*,*) 'Exit...'
         call exit(-1)
      endif
c
c     For NLO mapping type
      min_ab = MIN(a,b)
      max_ab = MAX(a,b)
      mapped_flavours = leg_pdgs
c
c     FaFb mapping : min_ab > 2,  max_ab > 2
c
      if(min_ab.gt.2)then
         do i=1,n
            if(i.lt.min_ab)then
               mapped_labels(i) = i
            elseif(i.eq.min_ab)then
               mapped_labels(i) = 0
            elseif(i.eq.max_ab)then
               mapped_labels(i) = n-1
            elseif((i.gt.min_ab).and.(i.ne.max_ab))then
               mapped_labels(i) = min_ab + j
               j = j + 1
            endif
         enddo
         mapped_flavours(min_ab) = 0
c        Notation: given (abc), [max_ab] > a + b
c        q(bq) > q(bq) + g
         if(leg_PDGs(a).ne.21.and.leg_PDGs(b).eq.21)then
            mapped_flavours(max_ab) = leg_PDGs(a)
c        q(bq) > g + q(bq)
         elseif(leg_PDGs(a).eq.21.and.leg_PDGs(b).ne.21)then
            mapped_flavours(max_ab) = leg_PDGs(b)
c        g > q(bq) bq(q)
         elseif(leg_PDGs(a).eq.(-leg_PDGs(b)))then
            mapped_flavours(max_ab) = 21
c        g > g + g
         elseif(leg_PDGs(a).eq.21.and.leg_PDGs(b).eq.21)then
            mapped_flavours(max_ab) = 21
         endif
c
c     FaIb mapping : min_ab <= 2, max_ab> 2
c
      elseif(min_ab.le.2)then
         do i=1,n
            if(i.lt.max_ab)then
               mapped_labels(i) = i
            elseif(i.eq.max_ab)then
               mapped_labels(i) = 0
            elseif(i.gt.max_ab)then
               mapped_labels(i) = max_ab + j
               j = j + 1
            endif
         enddo
         mapped_flavours(max_ab) = 0
c        Notation: given (a,b,c), b > [min_ab] + a
c        q(bq) > g + q(bq), g > g + g 
         if(leg_PDGs(a).eq.leg_PDGs(b))then
            mapped_flavours(min_ab) = 21
c        g > q(bq) bq(q)
         elseif(leg_PDGs(a).ne.leg_PDGs(b).and.leg_PDGs(b).eq.21)then
            mapped_flavours(min_ab) = - leg_PDGs(a)
c        q(bq) > q(bq) + g
         elseif(leg_PDGs(a).ne.leg_PDGs(b).and.leg_PDGs(a).eq.21)then
            mapped_flavours(min_ab) = leg_PDGs(b)
         endif
      endif
      end
