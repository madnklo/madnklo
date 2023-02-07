      subroutine get_soft_mapped_labels(a,b,c,n,mapped_labels)
      integer a,b,c,n
      integer mapped_labels(n)
      integer i,j
c
c     initialise
      j = 0
c
c     check on (a,b,c) PDGs
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
      do i=1,n
         if(i.lt.a)then
            mapped_labels(i) = i
         elseif(i.eq.a)then
            mapped_labels(i) = 0
         elseif(i.gt.a)then
            mapped_labels(i) = a + j
            j = j + 1
         endif
      enddo
      end


      subroutine get_collinear_mapped_labels(a,b,c,n,mapped_labels)
      integer a,b,c,n
      integer min_ab, max_ab
      integer mapped_labels(n)
      integer i,j
c
c     initialise
      min_ab = 0
      max_ab = 0
      j = 0
c
c     check on (a,b,c) PDGs
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
      endif
c
c     FaIb mapping : min_ab <= 2, max_ab> 2
c
      j = 0
      if(min_ab.le.2)then
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
      endif

      subroutine get_mapped_labels(a,b,c,n,mapped_labels)
      integer a,b,c,n
      integer min_ab
      integer mapped_labels(n)
      integer i,j
c
c     initialise
      min_ab = 0
      j = 0
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
      do i=1,n
         if(i.lt.min_ab)then
            mapped_labels(i) = i
         elseif(i.eq.a)then
            mapped_labels(i) = min_ab
         elseif(i.eq.b)then
            mapped_labels(i) = min_ab
         else
            j = j + 1
            mapped_labels(i) = min_ab + j
         endif
      enddo
      end
