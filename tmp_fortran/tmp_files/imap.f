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
