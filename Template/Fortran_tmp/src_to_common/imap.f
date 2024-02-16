      subroutine get_mapped_labels(maptype,a,b,c,n,leg_pdgs,
     $           mapped_labels,mapped_flavours,isLOQCDparton)
      implicit none
      integer a,b,c,n,isec,jsec
      integer leg_pdgs(n)
      integer mapped_labels(n),mapped_flavours(n)
      integer i,j
      logical isLOQCDparton(n-1)
      character*1 maptype
      
      if(maptype.eq.'S')then
         call get_soft_mapped_labels(a,b,c,n,leg_pdgs,mapped_labels,
     $           mapped_flavours,isLOQCDparton)
      elseif(maptype.eq.'C')then
         call get_collinear_mapped_labels(a,b,c,n,leg_pdgs,
     $           mapped_labels,mapped_flavours)
      else
         write(*,*) 'get_mapped_labels: '
         write(*,*) 'Invalid maptype, must be S or C!', maptype
         stop
      endif

      return
      end


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
         write(*,*) 'get_soft_mapped_labels: '
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
            if(abs(mapped_flavours(i)).le.6.or.mapped_flavours(i).eq.21) then
               isLOQCDparton(mapped_labels(i)) = .true.
            endif
         endif 
      enddo

      end


      subroutine get_collinear_mapped_labels(a,b,c,n,leg_pdgs,
     $           mapped_labels,mapped_flavours)
      implicit none
      include 'virtual_recoilers.inc'
      integer i,j
      integer a,b,c,n
      integer isec,jsec
      common/cnlosecindices/isec,jsec
      integer rm_leg,parent_leg,rec_leg
      integer leg_pdgs(n)
      integer mapped_labels(n),mapped_flavours(n)
      integer Born_leg_PDGs(n-1)
      integer parent_from_v,rec_from_v
c
c     initialise
      j = 0
      rm_leg = 0
      parent_leg = 0
      mapped_labels = 0
      mapped_flavours = 0

c     check mapping structure
      if(a.le.2)then
         write(*,*) 'get_collinear_mapped_labels: '
         write(*,*) 'The first particle must be in the final state!'
         write(*,*) a,b,c
         write(*,*) 'Exit...'
         call exit(-1)
      endif

c     Associate a,b,c to the collinear particles in the sector
      if(a.eq.isec) then
         rm_leg = a
         if(b.eq.jsec) then
            parent_leg = b
            rec_leg = c
         elseif(c.eq.jsec) then
            parent_leg = c
            rec_leg = b
         endif
      elseif(a.eq.jsec) then
         rm_leg = a
         if(b.eq.isec) then
            parent_leg = b
            rec_leg = c
         elseif(c.eq.isec) then
            parent_leg = c
            rec_leg = b
         endif
      elseif((rm_leg.eq.0).or.(parent_leg.eq.0).or.(rm_leg.eq.parent_leg)) then
         write(*,*) 'get_collinear_mapped_labels: '
         write(*,*) 'Mapping tuple (abc) does not include'
         write(*,*) 'particles defining the singular sector (ij) !'
         write(*,*) 'a, b, c = ', a, b, c
         write(*,*) 'i, j = ', isec, jsec
         stop
      endif
c
c     TODO: consistency check on (a,b,c) PDGs
c
      mapped_flavours = leg_pdgs
c
c     FaFb mapping : isec is always > 2,  jsec > 2
c     Notation: given (abc), [ab] > a + b
      
      if(jsec.gt.2)then
c        q(barq) > q(barq) + g
         if(leg_PDGs(rm_leg).ne.21.and.leg_PDGs(parent_leg).eq.21)then
            mapped_flavours(parent_leg) = leg_PDGs(rm_leg)
c        q(barq) > g + q(barq)
         elseif(leg_PDGs(rm_leg).eq.21.and.leg_PDGs(parent_leg).ne.21)then
            mapped_flavours(parent_leg) = leg_PDGs(parent_leg)
c        g > q(barq) barq(q)
         elseif(leg_PDGs(rm_leg).eq.(-leg_PDGs(parent_leg)))then
            mapped_flavours(parent_leg) = 21
c        g > g + g
         elseif(leg_PDGs(rm_leg).eq.21.and.leg_PDGs(parent_leg).eq.21)then
            mapped_flavours(parent_leg) = 21
         endif
c        remove the first particle in the mapping
         mapped_flavours(rm_leg) = 0

         call get_Born_PDGs(isec,jsec,n-1,Born_leg_PDGs)

         
c     Identify mapped_labels
c     TODO: check for more involved cases
c$$$         do i=1,n-1
c$$$            do j=1,n
c$$$               if(Born_leg_PDGs(i).eq.mapped_flavours(j)) then
c$$$                  if(mapped_labels(j).eq.0) mapped_labels(j) = i
c$$$                  exit
c$$$               endif
c$$$            enddo
c$$$  enddo
         do i=1,n-1
            do j=1,n
c               write(*,*) 'i,j', i, j
               if(Born_leg_PDGs(i).eq.mapped_flavours(j)) then
                  if(mapped_labels(j).eq.0) then
                     mapped_labels(j) = i
c                     write(*,*) 'mapped_labels', mapped_labels
                     exit
                  endif
               endif
            enddo
         enddo
               
c
c     FaIb mapping : isec is always > 2, jsec < 2
c
      elseif(jsec.le.2)then
c        Notation: given (a,b,c), b > [ab] + a
c        q(barq) > g + q(barq), g > g + g 
         if(leg_PDGs(rm_leg).eq.leg_PDGs(parent_leg))then
            mapped_flavours(parent_leg) = 21
c        g > q(barq) barq(q)
         elseif(leg_PDGs(rm_leg).ne.leg_PDGs(parent_leg).and.leg_PDGs(parent_leg).eq.21)then
c           check the correct order
            mapped_flavours(parent_leg) = - leg_PDGs(rm_leg)
c        q(barq) > q(barq) + g
         elseif(leg_PDGs(rm_leg).ne.leg_PDGs(parent_leg).and.leg_PDGs(rm_leg).eq.21)then
            mapped_flavours(parent_leg) = leg_PDGs(parent_leg)
         endif
         mapped_flavours(rm_leg) = 0
         
         call get_Born_PDGs(isec,jsec,n-1,Born_leg_PDGs)
         
c        rescaling of mapped_labels
         j = 1
         do i=1,n-1
            if(mapped_flavours(j).eq.0) then
               j = j + 1
            endif
            if(Born_leg_PDGs(i).eq.mapped_flavours(j)) then
               mapped_labels(j) = i
               j = j + 1
            endif
         enddo
      endif


c     Check that the recoiler mapped_flavour matches the recoiler particle
c     chosen for the virtual process (reported in virtual_recoilers.inc)
      do i=1,LEN_IREF
         parent_from_v = IREF(1,i)
         rec_from_v = IREF(2,i)
         if((mapped_labels(parent_leg).eq.parent_from_v)
     &        .and.(mapped_labels(rec_leg).eq.rec_from_v)) then
            if(mapped_flavours(rec_leg).ne.Born_leg_PDGs(rec_from_v))
     &           then
               write(*,*) 'get_collinear_mapped_labels: '
               write(*,*) 'Recoiler flavour from mapped_flavours'
               write(*,*) 'and recoiler flavour from virtual process does not match!'
               write(*,*) mapped_flavours(rec_leg), Born_leg_PDGs(rec_from_v)
               stop
            endif
         endif
      enddo
      
      end
