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
      if(leg_pdgs(a)+leg_pdgs(b).eq.0)mapped_flavours(b)=21
      if(leg_pdgs(b).eq.21)mapped_flavours(b)=leg_pdgs(a)
c
      return
      end

c$$$
c$$$      if(maptype.eq.'S')then
c$$$         call get_soft_mapped_labels(a,b,c,n,leg_pdgs,mapped_labels,
c$$$     $           mapped_flavours,isLOQCDparton)
c$$$      elseif(maptype.eq.'C')then
c$$$         call get_collinear_mapped_labels(a,b,c,n,leg_pdgs,
c$$$     $           mapped_labels,mapped_flavours)
c$$$      else
c$$$         write(*,*) 'get_mapped_labels: '
c$$$         write(*,*) 'Invalid maptype, must be S or C!', maptype
c$$$         stop
c$$$      endif
c$$$
c$$$      return
c$$$      end
c$$$
c$$$
c$$$      subroutine get_soft_mapped_labels(a,b,c,n,leg_pdgs,mapped_labels,
c$$$     $           mapped_flavours,isLOQCDparton)
c$$$      implicit none
c$$$      integer a,b,c,n
c$$$      integer leg_pdgs(n)
c$$$      integer mapped_labels(n),mapped_flavours(n)
c$$$      integer i,j
c$$$      logical isLOQCDparton(n-1)
c$$$c
c$$$c     initialise
c$$$      j = 0
c$$$      mapped_labels = 0
c$$$      mapped_flavours = 0
c$$$      isLOQCDparton = .false.
c$$$c
c$$$c     TODO: consistency check on (a,b,c) PDGs
c$$$c
c$$$c     check mapping structure
c$$$      if(a.le.2)then
c$$$         write(*,*) 'get_soft_mapped_labels: '
c$$$         write(*,*) 'The first particle must be in the final state!'
c$$$         write(*,*) a,b,c
c$$$         write(*,*) 'Exit...'
c$$$         call exit(-1)
c$$$      endif
c$$$c
c$$$c     For NLO mapping type
c$$$      mapped_flavours = leg_pdgs
c$$$      do i=1,n
c$$$         if(i.lt.a)then
c$$$            mapped_labels(i) = i
c$$$         elseif(i.eq.a)then
c$$$            mapped_labels(i) = 0
c$$$            mapped_flavours(i) = 0
c$$$         elseif(i.gt.a)then
c$$$            mapped_labels(i) = a + j
c$$$            j = j + 1
c$$$         endif
c$$$c     write isLOQCDparton
c$$$c     exclude the mapped_flavours=0 value of the removed gluon
c$$$         if(mapped_flavours(i).ne.0) then
c$$$            if(abs(mapped_flavours(i)).le.6.or.mapped_flavours(i).eq.21) then
c$$$               isLOQCDparton(mapped_labels(i)) = .true.
c$$$            endif
c$$$         endif 
c$$$      enddo
c$$$
c$$$      end
c$$$
c$$$
c$$$      subroutine get_collinear_mapped_labels(a,b,c,n,leg_pdgs,
c$$$     $           mapped_labels,mapped_flavours)
c$$$      implicit none
c$$$c      include 'virtual_recoilers.inc'
c$$$      integer i,j
c$$$      integer a,b,c,n
c$$$      integer isec,jsec
c$$$      integer ksec,lsec
c$$$c     common/cnlosecindices/isec,jsec
c$$$      common/csecindices/isec,jsec,ksec,lsec
c$$$      integer rm_leg,parent_leg,rec_leg
c$$$      integer leg_pdgs(n)
c$$$      integer mapped_labels(n),mapped_flavours(n)
c$$$      integer Born_leg_PDGs(n-1)
c$$$      integer UnderLying_leg_PDGs(n-1)
c$$$      integer parent_from_v,rec_from_v
c$$$c
c$$$c     initialise
c$$$      j = 0
c$$$      rm_leg = 0
c$$$      parent_leg = 0
c$$$      mapped_labels = 0
c$$$      mapped_flavours = 0
c$$$
c$$$
c$$$c     check mapping structure
c$$$      if(a.le.2)then
c$$$         write(*,*) 'get_collinear_mapped_labels: '
c$$$         write(*,*) 'The first particle must be in the final state!'
c$$$         write(*,*) a,b,c
c$$$         write(*,*) 'Exit...'
c$$$         call exit(-1)
c$$$      endif
c$$$
c$$$c     Associate a,b,c to the collinear particles in the sector
c$$$
c$$$c      write(*,*) 'a, b, c   = ', a, b, c
c$$$c      write(*,*) 'isec,jsec = ',isec,jsec
c$$$c     pause
c$$$
c$$$      
c$$$      rm_leg = a
c$$$
c$$$      if(ksec.eq.0 .and. lsec .eq.0) then ! This is NLO mapping
c$$$         parent_leg = b
c$$$         rec_leg = c
c$$$      elseif(ksec.ne.0 .and. lsec .eq. 0) then ! NNLO Mapping for 3 index sector (i,j,k)
c$$$c         if((a.eq.isec.and.b.eq.jsec).or.(a.eq.jsec .and. b.eq.ksec)) then
c$$$         parent_leg = b
c$$$         rec_leg = c
c$$$         ! TODO : implement correct checks on mapping
c$$$c$$$         else
c$$$c$$$            write(*,*) 'imap : NNLO mapping error...'
c$$$c$$$            write(*,*) 'It should be a = isec, b = jsec
c$$$c$$$     $       or a = jsec, b = ksec'
c$$$c$$$            write(*,*) 'a, b, isec, jsec, ksec = ', a, b, isec, jsec, ksec
c$$$c$$$            stop
c$$$c$$$         endif
c$$$      elseif(ksec.ne.0 .and. lsec .ne. 0) then ! NNLO Mapping for 4 index sector (i,j,k,l)
c$$$         write(*,*) 'NNLO mapping for 4 index sector to be implemented'
c$$$         stop
c$$$      endif
c$$$c$$$      elseif((a.eq.isec.and.c.eq.jsec).or.(a.eq.jsec.and.c.eq.isec)) then
c$$$c$$$         parent_leg = c
c$$$c$$$         rec_leg = b
c$$$c$$$      endif
c$$$
c$$$c$$$      
c$$$c$$$      
c$$$c$$$      if(a.eq.isec) then
c$$$c$$$         rm_leg = a
c$$$c$$$         if(b.eq.jsec) then
c$$$c$$$            parent_leg = b
c$$$c$$$            rec_leg = c
c$$$c$$$         elseif(c.eq.jsec) then
c$$$c$$$            parent_leg = c
c$$$c$$$            rec_leg = b
c$$$c$$$         endif
c$$$c$$$      elseif(a.eq.jsec) then
c$$$c$$$         rm_leg = a
c$$$c$$$         if(b.eq.isec) then
c$$$c$$$            parent_leg = b
c$$$c$$$            rec_leg = c
c$$$c$$$         elseif(c.eq.isec) then
c$$$c$$$            parent_leg = c
c$$$c$$$            rec_leg = b
c$$$c$$$         endif
c$$$c$$$      elseif((rm_leg.eq.0).or.(parent_leg.eq.0).or.(rm_leg.eq.parent_leg)) then
c$$$c$$$         write(*,*) 'get_collinear_mapped_labels: '
c$$$c$$$         write(*,*) 'Mapping tuple (abc) does not include'
c$$$c$$$         write(*,*) 'particles defining the singular sector (ij) !'
c$$$c$$$         write(*,*) 'a, b, c = ', a, b, c
c$$$c$$$         write(*,*) 'i, j = ', isec, jsec
c$$$c$$$         stop
c$$$c$$$      endif
c$$$c
c$$$c     TODO: consistency check on (a,b,c) PDGs
c$$$c
c$$$      mapped_flavours = leg_pdgs
c$$$
c$$$
c$$$
c$$$c
c$$$c     FaFb mapping : isec is always > 2,  jsec > 2
c$$$c     Notation: given (abc), [ab] > a + b
c$$$      
c$$$      if(jsec.gt.2)then
c$$$c        q(barq) > q(barq) + g
c$$$         if(leg_PDGs(rm_leg).ne.21.and.leg_PDGs(parent_leg).eq.21)then
c$$$            mapped_flavours(parent_leg) = leg_PDGs(rm_leg)
c$$$c        q(barq) > g + q(barq)
c$$$         elseif(leg_PDGs(rm_leg).eq.21.and.leg_PDGs(parent_leg).ne.21)then
c$$$            mapped_flavours(parent_leg) = leg_PDGs(parent_leg)
c$$$c        g > q(barq) barq(q)
c$$$         elseif(leg_PDGs(rm_leg).eq.(-leg_PDGs(parent_leg)))then
c$$$            mapped_flavours(parent_leg) = 21
c$$$c        g > g + g
c$$$         elseif(leg_PDGs(rm_leg).eq.21.and.leg_PDGs(parent_leg).eq.21)then
c$$$            mapped_flavours(parent_leg) = 21
c$$$         endif
c$$$c        remove the first particle in the mapping
c$$$         mapped_flavours(rm_leg) = 0
c$$$
c$$$c         call get_Born_PDGs(isec,jsec,n-1,Born_leg_PDGs)
c$$$         call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,N-1
c$$$     $        ,UNDERLYING_LEG_PDGS)
c$$$
c$$$
c$$$         
c$$$c     Identify mapped_labels
c$$$c     TODO: check for more involved cases
c$$$c$$$         do i=1,n-1
c$$$c$$$            do j=1,n
c$$$c$$$               if(Born_leg_PDGs(i).eq.mapped_flavours(j)) then
c$$$c$$$                  if(mapped_labels(j).eq.0) mapped_labels(j) = i
c$$$c$$$                  exit
c$$$c$$$               endif
c$$$c$$$            enddo
c$$$c$$$  enddo
c$$$         do i=1,n-1
c$$$            do j=1,n
c$$$c               write(*,*) 'i,j', i, j
c$$$               if(UnderLying_leg_PDGs(i).eq.mapped_flavours(j)) then
c$$$                  if(mapped_labels(j).eq.0) then
c$$$                     mapped_labels(j) = i
c$$$c                     write(*,*) 'mapped_labels', mapped_labels
c$$$                     exit
c$$$                  endif
c$$$               endif
c$$$            enddo
c$$$         enddo
c$$$c$$$         write(*,*) 'a,b,c,n', a,b,c,n
c$$$c$$$         write(*,*) 'mapped_flavours', mapped_flavours
c$$$c$$$         write(*,*) 'mapped_labels', mapped_labels
c$$$c$$$         pause
c$$$c
c$$$c     FaIb mapping : isec is always > 2, jsec < 2
c$$$c
c$$$      elseif(jsec.le.2)then
c$$$c        Notation: given (a,b,c), b > [ab] + a
c$$$c        q(barq) > g + q(barq), g > g + g 
c$$$         if(leg_PDGs(rm_leg).eq.leg_PDGs(parent_leg))then
c$$$            mapped_flavours(parent_leg) = 21
c$$$c        g > q(barq) barq(q)
c$$$         elseif(leg_PDGs(rm_leg).ne.leg_PDGs(parent_leg).and.leg_PDGs(parent_leg).eq.21)then
c$$$c           check the correct order
c$$$            mapped_flavours(parent_leg) = - leg_PDGs(rm_leg)
c$$$c        q(barq) > q(barq) + g
c$$$         elseif(leg_PDGs(rm_leg).ne.leg_PDGs(parent_leg).and.leg_PDGs(rm_leg).eq.21)then
c$$$            mapped_flavours(parent_leg) = leg_PDGs(parent_leg)
c$$$         endif
c$$$         mapped_flavours(rm_leg) = 0
c$$$         
c$$$c         call get_Born_PDGs(isec,jsec,n-1,Born_leg_PDGs)
c$$$         call GET_UNDERLYING_PDGS(ISEC,JSEC,0,0,N-1
c$$$     $ ,UnderLying_LEG_PDGS)
c$$$
c$$$c        rescaling of mapped_labels
c$$$         j = 1
c$$$         do i=1,n-1
c$$$            if(mapped_flavours(j).eq.0) then
c$$$               j = j + 1
c$$$            endif
c$$$            if(UnderLying_leg_PDGs(i).eq.mapped_flavours(j)) then
c$$$               mapped_labels(j) = i
c$$$               j = j + 1
c$$$            endif
c$$$         enddo
c$$$      endif
c$$$
c$$$
c$$$c     Check that the recoiler mapped_flavour matches the recoiler particle
c$$$c     chosen for the virtual process (reported in virtual_recoilers.inc)
c$$$c$$$      do i=1,LEN_IREF
c$$$c$$$         parent_from_v = IREF(1,i)
c$$$c$$$         rec_from_v = IREF(2,i)
c$$$c$$$         if((mapped_labels(parent_leg).eq.parent_from_v)
c$$$c$$$     &        .and.(mapped_labels(rec_leg).eq.rec_from_v)) then
c$$$c$$$            if(mapped_flavours(rec_leg).ne.UnderLying_leg_PDGs(rec_from_v))
c$$$c$$$     &           then
c$$$c$$$               write(*,*) 'get_collinear_mapped_labels: '
c$$$c$$$               write(*,*) 'Recoiler flavour from mapped_flavours'
c$$$c$$$               write(*,*) 'and recoiler flavour from virtual process does not match!'
c$$$c$$$               write(*,*) mapped_flavours(rec_leg), UnderLying_leg_PDGs(rec_from_v)
c$$$c$$$               stop
c$$$c$$$            endif
c$$$c$$$         endif
c$$$c$$$      enddo
c$$$      
c$$$      end
c$$$
c$$$
c$$$
c$$$c$$$      subroutine check_mapping(a,b,c,n)
c$$$c$$$      implicit none
c$$$c$$$      include 'nexternal.inc'
c$$$c$$$      include 'leg_PDGs.inc'
c$$$c$$$      integer isec,jsec,ksec,lsec
c$$$c$$$      common/secindices/isec,jsec,ksec,lsec
c$$$c$$$      integer a, b, c, n
c$$$c$$$      integer mapped_labels
c$$$c$$$
c$$$c$$$      if(ksec.eq.0 .and. lsec.eq.0) then
c$$$c$$$         if(a.ne.isec.or.b.ne.jsec) then
c$$$c$$$            write(*,*) 'imap : NLO mapping error...'
c$$$c$$$            write(*,*) 'It should be a = isec, b = jsec'
c$$$c$$$            write(*,*) 'a, b, isec, jsec = ', a, b, isec, jsec
c$$$c$$$            stop
c$$$c$$$         endif
c$$$c$$$      elseif(ksec.ne.0 .and. lsec .eq. 0) then
c$$$c$$$         call get_collinear_mapped_labels(isec,jsec,iref,nexternal,leg_pdgs,
c$$$c$$$     $     mapped_labels,mapped_flavours)
c$$$c$$$
c$$$c$$$         if((a.eq.isec.and.b.eq.jsec).or.(a.eq.jsec .and. b.eq.ksec))
c$$$c$$$         
c$$$c$$$      endif
c$$$c$$$      
c$$$c$$$      
c$$$c$$$      end
