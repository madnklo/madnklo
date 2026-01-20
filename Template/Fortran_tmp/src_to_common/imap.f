      subroutine get_soft_mapped_labels(npart,a,leg_pdgs,underlying_pdgs
     $     ,mapped_labels,mapped_flavours,mapped_indices_shuff)
c     assigns labels and flavours of particles after a mapping (a,x,y)
c     that removes gluon a from an npart-body final state
c     inputs : npart,a,leg_pdgs,underlying_pdgs
c     outputs: mapped_labels,mapped_flavours,is_mapped_qcd_parton,
c              mapped_indices_shuff
      implicit none
      integer npart,a,i,ii,jj
      integer leg_pdgs(npart),mapped_labels(npart)
      integer mapped_flavours(npart-1),underlying_pdgs(npart-1)
     $     ,mapped_indices_shuff(npart-1)
      logical looked(npart-1)
c
c     preliminary checks
      if(a.lt.3.or.leg_pdgs(a).ne.21)then
         write(*,*)'wrong a in get_soft_mapped_labels',a,leg_pdgs(a)
         stop
      endif
c
c     initial assignment of mapped_labels and mapped_flavours:
c     mapped_labels has npart entries, one per particle to be mapped
c     (including a), while mapped_flavours is the flavour of the
c     npart-1 underlying partons, hence it has npart-1 entries
c
c     example, npart=6:
c     i                            1  2  3  (4)  5  6
c     leg_pdgs(i)                  u  ub g   g   d  db
c     mapped_labels(i)             1  2  3   0   4  5
c     ii                           1  2  3       4  5
c     mapped_flavours(ii)          u  ub g       d  db
c
      ii=0
      mapped_labels=0
      do i=1,npart
         if(i.eq.a)cycle
         ii=ii+1
         mapped_labels(i)=ii
         mapped_flavours(ii)=leg_pdgs(i)
      enddo
c
c     now the underlying indices ii runs in 1,..,npart-1. Reshuffle them
c     to mapped_indices_shuff to align them to the conventions of the
c     underlying_pdgs
c
c     example:
c     ii                           1  2  3       4  5
c     mapped_flavours(ii)          u  ub g       d  db
c     underlying_pdgs(ii)          g  u  d       ub db
c     mapped_indices_shuff(ii)     2  4  1       3  5    
c
      looked=.false.
      do ii=1,npart-1
         do jj=1,npart-1
            if(mapped_flavours(ii).eq.underlying_pdgs(jj)
     &         .and..not.looked(jj))then
               looked(jj)=.true.
               mapped_indices_shuff(ii)=jj
            endif
         enddo
      enddo
c
      return
      end


      subroutine get_collinear_mapped_labels(npart,a,b,leg_pdgs,underlying_pdgs
     $     ,mapped_labels,mapped_flavours,mapped_indices_shuff)
c     assigns labels and flavours of particles after a mapping (a,b,y)
c     that clusters partons (a,b) in an npart-body final state
c     inputs : npart,a,b,leg_pdgs,underlying_pdgs
c     outputs: mapped_labels,mapped_flavours,mapped_indices_shuff
      implicit none
      integer npart,a,b,i,ii,jj,i1,i2
      integer leg_pdgs(npart),mapped_labels(npart)
      integer mapped_flavours(npart-1),underlying_pdgs(npart-1)
     $     ,mapped_indices_shuff(npart-1)
      logical isgg,isqq,isqg,isgq,looked(npart-1)
c
c     preliminary checks
      if(a.lt.3.or.b.lt.3)then
         write(*,*)'wrong a, b in get_collinear_mapped_labels',a,b
         stop
      endif
c
c     possible QCD pairs
      i1=min(a,b)
      i2=max(a,b)
      isqq=abs(leg_pdgs(i1)).le.6.and.leg_pdgs(i1)+leg_pdgs(i2).eq.0
      isgg=leg_pdgs(i1).eq.21.and.leg_pdgs(i2).eq.21
      isqg=abs(leg_pdgs(i1)).le.6.and.leg_pdgs(i2).eq.21
      isgq=leg_pdgs(i1).eq.21.and.abs(leg_pdgs(i2)).le.6
      if (.not.(isgg.or.isqq.or.isqg.or.isgq)) then
         write(*,*)'inconsistent pair in get_collinear_mapped_labels'
         write(*,*)i1,i2,a,b,leg_pdgs(a),leg_pdgs(b)
         stop
      endif
c
c     initial assignment of mapped_labels and mapped_flavours:
c     mapped_labels has npart entries, one per particle to be mapped
c     (including a), while mapped_flavours is the flavour of the
c     npart-1 underlying partons, hence it has npart-1 entries
c
c     example, npart=6:
c     i                            1  (2)  3  (4)  5  6
c     leg_pdgs(i)                  u  ub   g   g   d  db
c     mapped_labels(i)             1  3    2   3   4  5
c     ii                           1       2   3   4  5
c     mapped_flavours(ii)          u       g   ub  d  db
c
      ii=0
      mapped_labels=0
      do i=1,npart
         if(i.eq.i1)cycle
         ii=ii+1
         mapped_labels(i)=ii
         if(i.ne.i2)then
            mapped_flavours(ii)=leg_pdgs(i)
         else
            if(isqq.or.isgg)then
               mapped_flavours(ii)=21
            elseif(isqg)then
               mapped_flavours(ii)=leg_pdgs(i1)
            elseif(isgq)then
               mapped_flavours(ii)=leg_pdgs(i2)
            endif
         endif
      enddo
      mapped_labels(i1)=mapped_labels(i2)
c
c     now the underlying indices ii runs in 1,..,npart-1. Reshuffle them
c     to mapped_indices_shuff to align them to the conventions of the
c     underlying_pdgs
c
c     example:
c     ii                           1       2   3   4  5
c     mapped_flavours(ii)          u       g   ub  d  db
c     underlying_pdgs(ii)          g       u   d   ub db
c     mapped_indices_shuff(ii)     2       1   4   3  5    
c
      looked=.false.
      do ii=1,npart-1
         do jj=1,npart-1
            if(mapped_flavours(ii).eq.underlying_pdgs(jj)
     &         .and..not.looked(jj))then
               looked(jj)=.true.
               mapped_indices_shuff(ii)=jj
            endif
         enddo
      enddo
c
      return
      end

