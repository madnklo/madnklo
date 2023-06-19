      subroutine get_Z_NLO(xs,sCM,alpha,i1,i2,Z_NLO,sector_type,ierr)
c     i1 and i2 are the two sector indices.
c     This function is meant to be called with (i1,i2) = (isec,jsec)
c     or (i1,i2) = (jsec,isec), but in any case i1 must be in the
c     final state 
      implicit none
      include 'nexternal.inc'
      include 'all_sector_list.inc'
      include 'leg_PDGs.inc'
      integer i,k,l
      integer i1,i2,ierr
      double precision alpha
      double precision xs(nexternal,nexternal)
      double precision wgt
      double precision Z_NLO
      character*1 sector_type
      double precision num,sigma,sCM,sigma_kl
c
c     initialise
      Z_NLO = 0d0
      ierr = 0
c
c     safety checks
      if(i1.le.2) then
         write(*,*) 'First sector index must be in final state',i1
         stop
      endif
      if(abs(sCM-xs(1,2))/sCM.gt.1d-8)then
         write(77,*)'Wrong invariants in Z_NLO',sCM,xs(1,2)
         goto 999
      endif
c
c     build Z_NLO
      if(sector_type.eq.'F') then
         num = 0d0
         sigma = 0d0
         do i=1,lensectors
            k=all_sector_list(1,i)
            l=all_sector_list(2,i)
            if(l.eq.k)then
               write(*,*)'Wring indices in Z_NLO',k,l
               stop
            endif
c     here sigma_kl means sigma_kl + sigma_lk
            sigma_kl=(xs(l,1)+xs(l,2))/xs(k,l) +
     &               (xs(k,1)+xs(k,2))/xs(k,l)
            sigma = sigma + sigma_kl
            if( (k.eq.i1.and.l.eq.i2) .or.
     &          (l.eq.i1.and.k.eq.i2) ) num = num + sigma_kl
         enddo
         if(sigma.le.0d0)then
            write(*,*)'Wrong sigma in Z_NLO',sigma
            stop
         endif
         Z_NLO = num/sigma
c
c     build ZS_NLO
      elseif(sector_type.eq.'S') then
         num = 0d0
         sigma = 0d0
c$$$c check
c$$$c 3=g, 4=g, 5=q
c$$$c sector 34 -> S3.Z34 = (1/w34) / (1/w34 + 1/w35) (this function is called with i1=3,i2=4)
c$$$c sector 34 -> S4.Z34 = (1/w34) / (1/w34 + 1/w45) (this function is called with i1=4,i2=3)
c$$$c sector 35 -> S3.Z35 = (1/w35) / (1/w34 + 1/w35) (this function is called with i1=3,i2=5)
c$$$c sector 45 -> S4.Z45 = (1/w45) / (1/w34 + 1/w45) (this function is called with i1=4,i2=5)
c$$$c kl = 34 35 45
         do i=1,lensectors
            k=all_sector_list(1,i)
            l=all_sector_list(2,i)
            if(k.ne.i1.and.k.ne.i2.and.l.ne.i1.and.l.ne.i2)cycle
            sigma_kl=(xs(k,1)+xs(k,2))*(xs(l,1)+xs(l,2))/sCM/xs(k,l) ! here sigma_kl = 1/w_kl
            if(k.eq.i1)then
               sigma = sigma + sigma_kl
               if(l.eq.i2)num = num + sigma_kl
            endif
            if(l.eq.i1)then
               sigma = sigma + sigma_kl
               if(k.eq.i2)num = num + sigma_kl
            endif
         enddo
         if(sigma.le.0d0)then
            write(*,*)'Wrong sigma in ZS_NLO',sigma
            stop
         endif
         Z_NLO = num/sigma
      else
         write(*,*) 'Sector_type not recognised: ',sector_type
         stop
      endif
c
c     sanity check
      if(abs(Z_NLO).ge.huge(1d0).or.isnan(Z_NLO))then
        write(77,*)'Exception caught in Z_NLO',Z_NLO
        goto 999
      endif
c
      return
 999  ierr=1
      return
      end
