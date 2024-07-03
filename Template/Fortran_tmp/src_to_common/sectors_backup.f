      subroutine get_Z_NLO(xs,sCM,alpha,i1,i2,Z_NLO,ierr)
c     NLO sector functions Z
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
      double precision wgt,Z_NLO,ek,el,wkl
      double precision num,sigma,sCM,sigma_kl
c
c     initialise
      Z_NLO = 0d0
      ierr = 0
      wkl = 0d0
      ek = 0d0
      el = 0d0
      num = 0d0
      sigma = 0d0
c
c     safety checks
      if(sCM.le.0d0)then
         write(77,*)'Wrong sCM in Z_NLO',sCM
         stop
      endif
      if(alpha.lt.1d0)then
         write(77,*)'Wrong alpha in Z_NLO',alpha
         stop
      endif
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
      do i=1,lensectors
         k=all_sector_list(1,i)
         l=all_sector_list(2,i)
         if(l.eq.k)then
            write(*,*)'Wrong indices in Z_NLO',k,l
            stop
         endif
         if((xs(k,1)+xs(k,2))*(xs(l,1)+xs(l,2))*xs(k,l).ne.0d0)then
            ek=(xs(k,1)+xs(k,2))/sCM
            el=(xs(l,1)+xs(l,2))/sCM
            wkl=sCM*xs(k,l)/(xs(k,1)+xs(k,2))/(xs(l,1)+xs(l,2))
         else
            goto 999
         endif
c     here sigma_kl means sigma_kl + sigma_lk
         sigma_kl=(1d0/ek/wkl)**alpha+(1d0/el/wkl)**alpha
         sigma = sigma + sigma_kl
         if( (k.eq.i1.and.l.eq.i2) .or.
     &       (l.eq.i1.and.k.eq.i2) ) num = num + sigma_kl
      enddo
      if(sigma.le.0d0)then
         write(*,*)'Wrong sigma in Z_NLO',sigma
         stop
      endif
      Z_NLO = num/sigma
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




      subroutine get_ZS_NLO(xs,sCM,alpha,i1,i2,ZS_NLO,ierr)
c     NLO soft sector functions ZS
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
      double precision wgt,ZS_NLO,ek,el,wkl
      double precision num,sigma,sCM,sigma_kl
c
c     initialise
      ZS_NLO = 0d0
      ierr = 0
      wkl = 0d0
      ek = 0d0
      el = 0d0
      num = 0d0
      sigma = 0d0
c
c     safety checks
      if(sCM.le.0d0)then
         write(77,*)'Wrong sCM in ZS_NLO',sCM
         stop
      endif
      if(alpha.lt.1d0)then
         write(77,*)'Wrong alpha in ZS_NLO',alpha
         stop
      endif
      if(i1.le.2) then
         write(*,*) 'First sector index must be in final state',i1
         stop
      endif
      if(abs(sCM-xs(1,2))/sCM.gt.1d-8)then
         write(77,*)'Wrong invariants in ZS_NLO',sCM,xs(1,2)
         goto 999
      endif
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
         if((xs(k,1)+xs(k,2))*(xs(l,1)+xs(l,2))*xs(k,l).ne.0d0)then
            wkl=sCM*xs(k,l)/(xs(k,1)+xs(k,2))/(xs(l,1)+xs(l,2))
         else
            goto 999
         endif
c     here sigma_kl means 1/w_kl
         sigma_kl=(1d0/wkl)**alpha
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
      ZS_NLO = num/sigma
c
c     sanity check
      if(abs(ZS_NLO).ge.huge(1d0).or.isnan(ZS_NLO))then
        write(77,*)'Exception caught in ZS_NLO',ZS_NLO
        goto 999
      endif
c
      return
 999  ierr=1
      return
      end















      subroutine get_Z_NNLO(xs,sCM,alpha,i1,i2,i3,Z_NNLO,sector_type,ierr)
c     NNLO sector functions Z with three indices
c     i1, i2, and i3 are the three sector indices.
c     In any case i1 must be in the final state
      implicit none
      include 'nexternal.inc'
      include 'all_sector_list.inc'
      include 'leg_PDGs.inc'
      integer i,k,l,m
      integer i1,i2,i3,ierr
      double precision alpha
      double precision xs(nexternal,nexternal)
      double precision wgt,Z_NNLO,ek,el,em,wkl,wkm,wlm
      character*1 sector_type
      double precision num,sigma,sCM,sigma_klm
c
c     initialise
      Z_NNLO = 0d0
      ierr = 0
c
c     safety checks
      if(i1.le.2) then
         write(*,*) 'First sector index must be in final state',i1
         stop
      endif
      if(abs(sCM-xs(1,2))/sCM.gt.1d-8)then
         write(77,*)'Wrong invariants in Z_NNLO',sCM,xs(1,2)
         goto 999
      endif
      if(alpha.lt.1d0)then
         write(77,*)'Wrong alpha in Z_NNLO',alpha
         stop
      endif
c
c     build Z_NNLO
      if(sector_type.eq.'F') then
         num = 0d0
         sigma = 0d0
         do i=1,lensectors
            k=all_sector_list(1,i)
            l=all_sector_list(2,i)
            m=all_sector_list(3,i)
            if(l.eq.k.or.l.eq.m.or.k.eq.m)then
               write(*,*)'Wrong indices in Z_NNLO',k,l,m
               stop
            endif
            ek=(xs(k,1)+xs(k,2))/sCM
            el=(xs(l,1)+xs(l,2))/sCM
            em=(xs(m,1)+xs(m,2))/sCM
            wkl=sCM*xs(k,l)/(xs(k,1)+xs(k,2))/(xs(l,1)+xs(l,2))
            wkm=sCM*xs(k,m)/(xs(k,1)+xs(k,2))/(xs(m,1)+xs(m,2))
            wlm=sCM*xs(l,m)/(xs(l,1)+xs(l,2))/(xs(m,1)+xs(m,2))
c     here sigma_klm is the sum of 12 permutations
            sigma_klm=(1d0/ek/wkl)**alpha*(1d0/(ek+el)+1d0/em)*1d0/wlm
     &               +(1d0/el/wkl)**alpha*(1d0/(ek+el)+1d0/em)*1d0/wkm
     &               +(1d0/em/wkm)**alpha*(1d0/(ek+em)+1d0/el)*1d0/wkl
     &               +(1d0/ek/wkm)**alpha*(1d0/(ek+em)+1d0/el)*1d0/wlm
     &               +(1d0/el/wlm)**alpha*(1d0/(el+em)+1d0/ek)*1d0/wkm
     &               +(1d0/em/wlm)**alpha*(1d0/(el+em)+1d0/ek)*1d0/wkl
            sigma = sigma + sigma_klm
            if( (k.eq.i1.and.l.eq.i2.and.m.eq.i3) .or.
     &          (k.eq.i1.and.m.eq.i2.and.l.eq.i3) .or.
     &          (l.eq.i1.and.k.eq.i2.and.m.eq.i3) .or.
     &          (l.eq.i1.and.m.eq.i2.and.k.eq.i3) .or.
     &          (m.eq.i1.and.k.eq.i2.and.l.eq.i3) .or.
     &          (m.eq.i1.and.l.eq.i2.and.k.eq.i3) ) num = num + sigma_klm
         enddo
         if(sigma.le.0d0)then
            write(*,*)'Wrong sigma in Z_NNLO',sigma
            stop
         endif
         Z_NNLO = num/sigma
c$$$c
c$$$c     build ZSS_NNLO
c$$$      elseif(sector_type.eq.'S') then
c$$$         num = 0d0
c$$$         sigma = 0d0
c$$$c$$$c check
c$$$c$$$c 3=g, 4=g, 5=q
c$$$c$$$c sector 34 -> S3.Z34 = (1/w34) / (1/w34 + 1/w35) (this function is called with i1=3,i2=4)
c$$$c$$$c sector 34 -> S4.Z34 = (1/w34) / (1/w34 + 1/w45) (this function is called with i1=4,i2=3)
c$$$c$$$c sector 35 -> S3.Z35 = (1/w35) / (1/w34 + 1/w35) (this function is called with i1=3,i2=5)
c$$$c$$$c sector 45 -> S4.Z45 = (1/w45) / (1/w34 + 1/w45) (this function is called with i1=4,i2=5)
c$$$c$$$c kl = 34 35 45
c$$$         do i=1,lensectors
c$$$            k=all_sector_list(1,i)
c$$$            l=all_sector_list(2,i)
c$$$            if(k.ne.i1.and.k.ne.i2.and.l.ne.i1.and.l.ne.i2)cycle
c$$$c     here sigma_kl means 1/w_kl
c$$$            sigma_kl=((xs(k,1)+xs(k,2))*(xs(l,1)+xs(l,2))/sCM/xs(k,l))**alpha
c$$$            if(k.eq.i1)then
c$$$               sigma = sigma + sigma_kl
c$$$               if(l.eq.i2)num = num + sigma_kl
c$$$            endif
c$$$            if(l.eq.i1)then
c$$$               sigma = sigma + sigma_kl
c$$$               if(k.eq.i2)num = num + sigma_kl
c$$$            endif
c$$$         enddo
c$$$         if(sigma.le.0d0)then
c$$$            write(*,*)'Wrong sigma in ZS_NLO',sigma
c$$$            stop
c$$$         endif
c$$$         Z_NLO = num/sigma
      else
         write(*,*) 'Sector_type not recognised: ',sector_type
         stop
      endif
c
c     sanity check
      if(abs(Z_NNLO).ge.huge(1d0).or.isnan(Z_NNLO))then
        write(77,*)'Exception caught in Z_NNLO',Z_NNLO
        goto 999
      endif
c
      return
 999  ierr=1
      return
      end
