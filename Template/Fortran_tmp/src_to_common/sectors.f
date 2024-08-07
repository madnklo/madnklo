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
      integer i,a,b
      integer i1,i2,ierr
      double precision alpha
      double precision xs(nexternal,nexternal)
      double precision wgt,Z_NLO,ea,eb,wab
      double precision num,sigma,sCM,sigma_ab
c
c     initialise
      Z_NLO = 0d0
      ierr = 0
      wab = 0d0
      ea = 0d0
      eb = 0d0
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
         a=all_sector_list(1,i)
         b=all_sector_list(2,i)
         if(a.eq.b)then
            write(*,*)'Wrong indices in Z_NLO',a,b
            stop
         endif
         if((xs(a,1)+xs(a,2))*(xs(b,1)+xs(b,2))*xs(a,b).ne.0d0)then
            ea=(xs(a,1)+xs(a,2))/sCM
            eb=(xs(b,1)+xs(b,2))/sCM
            wab=sCM*xs(a,b)/(xs(a,1)+xs(a,2))/(xs(b,1)+xs(b,2))
         else
            goto 999
         endif
c     symmetrised sigma_ab
         sigma_ab=(1d0/ea/wab)**alpha+(1d0/eb/wab)**alpha
         sigma = sigma + sigma_ab
         if((a.eq.i1.and.b.eq.i2).or.(a.eq.i2.and.b.eq.i1)) num = num + sigma_ab
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
c     This function is meant to be called with (i1,i2) = perm(isec,jsec)
c     i1 must be in the final state as associated with the soft singularity
      implicit none
      include 'nexternal.inc'
      include 'all_sector_list.inc'
      include 'leg_PDGs.inc'
      integer i,a,b
      integer i1,i2,ierr
      double precision alpha
      double precision xs(nexternal,nexternal)
      double precision wgt,ZS_NLO,ea,eb,wab
      double precision num,sigma,sCM
c
c     initialise
      ZS_NLO = 0d0
      ierr = 0
      ea = 0d0
      eb = 0d0
      wab = 0d0
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
         write(*,*) 'First ZS_NLO index must be > 2',i1
         stop
      endif
      if(abs(sCM-xs(1,2))/sCM.gt.1d-8)then
         write(77,*)'Wrong invariants in ZS_NLO',sCM,xs(1,2)
         goto 999
      endif
c
c     build ZS_NLO
      do i=1,lensectors
         a=all_sector_list(1,i)
         b=all_sector_list(2,i)
         if(a.eq.b)then
            write(*,*)'Wrong indices in ZS_NLO',a,b
            stop
         endif
         if(a.ne.i1.and.a.ne.i2.and.b.ne.i1.and.b.ne.i2)cycle
         if((xs(a,1)+xs(a,2))*(xs(b,1)+xs(b,2))*xs(a,b).ne.0d0)then
            ea=(xs(a,1)+xs(a,2))/sCM
            eb=(xs(b,1)+xs(b,2))/sCM
            wab=sCM*xs(a,b)/(xs(a,1)+xs(a,2))/(xs(b,1)+xs(b,2))
         else
            goto 999
         endif
c$$$c check
c$$$c 3=g, 4=g, 5=q
c$$$c sector 34 -> S3.Z34 = (1/w34) / (1/w34 + 1/w35) (this function is called with i1=3,i2=4)
c$$$c sector 34 -> S4.Z34 = (1/w34) / (1/w34 + 1/w45) (this function is called with i1=4,i2=3)
c$$$c sector 35 -> S3.Z35 = (1/w35) / (1/w34 + 1/w35) (this function is called with i1=3,i2=5)
c$$$c sector 45 -> S4.Z45 = (1/w45) / (1/w34 + 1/w45) (this function is called with i1=4,i2=5)
c$$$c ab = 34 35 45
         if(a.eq.i1) sigma = sigma + (1d0/ea/wab)**alpha
         if(b.eq.i1) sigma = sigma + (1d0/eb/wab)**alpha
         if(a.eq.i1.and.b.eq.i2) num = num + (1d0/ea/wab)**alpha
         if(a.eq.i2.and.b.eq.i1) num = num + (1d0/eb/wab)**alpha
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


      subroutine get_Z_NNLO(xs,sCM,alpha,i1,i2,i3,i4,Z_NNLO,ierr)
c     NNLO sector functions Z with three (i4=0) or four (i4!=0) indices
c     i1, i2, i3, and i4 are the three (four) sector indices.
c     In any case i1 must be in the final state
      implicit none
      include 'nexternal.inc'
      include 'all_sector_list.inc'
      include 'leg_PDGs.inc'
      integer i,a,b,c,d
      integer i1,i2,i3,i4,ierr
      double precision alpha
      double precision xs(nexternal,nexternal)
      double precision wgt,Z_NNLO,ea,eb,ec,ed
      double precision wab,wac,wad,wbc,wbd,wcd
      double precision num,sigma,sCM,sigma_abcd
c
c     initialise
      Z_NNLO = 0d0
      ierr = 0
      wab = 0d0
      wac = 0d0
      wad = 0d0
      wbc = 0d0
      wbd = 0d0
      wcd = 0d0
      ea = 0d0
      eb = 0d0
      ec = 0d0
      ed = 0d0
      num = 0d0
      sigma = 0d0
c
c     safety checks
      if(sCM.le.0d0)then
         write(77,*)'Wrong sCM in Z_NNLO',sCM
         stop
      endif
      if(alpha.lt.1d0)then
         write(77,*)'Wrong alpha in Z_NNLO',alpha
         stop
      endif
      if(i1.le.2) then
         write(*,*) 'First sector index must be in final state',i1
         stop
      endif
      if(abs(sCM-xs(1,2))/sCM.gt.1d-8)then
         write(77,*)'Wrong invariants in Z_NNLO',sCM,xs(1,2)
         goto 999
      endif
c
c     build Z_NNLO
      do i=1,lensectors
         a=all_sector_list(1,i)
         b=all_sector_list(2,i)
         c=all_sector_list(3,i)
         d=all_sector_list(4,i)
         if(b.eq.a.or.c.eq.a.or.d.eq.a.or.d.eq.c)then
            write(*,*)'Wrong indices in Z_NNLO',a,b,c,d
            stop
         endif
         if(d.eq.0) then
            if((xs(a,1)+xs(a,2))*(xs(b,1)+xs(b,2))*(xs(c,1)+xs(c,2))*
     &           xs(a,b)*xs(a,c)*xs(b,c).ne.0d0)then
               ea=(xs(a,1)+xs(a,2))/sCM
               eb=(xs(b,1)+xs(b,2))/sCM
               ec=(xs(c,1)+xs(c,2))/sCM
               wab=sCM*xs(a,b)/(xs(a,1)+xs(a,2))/(xs(b,1)+xs(b,2))
               wac=sCM*xs(a,c)/(xs(a,1)+xs(a,2))/(xs(c,1)+xs(c,2))
               wbc=sCM*xs(b,c)/(xs(b,1)+xs(b,2))/(xs(c,1)+xs(c,2))
            else
               goto 999
            endif
         elseif(d.ne.0) then
            if((xs(a,1)+xs(a,2))*(xs(b,1)+xs(b,2))*(xs(c,1)+xs(c,2))*
     &           (xs(d,1)+xs(d,2))*xs(a,b)*xs(a,c)*xs(a,d)*xs(b,c)*
     &           xs(b,d)*xs(c,d).ne.0d0)then
               ea=(xs(a,1)+xs(a,2))/sCM
               eb=(xs(b,1)+xs(b,2))/sCM
               ec=(xs(c,1)+xs(c,2))/sCM
               wab=sCM*xs(a,b)/(xs(a,1)+xs(a,2))/(xs(b,1)+xs(b,2))
               wac=sCM*xs(a,c)/(xs(a,1)+xs(a,2))/(xs(c,1)+xs(c,2))
               wbc=sCM*xs(b,c)/(xs(b,1)+xs(b,2))/(xs(c,1)+xs(c,2))
               ed=(xs(d,1)+xs(d,2))/sCM
               wad=sCM*xs(a,d)/(xs(a,1)+xs(a,2))/(xs(d,1)+xs(d,2))
               wbd=sCM*xs(b,d)/(xs(b,1)+xs(b,2))/(xs(d,1)+xs(d,2))
               wcd=sCM*xs(c,d)/(xs(c,1)+xs(c,2))/(xs(d,1)+xs(d,2))
            else
               goto 999
            endif
         endif
         if(d.eq.0)then
            sigma_abcd=(1d0/ea/wab)**alpha*(1d0/(ea+eb)+1d0/ec)*1d0/wbc
     &                +(1d0/eb/wab)**alpha*(1d0/(ea+eb)+1d0/ec)*1d0/wac
     &                +(1d0/ea/wac)**alpha*(1d0/(ea+ec)+1d0/eb)*1d0/wbc
     &                +(1d0/ec/wac)**alpha*(1d0/(ea+ec)+1d0/eb)*1d0/wab
     &                +(1d0/eb/wbc)**alpha*(1d0/(eb+ec)+1d0/ea)*1d0/wac
     &                +(1d0/ec/wbc)**alpha*(1d0/(eb+ec)+1d0/ea)*1d0/wab
         else
            sigma_abcd=((1d0/ea/wab)**alpha+(1d0/eb/wab)**alpha)*(1d0/ec+1d0/ed)*1d0/wcd
     &                +((1d0/ec/wcd)**alpha+(1d0/ed/wcd)**alpha)*(1d0/ea+1d0/eb)*1d0/wab
         endif
         sigma = sigma + sigma_abcd
         if(i4.eq.0)then
            if( (a.eq.i1.and.b.eq.i2.and.c.eq.i3) .or.
     &          (a.eq.i1.and.b.eq.i3.and.c.eq.i2) .or.
     &          (a.eq.i2.and.b.eq.i1.and.c.eq.i3) .or.
     &          (a.eq.i2.and.b.eq.i3.and.c.eq.i1) .or.
     &          (a.eq.i3.and.b.eq.i1.and.c.eq.i2) .or.
     &          (a.eq.i3.and.b.eq.i2.and.c.eq.i1) ) num = num + sigma_abcd
         else
            if( (a.eq.i1.and.b.eq.i2.and.c.eq.i3.and.d.eq.i4) .or.
     &          (a.eq.i1.and.b.eq.i2.and.c.eq.i4.and.d.eq.i3) .or.
     &          (a.eq.i2.and.b.eq.i1.and.c.eq.i3.and.d.eq.i4) .or.
     &          (a.eq.i2.and.b.eq.i1.and.c.eq.i4.and.d.eq.i3) .or.
     &          (a.eq.i3.and.b.eq.i4.and.c.eq.i1.and.d.eq.i2) .or.
     &          (a.eq.i3.and.b.eq.i4.and.c.eq.i2.and.d.eq.i1) .or.
     &          (a.eq.i4.and.b.eq.i3.and.c.eq.i1.and.d.eq.i2) .or.
     &           (a.eq.i4.and.b.eq.i3.and.c.eq.i2.and.d.eq.i1) ) num = num + sigma_abcd
         endif
         enddo
         if(sigma.le.0d0)then
            write(*,*)'Wrong sigma in Z_NNLO',sigma
            stop
         endif
         Z_NNLO = num/sigma
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


      subroutine get_ZSS_NNLO(xs,sCM,alpha,i1,i2,i3,i4,ZSS_NNLO,ierr)
c     NNLO double-soft sector functions ZSS with three (i4=0) or four (i4!=0) indices
c     i1, i2, i3, and i4 are the three (four) sector indices.
c     This function is meant to be called with (i1,i2,i3,i4) = perm(isec,jsec,ksec,lsec)
c     i1 and i3 must be in the final state as associated with the double-soft singularity
      implicit none
      include 'nexternal.inc'
      include 'all_sector_list.inc'
      include 'leg_PDGs.inc'
      integer i,a,b,c,d
      integer i1,i2,i3,i4,ierr
      double precision alpha
      double precision xs(nexternal,nexternal)
      double precision wgt,ZSS_NNLO,ea,eb,ec,ed
      double precision wab,wac,wad,wbc,wbd,wcd
      double precision num,sigma,sCM
c
c     initialise
      ZSS_NNLO = 0d0
      ierr = 0
      wab = 0d0
      wac = 0d0
      wad = 0d0
      wbc = 0d0
      wbd = 0d0
      wcd = 0d0
      ea = 0d0
      eb = 0d0
      ec = 0d0
      ed = 0d0
      num = 0d0
      sigma = 0d0
c
c     safety checks
      if(sCM.le.0d0)then
         write(77,*)'Wrong sCM in ZSS_NNLO',sCM
         stop
      endif
      if(alpha.lt.1d0)then
         write(77,*)'Wrong alpha in ZSS_NNLO',alpha
         stop
      endif
      if(i1.le.2.or.i3.le.2) then
         write(*,*) 'First and third ZSS_NNLO indices must be > 2',i1,i3
         stop
      endif
      if(abs(sCM-xs(1,2))/sCM.gt.1d-8)then
         write(77,*)'Wrong invariants in ZSS_NNLO',sCM,xs(1,2)
         goto 999
      endif
c
c     build ZSS_NNLO
      do i=1,lensectors
         a=all_sector_list(1,i)
         b=all_sector_list(2,i)
         c=all_sector_list(3,i)
         d=all_sector_list(4,i)
         if(b.eq.a.or.c.eq.a.or.d.eq.a.or.d.eq.c)then
            write(*,*)'Wrong indices in ZSS_NNLO',a,b,c,d
            stop
         endif
         if(a.ne.i1.and.a.ne.i2.and.a.ne.i3.and.a.ne.i4.and.
     &      b.ne.i1.and.b.ne.i2.and.b.ne.i3.and.b.ne.i4.and.
     &      c.ne.i1.and.c.ne.i2.and.c.ne.i3.and.c.ne.i4.and.
     &      d.ne.i1.and.d.ne.i2.and.d.ne.i3.and.d.ne.i4)cycle
         if((xs(a,1)+xs(a,2))*(xs(b,1)+xs(b,2))*(xs(c,1)+xs(c,2))*(xs(d,1)+xs(d,2))*
     &      xs(a,b)*xs(a,c)*xs(a,d)*xs(b,c)*xs(b,d)*xs(c,d).ne.0d0)then
            ea=(xs(a,1)+xs(a,2))/sCM
            eb=(xs(b,1)+xs(b,2))/sCM
            ec=(xs(c,1)+xs(c,2))/sCM
            ed=(xs(d,1)+xs(d,2))/sCM
            wab=sCM*xs(a,b)/(xs(a,1)+xs(a,2))/(xs(b,1)+xs(b,2))
            wac=sCM*xs(a,c)/(xs(a,1)+xs(a,2))/(xs(c,1)+xs(c,2))
            wad=sCM*xs(a,d)/(xs(a,1)+xs(a,2))/(xs(d,1)+xs(d,2))
            wbc=sCM*xs(b,c)/(xs(b,1)+xs(b,2))/(xs(c,1)+xs(c,2))
            wbd=sCM*xs(b,d)/(xs(b,1)+xs(b,2))/(xs(d,1)+xs(d,2))
            wcd=sCM*xs(c,d)/(xs(c,1)+xs(c,2))/(xs(d,1)+xs(d,2))
         else
            goto 999
         endif
         if(i4.eq.0)then
c SHOULD IT BE d=0 OR i4=0???
            if(a.eq.i1.and.c.eq.i3) sigma = sigma + (1d0/ea/wab)**alpha*1d0/ec*1d0/wbc + (1d0/ea/wac)**alpha*1d0/(ea+ec)*1d0/wbc
            if(a.eq.i3.and.c.eq.i1) sigma = sigma + (1d0/ec/wbc)**alpha*1d0/ea*1d0/wab + (1d0/ec/wac)**alpha*1d0/(ea+ec)*1d0/wab
            if(a.eq.i1.and.b.eq.i3) sigma = sigma + (1d0/ea/wac)**alpha*1d0/eb*1d0/wbc + (1d0/ea/wab)**alpha*1d0/(ea+eb)*1d0/wbc
            if(a.eq.i3.and.b.eq.i1) sigma = sigma + (1d0/eb/wbc)**alpha*1d0/ea*1d0/wac + (1d0/eb/wab)**alpha*1d0/(ea+eb)*1d0/wac
            if(b.eq.i1.and.c.eq.i3) sigma = sigma + (1d0/eb/wab)**alpha*1d0/ec*1d0/wac + (1d0/eb/wbc)**alpha*1d0/(eb+ec)*1d0/wac
            if(b.eq.i3.and.c.eq.i1) sigma = sigma + (1d0/ec/wac)**alpha*1d0/eb*1d0/wab + (1d0/ec/wbc)**alpha*1d0/(eb+ec)*1d0/wab
            if(a.eq.i1.and.b.eq.i2.and.c.eq.i3) num = num + (1d0/ea/wab)**alpha*1d0/ec*1d0/wbc + (1d0/ea/wac)**alpha*1d0/(ea+ec)*1d0/wbc
            if(a.eq.i3.and.b.eq.i2.and.c.eq.i1) num = num + (1d0/ec/wbc)**alpha*1d0/ea*1d0/wab + (1d0/ec/wac)**alpha*1d0/(ea+ec)*1d0/wab
            if(a.eq.i1.and.b.eq.i3.and.c.eq.i2) num = num + (1d0/ea/wac)**alpha*1d0/eb*1d0/wbc + (1d0/ea/wab)**alpha*1d0/(ea+eb)*1d0/wbc
            if(a.eq.i3.and.b.eq.i1.and.c.eq.i2) num = num + (1d0/eb/wbc)**alpha*1d0/ea*1d0/wac + (1d0/eb/wab)**alpha*1d0/(ea+eb)*1d0/wac
            if(a.eq.i2.and.b.eq.i1.and.c.eq.i3) num = num + (1d0/eb/wab)**alpha*1d0/ec*1d0/wac + (1d0/eb/wbc)**alpha*1d0/(eb+ec)*1d0/wac
            if(a.eq.i2.and.b.eq.i3.and.c.eq.i1) num = num + (1d0/ec/wac)**alpha*1d0/eb*1d0/wab + (1d0/ec/wbc)**alpha*1d0/(eb+ec)*1d0/wab
         else
            if(a.eq.i1.and.c.eq.i3) sigma = sigma + (1d0/ea/wab)**alpha*1d0/ec*1d0/wcd
            if(a.eq.i3.and.c.eq.i1) sigma = sigma + (1d0/ec/wcd)**alpha*1d0/ea*1d0/wab
            if(b.eq.i1.and.c.eq.i3) sigma = sigma + (1d0/eb/wab)**alpha*1d0/ec*1d0/wcd
            if(b.eq.i3.and.c.eq.i1) sigma = sigma + (1d0/ec/wcd)**alpha*1d0/eb*1d0/wab
            if(a.eq.i1.and.d.eq.i3) sigma = sigma + (1d0/ea/wab)**alpha*1d0/ed*1d0/wcd
            if(a.eq.i3.and.d.eq.i1) sigma = sigma + (1d0/ed/wcd)**alpha*1d0/ea*1d0/wab
            if(b.eq.i1.and.d.eq.i3) sigma = sigma + (1d0/eb/wab)**alpha*1d0/ed*1d0/wcd
            if(b.eq.i3.and.d.eq.i1) sigma = sigma + (1d0/ed/wcd)**alpha*1d0/eb*1d0/wab
            if(a.eq.i1.and.b.eq.i2.and.c.eq.i3.and.d.eq.i4) num = num + (1d0/ea/wab)**alpha*1d0/ec*1d0/wcd
            if(a.eq.i3.and.b.eq.i4.and.c.eq.i1.and.d.eq.i2) num = num + (1d0/ec/wcd)**alpha*1d0/ea*1d0/wab
            if(a.eq.i2.and.b.eq.i1.and.c.eq.i3.and.d.eq.i4) num = num + (1d0/eb/wab)**alpha*1d0/ec*1d0/wcd
            if(a.eq.i4.and.b.eq.i3.and.c.eq.i1.and.d.eq.i2) num = num + (1d0/ec/wcd)**alpha*1d0/eb*1d0/wab
            if(a.eq.i1.and.b.eq.i2.and.c.eq.i4.and.d.eq.i3) num = num + (1d0/ea/wab)**alpha*1d0/ed*1d0/wcd
            if(a.eq.i3.and.b.eq.i4.and.c.eq.i2.and.d.eq.i1) num = num + (1d0/ed/wcd)**alpha*1d0/ea*1d0/wab
            if(a.eq.i2.and.b.eq.i1.and.c.eq.i4.and.d.eq.i3) num = num + (1d0/eb/wab)**alpha*1d0/ed*1d0/wcd
            if(a.eq.i4.and.b.eq.i3.and.d.eq.i2.and.d.eq.i1) num = num + (1d0/ed/wcd)**alpha*1d0/eb*1d0/wab
         endif
      enddo
      if(sigma.le.0d0)then
         write(*,*)'Wrong sigma in Z_NNLO',sigma
         stop
      endif
      ZSS_NNLO = num/sigma
c
c     sanity check
      if(abs(ZSS_NNLO).ge.huge(1d0).or.isnan(ZSS_NNLO))then
        write(77,*)'Exception caught in ZSS_NNLO',ZSS_NNLO
        goto 999
      endif
c
      return
 999  ierr=1
      return
      end
