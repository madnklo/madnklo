      subroutine invariants_from_p(p,nparticles,xs,ierr)
      implicit none
      include 'math.inc'
      integer nparticles
c     nparticles is the number of (inital+final)-state particles
      integer i,j,ii,ierr
      double precision p(0:3,nparticles)
      double precision xs(nparticles,nparticles)
      double precision xstmp(nparticles),xstmpsum
      double precision xsq(3:nparticles)
      double precision dot
c     TODO: set input value q()=(sqrt(s),0,0,0)
c     TODO: read sCM value 
      double precision q(0:3)
      logical do_checks
      parameter(do_checks=.true.)
c     
c     initialise
      xs=0d0
      xstmp=0d0
      xstmpsum=0d0
      xsq=0d0
      ierr=0
c
      sCM=q(0)**2
c
c     build invariants from p
      do i=1,nparticles-1
         do j=i+1,nparticles
            xs(i,j)=2d0*dot(p(0,i),p(0,j))
            xs(j,i)=xs(i,j)
         enddo
      enddo
c
      do i=3,nparticles
         xsq(i)=2d0*dot(q(0),p(0,j))
      enddo
c
c     output checks
      if(do_checks)then
         do i=1,nparticles
            do j=i+1,nparticles
               if(xs(i,j).le.0d0)then
                  write(77,*)'Inaccuracy 1 in invariants_from_p',i,j,xs(i,j)
                  goto 999
               endif
            enddo
         enddo
         do i=3,nparticles
            do j=3,nparticles
               xstmp(i)=xstmp(i)+xs(i,j)
            enddo
            if(abs(xsq(i)-xstmp(i))/abs(xsq(i)).gt.tiny1)then
               write(77,*)'Inaccuracy 2 in invariants_from_p',
     &         i,abs(xsq(i)-xstmp(i))/abs(xsq(i)),xsq(i),xstmp(i)
               goto 999
            endif
         enddo
         do i=3,nparticles
            if(abs(xsq(i)-xs(1,i)-xs(2,i))/abs(xsq(i)).gt.tiny1)then
               write(77,*)'Inaccuracy 3 in invariants_from_p',
     &         i,abs(xsq(i)-xs(1,i)-xs(2,i))/abs(xsq(i)),xsq(i),xs(1,i),xs(2,i)
               goto 999
            endif
         enddo
         do i=3,nparticles
            xstmpsum=xstmpsum+xsq(i)/2d0
         enddo
         if(abs(xstmpsum-sCM)/sCM.gt.tiny1)then
            write(77,*)'Inaccuracy 4 in invariants_from_p',
     &      abs(xstmpsum-sCM)/sCM,xstmpsum,sCM
            goto 999
         endif
         do i=1,nparticles
            if(xs(i,i).ne.0d0)then
               write(*,*)'Error 1 in invariants_from_p',i,xs(i,i)
               stop
            endif
         enddo
      endif
c
      return
 999  ierr=1
      return
      end


      subroutine mapping_p_to_pbar_FFF(iU,iS,iB,xp,xs,nparticles,xpbar,xsbar,ierr)
c     mapping from unbarred to barred momenta for U,S,B in the final state
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      integer iU,iS,iB,i,j,k,l,ii,imapsum1,imapsum2,nparticles,ierr
      double precision xs(-2:maxdim,-2:maxdim),xsbar(-2:maxdim,-2:maxdim)
      double precision xp(0:3,-2:nparticles),xpbar(0:3,-2:nparticles-1)
      double precision psq,y,stotbar
      logical :: is_first = .true.
c
c     initialise
      xsbar=0d0
      stotbar=0d0
      imapsum1=0
      imapsum2=0
      ierr=0
c
c     input checks
      if(is_first)then
         if(nparticles.ne.npartNLO.and.nparticles.ne.npartNNLO)then
            write(*,*)'Error 1 in mapping_p_to_pbar_FFF',nparticles
            stop
         endif
         if(imap(iU,iS,iU,0,0,nparticles).ne.imap(iS,iS,iU,0,0,nparticles))then
            write(*,*)'Error 2 in mapping_p_to_pbar_FFF'
            write(*,*)imap(iU,iS,iU,0,0,nparticles),imap(iS,iS,iU,0,0,nparticles)
            write(*,*)iS,iU,nparticles
            stop
         endif
         do j=-2,nparticles
            if(j.eq.0)cycle
            if(imap(j,iS,iU,0,0,nparticles).lt.-2.or.
     &         imap(j,iS,iU,0,0,nparticles).gt.nparticles-1.or.
     &         imap(j,iS,iU,0,0,nparticles).eq.0)then
               write(*,*)'Error 3 in mapping_p_to_pbar_FFF'
               write(*,*)imap(j,iS,iU,0,0,nparticles),j,iS,iU,nparticles
               stop
            endif
            if(imap(j,iS,iU,0,0,nparticles).ne.imap(j,iU,iS,0,0,nparticles))then
               write(*,*)'Error 4 in mapping_p_to_pbar_FFF'
               write(*,*)imap(j,iS,iU,0,0,nparticles),imap(j,iU,iS,0,0,nparticles)
               write(*,*)j,iS,iU,nparticles
               stop
            endif
            if(j.ne.iS)imapsum1=imapsum1+imap(j,iS,iU,0,0,nparticles)
            if(j.ne.iU)imapsum2=imapsum2+imap(j,iS,iU,0,0,nparticles)
         enddo
         if(imapsum1.ne.imapsum2)then
            write(*,*)'Error 5 in mapping_p_to_pbar_FFF'
            write(*,*)imapsum1,imapsum2,iS,iU,nparticles
            stop
         endif
         if(2*imapsum1.ne.(nparticles-1)*nparticles-6)then
            write(*,*)'Error 6 in mapping_p_to_pbar_FFF'
            write(*,*)imapsum1,iS,iU,nparticles
            stop
         endif
         is_first=.false.
      endif
c
      psq=xs(iB,iS)+xs(iB,iU)+xs(iS,iU)
      y=xs(iS,iU)/psq
      if(psq.eq.0d0)then
         write(*,*)'Wrong psq in mapping_p_to_pbar_FFF',psq
         stop
      endif
      if(y.ge.1d0.or.y.le.0d0)then
         write(77,*)'Inaccurate y in mapping_p_to_pbar_FFF',y
         goto 999
      endif
c
c     barred momenta
      do i=0,3
         xpbar(i,imap(iB,iS,iU,0,0,nparticles))=xp(i,iB)/(1d0-y)
         xpbar(i,imap(iS,iS,iU,0,0,nparticles))=xp(i,iU)+xp(i,iS)-xp(i,iB)*y/(1d0-y)
      enddo
      do j=-2,nparticles
         if(j.eq.iB.or.j.eq.iS.or.j.eq.iU.or.j.eq.0)cycle
         do i=0,3
            xpbar(i,imap(j,iS,iU,0,0,nparticles))=xp(i,j)
         enddo
      enddo
      do i=0,3
         xpbar(i,0)=xp(i,0)
      enddo
c
c     barred invariants
      call invariants_from_p(xpbar,nparticles-1,xsbar,ierr)
c
      return
 999  ierr=1
      return
      end


      function dot(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision dot,tmp
c
      tmp=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      dot=tmp
c
      return
      end


      subroutine rotate_to_z(n,p,pr)
c     performs on p = (p(1),p(2),p(3)) the rotation that brings
c     unit vector n = (n(1),n(2),n(3)) into (0,0,-1)
c     input: n, p
c     output pr
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      double precision n(3),nsqred,nmod
      double precision p(3),pr(3)
c
c     input check
      nsqred=n(1)**2+n(2)**2
      nmod=sqrt(nsqred+n(3)**2)
      if(abs(nmod-1d0).gt.tiny1)then
         write(*,*)'Wrong unit vector in rotate_to_z',n(1),n(2),n(3)
         stop
      endif
c
      pr(1)=(n(2)**2*p(1)-n(1)*n(2)*(1d0+n(3))*p(2)+n(1)*(-n(1)*n(3)*p(1)+nsqred*p(3)))/nsqred
      pr(2)=(-n(1)*n(2)*(1d0+n(3))*p(1)+n(1)**2*p(2)+n(2)*(-n(2)*n(3)*p(2)+nsqred*p(3)))/nsqred
      pr(3)=-n(1)*p(1)-n(2)*p(2)-n(3)*p(3)
c
      return
      end


      subroutine rotate_to_z_inv(n,p,pr)
c     performs on p = (p(1),p(2),p(3)) the rotation that brings
c     unit vector (0,0,-1) into n = (n(1),n(2),n(3)) 
c     input: n, p
c     output pr
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      double precision n(3),nsqred,nmod
      double precision p(3),pr(3)
c
c     input check
      nsqred=n(1)**2+n(2)**2
      nmod=sqrt(nsqred+n(3)**2)
      if(abs(nmod-1d0).gt.tiny1)then
         write(*,*)'Wrong unit vector in rotate_to_z',n(1),n(2),n(3)
         stop
      endif
c
      pr(1)=(n(2)**2*p(1)-n(1)*n(2)*(1d0+n(3))*p(2)-n(1)*(n(1)*n(3)*p(1)+nsqred*p(3)))/nsqred
      pr(2)=(-n(1)*n(2)*(1d0+n(3))*p(1)+n(1)**2*p(2)-n(2)*(n(2)*n(3)*p(2)+nsqred*p(3)))/nsqred
      pr(3)=n(1)*p(1)+n(2)*p(2)-n(3)*p(3)
c
      return
      end


      subroutine get_kt(ia,ib,ir,p,pb,nparticles,wa,wb,wr,pkt,ktkt,ierr)
      implicit none
      include 'math.inc'
      integer nparticles
c     nparticles is the number of final-state particles
      integer i,j,ia,ib,ir,ierr
      double precision p(0:3,nparticles),pb(0:3,nparticles-1)
      double precision pkt(nparticles-1),ktkt,ktkt2
      double precision kt(0:3),wa,wb,wr,dot
c     
c     initialise
      pkt=0d0
      ktkt=0d0
      ierr=0
c
c     define kperp (kt)
      do i=0,3
         kt(i)=wa*p(i,ia)+wb*p(i,ib)+wr*p(i,ir)
      enddo
c
c     build invariants p(i).eps
      do i=1,nparticles-1
         pkt(i)=dot(kt,pb(0,i))
      enddo
      ktkt=wa*wb*2d0*dot(p(0,ia),p(0,ib))
c
c     check
      ktkt2=dot(kt,kt)
c     TODO: read sCM
      if(abs(ktkt-ktkt2)/sCM.gt.tiny1)then
         write(77,*)'Inaccuracy 1 in get_eps',ktkt,ktkt2
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
