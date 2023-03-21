      subroutine invariants_from_p(p,nparticles,xs,ierr)
      implicit none
      include 'math.inc'
      integer nparticles
c     nparticles is the number of (inital+final)-state particles
      integer i,j,ierr
      double precision p(0:3,nparticles)
      double precision xs(nparticles,nparticles)
      double precision dot
c     
c     initialise
      xs=0d0
      ierr=0
c
c     build invariants from p
      do i=1,nparticles-1
         do j=i+1,nparticles
            xs(i,j)=2d0*dot(p(0,i),p(0,j))
            xs(j,i)=xs(i,j)
         enddo
      enddo
c
      return
 999  ierr=1
      return
      end


      double precision function dot(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision tmp
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
      double precision n(3),nsqred,nmod
      double precision p(3),pr(3)
c
c     input check
      nsqred=n(1)**2+n(2)**2
      nmod=sqrt(nsqred+n(3)**2)
      if(abs(nmod-1d0).gt.1d-8)then
         write(*,*)'Wrong unit vector in rotate_to_z',n(1),n(2),n(3)
         stop
      endif
c
      pr(1)=(n(2)**2*p(1)-n(1)*n(2)*(1d0+n(3))*p(2)+n(1)*
     & (-n(1)*n(3)*p(1)+nsqred*p(3)))/nsqred
      pr(2)=(-n(1)*n(2)*(1d0+n(3))*p(1)+n(1)**2*p(2)+n(2)*
     & (-n(2)*n(3)*p(2)+nsqred*p(3)))/nsqred
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
      double precision n(3),nsqred,nmod
      double precision p(3),pr(3)
c
c     input check
      nsqred=n(1)**2+n(2)**2
      nmod=sqrt(nsqred+n(3)**2)
      if(abs(nmod-1d0).gt.1d-8)then
         write(*,*)'Wrong unit vector in rotate_to_z',n(1),n(2),n(3)
         stop
      endif
c
      pr(1)=(n(2)**2*p(1)-n(1)*n(2)*(1d0+n(3))*p(2)-n(1)*
     & (n(1)*n(3)*p(1)+nsqred*p(3)))/nsqred
      pr(2)=(-n(1)*n(2)*(1d0+n(3))*p(1)+n(1)**2*p(2)-n(2)*
     & (n(2)*n(3)*p(2)+nsqred*p(3)))/nsqred
      pr(3)=n(1)*p(1)+n(2)*p(2)-n(3)*p(3)
c
      return
      end


