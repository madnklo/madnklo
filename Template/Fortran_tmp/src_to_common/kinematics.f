      subroutine invariants_from_p(p,nparticles,xs,ierr)
      implicit none
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
c     safety measure
            if(xs(i,j).lt.0d0)then
               write(77,*)'negative invariants in invariants_from_p'
               write(77,*)i,j,xs(i,j)
               goto 999
            endif
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


      subroutine boost(q,pboost,qprime,ierr)
c     this subroutine performs the boost according to the vector pboost of 
c     the fourvector q to the fourvector qprime. The vector q and the 
c     resulting qprime might also have the same name in the call       
c     Let's call rmboost the mass of pboost. With this transformation
c     if q= (rmboost,0,0,0), qprime=pboost
      implicit none
      integer ierr
      double precision q(0:3),pboost(0:3),qprime(0:3)
      double precision rmboost,aux,aaux
c
      ierr=0
      rmboost=pboost(0)**2-pboost(1)**2-pboost(2)**2-pboost(3)**2
      if(rmboost.gt.0d0)then
         rmboost=sqrt(rmboost)
      else
         goto 999
      endif
      aux=(q(0)*pboost(0)+q(1)*pboost(1)+q(2)*pboost(2)+q(3)*pboost(3))
     &      /rmboost
      if(pboost(0)+rmboost.gt.0d0)then
         aaux=(aux+q(0))/(pboost(0)+rmboost)
      else
         goto 999
      endif
c
      qprime(0)=aux
      qprime(1)=q(1)+aaux*pboost(1)
      qprime(2)=q(2)+aaux*pboost(2)
      qprime(3)=q(3)+aaux*pboost(3)

      return
 999  ierr=1
      return
      end


      subroutine boostinv(q,pboost,qprime,ierr)
c     this subroutine performs the inverse boost according to the vector 
c     pboost  ( which corresponds to the boost according to pboost(0), 
c     -pboost(i) ) of the fourvector q to the fourvector qprime. 
c     The vector q and the resulting qprime might also have the same 
c     name in the call       
c     Let's call rmboost the mass of pboost. With this transformation
c     if q=pboost, qprime=(rmboost,0,0,0).   
      implicit none
      integer ierr
      double precision q(0:3),pboost(0:3),qprime(0:3)
      double precision rmboost,aux,aaux
c
      ierr=0
      rmboost=pboost(0)**2-pboost(1)**2-pboost(2)**2-pboost(3)**2
      if(rmboost.gt.0d0)then
         rmboost=sqrt(rmboost)
      else
         goto 999
      endif
      aux=(q(0)*pboost(0)-q(1)*pboost(1)-q(2)*pboost(2)-q(3)*pboost(3))
     &      /rmboost
      if(pboost(0)+rmboost.gt.0d0)then
         aaux=(aux+q(0))/(pboost(0)+rmboost)
      else
         goto 999
      endif
c
      qprime(0)=aux
      qprime(1)=q(1)-aaux*pboost(1)
      qprime(2)=q(2)-aaux*pboost(2)
      qprime(3)=q(3)-aaux*pboost(3)

      return
 999  ierr=1
      return
      end




      subroutine rot(p,q,pp)

* the subroutine performs a rotation such that the 3-vector p goes into
* the 3-vector pp. The rotation is such that the vector q', of the same
*  lenght as q but along the z axis, becomes q after the rotation
      
      implicit double precision (a-h,o-z)

      dimension pp(3), p(3), q(3) 


      qmodt=q(1)**2+q(2)**2
      qmod=qmodt+q(3)**2
      qmodt=sqrt(qmodt)
      qmod=sqrt(qmod)
      if (qmod.eq.0.d0) then
        print*, ' ERROR in subroutine rot '
        print*, ' spatial q components are 0.d0 ! '
        stop
      endif

      cth=q(3)/qmod
      sth=1.d0-cth**2
      if (sth.eq.0.d0) then
        pp(1)=p(1)
        pp(2)=p(2)
        pp(3)=p(3)
        return
      endif
      sth=sqrt(sth)

      if (qmodt.eq.0.d0) then
        pp(1)=p(1)
        pp(2)=p(2)
        pp(3)=p(3)
        return
      endif
        
      cfi=q(1)/qmodt
      sfi=q(2)/qmodt

c  avoid possible problems if p and pp are the same vector:
      
      p1=p(1)
      p2=p(2)
      p3=p(3)

c  perform the rotation

      pp(1)=cth*cfi*p1-sfi*p2+sth*cfi*p3
      pp(2)=cth*sfi*p1+cfi*p2+sth*sfi*p3
      pp(3)=-sth   *p1       +cth*    p3


      return
      end



      subroutine rotinv(p,q,pp)

* the subroutine performs a rotation such that the 3-vector p goes into
* the 3-vector pp. The rotation is such that the vector q goes to the z axis.
* It performs first a rotation around z which brings the vectot in the xz-plane and
* then a rotaion around y which brings the vector along the z axis
      
      implicit real*8 (a-h,o-z)

      dimension pp(3), p(3), q(3) 


      qmodt=q(1)**2+q(2)**2
      qmod=qmodt+q(3)**2
      qmodt=sqrt(qmodt)
      qmod=sqrt(qmod)

      if (qmod.eq.0.d0) then
        print*, ' ERROR in subroutine rot '
        print*, ' spatial q components are 0.d0 ! '
        stop
      endif

      cth=q(3)/qmod
      sth=1.d0-cth**2
      if (sth.eq.0.d0) then
        pp(1)=p(1)
        pp(2)=p(2)
        pp(3)=p(3)
        return
      endif
      sth=sqrt(sth)

      if (qmodt.eq.0.d0) then
        pp(1)=p(1)
        pp(2)=p(2)
        pp(3)=p(3)
        return
      endif
        
      cfi=q(1)/qmodt
      sfi=q(2)/qmodt

c  avoid possible problems if p and pp are the same vector:
      
      p1=p(1)
      p2=p(2)
      p3=p(3)

c  perform the rotation

c      pp(1)=cth*cfi*p1-sfi*p2+sth*cfi*p3
c      pp(2)=cth*sfi*p1+cfi*p2+sth*sfi*p3
c      pp(3)=-sth   *p1       +cth*    p3

      pp(1)=cth*cfi*p1+cth*sfi*p2-sth*p3
      pp(2)=-sfi*p1   +cfi*p2
      pp(3)=+sth*cfi*p1+sth*sfi*p2+cth*p3


      return
      end
      
*****************************************************************************
* Auxiliary subroutine
* INPUT:
*       pbst        : intermediate particle momentum
*       shat        : the corresponding mass (assumed positive)
*       pmod        : the absolute value of the 3-momentum of the two
*                     daughter particles in the decay CM
*       xm1sq,xm2sq : the squared masses of the daughter particles
*       x1,x2       : two random numbers uniformly distributed in [0,1]
* OUTPUT:
*       p1,p2       : the momenta of the two daughters in the frame in
*                     which the intermediate particle has momentum pbst 
*                             p1 + p2 = pbst
*
* If pbst=(sqrt(shat),0,0,0) the boost reduces to the identity
* Ezio May 2 2007
* A rotation of the generated momentum has been added in such a way to assume as
* z-axis the boost momentum
* 
*****************************************************************************
      SUBROUTINE TWOBODY(pbst,shat,pmod,xm1sq,xm2sq,x1,x2,p1,p2)
      IMPLICIT NONE
      DOUBLE PRECISION pbst(0:3),shat,pmod,xm1sq,xm2sq,x1,x2,p1(0:3),p2(0:3)
* Local variables
      DOUBLE PRECISION app,costh,sinth,ph,rmboost,aapp,pi,twopi
      INTEGER i
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi)
       
      costh=min(x1*2.d0-1.d0,1.d0)
      costh=max(costh,-1.d0)
      sinth=sqrt(1.d0-costh*costh)
      ph=x2*twopi
* p1
      p1(0)=sqrt(pmod*pmod+xm1sq)
      p1(1)=pmod*sinth*cos(ph)
      p1(2)=pmod*sinth*sin(ph)
      p1(3)=pmod*costh
* rotate in such a way that angles are relative to the direction of pbst
      call rot(p1(1),pbst(1),p1(1))
* boost back
      rmboost=sqrt(shat)
      app=(p1(0)*pbst(0)+p1(1)*pbst(1)+p1(2)*pbst(2)+p1(3)*pbst(3))
     &      /rmboost
      aapp=(app+p1(0))/(pbst(0)+rmboost)
      p1(0)=app
      p1(1)=p1(1)+aapp*pbst(1)
      p1(2)=p1(2)+aapp*pbst(2)
      p1(3)=p1(3)+aapp*pbst(3)
      DO i=0,3
        p2(i)=pbst(i)-p1(i)
      ENDDO
      
      RETURN
      END





      subroutine getpt(p,pt)
      implicit none
      real * 8 p(0:3), pt
      pt = dsqrt(p(1)**2+p(2)**2)
      end
