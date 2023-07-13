      subroutine phase_space_n(x,shat,p,npart,xjac)
      implicit none
      include 'math.inc'
      integer npart
      double precision x(3*(npart-2)-4-1+2),shat
      double precision p(0:3,npart),xjac

      if(npart.eq.4) then
         call born_phase_space_2j(x,shat,p,npart,xjac)
      elseif(npart.eq.5) then
         call born_phase_space_3j(x,shat,p,npart,xjac)
      else
         write(*,*) 'The number of external particles must be either
     $               4 or 5'
         write(*,*) 'nparticles = ', npart
         write(*,*) 'Exit...'
         stop
      endif

      return
      end

      subroutine born_phase_space_2j(x,shat,p,npart,xjac)
      implicit none
      include 'math.inc'
      integer npart
      double precision x(1),shat
      double precision p(0:3,npart),xjac
c
c     local variables
      double precision app,pmod,costh,sinth,ph,ran2
      integer idum
      common/rand/idum
c
c     initialise
      p=0d0
      xjac=GeVtopb
c
      if(npart.ne.4)then
         write(*,*) 'Four particles expected in Born phase space', npart
         stop
      endif
c
c     app is modulus squared of three-momentum
      app=shat/4d0
      if(app.le.0d0)then
        xjac=0d0
        return
      endif
      pmod=sqrt(app)
      costh=max(-1d0,min(x(1)*2d0-1d0,1d0))
      sinth=sqrt(1d0-costh*costh)
      ph=RAN2(idum)*twopi
c
c     final-state momenta
      p(0,3)=pmod
      p(1,3)=pmod*sinth*cos(ph)
      p(2,3)=pmod*sinth*sin(ph)
      p(3,3)=pmod*costh
      p(0,4)=pmod
      p(1,4)=-p(1,3)
      p(2,4)=-p(2,3)
      p(3,4)=-p(3,3)
c
c     initial-state and center-of-mass momenta
      p(0,1)=sqrt(shat)/2d0
      p(3,1)=sqrt(shat)/2d0
      p(0,2)=sqrt(shat)/2d0
      p(3,2)=-sqrt(shat)/2d0
c
c     phase-space weight
      xjac=xjac*pmod/sqrt(shat)/(fourpi)**2
      xjac=xjac*fourpi
c
      return
      end


      subroutine born_phase_space_3j(x,shat,p,npart,xjac)
      implicit none
      include 'math.inc'
      integer npart
      double precision x(4),shat
      double precision p(0:3,npart),xjac
c
c     local variables
      double precision app,pmod,costh,sinth,ph,ran2
      double precision e_sup,e_inf,s1,q1(0:3),e1,p1mod
      integer idum
      common/rand/idum
      double precision xm(7),xmsq(7)
      common/c_masses/xm,xmsq
c
c     initialise
      p=0d0
      xm=0d0
      xmsq=0d0
      xjac=GeVtopb
c
c     s1 = (p1+p2)**2
      e_sup=sqrt(shat)-xm(3)
      e_inf=xm(1)+xm(2)      
      if(e_inf.ge.e_sup)then
        xjac=0d0
        return
      endif     
      e1=(e_sup-e_inf)*x(1)+e_inf
      s1=e1*e1
      xjac=xjac*(e_sup-e_inf)*2d0*e1
      xjac=xjac/twopi
c
c     q1 and p3
      app=((shat-s1-xmsq(3))**2-4d0*s1*xmsq(3))/(4d0*shat)
      if(app.le.0d0)then
         xjac=0d0
         return
      endif
      pmod=sqrt(app)
      costh=max(-1d0,min(x(2)*2d0-1d0,1d0))
      sinth=sqrt(1d0-costh*costh)
      ph=RAN2(idum)*twopi
c
c     final-state momenta
c     p3
      q1(0)=sqrt(pmod*pmod+s1)
      q1(1)=pmod*sinth*cos(ph)
      q1(2)=pmod*sinth*sin(ph)
      q1(3)=pmod*costh
      p(0,npart)=sqrt(pmod*pmod+xmsq(3))
      p(1,npart)=-q1(1)
      p(2,npart)=-q1(2)
      p(3,npart)=-q1(3)
      xjac=xjac*pmod/sqrt(shat)/(fourpi)**2
      xjac=xjac*fourpi
c
c     p1 and p2      
      app=((s1-xmsq(1)-xmsq(2))**2-4d0*xmsq(1)*xmsq(2))/(4d0*s1)
      if(app.le.0d0)then
        xjac=0d0
        return
      endif
      p1mod=sqrt(app)
      call twobody(q1,s1,p1mod,xmsq(1),xmsq(2),x(3),x(4),p(0,3),p(0,4))
c
c     initial-state and center-of-mass momenta
!      p(0,0)=sqrt(shat)
      p(0,2)=dsqrt(shat)/2d0
      p(3,2)=dsqrt(shat)/2d0
      p(0,1)=dsqrt(shat)/2d0
      p(3,1)=-dsqrt(shat)/2d0
c
c     phase-space weight
      xjac=xjac*p1mod/sqrt(s1)/(fourpi)**2
      xjac=xjac*fourpi
c
      return
      end
