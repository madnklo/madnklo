      subroutine phase_space_n(x,shat,p,npart,xjac)
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
