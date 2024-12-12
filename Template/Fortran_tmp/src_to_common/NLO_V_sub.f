      subroutine nlo_v_sub(s,v,mk2,ml2,mu,ccBLO,INLO)
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      double precision s,v,yp
      double precision mk2,ml2,mu
      double precision INLO,ddilog,ccBLO,Q2
      double precision lam, eta, alpha, beta
      double precision Q, mk, ml

      
      Q2=s+ml2+mk2
      lam=(s*v)**2d0
      Q = dsqrt(Q2)
      mk = dsqrt(mk2)
      ml = dsqrt(ml2)
      eta = (1d0-v)/(1d0+v) !(s-dsqrt(lam))/(s+dsqrt(lam))
      alpha = (Q-ml)/(Q+ml)
      beta = dsqrt( ((Q-ml)**2-mk2)/((Q+ml)**2-mk2) )

      if(alpha.lt.beta.or.alpha.gt.1d0.or.beta.lt.0d0)then
         write(*,*)'Wrong alpha, beta, values in NLO_V_sub'
         write(*,*)alpha,beta
         stop
      endif

c$$$      INLO = - CCBLO * (
c$$$     $ s/dsqrt(lam)*(
c$$$     $ 1d0/2d0*dlog(eta)*(-dlog(s/mu**2)-dlog(lam**2/(mk2*ml2*Q2*s)))-
c$$$     $ 1d0/4d0*dlog(eta)**2+dlog(eta)+2d0*dsqrt(lam)/s+
c$$$     $ dlog(((Q+ml)**2-mk2)/(ml*(Q+ml)))*dlog((1d0-beta)/(1d0+beta))+
c$$$     $ ddilog(-2d0*beta/(1d0-beta))-ddilog(2d0*beta/(1d0+beta))+
c$$$     $ ddilog((1d0+beta)/2d0)-ddilog((1d0-beta)/2d0)-
c$$$     $ ddilog((1d0-alpha)/(1d0-beta))+ddilog((1d0-alpha)/(1d0+beta))-
c$$$     $ ddilog((1d0-beta)/(1d0+alpha))+ddilog((1d0+beta)/(1d0+alpha))+
c$$$     $ ddilog(-2d0*beta/(alpha-beta))-ddilog(2d0*beta/(alpha+beta))-
c$$$     $ ddilog(1d0-eta))+
c$$$     $ 2d0-(2d0*ml2-Q2+s)*dlog((1d0-beta)/(1d0+beta))/dsqrt(lam)+
c$$$     $ (2d0*ml2-2d0*Q2-s)*dlog(eta)/(2d0*dsqrt(lam))-
c$$$     $     dlog(lam**2/(mk2*ml2*Q2*s)) - dlog(s/mu**2) )



      INLO = - 0.5d0*CCBLO*(8d0+(4d0*mk2*dlog((1d0-beta)/(1d0+beta)))/dsqrt(lam)
     $     -(4d0*ml2*dlog((1d0-beta)/(1d0+beta)))/dsqrt(lam) - 
     $     (2d0*s*dlog((1d0-beta)/(1d0+beta))*dlog((1d0-beta**2)/
     $     (2d0*(1d0+alpha))))/dsqrt(lam)+(2d0*ml2*dlog(eta))/dsqrt(lam) - 
     $     (2d0*Q2*dlog(eta))/dsqrt(lam)+(s*dlog(eta))/dsqrt(lam)-
     $     (s*dlog(eta)**2)/(2d0*dsqrt(lam)) - 
     $     2d0*dlog(lam**2/(mk2*ml2*Q2*s))-(s*dlog(eta)*
     $     dlog(lam**2/(mk2*ml2*Q2*s)))/dsqrt(lam) - 
     $     2d0*dlog(s/mu**2)-(s*dlog(eta)*dlog(s/mu**2))/dsqrt(lam)-
     $     (2d0*s*ddilog((1d0-alpha)/(1d0-beta)))/dsqrt(lam) - 
     $     (2d0*s*ddilog((1d0-beta)/2d0))/dsqrt(lam) -
     $     (2d0*s*ddilog((1d0-beta)/(1d0+alpha)))/dsqrt(lam) + 
     $     (2d0*s*ddilog((-2d0*beta)/(1d0-beta)))/dsqrt(lam) +
     $     (2d0*s*ddilog((-2d0*beta)/(alpha-beta)))/dsqrt(lam) + 
     $     (2d0*s*ddilog((1d0-alpha)/(1d0+beta)))/dsqrt(lam) -
     $     (2d0*s*ddilog((2d0*beta)/(1d0+beta)))/dsqrt(lam) + 
     $     (2d0*s*ddilog((1d0+beta)/2d0))/dsqrt(lam) +
     $     (2d0*s*ddilog((1d0+beta)/(1d0+alpha)))/dsqrt(lam) - 
     $     (2d0*s*ddilog((2d0*beta)/(alpha+beta)))/dsqrt(lam) -
     $     (2d0*s*ddilog(1d0-eta))/dsqrt(lam))

      end

