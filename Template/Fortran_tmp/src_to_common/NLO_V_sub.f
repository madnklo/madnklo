      subroutine nlo_v_sub(s,v,mk2,ml2,mu,ccBLO,INLO)
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      double precision s,v,yp
      double precision mk2,ml2,mu
      double precision INLO,ddilog,ccBLO,Q2

      Q2=s+ml2+mk2
      YP=1D0+(DSQRT(ML2)-DSQRT(Q2))*2D0*DSQRT(ML2)/S
!Jkl1
      INLO = INLO - CCBLO*(
     $     (dlog((1 + v)/(1 - v))*(-2*dlog(mk2) - 4*dlog(mu) + dlog((1 + v)/(1 - v)) + 
     $     4*dlog(s*yp)) + 4*ddilog((2*v)/(-1 + v)) + 
     $     4*(-0.5*dlog(((2*mk2 + s)*(1 + v))/(2*mk2 + s - s*v**2))**2 + 
     $     (s**2*(-1 + v)*v**3*
     $     (-1 + dlog((s*(-1 + v)*v)/((-2*mk2 + s*(-1 + v))*yp))))/
     $     ((2*mk2 + s)*(2*mk2 + s - s*v)) - 
     $     dlog((v*(2*mk2 + s - s*v))/(2*mk2 + s - s*v**2))*
     $     dlog((s*(-1 + v)*v)/((-2*mk2 + s*(-1 + v))*yp)) - 
     $     (s**2*v**3*(1 + v)*
     $     (-1 + dlog((s*v*(1 + v))/((2*mk2 + s + s*v)*yp))))/
     $     ((2*mk2 + s)*(2*mk2 + s + s*v)) + 
     $     (dlog(v) + dlog((2*mk2 + s + s*v)/(2*mk2 + s - s*v**2)))*
     $     dlog((s*v*(1 + v))/((2*mk2 + s + s*v)*yp)) + 
     $     (s*v**2*(s*(1 + v)*(v - yp) - 2*mk2*yp)*
     $     (-1 + dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/((2*mk2 + s + s*v)*yp))))
     $     /((2*mk2 + s)*(2*mk2 + s + s*v)) + 
     $     dlog(((2*mk2 + s)*(-(s*(1 + v)*(v - yp)) + 2*mk2*yp))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2))))**2/2. + 
     $     (2*s*v**2*(dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/
     $     ((2*mk2 + s - s*v)*yp)) - 
     $     dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/((2*mk2 + s + s*v)*yp)))*
     $     dlog(1 - ((2*mk2 + s)*yp)/(s*v**2)))/(2*mk2 + s) + 
     $     (-dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/((2*mk2 + s - s*v)*yp)) + 
     $     dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/((2*mk2 + s + s*v)*yp)))*
     $     dlog(1 - ((2*mk2 + s)*yp)/(s*v**2)) + 
     $     (s*v**2*(-dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/
     $     ((2*mk2 + s - s*v)*yp)) + 
     $     dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/((2*mk2 + s + s*v)*yp)))*
     $     ((2*mk2 + s)*yp + s*v**2*dlog(1 - ((2*mk2 + s)*yp)/(s*v**2))))/
     $     (2*mk2 + s)**2 - (s*v**2*(-2*mk2*yp + s*(-1 + v)*(v + yp))*
     $     (-1 + dlog((-2*mk2*yp + s*(-1 + v)*(v + yp))/
     $     ((-2*mk2 + s*(-1 + v))*yp))))/((2*mk2 + s)*(2*mk2 + s - s*v)) + 
     $     dlog((-2*mk2*yp + s*(-1 + v)*(v + yp))/((-2*mk2 + s*(-1 + v))*yp))*
     $     dlog(((2*mk2 + s - s*v)*(2*mk2*yp + s*(-v**2 + yp)))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2)))) - 
     $     dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/((2*mk2 + s + s*v)*yp))*
     $     dlog(((2*mk2 + s + s*v)*(2*mk2*yp + s*(-v**2 + yp)))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2)))) - 
     $     ddilog(-(((2*mk2 + s)*(-1 + v))/(2*mk2 + s - s*v**2))) + 
     $     (2*s*v**2*(dlog((v*(2*mk2 + s - s*v))/(2*mk2 + s - s*v**2))*
     $     dlog((s*(-1 + v)*v)/((-2*mk2 + s*(-1 + v))*yp)) + 
     $     ddilog(-(((2*mk2 + s)*(-1 + v))/(2*mk2 + s - s*v**2)))))/
     $     (2*mk2 + s) - (s**2*v**4*
     $     (dlog((v*(2*mk2 + s - s*v))/(2*mk2 + s - s*v**2))*
     $     dlog((s*(-1 + v)*v)/((-2*mk2 + s*(-1 + v))*yp)) + 
     $     ddilog(-(((2*mk2 + s)*(-1 + v))/(2*mk2 + s - s*v**2)))))/
     $     (2*mk2 + s)**2 - ddilog(
     $     (2*mk2 + s - s*v**2)/((2*mk2 + s)*(1 + v))) + 
     $     (s*v**2*(Pi**2 + 3*dlog(((2*mk2 + s)*(1 + v))/(2*mk2 + s - s*v**2))**
     $     2 - 6*(dlog(v) + dlog((2*mk2 + s + s*v)/(2*mk2 + s - s*v**2)))*
     $     dlog((s*v*(1 + v))/((2*mk2 + s + s*v)*yp)) + 
     $     6*ddilog((2*mk2 + s - s*v**2)/((2*mk2 + s)*(1 + v)))))/
     $     (3.*(2*mk2 + s)) - (s**2*v**4*
     $     (Pi**2 + 3*dlog(((2*mk2 + s)*(1 + v))/(2*mk2 + s - s*v**2))**2 - 
     $     6*(dlog(v) + dlog((2*mk2 + s + s*v)/(2*mk2 + s - s*v**2)))*
     $     dlog((s*v*(1 + v))/((2*mk2 + s + s*v)*yp)) + 
     $     6*ddilog((2*mk2 + s - s*v**2)/((2*mk2 + s)*(1 + v)))))/
     $     (6.*(2*mk2 + s)**2) + 
     $     ddilog(((2*mk2 + s - s*v)*yp)/(s*(-1 + v)*v)) - 
     $     ddilog(((2*mk2 + s + s*v)*yp)/(s*v*(1 + v))) + 
     $     ddilog((s*v*(-2*mk2 + s*(-1 + v**2)))/
     $     ((2*mk2 + s)*(-(s*(1 + v)*(v - yp)) + 2*mk2*yp))) - 
     $     (s*v**2*(Pi**2 + 3*dlog(((2*mk2 + s)*
     $     (-(s*(1 + v)*(v - yp)) + 2*mk2*yp))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2))))**2 - 
     $     6*dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/((2*mk2 + s + s*v)*yp))*
     $     dlog(((2*mk2 + s + s*v)*(-(s*v**2) + 2*mk2*yp + s*yp))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2)))) + 
     $     6*ddilog((s*v*(-2*mk2 + s*(-1 + v**2)))/
     $     ((2*mk2 + s)*(-(s*(1 + v)*(v - yp)) + 2*mk2*yp)))))/
     $     (3.*(2*mk2 + s)) + (s**2*v**4*
     $     (Pi**2 + 3*dlog(((2*mk2 + s)*(-(s*(1 + v)*(v - yp)) + 2*mk2*yp))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2))))**2 - 
     $     6*dlog((s*(1 + v)*(v - yp) - 2*mk2*yp)/((2*mk2 + s + s*v)*yp))*
     $     dlog(((2*mk2 + s + s*v)*(-(s*v**2) + 2*mk2*yp + s*yp))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2)))) + 
     $     6*ddilog((s*v*(-2*mk2 + s*(-1 + v**2)))/
     $     ((2*mk2 + s)*(-(s*(1 + v)*(v - yp)) + 2*mk2*yp)))))/
     $     (6.*(2*mk2 + s)**2) + 
     $     ddilog(((2*mk2 + s)*(-2*mk2*yp + s*(-1 + v)*(v + yp)))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2)))) - 
     $     (2*s*v**2*(dlog((-2*mk2*yp + s*(-1 + v)*(v + yp))/
     $     ((-2*mk2 + s*(-1 + v))*yp))*
     $     dlog(((2*mk2 + s - s*v)*(2*mk2*yp + s*(-v**2 + yp)))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2)))) + 
     $     ddilog(((2*mk2 + s)*(-2*mk2*yp + s*(-1 + v)*(v + yp)))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2))))))/(2*mk2 + s) + 
     $     (s**2*v**4*(dlog((-2*mk2*yp + s*(-1 + v)*(v + yp))/
     $     ((-2*mk2 + s*(-1 + v))*yp))*
     $     dlog(((2*mk2 + s - s*v)*(2*mk2*yp + s*(-v**2 + yp)))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2)))) + 
     $     ddilog(((2*mk2 + s)*(-2*mk2*yp + s*(-1 + v)*(v + yp)))/
     $     (s*v*(-2*mk2 + s*(-1 + v**2))))))/(2*mk2 + s)**2))/(4d0*v)
     $     )

            
!     Jkl2
            INLO = INLO - CCBLO*(
     $     yp-((mk2+s)*dlog(1d0+(s*yp)/mk2))/s
     $     )
!     Jkk
            INLO = INLO - CCBLO*(
     $           1d0+dlog((mk2*mu**2)/(s**2*yp**2))/2d0 +
     $           ((mk2+s)*dlog(1d0+(s*yp)/mk2))/s
     $           )
            
            !Jll
            INLO = INLO - CCBLO*(
     $           -0.5d0*((4*mk2*ml2*(v*dlog(mk2) + dlog((1 + v)/(1 - v)) - 
     $           2*v*dlog((s*yp)/mu)))/(v*(-1 + v**2)) + 
     $           (s*(8*ml2**2*v + s**2*(-1 + v)**2*(1 + v) - 2*ml2*s*(-1 + v**2))*
     $           dlog((2*ml2*(v - yp) + s*(-1 + v)*yp)/(2.*ml2*v)))/
     $           ((-2*ml2 + s*(-1 + v))*(1 + v)) + 
     $           (s*(8*ml2**2*v + s**2*(-1 + v)*(1 + v)**2 + 2*ml2*s*(-1 + v**2))*
     $           dlog((2*ml2*v)/(s*(1 + v)*yp + 2*ml2*(v + yp))))/
     $           ((-1 + v)*(2*ml2 + s + s*v)))/s**2
     $           )

      
      end
