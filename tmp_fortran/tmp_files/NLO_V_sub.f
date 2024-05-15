      subroutine nlo_v_sub(s,v,yp,mk2,ml2,mu,ccBLO,INLO)
      implicit none
      include 'math.inc'
      double precision s,v,yp
      double precision mk2,ml2,mu
      double precision INLO,ddilog
      INLO=0d0


!Jkl1
      INLO = INLO - CCBLO*(
     $     (dlog((1 + vv)/(1 - vv))*(-2*dlog(mk2) + dlog((1 + vv)/(1 - vv)) + 
     $     4*dlog((s*yp)/mu)) + 4*ddilog((2*vv)/(-1 + vv)) + 
     $     4*(-0.5*dlog(((2*mk2 + s)*(1 + vv))/(2*mk2 + s - s*vv**2))**2 + 
     $     (s**2*(-1 + vv)*vv**3*
     $     (-1 + dlog((s*(-1 + vv)*vv)/((-2*mk2 + s*(-1 + vv))*yp))))/
     $     ((2*mk2 + s)*(2*mk2 + s - s*vv)) - 
     $     dlog((vv*(2*mk2 + s - s*vv))/(2*mk2 + s - s*vv**2))*
     $     dlog((s*(-1 + vv)*vv)/((-2*mk2 + s*(-1 + vv))*yp)) - 
     $     (s**2*vv**3*(1 + vv)*
     $     (-1 + dlog((s*vv*(1 + vv))/((2*mk2 + s + s*vv)*yp))))/
     $     ((2*mk2 + s)*(2*mk2 + s + s*vv)) + 
     $     (dlog(vv) + dlog((2*mk2 + s + s*vv)/(2*mk2 + s - s*vv**2)))*
     $     dlog((s*vv*(1 + vv))/((2*mk2 + s + s*vv)*yp)) + 
     $     (s*vv**2*(s*(1 + vv)*(vv - yp) - 2*mk2*yp)*
     $     (-1 + dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/
     $     ((2*mk2 + s + s*vv)*yp))))/((2*mk2 + s)*(2*mk2 + s + s*vv))
     $     + dlog(((2*mk2 + s)*(-(s*(1 + vv)*(vv - yp)) + 2*mk2*yp))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2))))**2/2. + 
     $     (2*s*vv**2*(dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/
     $     ((2*mk2 + s - s*vv)*yp)) - 
     $     dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/((2*mk2 + s + s*vv)*yp)))*
     $     dlog(1 - ((2*mk2 + s)*yp)/(s*vv**2)))/(2*mk2 + s) + 
     $     (-dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/((2*mk2 + s - s*vv)*yp)) + 
     $     dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/((2*mk2 + s + s*vv)*yp)))*
     $     dlog(1 - ((2*mk2 + s)*yp)/(s*vv**2)) + 
     $     (s*vv**2*(-dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/
     $     ((2*mk2 + s - s*vv)*yp)) + 
     $     dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/((2*mk2 + s + s*vv)*yp)))*
     $     ((2*mk2 + s)*yp + s*vv**2*dlog(1 - ((2*mk2 + s)*yp)/(s*vv**2))))/
     $     (2*mk2 + s)**2 - (s*vv**2*(-2*mk2*yp + s*(-1 + vv)*(vv + yp))*
     $     (-1 + dlog((-2*mk2*yp + s*(-1 + vv)*(vv + yp))/
     $     ((-2*mk2 + s*(-1 + vv))*yp))))/
     $     ((2*mk2 + s)*(2*mk2 + s - s*vv)) + 
     $     dlog((-2*mk2*yp + s*(-1 + vv)*(vv + yp))/((-2*mk2 + s*(-1 + vv))*yp))*
     $     dlog(((2*mk2 + s - s*vv)*(2*mk2*yp + s*(-vv**2 + yp)))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2)))) - 
     $     dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/((2*mk2 + s + s*vv)*yp))*
     $     dlog(((2*mk2 + s + s*vv)*(2*mk2*yp + s*(-vv**2 + yp)))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2)))) - 
     $     ddilog(-(((2*mk2 + s)*(-1 + vv))/(2*mk2 + s - s*vv**2))) + 
     $     (2*s*vv**2*(dlog((vv*(2*mk2 + s - s*vv))/(2*mk2 + s - s*vv**2))*
     $     dlog((s*(-1 + vv)*vv)/((-2*mk2 + s*(-1 + vv))*yp)) + 
     $     ddilog(-(((2*mk2 + s)*(-1 + vv))/(2*mk2 + s - s*vv**2)))))/
     $     (2*mk2 + s) - (s**2*vv**4*
     $     (dlog((vv*(2*mk2 + s - s*vv))/(2*mk2 + s - s*vv**2))*
     $     dlog((s*(-1 + vv)*vv)/((-2*mk2 + s*(-1 + vv))*yp)) + 
     $     ddilog(-(((2*mk2 + s)*(-1 + vv))/(2*mk2 + s - s*vv**2)))))/
     $     (2*mk2 + s)**2 - ddilog(
     $     (2*mk2 + s - s*vv**2)/((2*mk2 + s)*(1 + vv))) + 
     $     (s*vv**2*(Pi**2 + 3*dlog(((2*mk2 + s)*(1 + vv))/
     $     (2*mk2 + s - s*vv**2))**2 - 
     $     6*(dlog(vv) + dlog((2*mk2 + s + s*vv)/(2*mk2 + s - s*vv**2)))*
     $     dlog((s*vv*(1 + vv))/((2*mk2 + s + s*vv)*yp)) + 
     $     6*ddilog((2*mk2 + s - s*vv**2)/((2*mk2 + s)*(1 + vv)))))/
     $     (3.*(2*mk2 + s)) - (s**2*vv**4*
     $     (Pi**2 + 3*dlog(((2*mk2 + s)*(1 + vv))/(2*mk2 + s - s*vv**2))**2 - 
     $     6*(dlog(vv) + dlog((2*mk2 + s + s*vv)/(2*mk2 + s - s*vv**2)))*
     $     dlog((s*vv*(1 + vv))/((2*mk2 + s + s*vv)*yp)) + 
     $     6*ddilog((2*mk2 + s - s*vv**2)/((2*mk2 + s)*(1 + vv)))))/
     $     (6.*(2*mk2 + s)**2) + 
     $     ddilog(((2*mk2 + s - s*vv)*yp)/(s*(-1 + vv)*vv)) - 
     $     ddilog(((2*mk2 + s + s*vv)*yp)/(s*vv*(1 + vv))) + 
     $     ddilog((s*vv*(-2*mk2 + s*(-1 + vv**2)))/
     $     ((2*mk2 + s)*(-(s*(1 + vv)*(vv - yp)) + 2*mk2*yp))) - 
     $     (s*vv**2*(Pi**2 + 3*dlog(((2*mk2 + s)*
     $     (-(s*(1 + vv)*(vv - yp)) + 2*mk2*yp))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2))))**2 - 
     $     6*dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/
     $     ((2*mk2 + s + s*vv)*yp))*
     $     dlog(((2*mk2 + s + s*vv)*(-(s*vv**2) + 2*mk2*yp + s*yp))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2)))) + 
     $     6*ddilog((s*vv*(-2*mk2 + s*(-1 + vv**2)))/
     $     ((2*mk2 + s)*(-(s*(1 + vv)*(vv - yp)) + 2*mk2*yp)))))/
     $     (3.*(2*mk2 + s)) + (s**2*vv**4*
     $     (Pi**2 + 3*dlog(((2*mk2 + s)*(-(s*(1 + vv)*(vv - yp)) + 2*mk2*yp))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2))))**2 - 
     $     6*dlog((s*(1 + vv)*(vv - yp) - 2*mk2*yp)/
     $     ((2*mk2 + s + s*vv)*yp))*
     $     dlog(((2*mk2 + s + s*vv)*(-(s*vv**2) + 2*mk2*yp + s*yp))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2)))) + 
     $     6*ddilog((s*vv*(-2*mk2 + s*(-1 + vv**2)))/
     $     ((2*mk2 + s)*(-(s*(1 + vv)*(vv - yp)) + 2*mk2*yp)))))/
     $     (6.*(2*mk2 + s)**2) + 
     $     ddilog(((2*mk2 + s)*(-2*mk2*yp + s*(-1 + vv)*(vv + yp)))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2)))) - 
     $     (2*s*vv**2*(dlog((-2*mk2*yp + s*(-1 + vv)*(vv + yp))/
     $     ((-2*mk2 + s*(-1 + vv))*yp))*
     $     dlog(((2*mk2 + s - s*vv)*(2*mk2*yp + s*(-vv**2 + yp)))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2)))) + 
     $     ddilog(((2*mk2 + s)*(-2*mk2*yp + s*(-1 + vv)*(vv + yp)))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2))))))/(2*mk2 + s) + 
     $     (s**2*vv**4*(dlog((-2*mk2*yp + s*(-1 + vv)*(vv + yp))/
     $     ((-2*mk2 + s*(-1 + vv))*yp))*
     $     dlog(((2*mk2 + s - s*vv)*(2*mk2*yp + s*(-vv**2 + yp)))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2)))) + 
     $     ddilog(((2*mk2 + s)*(-2*mk2*yp + s*(-1 + vv)*(vv + yp)))/
     $     (s*vv*(-2*mk2 + s*(-1 + vv**2))))))/(2*mk2 + s)**2))/(4d0*vv)
     $     )

            
!     Jkl2
            INLO = INLO - CCBLO*(
     $           (s*yp-(mk2 + s)*dlog((mk2+s*yp)/mk**2))/s
            )
!     Jkk
            INLO = INLO - CCBLO*(
     $           1d0+dlog(mu**2/s)+dlog(mk2/(mu**2*yp**2))/2d0 + 
     $           ((mk2+s)*dlog((mk2+s*yp)/mk2))/s
     $           )
            
            !Jll
            INLO = INLO - CCBLO*(
     $           ((s*(-2*ml2*s**2*(-1 + vv**2) - s**3*(-1 + vv**2) + 
     $           Sqrt(s**2)*vv*(8*ml**4 + s**2*(-1 + vv**2)))*
     $           Log((2*ml2 + s*ypl - (s*(2*ml2 + s)*ypl)/Sqrt(s**2*vv**2))/(2.*ml**2)))/
     $           ((2*ml2 + s - Sqrt(s**2)*vv)*(s + Sqrt(s**2)*vv)) + 
     $           (s*(2*ml2*s**2*(-1 + vv**2) + s**3*(-1 + vv**2) + 
     $           Sqrt(s**2)*vv*(8*ml**4 + s**2*(-1 + vv**2)))*
     $           Log((2*ml2*Sqrt(s**2)*vv)/
     $           (s*(s + Sqrt(s**2)*vv)*ypl + 2*ml2*(Sqrt(s**2)*vv + s*ypl))))/
     $           ((s - Sqrt(s**2)*vv)*(2*ml2 + s + Sqrt(s**2)*vv)) + 
     $           (4*mk**2*ml2*(Log(mk**(-2)) - 
     $           (s*Log(1 + (2*s*vv)/(Sqrt(s**2) - s*vv)))/(Sqrt(s**2)*vv) + 
     $           2*Log((s*ypl)/mu)))/(-1 + vv**2))/(2d0*s**2)
     $           )

      
      end
