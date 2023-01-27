      subroutine trn_ffn_variables(pp,children,xis,kts,ss_ij,y,x)
      implicit none
      include 'nexternal.inc'
      integer i,j
      real * 8 pp(0:3,NEXTERNAL),allpfs(0:3,2)
      integer children(2)
      real * 8 ss_ij(1),ss_ijtot
      real * 8 xi,xj,kitil(0:3),kjtil(0:3),y,x
      real * 8 xis(2),kts(0:3,2)
      real * 8 prec(0:3)
      real * 8 dot
      integer rec_legnumber
      
      ss_ijtot=0d0

      !!!!!! This only works for NLO computation with 3 particles 
      !!!!!! in the final state.
      
      
      do i=1,2
         allpfs(:,i) = pp(:,children(i))
      enddo

!     TODO
!     call get_recoiler(...,prec,rec_legnumber)



      ss_ij(1) = 2d0*dot(allpfs(:,1),allpfs(:,2))
      
      do i=1,1 !!!!! TO BE GENERALIZED LATER
         ss_ijtot = ss_ijtot + ss_ij(i)
      enddo


      call GET_SUDAKOV_DECOMP_FF(allpfs(:,1),allpfs(:,2),prec,xi,xj,
     $     kitil,kjtil,y,x)
      
      xis = [xi,xj]
      kts(:,1) = kitil
      kts(:,2) = kjtil
      
      
      end




      subroutine trn_ifn_variables(pp,children,xis,kts,ss_ij,z,v)
      implicit none
      include 'nexternal.inc'
      integer i,j
      real * 8 pp(0:3,NEXTERNAL),allpfs(0:3,2)
      integer children(2)
      real * 8 ss_ij(1),ss_ijtot
      real * 8 xi,xj,kitil(0:3),kjtil(0:3),z,v,y,x
      real * 8 xis(2),kts(0:3,2)
      real * 8 prec(0:3)
      integer rec_legnumber
      real * 8 dot

      ss_ijtot=0d0

      !!!!!! This only works for NLO computation with 3 particles 
      !!!!!! in the final state.
      
      
      do i=1,2
         allpfs(:,i) = pp(:,children(i))
      enddo

!     TODO
!     call get_recoiler(...,prec,rec_legnumber)



      ss_ij(1) = 2d0*dot(allpfs(:,1),allpfs(:,2))
      
      do i=1,1 !!!!! TO BE GENERALIZED LATER
         ss_ijtot = ss_ijtot + ss_ij(i)
      enddo


      call GET_SUDAKOV_DECOMP_IF(allpfs(:,1),allpfs(:,2),prec,
     $     rec_legnumber,xi,xj,kitil,kjtil,y,x)
      
      xis = [xi,xj]
      kts(:,1) = kitil
      kts(:,2) = kjtil
      
      
      end



      SUBROUTINE GET_SUDAKOV_DECOMP_FF(ki,kj,kr,xi,xj,kitil,kjtil,y,x)
      IMPLICIT NONE
C
C---- returns the sudakov decomposition (x and kt) of ki and kj, with kr being the other
C---- reference momentum
C
C     
C     ARGUMENTS 
C
      double precision ki(0:3),kj(0:3),kr(0:3),kitil(0:3),kjtil(0:3)
      double precision xi,xj,y,x
CF2PY INTENT(IN)  :: ki
CF2PY INTENT(IN)  :: kj
CF2PY INTENT(IN)  :: kr
CF2PY INTENT(OUT) :: xi
CF2PY INTENT(OUT) :: xj
CF2PY INTENT(OUT) :: kitil
CF2PY INTENT(OUT) :: kjtil
CF2PY INTENT(OUT) :: y
CF2PY INTENT(OUT) :: x

C     LOCAL
      double precision sij,sir,sjr,k2
      double precision k(0:3)
C     
C     EXTERNAL
C     
      real * 8 dot
C      EXTERNAL dot


      sij=dot(ki,kj)
      sir=dot(ki,kr)
      sjr=dot(kj,kr)

      k=ki(:)+kj(:)
      k2=dot(k,k)
      xi=sir/(sir+sjr)
      xj=1d0-xi

      kitil=ki-k*xi-kr*((dot(k,ki)/k2)-xi)*(k2/dot(k,kr))
      kjtil=kj-k*xj-kr*((dot(k,kj)/k2)-xj)*(k2/dot(k,kr))

C---- for final kr
      y=sij/(sij + sir + sjr)
C---- for initial kr
      x=1d0-sij/(sir + sjr)

      END


      SUBROUTINE GET_SUDAKOV_DECOMP_IF(kj,ki,kr,kr_leg_number,
     $      xj,xi,kjtil,kitil,z,v)
      IMPLICIT NONE
C
C---- returns the sudakov decomposition (x and kt) of ki and kj, with kr being the other
C---- reference momentum for the IF case. kj=initial, ki=final
C
C     
C     ARGUMENTS 
C
      integer kr_leg_number
      double precision kj(0:3),ki(0:3),kr(0:3),kitil(0:3),kjtil(0:3)
      double precision xi,xj,z,v
CF2PY INTENT(IN)  :: ki
CF2PY INTENT(IN)  :: kj
CF2PY INTENT(IN)  :: kr
CF2PY INTENT(IN)  :: kr_leg_number
CF2PY INTENT(OUT) :: xi
CF2PY INTENT(OUT) :: xj
CF2PY INTENT(OUT) :: kitil
CF2PY INTENT(OUT) :: kjtil
CF2PY INTENT(OUT) :: z
CF2PY INTENT(OUT) :: v
C
C     LOCAL
C
      double precision sij,sir,sjr
      double precision pa(0:3),pAA(0:3),ktA(0:3),kti(0:3),pi(0:3)
      real*8 napAA, nbpAA, nanb, nbpa, napi, nbpi
C     
C     EXTERNAL
C     
      real * 8 dot
c      EXTERNAL dot
      double precision na(0:3), nb(0:3)



      sij=dot(ki,kj)
      sir=dot(ki,kr)
      sjr=dot(kj,kr)

      call collinear_and_reference(kj,na,nb)

      pa = kj
C---- Compute the sum of momenta
      pAA = pa - ki
C---- Pre-compute variables
      napAA = dot(na,pAA)
      nbpAA = dot(nb,pAA)
      nanb = dot(na,nb)
      nbpa = dot(nb,pa)
C---- zA = nbpA / nbpa
      ktA = pAA - (nbpAA * na + napAA * nb) / nanb
C---- Compute all kinematic variables
      pi = ki
      napi = dot(na,pi)
      nbpi = dot(nb,pi)
C---- zi = nbpi / nbpa
      kti = pi - (nbpi*na + napi*nb) / nanb
C---- kts
      kjtil = ktA
      kitil = kti
  
      IF (kr_leg_number > 2) THEN
          xj = 1d0 - sir / (sij+sjr)
          xi = 1d0 - xj
      ELSE
          xj = 1d0 - ((sij + sir) / sjr)
          xi = 1d0 - xj
      END IF
  
C---- Variable for initial splitting - final rec
      z = sij / (sij + sjr)
C---- Variables for initial splitting - initial rec
      v = sij / (sij + sir)

      END




      SUBROUTINE collinear_and_reference(p,na,nb)
      IMPLICIT NONE
C****************************************************************************
C     Given a momentum, return normalized vectors on the light-cone
C****************************************************************************
      integer I
      double precision p(0:3)
      double precision n(0:2),na(0:3),nb(0:3)
CF2PY INTENT(IN)  :: p
CF2PY INTENT(OUT) :: na
CF2PY INTENT(OUT) :: nb

      !!!! This subroutine takes initial state momenta as input
      !!!! They are taken in the center of mass frame so they 
      !!!! do not have a transverse compoment wrt to the collision axis


      if(dsqrt(p(1)**2+p(2)**2).gt.1d-8) then
         write(*,*) 'The initial state momenta cannot have '
         write(*,*) 'a non-vanishing transverse momentum!'
         write(*,*) 'EXIT...'
         call exit(-1)
      endif

      n(0:1) = p(1:2) 
      n(2) = p(3)/abs(p(3))
      

C---- For good phase space points, p(0) >= 0, but test PS points might not be valid
      IF (p(0).ge.0d0) THEN
          na = (/1d0, n(0), n(1), n(2)/)
          nb = (/1d0, -n(0), -n(1), -n(2)/)
      ELSE
          na = (/1d0, -n(0), -n(1), -n(2)/)
          nb = (/1d0, n(0), n(1), n(2)/)
      END IF

      END
