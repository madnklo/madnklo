CCCCCCCCCCCCC Fortran Functions for implementing universal_kernels.py CCCCCCCCCC

      function eikonaldip(pi,pj,ps) 
      implicit none
      real * 8 pi(0:3), pj(0:3), ps(0:3)
      real * 8 pipj, pips, pjps
      real * 8 dot,eikonaldip

      pipj = dot(pi,pj)
      pips = dot(pi,ps)
      pjps = dot(pj,ps)

      eikonaldip = pipj/(pips * pjps)

      end


CCCCCCCCCCCCC Splitting kernels CCCCCCCC
      
      
      subroutine Pgg(z,kt,res)
      implicit none
      include 'colfac.inc'
      real * 8 z
      real * 8 kt(0:3)
      real * 8 dot
      real * 8  colfac
      real * 8 res(2)

      colfac = CA

      
      res(1) = 2d0 * (z/(1d0-z) + (1d0-z)/z)
      res(2) = -2d0 * 2d0 * z * (1d0-z)/dot(kt,kt)
      
      res(1) = colfac * res(1)
      res(2) = colfac * res(2)
      
      end



      subroutine Pqq(z,kt,res)
      implicit none
      include 'colfac.inc'
      real * 8 z
      real * 8 kt(0:3)
      real * 8 dot
      real * 8  colfac
      real * 8 res(2)

      colfac = TR

      
      res(1) = 1d0
      res(2) = 4d0 * z * (1d0-z)/dot(kt,kt) 
      
      res(1) = colfac * res(1)
      res(2) = colfac * res(2)
      
      end



      subroutine Pqg(zz,res)
      implicit none
      real * 8 zz,res(2)
      
      call Pqgav(zz,res) 
      
      end


      subroutine Pgq(zz,res)
      implicit none
      real * 8 zz,res(2)
      real * 8 tmp

      tmp = 1d0-zz
      
      call Pqgav(tmp,res) 
      
      end

      

      subroutine Pqgav(z,res)
!     This subroutine fills the variable res with the coefficients
!     of the epsilon expansion. Namely:
!     Pqgav = eps^0 * res(1) + eps * res(2)
      implicit none
      include 'colfac.inc'
      real * 8 z
      real * 8 kt(0:3)
      real * 8 dot
      real * 8  colfac
      real * 8 res(2)

      colfac = CF

      
      res(1) = (1d0+z**2)/(1d0-z)
      res(2) = -(1d0-z)
      
      res(1) = colfac * res(1)
      res(2) = colfac * res(2)
      
      end





      
