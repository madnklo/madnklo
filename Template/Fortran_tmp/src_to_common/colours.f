      subroutine get_casimirs
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'leg_PDGs.inc'
      integer i
      double precision C_i(nexternal)
      common/ccasimirs/C_i
c
      C_i = 0d0
c
      do i=1,nexternal
        if(abs(leg_pdgs(i)).le.6)then
          C_i(i) = Cf
        elseif(abs(leg_pdgs(i)).eq.21)then
          C_i(i) = Ca
        endif
      enddo
c
      return
      end

      subroutine get_rhos(a,b,parent,ulegs, rho_ab,rho_ba,rhob)
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      integer a,b,parent,ulegs(nexternal-1)
      double precision C_i(nexternal),Cb_i
      common/ccasimirs/C_i
      double precision rho_ab,rho_ba,rhob
c
      rho_ab = 0d0
      rho_ba = 0d0
      rhob   = 0d0
c
      if(abs(ulegs(parent)).le.6)then
        Cb_i = Cf
        rho_ab = (Cb_i + C_i(a) - C_i(b)) / Cb_i
        rho_ba = (Cb_i + C_i(b) - C_i(a)) / Cb_i
        rhob   = (Cb_i - C_i(a) - C_i(b)) / Cb_i
      elseif(abs(ulegs(parent)).eq.21)then
        Cb_i = Ca
        rho_ab = (Cb_i + C_i(a) - C_i(b)) / Cb_i
        rho_ba = (Cb_i + C_i(b) - C_i(a)) / Cb_i
        rhob   = (Cb_i - C_i(a) - C_i(b)) / Cb_i
      endif
c
      return
      end
