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
