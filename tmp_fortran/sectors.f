!!!!! Fortran version of the sectors.py !!!!!!

      subroutine getsectorwgt(q,p_sec,wgt)
      implicit none
c$$$          """ the sector function sigma_ij of the Torino paper (eq. 2.13-2.14)
c$$$     without the normalisation.
c$$$     - q is the total momentum of the incoming particles
c$$$     - p_sector is a list of the momenta of the particles defining the sector
c$$$    """
c      integer npsec
      real * 8 q(0:3),p_sec(0:3,2)
      real * 8 s,sqi,sqj,sij,ei,ej,wij
      real * 8 dot,momsq
      real * 8 wgt


      call getsecinv(q, p_sec,s,sqi,sqj,sij)
      
      ei = sqi/s
      ej = sqj/s

      wij = s*sij/sqi/sqj

      wgt = 1d0/ei/wij
      
      end
      

      subroutine getsectorwgtS(q, p_sec,wgt)
      implicit none
      real * 8 q(0:3),p_sec(0:3,2)
      real * 8 wgt
      real * 8 s,sij,sqi,sqj,wij


      call getsecinv(q, p_sec,s,sqi,sqj,sij)
      
      wij =  s*sij/sqi/sqj

      wgt = 1d0/wij
      
      end


      subroutine getsecinv(q, p_sec,s,sqi,sqj,sij)
      implicit none
      real * 8 q(0:3),p_sec(0:3,2)
      real * 8 s,sij,sqi,sqj
      real * 8 dot,momsq

      s = momsq(q)
      sqi = 2d0 * dot(q(:),p_sec(:,1))
      sqj = 2d0 * dot(q(:),p_sec(:,2))
      sij = 2d0 * dot(p_sec(:,1),p_sec(:,2))
      
      end

