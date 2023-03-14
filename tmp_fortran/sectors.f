!!!!! Fortran version of the sectors.py !!!!!!
      subroutine get_Z_NLO(snlo,alpha,isec,jsec,Z_NLO,sector_type,ierr)
      implicit none
      include 'nexternal.inc'
      include 'all_sector_list.inc'
      integer isec,jsec,ierr
      double precision alpha
      double precision snlo(nexternal,nexternal) 
      double precision Z_NLO
      character*1 sector_type
      real * 8 sigma

      
!     build the total sigma (we need to know the number of sectors)

      sigma = 0d0
      
      if(sector_type.eq.'F') then
         do i=1,lensectors
            call getsectorwgt(xs,all_sector_list(i,1),all_sector_list(i,2),wgt)
            sigma = sigma + wgt
            if(all_sector_list(i,2).gt.2) then
               call getsectorwgt(xs,all_sector_list(i,2),all_sector_list(i,1),wgt)
               sigma = sigma + wgt
            endif
         enddo
         call getsectorwgt(xs,isec,jsec,wgt)
         Z_NLO = wgt/sigma
         if(jsec.gt.2) then
            call getsectorwgt(xs,jsec,isec,wgt)
            Z_NLO = Z_NLO + wgt/sigma ! Symmetrization
         endif
      elseif(sector_type.eq.'S') then
         
         
      else
         write(*,*) 'Not allowed value for sector_type: ', sector_type
         write(*,*) 'Exit...'
         stop
      endif
 
      end


      subroutine getsectorwgt(xs,sCM,isec,jsec,wgt)
      implicit none
c$$$          """ the sector function sigma_ij of the Torino paper (eq. 2.13-2.14)
c$$$     without the normalisation.
c$$$     - q is the total momentum of the incoming particles
c$$$     - p_sector is a list of the momenta of the particles defining the sector
c$$$    """
c     integer npsec
      include 'nexternal.inc'
      integer isec,jsec
      real * 8 xs(nexternal,nexternal)
      real * 8 wij,wgt
      real * 8 sCM,sqisec,sqjsec,eisec

      sqisec=0d0
      sqjsec=0d0
      
!build s,sqi,sqj 

      do j=2,nexternal
         sqisec = sqisec + xs(isec,j)
         sqjsec = sqjsec + xs(jsec,j)
      enddo

      eisec = sqisec/sCM
      wij = sCM*xs(isec,jsec)/sqisec/sqjsec
      
      wgt = 1d0/eisec/wij
      
      end
      

      subroutine getsectorwgtS(q, p_sec,wgt)
      implicit none
      real * 8 q(0:3),p_sec(0:3,2)
      real * 8 wgt
      real * 8 s,sij,sqi,sqj,wij
!      real * 8 getsecinv

      call getsecinv(q, p_sec,s,sqi,sqj,sij)
      
      wij =  s*sij/sqi/sqj

      wgt = 1d0/wij
      
      end

      subroutine getsecwgtC(q,p_sec,wgt)
      implicit none
      real * 8 q(0:3),p_sec(0:3,2)
      real * 8 wgt
      real * 8 s,sij,sqi,sqj,wij,ei
!      real * 8 getsecinv

      call getsecinv(q, p_sec,s,sqi,sqj,sij)
      
      ei = sqi/s
      
      wgt = 1d0/ei
      
      
      end

c$$$      subroutine getsecinv(q, p_sec,s,sqi,sqj,sij)
c$$$      implicit none
c$$$      real * 8 q(0:3),p_sec(0:3,2)
c$$$      real * 8 s,sij,sqi,sqj
c$$$      real * 8 dot,momsq
c$$$
c$$$      s = momsq(q)
c$$$      sqi = 2d0 * dot(q(:),p_sec(:,1))
c$$$      sqj = 2d0 * dot(q(:),p_sec(:,2))
c$$$      sij = 2d0 * dot(p_sec(:,1),p_sec(:,2))
c$$$      
c$$$      end

