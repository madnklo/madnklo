!!!!! Fortran version of the sectors.py !!!!!!
      subroutine get_Z_NLO(snlo,sCM,alpha,isec,jsec,Z_NLO,sector_type,ierr)
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

      if(isec.le.2) then
         write(*,*) 'The first particle sector must be in 
     $           the final state!'
         write(*,*) 'isec= ', isec
         write(*,*) 'Exit...'
         stop
      endif

      if(sector_type.eq.'F') then
!     build the total sigma 
         do i=1,lensectors
            call getsectorwgt(xs,sCM,all_sector_list(i,1),all_sector_list(i,2),wgt)
            sigma = sigma + wgt
            if(all_sector_list(i,2).gt.2) then
               call getsectorwgt(xs,sCM,all_sector_list(i,2),all_sector_list(i,1),wgt)
               sigma = sigma + wgt
            endif
         enddo
         call getsectorwgt(xs,sCM,isec,jsec,wgt)
         if(jsec.gt.2) then
            call getsectorwgt(xs,sCM,jsec,isec,wgt)
            Z_NLO = Z_NLO + wgt ! Symmetrization
         endif
      elseif(sector_type.eq.'S') then
         do i=1,lensectors
            if(isec.lt.jsec) then
               call getsectorwgts(xs,sCM,all_sector_list(i,1),all_sector_list(i,2),wgt)
               sigma=sigma+wgt
            else
               call getsectorwgts(xs,sCM,all_sector_list(i,2),all_sector_list(i,1),wgt)
               sigma=sigma+wgt
         enddo
         call getsectorwgts(xs,sCM,isec,jsec,wgt)
      else
         write(*,*) 'Not allowed value for sector_type: ', sector_type
         write(*,*) 'Exit...'
         stop
      endif

      Z_NLO = wgt/sigma
      
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
      

      subroutine getsectorwgts(xs,sCM,isec,jsec,wgt)
      implicit none
      include 'nexternal.inc'
      integer isec,jsec
      real * 8 xs(nexternal,nexternal)
      real * 8 wij,wgt
      real * 8 sCM,sqisec,sqjsec

      sqisec=0d0
      sqjsec=0d0
      
      do j=2,nexternal
         sqisec = sqisec + xs(isec,j)
         sqjsec = sqjsec + xs(jsec,j)
      enddo

      wij = sCM*xs(isec,jsec)/sqisec/sqjsec

      wgt = 1d0/wij
      
      end

c$$$      subroutine getsecwgtC(q,p_sec,wgt)
c$$$      implicit none
c$$$      real * 8 q(0:3),p_sec(0:3,2)
c$$$      real * 8 wgt
c$$$      real * 8 s,sij,sqi,sqj,wij,ei
c$$$!      real * 8 getsecinv
c$$$
c$$$      call getsecinv(q, p_sec,s,sqi,sqj,sij)
c$$$      
c$$$      ei = sqi/s
c$$$      
c$$$      wgt = 1d0/ei
c$$$      
c$$$      
c$$$      end

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

