      subroutine get_Z_NLO(xs,sCM,alpha,isec,jsec,Z_NLO,sector_type,ierr)
      implicit none
      include 'nexternal.inc'
      include 'all_sector_list.inc'
      include 'leg_PDGs.inc'
      integer i
      integer isec,jsec,ierr
      double precision alpha
      double precision xs(nexternal,nexternal)
      double precision wgt
      double precision Z_NLO
      character*1 sector_type
      real * 8 sigma,sCM

      Z_NLO = 0d0
!     build the total sigma (we need to know the number of sectors)
      sigma = 0d0
      ierr=0d0

      if(isec.le.2) then
         write(*,*) 'get_Z_NLO:'
         write(*,*) 'The first particle sector must be in 
     $           the final state!'
         write(*,*) 'isec= ', isec
         write(*,*) 'Exit...'
         stop
      endif

      if(sector_type.eq.'F') then
!     build the total sigma
         do i=1,lensectors
            call getsectorwgt(xs,sCM,all_sector_list(1,i),all_sector_list(2,i),wgt)
            sigma = sigma + wgt
            if(all_sector_list(2,i).gt.2) then
               call getsectorwgt(xs,sCM,all_sector_list(2,i),all_sector_list(1,i),wgt)
               sigma = sigma + wgt
            endif
         enddo
         call getsectorwgt(xs,sCM,isec,jsec,wgt)
         Z_NLO = wgt
         if(jsec.gt.2) then
            call getsectorwgt(xs,sCM,jsec,isec,wgt)
            Z_NLO = Z_NLO + wgt ! Symmetrization
         endif
      elseif(sector_type.eq.'S') then
         do i=1,lensectors
            if(all_sector_list(1,i).eq.isec) then
               call getsectorwgt_S(xs,sCM,all_sector_list(1,i),all_sector_list(2,i),wgt)
               sigma=sigma+wgt
            elseif(all_sector_list(2,i).eq.isec) then
               call getsectorwgt_S(xs,sCM,all_sector_list(2,i),all_sector_list(1,i),wgt)
               sigma=sigma+wgt
            endif
         enddo
         call getsectorwgt_S(xs,sCM,isec,jsec,wgt)
         Z_NLO = wgt
      else
         write(*,*) 'Not allowed value for sector_type: ', sector_type
         write(*,*) 'Exit...'
         stop
      endif

      Z_NLO = Z_NLO/sigma
      
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
      integer isec,jsec,j
      real * 8 xs(nexternal,nexternal)
      real * 8 wij,wgt
      real * 8 sCM,sqisec,sqjsec,eisec

      sqisec=0d0
      sqjsec=0d0
      wij=0d0
      wgt=0d0
      eisec=0d0
        
!build s,sqi,sqj 

!      do j=3,nexternal
!         sqisec = sqisec + xs(isec,j)
!         sqjsec = sqjsec + xs(jsec,j)
!      enddo


      sqisec = xs(isec,1)+xs(isec,2)
      sqjsec = xs(jsec,1)+xs(jsec,2)
c$$$      write(*,*) sqisec
c$$$      write(*,*) sqjsec
c$$$      write(*,*) sCM
c$$$      write(*,*) xs(isec,jsec)
c$$$      write(*,*)

      eisec = sqisec/sCM
      wij = sCM*xs(isec,jsec)/sqisec/sqjsec
c      if(xs(isec,jsec).lt.1d-10)then
c         write(*,*) 'isec - jsec', isec, jsec
c         write(*,*) 'xs', xs(isec,jsec)
c      endif
c      write(*,*) wij, sCM, sqisec, sqjsec   
      if(wij.lt.1d-15) then
c         write(*,*) 'wij too small in getsectorwgt', isec,jsec, wij, sqisec,sqjsec,sCM,xs(isec,jsec)
         wij=1d-15
      endif
      wgt = 1d0/eisec/wij
c      if(eisec.eq.0d0.or.wij.eq.0d0) write(*,*) eisec,wij
      end
      

      subroutine getsectorwgt_S(xs,sCM,isec,jsec,wgt)
      implicit none
      include 'nexternal.inc'
      integer isec,jsec,j
      real * 8 xs(nexternal,nexternal)
      real * 8 wij,wgt
      real * 8 sCM,sqisec,sqjsec

      sqisec=0d0
      sqjsec=0d0
      wij=0d0
      wgt=0d0
      
      do j=3,nexternal
         sqisec = sqisec + xs(isec,j)
         sqjsec = sqjsec + xs(jsec,j)
      enddo

      wij = sCM*xs(isec,jsec)/sqisec/sqjsec
      if(wij.lt.1d-15) then
c     write(*,*) 'wij too small in getsectorwgt', isec,jsec, wij, sqisec,sqjsec,sCM,xs(isec,jsec)
         wij=1d-15
      endif
      wgt = 1d0/wij
      
      end

c$$$      subroutine getsecwgt_C(q,p_sec,wgt)
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



